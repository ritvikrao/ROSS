/*
 * dragonfly: Charades' models/dragonfly (dragonfly.h, dragonfly.C,
 * terminal.C, router.C; itself derived from the CODES dragonfly model) for
 * ROSS, written to match it draw for draw, so both engines commit the same
 * events. Charades RESULTS.md section 171.
 *
 * Topology: groups of a routers (--routers-per-group), fully connected by
 * local channels; h (--global-channels, even) global channels per router and
 * p (--terminals-per-router) terminals per router; a*h + 1 groups. LP gids:
 * routers [0, R), terminals [R, R + T), MPI processes [R + T, R + 2T), with
 * R = a * groups, T = p * R. Each MPI process generates messages
 * (--traffic, --arrive-rate, and --gen-time or a fixed count), a terminal
 * breaks each packet into 16 chunks, and routers forward them with credit
 * flow control, a waiting list for a full output channel, and minimal or
 * non-minimal (Valiant) routing.
 *
 * Time is integer picoseconds (1 unit = 1 ps): every delay is
 * llround(1000 * d ns), at least 100 (0.1 ns). ROSS's double timestamps hold
 * these integers exactly.
 *
 * Equal timestamps: as in models/traffic, each event's ROSS user priority is
 * (sender gid * 2^32 + the sender's creation count) / 2^52, which orders
 * equal-timestamp events at an LP the way Charades does (sending LP, then the
 * order the sender created them; Charades numbers an event when it is
 * created, including one past the end time).
 *
 * Mapping (--mapping=1, the default): each rank holds a block of consecutive
 * routers with their terminals and MPI processes, as Charades' mapper puts a
 * router with its terminals and processes on one chare. With R routers over n
 * ranks the first R mod n ranks hold floor(R/n) + 1 routers and the rest
 * floor(R/n): Charm++'s default block map of Charades' chares at one router
 * per chare (--routers-per-chare=1). --mapping=0 is ROSS's linear mapping
 * (the LP count must divide).
 *
 * Statistics: the same totals as Charades' report_stats(), summed over ranks.
 */
#include <ross.h>
#include "check-revent/crv-state.h"
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Fixed parameters of the original model (ns, bytes, chunks) */
#define MESSAGE_SIZE 512.0
#define PACKET_SIZE 512.0
#define CHUNK_SIZE 32.0
#define GLOBAL_BANDWIDTH 4.7
#define LOCAL_BANDWIDTH 5.25
#define NODE_BANDWIDTH 5.25
#define CREDIT_SIZE 8
#define INJECTION_INTERVAL 3000
#define RESCHEDULE_DELAY 1.0
#define MEAN_PROCESS 1.0
#define NUM_VC 1
#define ROUTER_WAITING_PACK_COUNT (1 << 8)
#define TIME_PER_NS 1000.0
#define MIN_DELAY ((uint64_t)100)

enum { GLOBAL, LOCAL, TERMINAL };
enum { MINIMAL, NON_MINIMAL, ADAPTIVE };
enum { UNIFORM_RANDOM = 1, WORST_CASE, TRANSPOSE, NEAREST_NEIGHBOR, BISECTION, NEAREST_ROUTER };
enum { T_GENERATE, T_SEND, R_ARRIVE, R_SEND, R_WAIT, T_ARRIVE, CREDIT, FINISH, MPI_SEND };

/* Configuration (command line); defaults are the original #defines */
static unsigned int a_routers = 24, h_global = 12, p_terms = 12;
static unsigned int routing = MINIMAL, traffic = NEAREST_NEIGHBOR;
static double mean_interval = 200.0;
static unsigned int local_vc_size = 32768, global_vc_size = 65536, terminal_vc_size = 16384;
static double gen_time = 0.0;
static double lookahead = 100.0;
static unsigned int use_mapping = 1;
static double mult = 20.0;
static unsigned int optimistic_memory = 100000;

/* Derived */
static uint32_t num_groups, total_routers, total_terminals, total_procs, radix;
static int num_packets, num_chunks, max_packets;
static uint64_t gen_end;
static uint64_t rt_q, rt_rem; /* routers per rank: rt_q, plus one on the first rt_rem ranks */

static inline uint64_t to_time(double ns) {
  uint64_t t = (uint64_t)llround(ns * TIME_PER_NS);
  return t < MIN_DELAY ? MIN_DELAY : t;
}
static inline uint64_t now_of(tw_lp *lp) { return (uint64_t)tw_now(lp); }

static inline int is_router(uint64_t gid) { return gid < total_routers; }
static inline uint32_t terminal_index(uint64_t gid) { return (uint32_t)(gid - total_routers); }
static inline uint32_t terminal_router(uint64_t gid) { return terminal_index(gid) / p_terms; }
static inline uint32_t proc_index(uint64_t gid) { return (uint32_t)(gid - total_routers - total_terminals); }

/* One message type for all events; saved_* are written for reverse computation */
typedef struct {
  int type;
  uint32_t dest_terminal, travel_start, prev_lp, saved_served;
  int16_t intm_group;
  uint16_t chan, saved_chan;   /* chan: a WAIT's channel, or a CREDIT's vc */
  uint8_t hops, chunk_id, last_hop;
  uint64_t saved_time;         /* saved_available_time / _credit_time / _max_latency */
} df_msg;

typedef struct {
  uint32_t travel_start, dest_terminal;
  int16_t intm_group;
  uint16_t chan;
  uint8_t hops, chunk_id, last_hop, served, retired;
} wait_entry;

typedef struct {
  uint64_t sends;              /* events created: the tie-break sequence */
  /* router */
  uint32_t group_id;
  uint32_t *global_channels;
  uint64_t *next_output_available_time, *next_credit_available_time;
  int32_t *vc_occupancy;
  wait_entry *waiting;         /* ring of capacity wait_cap: entry i is waiting[(wait_head + i) % wait_cap] */
  uint32_t wait_cap, wait_head, wait_len, wait_base, wait_count;
  uint64_t dropped_packets, wait_drops, waits;
  uint8_t drop_reported, wait_drop_reported;
  /* terminal */
  uint32_t router_id;
  int32_t t_vc_occupancy[NUM_VC];
  uint64_t terminal_available_time, t_next_credit_available_time;
  uint64_t packets_injected, packets_finished, total_latency, total_hops, max_latency;
  /* MPI process */
  int32_t message_counter;
  uint64_t available_time;
  uint32_t p_router_id, p_group_id;
  int32_t row, col;
} df_state;

/* Creates an event with Charades' tie order; the caller sends it */
static tw_event *df_new(df_state *s, tw_lp *lp, tw_lpid dest, uint64_t offset, int type) {
  double prio = ldexp((double)lp->gid * 4294967296.0 + (double)(s->sends++), -52);
  tw_event *e = tw_event_new_user_prio(dest, (double)offset, lp, prio);
  df_msg *m = (df_msg *)tw_event_data(e);
  memset(m, 0, sizeof *m);
  m->type = type;
  return e;
}

/* ---------------------------------------------------------------- topology */

static uint32_t global_channel(uint32_t r, uint32_t i) {
  int64_t a = a_routers;
  int64_t offset = (int64_t)(r % a) * (h_global / 2) + 1 + i / 2;
  int64_t c;
  if (i % 2 != 0) c = ((int64_t)r + offset * a) % total_routers;
  else c = (int64_t)r - offset * a;
  if (c < 0) c += total_routers;
  return (uint32_t)c;
}

static uint32_t router_from_group_id(int64_t gid, int64_t group_id) {
  int64_t a = a_routers;
  int64_t group_begin = group_id * a;
  int64_t group_end = group_id * a + a - 1;
  int64_t offset = (gid * a - group_begin) / a;
  if ((gid * a) < group_begin) offset = (group_begin - gid * a) / a;
  int64_t half_channel = h_global / 2;
  int64_t index = (offset - 1) / (half_channel * a);
  offset = (offset - 1) % (half_channel * a);
  if (index % 2 != 0) return (uint32_t)(group_end - (offset / half_channel));
  return (uint32_t)(group_begin + (offset / half_channel));
}

/* ------------------------------------------------------------------ mapping */

static uint64_t rt_first(tw_peid pe) { return pe * rt_q + (pe < rt_rem ? pe : rt_rem); }
static uint64_t rt_count(tw_peid pe) { return rt_q + (pe < rt_rem ? 1 : 0); }
static tw_peid rt_owner(uint64_t router) {
  uint64_t b = rt_rem * (rt_q + 1);
  return (tw_peid)(router < b ? router / (rt_q + 1) : rt_rem + (router - b) / rt_q);
}

static tw_peid df_map(tw_lpid gid) {
  if (!use_mapping) return (tw_peid)(gid / g_tw_nlp);
  if (gid < total_routers) return rt_owner(gid);
  gid -= total_routers;
  if (gid < total_terminals) return rt_owner(gid / p_terms);
  return rt_owner((gid - total_terminals) / p_terms);
}

/* Local index on its rank of a gid: routers, then their terminals, then their
 * processes */
static tw_lpid df_local(tw_lpid gid) {
  if (gid < total_routers) return gid - rt_first(rt_owner(gid));
  gid -= total_routers;
  int term = gid < total_terminals;
  if (!term) gid -= total_terminals;
  tw_peid k = rt_owner(gid / p_terms);
  uint64_t rp = rt_count(k), off = gid - rt_first(k) * p_terms;
  return term ? rp + off : rp + rp * p_terms + off;
}

static tw_lpid df_global(tw_peid pe, tw_lpid local) {
  uint64_t rp = rt_count(pe), tp = rp * p_terms, first = rt_first(pe);
  if (local < rp) return first + local;
  if (local < rp + tp) return total_routers + first * p_terms + (local - rp);
  return (tw_lpid)total_routers + total_terminals + first * p_terms + (local - rp - tp);
}

static void df_custom_mapping(void) {
  unsigned int nlp_per_kp = (unsigned int)ceil((double)g_tw_nlp / (double)g_tw_nkp);
  tw_lpid lpid = 0;
  for (tw_kpid kpid = 0; kpid < g_tw_nkp; kpid++) {
    tw_kp_onpe(kpid, g_tw_pe);
    for (unsigned int j = 0; j < nlp_per_kp && lpid < g_tw_nlp; j++, lpid++) {
      tw_lp_onpe(lpid, g_tw_pe, df_global(g_tw_mynode, lpid));
      tw_lp_onkp(g_tw_lp[lpid], g_tw_kp[kpid]);
    }
  }
}

static tw_lp *df_global_to_local(tw_lpid gid) { return g_tw_lp[df_local(gid)]; }

/* ---------------------------------------------------------------- statistics */

static uint64_t st_messages, st_injected, st_finished, st_latency, st_hops, st_max_latency, st_dropped, st_waits;

/* ---------------------------------------------------------------- terminal */

static void term_init(df_state *s, tw_lp *lp) {
  memset(s, 0, sizeof *s);
  s->router_id = terminal_router(lp->gid);
}

static void term_event(df_state *s, tw_bf *bf, df_msg *msg, tw_lp *lp) {
  uint64_t now = now_of(lp);
  switch (msg->type) {
    case T_GENERATE:
      for (int i = 0; i < num_chunks; i++) {
        double ts = 0.1 + tw_rand_exponential(lp->rng, mean_interval / 200);
        int chan = -1;
        for (int j = 0; j < NUM_VC; j++) {
          if (s->t_vc_occupancy[j] < (int64_t)terminal_vc_size * num_chunks) { chan = j; break; }
        }
        if (chan == -1) tw_error(TW_LOC, "Dragonfly: terminal injection queue size exceeded");
        tw_event *e = df_new(s, lp, lp->gid, to_time(i + ts), T_SEND);
        df_msg *m = (df_msg *)tw_event_data(e);
        m->dest_terminal = msg->dest_terminal;
        m->travel_start = msg->travel_start;
        m->chunk_id = (uint8_t)i;
        tw_event_send(e);
      }
      break;
    case T_SEND: {
      bf->c3 = 0;
      if (msg->chunk_id == num_chunks - 1) { bf->c3 = 1; s->packets_injected++; }
      msg->saved_time = s->terminal_available_time;
      int head_delay = (1 / NODE_BANDWIDTH) * CHUNK_SIZE;
      double ts = head_delay + tw_rand_exponential(lp->rng, (double)head_delay / 200);
      if (s->terminal_available_time < now) s->terminal_available_time = now;
      s->terminal_available_time += to_time(ts);
      tw_event *e = df_new(s, lp, s->router_id, s->terminal_available_time - now, R_ARRIVE);
      df_msg *m = (df_msg *)tw_event_data(e);
      m->travel_start = msg->travel_start;
      m->dest_terminal = msg->dest_terminal;
      m->prev_lp = (uint32_t)lp->gid;
      m->intm_group = -1;
      m->hops = 0;
      m->chunk_id = msg->chunk_id;
      tw_event_send(e);
      s->t_vc_occupancy[0]++;
    } break;
    case T_ARRIVE: {
      bf->c1 = 0;
      if (msg->chunk_id == num_chunks - 1) {
        bf->c1 = 1;
        double ts = 0.1 + tw_rand_exponential(lp->rng, MEAN_PROCESS / 200);
        tw_event *e = df_new(s, lp, lp->gid, to_time(ts), FINISH);
        df_msg *m = (df_msg *)tw_event_data(e);
        m->travel_start = msg->travel_start;
        m->hops = msg->hops;
        tw_event_send(e);
      }
      int credit_delay = (1 / NODE_BANDWIDTH) * CREDIT_SIZE;
      double ts = credit_delay + tw_rand_exponential(lp->rng, credit_delay / 1000);  /* integer division: 0 */
      msg->saved_time = s->t_next_credit_available_time;
      if (s->t_next_credit_available_time < now) s->t_next_credit_available_time = now;
      s->t_next_credit_available_time += to_time(ts);
      tw_event *b = df_new(s, lp, s->router_id, s->t_next_credit_available_time - now, CREDIT);
      ((df_msg *)tw_event_data(b))->chan =
          (uint16_t)((a_routers + h_global + terminal_index(lp->gid) % p_terms) * NUM_VC);
      tw_event_send(b);
    } break;
    case CREDIT:
      if (msg->chan != 0) tw_error(TW_LOC, "Dragonfly: invalid terminal VC index");
      s->t_vc_occupancy[msg->chan]--;
      break;
    case FINISH: {
      bf->c3 = 0;
      uint64_t latency = now - msg->travel_start;
      s->packets_finished++;
      s->total_latency += latency;
      s->total_hops += msg->hops;
      if (s->max_latency < latency) {
        bf->c3 = 1;
        msg->saved_time = s->max_latency;
        s->max_latency = latency;
      }
    } break;
    default: tw_error(TW_LOC, "Dragonfly: terminal got event type %d", msg->type);
  }
}

static void term_event_rc(df_state *s, tw_bf *bf, df_msg *msg, tw_lp *lp) {
  switch (msg->type) {
    case T_GENERATE:
      for (int i = 0; i < num_chunks; i++) { tw_rand_reverse_unif(lp->rng); s->sends--; }
      break;
    case T_SEND:
      if (bf->c3) s->packets_injected--;
      s->terminal_available_time = msg->saved_time;
      tw_rand_reverse_unif(lp->rng);
      s->t_vc_occupancy[0]--;
      s->sends--;
      break;
    case T_ARRIVE:
      if (bf->c1) { tw_rand_reverse_unif(lp->rng); s->sends--; }
      tw_rand_reverse_unif(lp->rng);
      s->t_next_credit_available_time = msg->saved_time;
      s->sends--;
      break;
    case CREDIT:
      s->t_vc_occupancy[msg->chan]++;
      break;
    case FINISH: {
      uint64_t latency = now_of(lp) - msg->travel_start;
      s->packets_finished--;
      s->total_latency -= latency;
      s->total_hops -= msg->hops;
      if (bf->c3) s->max_latency = msg->saved_time;
    } break;
  }
}

static void term_final(df_state *s, tw_lp *lp) {
  (void)lp;
  st_injected += s->packets_injected;
  st_finished += s->packets_finished;
  st_latency += s->total_latency;
  st_hops += s->total_hops;
  if (s->max_latency > st_max_latency) st_max_latency = s->max_latency;
}

/* ------------------------------------------------------------- MPI process */

static void proc_init(df_state *s, tw_lp *lp) {
  memset(s, 0, sizeof *s);
  uint32_t a = a_routers, p = p_terms;
  s->message_counter = 0;
  s->p_router_id = proc_index(lp->gid) / p;
  s->p_group_id = proc_index(lp->gid) / (p * a);
  s->row = proc_index(lp->gid) / (a * p);
  s->col = proc_index(lp->gid) % (a * p + 1);
  double ts = tw_rand_exponential(lp->rng, mean_interval);
  tw_event_send(df_new(s, lp, lp->gid, to_time(ts), MPI_SEND));
}

static inline uint64_t wrap_terminal(int64_t dst) {
  int64_t R = total_routers, T = total_terminals;
  return (uint64_t)(R + (((dst - R) % T) + T) % T);
}

static void proc_event(df_state *s, tw_bf *bf, df_msg *msg, tw_lp *lp) {
  int64_t a = a_routers, p = p_terms;
  int64_t R = total_routers, T = total_terminals, G = num_groups;
  int64_t offset = p * a - 1;
  int64_t dst_lp = -1;
  uint64_t now = now_of(lp);
  bf->c1 = bf->c2 = bf->c3 = 0;
  if (gen_end ? now >= gen_end : s->message_counter >= max_packets) { bf->c1 = 1; return; }
  switch (traffic) {
    case BISECTION: {
      bf->c3 = 1; bf->c2 = 1;
      int64_t grp = tw_rand_integer(lp->rng, 0, G - 1);
      if (grp == s->p_group_id) grp = (s->p_group_id + (G / 2)) % G;
      int64_t next_group_begin = grp * p * a;
      dst_lp = tw_rand_integer(lp->rng, R + next_group_begin, R + next_group_begin + offset);
    } break;
    case WORST_CASE: {
      bf->c3 = 1;
      int64_t next_group_begin = ((s->p_group_id + 1) % G) * p * a;
      dst_lp = tw_rand_integer(lp->rng, R + next_group_begin, R + next_group_begin + offset);
    } break;
    case UNIFORM_RANDOM:
      bf->c3 = 1;
      dst_lp = tw_rand_integer(lp->rng, R, R + T - 1);
      break;
    case TRANSPOSE:
      if (s->col == s->row) { bf->c1 = 1; return; }
      dst_lp = wrap_terminal(R + ((int64_t)s->col * (a * p) + s->row));
      break;
    case NEAREST_NEIGHBOR: {
      int64_t local_offset = proc_index(lp->gid) % p;
      int64_t group_begin = R + ((int64_t)s->p_router_id * p);
      dst_lp = wrap_terminal(group_begin + (p - local_offset));
    } break;
    case NEAREST_ROUTER: {
      int64_t local_offset = proc_index(lp->gid) % p;
      int64_t next_router_begin = R + (((s->p_router_id + 1) % R) * p);
      dst_lp = wrap_terminal(next_router_begin + (p - local_offset));
    } break;
    default: tw_error(TW_LOC, "Dragonfly: undefined traffic pattern");
  }
  if (now > UINT32_MAX) tw_error(TW_LOC, "Dragonfly: packet generated after 2^32 ps");
  msg->saved_time = s->available_time;
  for (int i = 0; i < num_packets; i++) {
    if (s->available_time < now) s->available_time = now;
    double ts = 0.1 + tw_rand_exponential(lp->rng, mean_interval / 100);
    s->available_time += to_time(ts);
    tw_event *e = df_new(s, lp, lp->gid - total_procs, s->available_time - now, T_GENERATE);
    df_msg *m = (df_msg *)tw_event_data(e);
    m->dest_terminal = (uint32_t)dst_lp;
    m->travel_start = (uint32_t)now;
    s->message_counter++;
    tw_event_send(e);
  }
  double ts = mean_interval + tw_rand_exponential(lp->rng, mean_interval / 100);
  tw_event_send(df_new(s, lp, lp->gid, to_time(ts), MPI_SEND));
}

static void proc_event_rc(df_state *s, tw_bf *bf, df_msg *msg, tw_lp *lp) {
  if (bf->c1) return;
  if (bf->c3) tw_rand_reverse_unif(lp->rng);
  if (bf->c2) tw_rand_reverse_unif(lp->rng);
  for (int i = 0; i < num_packets; i++) {
    tw_rand_reverse_unif(lp->rng);
    s->message_counter--;
    s->sends--;
  }
  s->available_time = msg->saved_time;
  tw_rand_reverse_unif(lp->rng);
  s->sends--;
}

static void proc_final(df_state *s, tw_lp *lp) { (void)lp; st_messages += (uint64_t)s->message_counter; }

/* ------------------------------------------------------------------ router */

static wait_entry *wait_at(df_state *s, uint32_t i) { return &s->waiting[(s->wait_head + i) % s->wait_cap]; }

static void wait_push_back(df_state *s, const wait_entry *w) {
  if (s->wait_len == s->wait_cap) {  /* grow: unroll the ring into a bigger array */
    uint32_t cap = s->wait_cap ? 2 * s->wait_cap : 16;
    wait_entry *n = (wait_entry *)malloc(sizeof(wait_entry) * cap);
    for (uint32_t i = 0; i < s->wait_len; i++) n[i] = *wait_at(s, i);
    free(s->waiting);
    s->waiting = n; s->wait_cap = cap; s->wait_head = 0;
  }
  *wait_at(s, s->wait_len++) = *w;
}

static void router_init(df_state *s, tw_lp *lp) {
  memset(s, 0, sizeof *s);
  s->group_id = (uint32_t)(lp->gid / a_routers);
  s->global_channels = (uint32_t *)malloc(sizeof(uint32_t) * h_global);
  for (uint32_t i = 0; i < h_global; i++) s->global_channels[i] = global_channel((uint32_t)lp->gid, i);
  s->next_output_available_time = (uint64_t *)calloc(radix, sizeof(uint64_t));
  s->next_credit_available_time = (uint64_t *)calloc(radix, sizeof(uint64_t));
  s->vc_occupancy = (int32_t *)calloc(radix * NUM_VC, sizeof(int32_t));
}

static uint8_t last_hop_from(df_state *s, uint32_t prev_lp) {
  if (!is_router(prev_lp)) return TERMINAL;
  return prev_lp / a_routers == s->group_id ? LOCAL : GLOBAL;
}

static uint32_t sender_chan(df_state *s, tw_lp *lp, uint32_t prev_lp) {
  if (!is_router(prev_lp)) return 0;
  if (prev_lp / a_routers == s->group_id) return (uint32_t)((lp->gid % a_routers) * NUM_VC);
  for (uint32_t i = 0; i < h_global; i++) {
    if (global_channel(prev_lp, i) == lp->gid) return (a_routers + i) * NUM_VC;
  }
  tw_error(TW_LOC, "Dragonfly: chunk from a router with no channel to this one");
  return 0;
}

static uint32_t next_stop(df_state *s, tw_lp *lp, df_msg *msg, tw_bf *bf, int16_t *intm) {
  uint32_t a = a_routers;
  uint32_t dest_router_id = terminal_router(msg->dest_terminal);
  *intm = msg->intm_group;
  if (dest_router_id == lp->gid) return msg->dest_terminal;
  if (msg->last_hop == TERMINAL && routing == NON_MINIMAL) {
    if (dest_router_id / a != s->group_id) {
      bf->c2 = 1;
      *intm = (int16_t)tw_rand_integer(lp->rng, 0, num_groups - 1);
    }
  }
  if (*intm == (int16_t)s->group_id) *intm = -1;
  uint32_t dest_group_id = *intm >= 0 ? (uint32_t)*intm : dest_router_id / a;
  if (s->group_id == dest_group_id) return dest_router_id;
  uint32_t dest_lp = router_from_group_id(dest_group_id, s->group_id);
  if (dest_lp == lp->gid) {
    for (uint32_t i = 0; i < h_global; i++) {
      if (s->global_channels[i] / a == dest_group_id) dest_lp = s->global_channels[i];
    }
  }
  return dest_lp;
}

static uint32_t output_port(df_state *s, uint32_t next, uint32_t dest_terminal) {
  uint32_t a = a_routers;
  if (next == dest_terminal) return a + h_global + (terminal_index(dest_terminal) % p_terms);
  if (next / a != s->group_id) {
    int port = -1;
    for (uint32_t i = 0; i < h_global; i++) if (s->global_channels[i] == next) port = (int)(a + i);
    if (port < 0) tw_error(TW_LOC, "Dragonfly: no global channel to the next stop");
    return (uint32_t)port;
  }
  return next % a;
}

static void router_event(df_state *s, tw_bf *bf, df_msg *msg, tw_lp *lp) {
  uint64_t now = now_of(lp);
  uint32_t a = a_routers;
  switch (msg->type) {
    case R_ARRIVE: {
      uint8_t last_hop = last_hop_from(s, msg->prev_lp);
      double ts = 0.1 + tw_rand_exponential(lp->rng, mean_interval / 200);
      tw_event *e = df_new(s, lp, lp->gid, to_time(ts), R_SEND);
      df_msg *m = (df_msg *)tw_event_data(e);
      m->travel_start = msg->travel_start;
      m->dest_terminal = msg->dest_terminal;
      m->intm_group = msg->intm_group;
      m->hops = msg->hops + 1;
      m->chunk_id = msg->chunk_id;
      m->last_hop = last_hop;
      int credit_delay;
      if (last_hop == TERMINAL) credit_delay = (1 / NODE_BANDWIDTH) * CREDIT_SIZE;
      else if (last_hop == GLOBAL) credit_delay = (1 / GLOBAL_BANDWIDTH) * CREDIT_SIZE;
      else credit_delay = (1 / LOCAL_BANDWIDTH) * CREDIT_SIZE;
      uint32_t chan = sender_chan(s, lp, msg->prev_lp);
      uint32_t port = chan / NUM_VC;
      msg->saved_time = s->next_credit_available_time[port];
      if (s->next_credit_available_time[port] < now) s->next_credit_available_time[port] = now;
      double cts = credit_delay + tw_rand_exponential(lp->rng, (double)credit_delay / 1000);
      s->next_credit_available_time[port] += to_time(cts);
      tw_event *b = df_new(s, lp, msg->prev_lp, s->next_credit_available_time[port] - now, CREDIT);
      ((df_msg *)tw_event_data(b))->chan = (uint16_t)chan;
      tw_event_send(b);
      tw_event_send(e);
    } break;
    case R_SEND: {
      bf->c2 = bf->c3 = bf->c4 = 0;
      int16_t intm;
      uint32_t next = next_stop(s, lp, msg, bf, &intm);
      uint32_t port = output_port(s, next, msg->dest_terminal);
      uint32_t output_chan = port * NUM_VC;
      float bandwidth = LOCAL_BANDWIDTH;
      uint32_t buf_size = local_vc_size;
      if (port >= a && port < a + h_global) { bandwidth = GLOBAL_BANDWIDTH; buf_size = global_vc_size; }
      if (port >= a + h_global) buf_size = terminal_vc_size;
      if (s->vc_occupancy[output_chan] >= (int32_t)buf_size) {
        if (msg->last_hop == TERMINAL) {
          bf->c4 = 1;
          s->dropped_packets++;
        } else {
          bf->c3 = 1;
          double ts = 0.1 + tw_rand_exponential(lp->rng, RESCHEDULE_DELAY / 10);
          tw_event *e = df_new(s, lp, lp->gid, to_time(ts), R_WAIT);
          df_msg *m = (df_msg *)tw_event_data(e);
          m->travel_start = msg->travel_start;
          m->dest_terminal = msg->dest_terminal;
          m->intm_group = intm;
          m->chan = (uint16_t)output_chan;
          m->hops = msg->hops;
          m->chunk_id = msg->chunk_id;
          m->last_hop = msg->last_hop;
          tw_event_send(e);
        }
        return;
      }
      msg->saved_time = s->next_output_available_time[port];
      msg->saved_chan = (uint16_t)output_chan;
      double ts = ((1 / bandwidth) * CHUNK_SIZE) + tw_rand_exponential(lp->rng, (double)CHUNK_SIZE / 200);
      uint64_t base = s->next_output_available_time[port] > now ? s->next_output_available_time[port] : now;
      s->next_output_available_time[port] = base + to_time(ts);
      uint64_t off = s->next_output_available_time[port] - now;
      if (next == msg->dest_terminal) {
        tw_event *e = df_new(s, lp, next, off, T_ARRIVE);
        df_msg *m = (df_msg *)tw_event_data(e);
        m->travel_start = msg->travel_start;
        m->hops = msg->hops;
        m->chunk_id = msg->chunk_id;
        tw_event_send(e);
      } else {
        tw_event *e = df_new(s, lp, next, off, R_ARRIVE);
        df_msg *m = (df_msg *)tw_event_data(e);
        m->travel_start = msg->travel_start;
        m->dest_terminal = msg->dest_terminal;
        m->prev_lp = (uint32_t)lp->gid;
        m->intm_group = intm;
        m->hops = msg->hops;
        m->chunk_id = msg->chunk_id;
        tw_event_send(e);
      }
      s->vc_occupancy[output_chan]++;
    } break;
    case R_WAIT: {
      bf->c3 = 0;
      if (s->wait_count >= ROUTER_WAITING_PACK_COUNT) { bf->c3 = 1; s->wait_drops++; return; }
      wait_entry w;
      memset(&w, 0, sizeof w);
      w.travel_start = msg->travel_start;
      w.dest_terminal = msg->dest_terminal;
      w.intm_group = msg->intm_group;
      w.chan = msg->chan;
      w.hops = msg->hops;
      w.chunk_id = msg->chunk_id;
      w.last_hop = msg->last_hop;
      wait_push_back(s, &w);
      s->wait_count++;
      s->waits++;
    } break;
    case CREDIT: {
      bf->c3 = 0;
      uint32_t chan = msg->chan;
      s->vc_occupancy[chan]--;
      if (s->wait_count == 0) return;
      for (uint32_t j = 0; j < s->wait_len; j++) {
        wait_entry *w = wait_at(s, j);
        if (w->served || w->chan != chan) continue;
        bf->c3 = 1;
        double ts = tw_rand_exponential(lp->rng, mean_interval / 1000);
        tw_event *e = df_new(s, lp, lp->gid, to_time(ts), R_SEND);
        df_msg *m = (df_msg *)tw_event_data(e);
        m->travel_start = w->travel_start;
        m->dest_terminal = w->dest_terminal;
        m->intm_group = w->intm_group;
        m->hops = w->hops;
        m->chunk_id = w->chunk_id;
        m->last_hop = w->last_hop;
        tw_event_send(e);
        w->served = 1;
        s->wait_count--;
        msg->saved_served = s->wait_base + j;
        break;
      }
    } break;
    default: tw_error(TW_LOC, "Dragonfly: router got event type %d", msg->type);
  }
}

static void router_event_rc(df_state *s, tw_bf *bf, df_msg *msg, tw_lp *lp) {
  switch (msg->type) {
    case R_ARRIVE: {
      tw_rand_reverse_unif(lp->rng);
      tw_rand_reverse_unif(lp->rng);
      uint32_t port = sender_chan(s, lp, msg->prev_lp) / NUM_VC;
      s->next_credit_available_time[port] = msg->saved_time;
      s->sends -= 2;
    } break;
    case R_SEND:
      if (bf->c2) tw_rand_reverse_unif(lp->rng);
      if (bf->c4) { s->dropped_packets--; return; }
      tw_rand_reverse_unif(lp->rng);
      s->sends--;
      if (bf->c3) return;
      s->next_output_available_time[msg->saved_chan / NUM_VC] = msg->saved_time;
      s->vc_occupancy[msg->saved_chan]--;
      break;
    case R_WAIT:
      if (bf->c3) { s->wait_drops--; return; }
      s->wait_len--;
      s->wait_count--;
      s->waits--;
      break;
    case CREDIT:
      s->vc_occupancy[msg->chan]++;
      if (bf->c3) {
        tw_rand_reverse_unif(lp->rng);
        wait_at(s, msg->saved_served - s->wait_base)->served = 0;
        s->wait_count++;
        s->sends--;
      }
      break;
  }
}

static void router_commit(df_state *s, tw_bf *bf, df_msg *msg, tw_lp *lp) {
  switch (msg->type) {
    case R_SEND:
      if (bf->c4 && !s->drop_reported) {
        s->drop_reported = 1;
        printf("Router %llu buffers overflowed from incoming terminals: chunk dropped (further drops not reported)\n",
               (unsigned long long)lp->gid);
      }
      break;
    case R_WAIT:
      if (bf->c3 && !s->wait_drop_reported) {
        s->wait_drop_reported = 1;
        printf("Router %llu reached maximum count of waiting packets %d: chunk dropped (further drops not reported)\n",
               (unsigned long long)lp->gid, ROUTER_WAITING_PACK_COUNT);
      }
      break;
    case CREDIT:
      if (bf->c3) {
        wait_at(s, msg->saved_served - s->wait_base)->retired = 1;
        while (s->wait_len && wait_at(s, 0)->retired) {
          s->wait_head = (s->wait_head + 1) % s->wait_cap;
          s->wait_len--;
          s->wait_base++;
        }
      }
      break;
  }
}

static void router_final(df_state *s, tw_lp *lp) {
  (void)lp;
  st_dropped += s->dropped_packets + s->wait_drops;
  st_waits += s->waits;
}

/*
 * --synch=6 (ROSS's rollback check) compares an LP's state before an event
 * and after the event and its reverse. A router's arrays live outside the
 * state struct, and its waiting list may grow (a new array) inside a forward
 * handler, so the default byte comparison would check pointers. This
 * checkpointer copies and compares the router's logical state instead:
 * scalars, the per-port arrays, and the waiting entries in order.
 */
/* The checkpoint is a df_state whose pointers hold deep copies (the waiting
 * list unrolled: head 0, capacity = length). It must not be larger than the
 * state: ROSS's rollback-check scheduler sizes its checkpoint buffers from
 * the LP state sizes before it reads sz_storage. */
static void router_save(void *into, void const *from) {
  df_state *c = (df_state *)into;
  df_state const *s = (df_state const *)from;
  *c = *s;
  c->next_output_available_time = (uint64_t *)malloc(sizeof(uint64_t) * radix);
  c->next_credit_available_time = (uint64_t *)malloc(sizeof(uint64_t) * radix);
  c->vc_occupancy = (int32_t *)malloc(sizeof(int32_t) * radix * NUM_VC);
  memcpy(c->next_output_available_time, s->next_output_available_time, sizeof(uint64_t) * radix);
  memcpy(c->next_credit_available_time, s->next_credit_available_time, sizeof(uint64_t) * radix);
  memcpy(c->vc_occupancy, s->vc_occupancy, sizeof(int32_t) * radix * NUM_VC);
  c->waiting = (wait_entry *)malloc(sizeof(wait_entry) * (s->wait_len ? s->wait_len : 1));
  for (uint32_t i = 0; i < s->wait_len; i++) c->waiting[i] = s->waiting[(s->wait_head + i) % s->wait_cap];
  c->wait_head = 0;
  c->wait_cap = s->wait_len;
}

static void router_clean(void *state) {
  df_state *c = (df_state *)state;
  free(c->next_output_available_time); free(c->next_credit_available_time);
  free(c->vc_occupancy); free(c->waiting);
}

/*
 * ROSS passes the checkpoint first, then the live state. It checks twice:
 * before an event against after it and its reverse, and after the event
 * against after the reverse, the event again and its commit. Commit handlers
 * here change state on purpose: the print-once drop flags, and a committed
 * credit retires the entry it released and pops retired entries off the
 * front of the waiting list (as in Charades). So those flags, the retired
 * marks and the list's front are not compared; entries are compared by
 * absolute index from the later of the two fronts.
 */
static bool router_check(void *ckpt, void const *live) {
  df_state const *b = (df_state const *)ckpt, *s = (df_state const *)live;
  if (b->sends != s->sends || b->group_id != s->group_id || b->wait_count != s->wait_count ||
      b->wait_base + b->wait_len != s->wait_base + s->wait_len ||
      b->dropped_packets != s->dropped_packets || b->wait_drops != s->wait_drops || b->waits != s->waits)
    return false;
  if (memcmp(b->next_output_available_time, s->next_output_available_time, sizeof(uint64_t) * radix) ||
      memcmp(b->next_credit_available_time, s->next_credit_available_time, sizeof(uint64_t) * radix) ||
      memcmp(b->vc_occupancy, s->vc_occupancy, sizeof(int32_t) * radix * NUM_VC))
    return false;
  uint32_t from = b->wait_base > s->wait_base ? b->wait_base : s->wait_base;
  for (uint32_t k = from; k < s->wait_base + s->wait_len; k++) {
    wait_entry x = b->waiting[(b->wait_head + (k - b->wait_base)) % (b->wait_cap ? b->wait_cap : 1)];
    wait_entry y = s->waiting[(s->wait_head + (k - s->wait_base)) % s->wait_cap];
    x.retired = y.retired = 0;
    if (memcmp(&x, &y, sizeof(wait_entry))) return false;
  }
  return true;
}

static void no_commit(df_state *s, tw_bf *bf, df_msg *m, tw_lp *lp) { (void)s; (void)bf; (void)m; (void)lp; }

static tw_lptype df_lps[] = {
  {(init_f)router_init, (pre_run_f)NULL, (event_f)router_event, (revent_f)router_event_rc,
   (commit_f)router_commit, (final_f)router_final, (map_f)df_map, sizeof(df_state)},
  {(init_f)term_init, (pre_run_f)NULL, (event_f)term_event, (revent_f)term_event_rc,
   (commit_f)no_commit, (final_f)term_final, (map_f)df_map, sizeof(df_state)},
  {(init_f)proc_init, (pre_run_f)NULL, (event_f)proc_event, (revent_f)proc_event_rc,
   (commit_f)no_commit, (final_f)proc_final, (map_f)df_map, sizeof(df_state)},
  {0},
};

static const tw_optdef app_opt[] = {
  TWOPT_GROUP("Dragonfly (Charades-compatible)"),
  TWOPT_UINT("routers-per-group", a_routers, "routers per group (a)"),
  TWOPT_UINT("global-channels", h_global, "global channels per router (h, even)"),
  TWOPT_UINT("terminals-per-router", p_terms, "terminals per router (p)"),
  TWOPT_UINT("traffic", traffic, "UNIFORM RANDOM=1, WORST CASE=2, TRANSPOSE=3, NEAREST NEIGHBOR=4, BISECTION=5, NEAREST ROUTER=6"),
  TWOPT_UINT("routing", routing, "MINIMAL=0, NON_MINIMAL=1, ADAPTIVE=2 (routes as MINIMAL)"),
  TWOPT_DOUBLE("arrive-rate", mean_interval, "mean packet inter-arrival time (ns)"),
  TWOPT_UINT("local-vc-size", local_vc_size, "local channel buffer (chunks)"),
  TWOPT_UINT("global-vc-size", global_vc_size, "global channel buffer (chunks)"),
  TWOPT_UINT("terminal-vc-size", terminal_vc_size, "terminal channel buffer (chunks)"),
  TWOPT_DOUBLE("lookahead", lookahead, "minimum delay (ps), at most 100; sets g_tw_lookahead"),
  TWOPT_DOUBLE("gen-time", gen_time, "generate until this virtual time (ns) instead of 3000 / arrive-rate messages (0: off)"),
  TWOPT_UINT("mapping", use_mapping, "1: each rank holds whole routers with their terminals and processes; 0: ROSS linear"),
  TWOPT_DOUBLE("mult", mult, "event pool = mult * LPs per rank + memory (+ --extramem)"),
  TWOPT_UINT("memory", optimistic_memory, "additional event buffers per rank"),
  TWOPT_END()
};

int main(int argc, char **argv) {
  tw_opt_add(app_opt);
  tw_init(&argc, &argv);

  if (a_routers == 0 || p_terms == 0 || h_global == 0 || h_global % 2 != 0)
    tw_error(TW_LOC, "Dragonfly: need a, p > 0 and an even h > 0");
  if (routing > ADAPTIVE) tw_error(TW_LOC, "Dragonfly: unknown routing algorithm");
  if (traffic < UNIFORM_RANDOM || traffic > NEAREST_ROUTER) tw_error(TW_LOC, "Dragonfly: unknown traffic pattern");
  if (mean_interval <= 0) tw_error(TW_LOC, "Dragonfly: --arrive-rate must be positive");
  if (lookahead > (double)MIN_DELAY) tw_error(TW_LOC, "Dragonfly: lookahead must be at most 100");
  g_tw_lookahead = lookahead;
  if (gen_time < 0 || gen_time * TIME_PER_NS >= (double)UINT32_MAX) tw_error(TW_LOC, "Dragonfly: --gen-time must be in [0, 4.29e6) ns");
  gen_end = (uint64_t)llround(gen_time * TIME_PER_NS);

  num_packets = (int)(MESSAGE_SIZE / PACKET_SIZE);
  num_chunks = (int)(PACKET_SIZE / CHUNK_SIZE);
  max_packets = (int)(INJECTION_INTERVAL / mean_interval);
  num_groups = a_routers * h_global + 1;
  total_routers = a_routers * num_groups;
  total_terminals = total_routers * p_terms;
  total_procs = total_terminals;
  radix = a_routers + h_global + p_terms;
  uint64_t total_lps = (uint64_t)total_routers + 2ull * total_terminals;
  if (total_lps >= (1ull << 20)) tw_error(TW_LOC, "Dragonfly: 2^20 LPs or more: the tie-break priority would not be exact");

  unsigned int nlp;
  if (use_mapping) {
    if (total_routers < tw_nnodes()) tw_error(TW_LOC, "Dragonfly: %u routers, fewer than %u ranks", total_routers, tw_nnodes());
    rt_q = total_routers / tw_nnodes();
    rt_rem = total_routers % tw_nnodes();
    nlp = (unsigned int)(rt_count(g_tw_mynode) * (1 + 2ull * p_terms));
    g_tw_mapping = CUSTOM;
    g_tw_custom_initial_mapping = &df_custom_mapping;
    g_tw_custom_lp_global_to_local_map = &df_global_to_local;
  } else {
    if (total_lps % tw_nnodes()) tw_error(TW_LOC, "Dragonfly: %llu LPs not divisible by %u ranks", (unsigned long long)total_lps, tw_nnodes());
    nlp = (unsigned int)(total_lps / tw_nnodes());
  }
  g_tw_events_per_pe = (unsigned int)(mult * nlp) + optimistic_memory;
  tw_define_lps(nlp, sizeof(df_msg));
  if (g_tw_synchronization_protocol == SEQUENTIAL_ROLLBACK_CHECK) {
    static crv_checkpointer router_chk;
    memset(&router_chk, 0, sizeof router_chk);
    router_chk.lptype = &df_lps[0];
    router_chk.sz_storage = sizeof(df_state);
    router_chk.save_lp = router_save;
    router_chk.clean_lp = router_clean;
    router_chk.check_lps = router_check;
    crv_add_custom_state_checkpoint(&router_chk);
  }
  for (tw_lpid i = 0; i < g_tw_nlp; i++) {
    tw_lpid gid = use_mapping ? df_global(g_tw_mynode, i) : g_tw_mynode * (tw_lpid)nlp + i;
    int t = gid < total_routers ? 0 : gid < (tw_lpid)total_routers + total_terminals ? 1 : 2;
    tw_lp_settype(i, &df_lps[t]);
  }

  if (g_tw_mynode == 0) {
    printf("========================================\n");
    printf("Dragonfly (Charades-compatible)........\n");
    printf("   Groups.................%u\n", num_groups);
    printf("   Routers................%u (%u per group)\n", total_routers, a_routers);
    printf("   Global channels........%u per router\n", h_global);
    printf("   Terminals..............%u (%u per router)\n", total_terminals, p_terms);
    printf("   Routing / traffic......%u / %u\n", routing, traffic);
    if (gen_end) printf("   Arrive rate............%lf ns, generating until %.0f ns\n", mean_interval, gen_time);
    else printf("   Arrive rate............%lf ns (%d messages per process)\n", mean_interval, max_packets);
    printf("   VC sizes...............local %u, global %u, terminal %u chunks\n", local_vc_size, global_vc_size, terminal_vc_size);
    if (use_mapping) printf("   Mapping................routers with their terminals and processes, %llu-%llu routers per rank\n", (unsigned long long)rt_q, (unsigned long long)(rt_q + (rt_rem ? 1 : 0)));
    else printf("   Mapping................linear\n");
    printf("========================================\n\n");
  }

  tw_run();

  uint64_t mine[7] = {st_messages, st_injected, st_finished, st_latency, st_hops, st_dropped, st_waits}, sum[7];
  uint64_t mx;
  MPI_Reduce(mine, sum, 7, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_ROSS);
  MPI_Reduce(&st_max_latency, &mx, 1, MPI_UINT64_T, MPI_MAX, 0, MPI_COMM_ROSS);
  if (g_tw_mynode == 0) {
    double n = sum[2] ? (double)sum[2] : 1.0;
    printf("========================================\n");
    printf("Dragonfly Statistics....................\n");
    printf("   MPI messages sent......%llu\n", (unsigned long long)sum[0]);
    printf("   Packets injected.......%llu\n", (unsigned long long)sum[1]);
    printf("   Packets finished.......%llu\n", (unsigned long long)sum[2]);
    printf("   Total latency..........%llu ps\n", (unsigned long long)sum[3]);
    printf("   Total hops.............%llu\n", (unsigned long long)sum[4]);
    printf("   Average latency........%.6f ns\n", sum[3] / n / TIME_PER_NS);
    printf("   Maximum latency........%.6f ns\n", mx / TIME_PER_NS);
    printf("   Average hops...........%.6f\n", sum[4] / n);
    printf("   Chunks dropped.........%llu\n", (unsigned long long)sum[5]);
    printf("   Chunks that waited.....%llu\n", (unsigned long long)sum[6]);
    printf("========================================\n");
  }
  tw_end();
  return 0;
}
