/*
 * Charades' VLSI model (charades models/vlsi, RESULTS.md section 190) in ROSS:
 * gate-level simulation of an ITC'99 .bench netlist, copies joined in a ring.
 * The netlist code is Charades' own netlist.C (VLSI_NETLIST_DIR at configure
 * time), so both engines build the same circuit, layout and gids.
 *
 * Events: a pin change (even times), the gate's single evaluation 1 ps later
 * (odd times), and clock edges (multiples of the period) for flip-flops and
 * random inputs. Equal-timestamp events at an LP commute, so no tie order is
 * needed. Times are integer ps carried in ROSS's double.
 *
 * Mapping: rank k holds gids [k q + min(k, r), (k+1) q + min(k+1, r)), q and
 * r from the LP count over the ranks.
 */
#include "ross.h"
#include "netlist.h"

#include <cinttypes>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>

using namespace vlsi;

static Netlist G;
static char netlist_path[1024] = "";
static unsigned int n_copies = 1, n_span = 1, n_chunk = 512, n_layout_bins = 0;
static unsigned long long clock_period = 0, input_seed = 1;
static double mult = 2.0;
static unsigned int optimistic_memory = 100000;

enum { M_PIN, M_EVAL, M_CLOCK };

struct vlsi_msg {
  uint64_t saved_pending;
  uint8_t kind, pin, value;
};

struct vlsi_state {
  uint32_t copy, local, pins;
  uint8_t type, nin, out;
  uint64_t pending;
};

static uint64_t st_gates, st_flops, st_inputs, st_pins, st_evals;

/* ------------------------------------------------------------------ mapping */

static uint64_t total_lps, map_q, map_r;
static uint64_t first_of(tw_peid k) { return k * map_q + (k < map_r ? k : map_r); }
static tw_peid vlsi_map(tw_lpid gid) {
  uint64_t b = map_r * (map_q + 1);
  return (tw_peid)(gid < b ? gid / (map_q + 1) : map_r + (gid - b) / map_q);
}
static void vlsi_custom_mapping(void) {
  unsigned int nlp_per_kp = (unsigned int)((g_tw_nlp + g_tw_nkp - 1) / g_tw_nkp);
  tw_lpid lpid = 0;
  uint64_t first = first_of(g_tw_mynode);
  for (tw_kpid kpid = 0; kpid < g_tw_nkp; kpid++) {
    tw_kp_onpe(kpid, g_tw_pe);
    for (unsigned int j = 0; j < nlp_per_kp && lpid < g_tw_nlp; j++, lpid++) {
      tw_lp_onpe(lpid, g_tw_pe, first + lpid);
      tw_lp_onkp(g_tw_lp[lpid], g_tw_kp[kpid]);
    }
  }
}
static tw_lp* vlsi_global_to_local(tw_lpid gid) { return g_tw_lp[gid - first_of(g_tw_mynode)]; }

/* ---------------------------------------------------------------- the gates */

static uint64_t now_of(tw_lp* lp) { return (uint64_t)tw_now(lp); }

static void send_output(vlsi_state* s, tw_lp* lp, uint64_t offset) {
  for (uint32_t k = G.fan_start[s->local]; k < G.fan_start[s->local + 1]; k++) {
    tw_event* e = tw_event_new(G.gid(s->copy, G.fan[k].dest), (tw_stime)offset, lp);
    vlsi_msg* m = (vlsi_msg*)tw_event_data(e);
    m->kind = M_PIN; m->pin = G.fan[k].pin; m->value = s->out;
    tw_event_send(e);
  }
  uint32_t next = (s->copy + 1) % G.copies;
  for (uint32_t k = G.ring_start[s->local]; k < G.ring_start[s->local + 1]; k++) {
    tw_event* e = tw_event_new(G.gid(next, G.ring[k]), (tw_stime)offset, lp);
    vlsi_msg* m = (vlsi_msg*)tw_event_data(e);
    m->kind = M_PIN; m->pin = 0; m->value = s->out;
    tw_event_send(e);
  }
}

static void gate_init(vlsi_state* s, tw_lp* lp) {
  memset(s, 0, sizeof *s);
  if (!G.locate(lp->gid, &s->copy, &s->local)) return;   // padding: inert (type G_CONST0), never sends
  s->type = G.type[s->local];
  s->nin = G.nin[s->local];
  s->out = G.init[s->local];
  s->pins = G.init_pins(s->local);
  if (s->type == G_DFF || s->type == G_INPUT) {
    tw_event* e = tw_event_new(lp->gid, (tw_stime)G.period, lp);
    ((vlsi_msg*)tw_event_data(e))->kind = M_CLOCK;
    tw_event_send(e);
  }
}

static void gate_event(vlsi_state* s, tw_bf* bf, vlsi_msg* m, tw_lp* lp) {
  bf->c0 = bf->c1 = 0;
  if (m->kind == M_PIN) {
    uint32_t bit = 1u << m->pin;
    bf->c0 = (((s->pins & bit) != 0) != (m->value != 0));
    if (!bf->c0) return;
    s->pins ^= bit;
    if (s->type >= G_NOT) {
      uint64_t t = now_of(lp) + EVAL_STEP;
      if (s->pending != t) {
        bf->c1 = 1;
        m->saved_pending = s->pending;
        s->pending = t;
        tw_event* e = tw_event_new(lp->gid, (tw_stime)EVAL_STEP, lp);
        ((vlsi_msg*)tw_event_data(e))->kind = M_EVAL;
        tw_event_send(e);
      }
    }
  } else if (m->kind == M_EVAL) {
    int v = eval(s->type, s->pins, s->nin);
    bf->c0 = (v != s->out);
    if (!bf->c0) return;
    s->out = (uint8_t)v;
    send_output(s, lp, G.delay[s->local] - EVAL_STEP);
  } else {
    int v = (s->type == G_DFF) ? (int)(s->pins & 1)
                               : G.input_value(s->copy, s->local, now_of(lp) / G.period);
    bf->c0 = (v != s->out);
    if (bf->c0) {
      s->out = (uint8_t)v;
      send_output(s, lp, s->type == G_DFF ? CLK_TO_Q : INPUT_DELAY);
    }
    tw_event* e = tw_event_new(lp->gid, (tw_stime)G.period, lp);
    ((vlsi_msg*)tw_event_data(e))->kind = M_CLOCK;
    tw_event_send(e);
  }
}

static void gate_event_rc(vlsi_state* s, tw_bf* bf, vlsi_msg* m, tw_lp* lp) {
  if (m->kind == M_PIN) {
    if (bf->c1) s->pending = m->saved_pending;
    if (bf->c0) s->pins ^= 1u << m->pin;
  } else if (bf->c0) {
    s->out ^= 1;
  }
}

static void gate_commit(vlsi_state* s, tw_bf* bf, vlsi_msg* m, tw_lp* lp) {
  if (m->kind == M_PIN) st_pins++;
  else if (m->kind == M_EVAL) { st_evals++; if (bf->c0) st_gates++; }
  else if (bf->c0) { if (s->type == G_DFF) st_flops++; else st_inputs++; }
}

static void gate_final(vlsi_state* s, tw_lp* lp) {}

static tw_lptype vlsi_lps[] = {
  {(init_f)gate_init, (pre_run_f)NULL, (event_f)gate_event, (revent_f)gate_event_rc,
   (commit_f)gate_commit, (final_f)gate_final, (map_f)vlsi_map, sizeof(vlsi_state)},
  {0},
};

static const tw_optdef app_opt[] = {
  TWOPT_GROUP("VLSI (Charades-compatible)"),
  TWOPT_CHAR("netlist", netlist_path, "ITC'99 .bench netlist"),
  TWOPT_UINT("copies", n_copies, "copies joined in a ring"),
  TWOPT_UINT("span", n_span, "spread each copy over this many groups of LPs; divides --copies"),
  TWOPT_ULONGLONG("clock-period", clock_period, "ps, even, longer than the longest path (0: longest path + 25%, rounded up to 100)"),
  TWOPT_ULONGLONG("input-seed", input_seed, "seed of the random inputs"),
  TWOPT_UINT("chunk", n_chunk, "layout: chunks of this many nets balanced by pre-simulated activity (0: file order)"),
  TWOPT_UINT("layout-bins", n_layout_bins, "bins per group of span copies (0: ranks x span / copies, one bin per rank)"),
  TWOPT_DOUBLE("mult", mult, "event pool = mult * LPs per rank + memory (+ --extramem)"),
  TWOPT_UINT("memory", optimistic_memory, "additional event buffers per rank"),
  TWOPT_END()
};

int main(int argc, char** argv) {
  tw_opt_add(app_opt);
  tw_init(&argc, &argv);
  if (!netlist_path[0]) tw_error(TW_LOC, "VLSI: give a netlist (--netlist=FILE.bench)");
  try {
    G.load(netlist_path, n_copies, n_span, clock_period, input_seed);
  } catch (const std::exception& e) {
    tw_error(TW_LOC, "%s", e.what());
  }
  {
    unsigned long long b = (unsigned long long)tw_nnodes() * n_span / n_copies;
    G.balance(n_chunk, n_layout_bins ? n_layout_bins : (unsigned int)(b ? b : 1));
  }
  total_lps = G.total();
  if (total_lps < tw_nnodes()) tw_error(TW_LOC, "VLSI: fewer LPs than ranks");
  map_q = total_lps / tw_nnodes();
  map_r = total_lps % tw_nnodes();
  unsigned int nlp = (unsigned int)(map_q + (g_tw_mynode < map_r ? 1 : 0));
  g_tw_mapping = CUSTOM;
  g_tw_custom_initial_mapping = &vlsi_custom_mapping;
  g_tw_custom_lp_global_to_local_map = &vlsi_global_to_local;
  g_tw_events_per_pe = (unsigned int)(mult * nlp) + optimistic_memory;
  g_tw_lookahead = (double)EVAL_STEP;
  tw_define_lps(nlp, sizeof(vlsi_msg));
  for (tw_lpid i = 0; i < g_tw_nlp; i++) tw_lp_settype(i, &vlsi_lps[0]);

  if (g_tw_mynode == 0) {
    printf("========================================\n");
    printf("VLSI (Charades-compatible)..............\n");
    printf("   Netlist................%s\n", G.name.c_str());
    printf("   Per copy...............%u nets: %u gates, %u flip-flops, %u random inputs, %zu bus bits\n",
           G.nets, G.gates, G.flops, G.random_inputs, G.bus.size());
    printf("   Copies.................%u (ring), span %u\n", G.copies, G.span);
    printf("   LPs....................%llu (%llu-%llu per rank)\n", (unsigned long long)total_lps,
           (unsigned long long)map_q, (unsigned long long)(map_q + (map_r ? 1 : 0)));
    if (G.chunk) printf("   Layout.................chunks of %u nets, %u per bin, %u bins per group (busiest bin %.2fx the mean, estimated)\n",
                        G.chunk, G.slots, G.bins, G.busiest);
    else printf("   Layout.................file order\n");
    printf("   Longest path...........%llu ps\n", (unsigned long long)G.longest_path);
    printf("   Clock period...........%llu ps\n", (unsigned long long)G.period);
    printf("========================================\n");
  }
  tw_run();

  uint64_t mine[5] = {st_gates, st_flops, st_inputs, st_pins, st_evals}, sum[5];
  MPI_Reduce(mine, sum, 5, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_ROSS);
  if (g_tw_mynode == 0) {
    printf("========================================\n");
    printf("VLSI Statistics.........................\n");
    printf("   Gate output changes....%" PRIu64 "\n", sum[0]);
    printf("   Flip-flop changes......%" PRIu64 "\n", sum[1]);
    printf("   Input changes..........%" PRIu64 "\n", sum[2]);
    printf("   Pin changes (messages).%" PRIu64 "\n", sum[3]);
    printf("   Evaluations............%" PRIu64 "\n", sum[4]);
    printf("========================================\n");
  }
  tw_end();
  return 0;
}
