/*
 * traffic: Charades' models/traffic (traffic.C, intersection.C/.h) for ROSS,
 * written to match it draw for draw, so both engines commit the same events.
 *
 * A num_x x num_y grid of intersections (one LP each, gid = y * num_x + x,
 * neighbours wrapping around: a torus). Cars arrive (ARRIVAL), choose a
 * direction (DIRECTION), and leave for the next intersection (DEPARTURE ->
 * the neighbour's ARRIVAL). Every delay is lookahead + Exp(mean - lookahead),
 * truncated to an integer as Charades' uint64 Time does.
 *
 * --trip-radius=R > 0 (Charades' weak-scaling form): each destination is
 * within R intersections in x and y, on the torus, never the car's own
 * intersection (one draw over the (2R+1)^2 - 1 other cells); a car that
 * reaches it picks the next one and keeps driving; routing takes the short
 * way around. R = 0 is the original model: destinations anywhere, routing
 * by coordinate comparison, and a car stops at its destination.
 *
 * Equal timestamps: the car an event carries gets that intersection's next
 * random draws, so the order of equal-timestamp events at one LP changes the
 * run. Charades orders them by sending LP, then by the order the sender
 * created them; here each event's ROSS user priority is
 * (sender gid * 2^32 + the sender's send count) / 2^52, which gives the same
 * order (ROSS compares the priority, in [0, 1], before its random
 * tiebreaker).
 *
 * RNG: each LP's one clcg4 stream, as in Charades; the reverse handlers undo
 * exactly the draws the forward handlers made. Not ported: the start cluster
 * (--percent-start), which none of the configs use.
 */
#include <ross.h>
#include <math.h>
#include <stdint.h>

enum { NORTH = 0, EAST = 1, SOUTH = 2, WEST = 3 };
enum { ARRIVAL = 0, DEPARTURE = 1, DIRECTION = 2 };

typedef struct {
  uint32_t src, id;
  uint32_t dest_x, dest_y;
  uint8_t current_lane;
} car_t;

typedef struct {
  int type;
  car_t car;
} tr_msg;

typedef struct {
  uint32_t this_x, this_y;
  uint32_t initial_cars, total_cars_arrived, total_cars_finished;
  uint32_t in_lane[4], out_lane[4];
  tw_lpid neighbors[4];
  uint64_t sends;  /* events this LP created: its tie-break sequence */
} tr_state;

static unsigned int num_x = 4, num_y = 4;
static unsigned int total_cars = 256;
static double mean_total = 1000.0;  /* Charades --mean; the exponential mean is mean - lookahead */
static double lookahead = 5.0;      /* Charades' default */
static double percent_end = 0.0;
static unsigned int end_size = 1, end_x = 0, end_y = 0;
static unsigned int trip_radius = 0;
static double mult = 4.0;
static unsigned int optimistic_memory = 10000;

static double exp_mean;
static tw_lpid ttl_lps;

static tw_peid tr_map(tw_lpid gid) { return (tw_peid)(gid / g_tw_nlp); }

/* Charades' Time offset = lookahead + Exp(mean), truncated */
static inline double tr_offset(tw_lp *lp) {
  return trunc(lookahead + tw_rand_exponential(lp->rng, exp_mean));
}

/* ROSS takes priorities in [0, 1]: (gid * 2^32 + sends) / 2^52 is exact for
 * gid < 2^20 (checked in main) and sends < 2^32 */
static void tr_send(tr_state *s, tw_lp *lp, tw_lpid dest, double off, int type, const car_t *car) {
  double prio = ldexp((double)lp->gid * 4294967296.0 + (double)(s->sends++), -52);
  tw_event *e = tw_event_new_user_prio(dest, off, lp, prio);
  tr_msg *m = (tr_msg *)tw_event_data(e);
  m->type = type;
  m->car = *car;
  tw_event_send(e);
}

static inline int64_t ring_delta(uint32_t a, uint32_t b, uint32_t n) {
  int64_t d = ((int64_t)b - a + n) % n;
  return d <= (int64_t)n / 2 ? d : d - n;
}

static void near_destination(tw_lp *lp, uint32_t x, uint32_t y, car_t *car) {
  int64_t R = trip_radius, w = 2 * R + 1, centre = (w * w) / 2;
  int64_t k = tw_rand_integer(lp->rng, 0, w * w - 2);
  if (k >= centre) k++;
  car->dest_x = (uint32_t)(((int64_t)x + k % w - R + num_x) % num_x);
  car->dest_y = (uint32_t)(((int64_t)y + k / w - R + num_y) % num_y);
}

static void tr_init(tr_state *s, tw_lp *lp) {
  uint64_t gid = lp->gid;
  s->initial_cars = total_cars / ttl_lps + (gid < total_cars % ttl_lps ? 1 : 0);
  s->total_cars_arrived = s->total_cars_finished = 0;
  for (int i = 0; i < 4; i++) s->in_lane[i] = s->out_lane[i] = 0;
  s->sends = 0;
  s->this_x = gid % num_x;
  s->this_y = gid / num_x;
  uint32_t x = s->this_x, y = s->this_y;
  s->neighbors[NORTH] = ((y - 1 + num_y) % num_y) * num_x + x;
  s->neighbors[SOUTH] = ((y + 1) % num_y) * num_x + x;
  s->neighbors[EAST] = y * num_x + ((x + 1) % num_x);
  s->neighbors[WEST] = y * num_x + ((x - 1 + num_x) % num_x);

  for (uint32_t i = 0; i < s->initial_cars; i++) {
    double off = tr_offset(lp);
    car_t car;
    car.src = (uint32_t)gid;
    car.id = i;
    car.current_lane = (uint8_t)tw_rand_integer(lp->rng, 0, 3);
    if (trip_radius) {
      near_destination(lp, x, y, &car);
    } else if (percent_end > tw_rand_unif(lp->rng)) {
      car.dest_x = (uint32_t)tw_rand_integer(lp->rng, 0, end_size - 1) + end_x;
      car.dest_y = (uint32_t)tw_rand_integer(lp->rng, 0, end_size - 1) + end_y;
    } else {
      car.dest_x = (uint32_t)tw_rand_integer(lp->rng, 0, num_x - 1);
      car.dest_y = (uint32_t)tw_rand_integer(lp->rng, 0, num_y - 1);
    }
    tr_send(s, lp, gid, off, ARRIVAL, &car);
  }
}

static void arrival(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  int here = s->this_x == m->car.dest_x && s->this_y == m->car.dest_y;
  if (here && !trip_radius) {
    s->total_cars_finished++;
    return;
  }
  car_t car = m->car;
  bf->c2 = here;
  if (here) {
    s->total_cars_finished++;
    near_destination(lp, s->this_x, s->this_y, &car);
  }
  uint8_t new_lane = (car.current_lane + 2) % 4;
  double off = tr_offset(lp);
  car.current_lane = new_lane;
  tr_send(s, lp, lp->gid, off, DIRECTION, &car);
  s->total_cars_arrived++;
  s->in_lane[new_lane]++;
}

static void arrival_rc(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  if (s->this_x == m->car.dest_x && s->this_y == m->car.dest_y && !trip_radius) {
    s->total_cars_finished--;
    return;
  }
  uint8_t new_lane = (m->car.current_lane + 2) % 4;
  tw_rand_reverse_unif(lp->rng);
  if (bf->c2) {
    tw_rand_reverse_unif(lp->rng);
    s->total_cars_finished--;
  }
  s->sends--;
  s->total_cars_arrived--;
  s->in_lane[new_lane]--;
}

static void departure(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  (void)bf;
  tw_lpid dest = s->neighbors[m->car.current_lane];
  double off = tr_offset(lp);
  tr_send(s, lp, dest, off, ARRIVAL, &m->car);
  s->out_lane[m->car.current_lane]--;
}

static void departure_rc(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  (void)bf;
  tw_rand_reverse_unif(lp->rng);
  s->sends--;
  s->out_lane[m->car.current_lane]++;
}

static void direction(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  uint32_t x = s->this_x, y = s->this_y;
  uint8_t new_lane, decision;
  if (x != m->car.dest_x && y != m->car.dest_y) {
    decision = (uint8_t)tw_rand_integer(lp->rng, 0, 1);
  } else if (x != m->car.dest_x) {
    decision = 0;
  } else if (y != m->car.dest_y) {
    decision = 1;
  } else {
    tw_error(TW_LOC, "Choosing a direction from our destination");
  }
  int east = trip_radius ? ring_delta(x, m->car.dest_x, num_x) > 0 : x < m->car.dest_x;
  int south = trip_radius ? ring_delta(y, m->car.dest_y, num_y) > 0 : y < m->car.dest_y;
  if (decision == 0) {
    bf->c0 = 0;
    if (east) { bf->c1 = 0; new_lane = EAST; } else { bf->c1 = 1; new_lane = WEST; }
  } else {
    bf->c0 = 1;
    if (south) { bf->c1 = 0; new_lane = SOUTH; } else { bf->c1 = 1; new_lane = NORTH; }
  }
  double off = tr_offset(lp);
  car_t car = m->car;
  car.current_lane = new_lane;
  tr_send(s, lp, lp->gid, off, DEPARTURE, &car);
  s->in_lane[m->car.current_lane]--;
  s->out_lane[new_lane]++;
}

static void direction_rc(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  uint8_t new_lane;
  if (s->this_x != m->car.dest_x && s->this_y != m->car.dest_y) tw_rand_reverse_unif(lp->rng);
  if (bf->c0 == 0) new_lane = bf->c1 == 0 ? EAST : WEST;
  else new_lane = bf->c1 == 0 ? SOUTH : NORTH;
  tw_rand_reverse_unif(lp->rng);
  s->sends--;
  s->in_lane[m->car.current_lane]++;
  s->out_lane[new_lane]--;
}

static void tr_event(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  switch (m->type) {
    case ARRIVAL: arrival(s, bf, m, lp); break;
    case DEPARTURE: departure(s, bf, m, lp); break;
    case DIRECTION: direction(s, bf, m, lp); break;
    default: tw_error(TW_LOC, "unknown event type %d", m->type);
  }
}

static void tr_event_rc(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  switch (m->type) {
    case ARRIVAL: arrival_rc(s, bf, m, lp); break;
    case DEPARTURE: departure_rc(s, bf, m, lp); break;
    case DIRECTION: direction_rc(s, bf, m, lp); break;
  }
}

static void tr_commit(tr_state *s, tw_bf *bf, tr_msg *m, tw_lp *lp) {
  (void)s; (void)bf; (void)m; (void)lp;
}

static unsigned long long finished_trips = 0;
static void tr_final(tr_state *s, tw_lp *lp) {
  (void)lp;
  finished_trips += s->total_cars_finished;
}

static tw_lptype tr_lps[] = {
  {(init_f)tr_init, (pre_run_f)NULL, (event_f)tr_event, (revent_f)tr_event_rc,
   (commit_f)tr_commit, (final_f)tr_final, (map_f)tr_map, sizeof(tr_state)},
  {0},
};

static const tw_optdef app_opt[] = {
  TWOPT_GROUP("Traffic (Charades-compatible)"),
  TWOPT_UINT("num-x", num_x, "width of the grid of intersections"),
  TWOPT_UINT("num-y", num_y, "height of the grid of intersections"),
  TWOPT_UINT("total-cars", total_cars, "cars in the simulation"),
  TWOPT_DOUBLE("mean", mean_total, "TOTAL mean delay incl. lookahead (Charades --mean)"),
  TWOPT_DOUBLE("lookahead", lookahead, "minimum delay; also sets g_tw_lookahead"),
  TWOPT_DOUBLE("percent-end", percent_end, "share of cars bound for the end cluster [0-1]"),
  /* Charades --end-size/-x/-y; renamed because ROSS matches option names by
     prefix, so "end-..." would capture --end */
  TWOPT_UINT("dest-size", end_size, "end cluster (Charades --end-size): X by X block"),
  TWOPT_UINT("dest-x", end_x, "end cluster (Charades --end-x): x of its upper left corner"),
  TWOPT_UINT("dest-y", end_y, "end cluster (Charades --end-y): y of its upper left corner"),
  TWOPT_UINT("trip-radius", trip_radius, "cars keep driving between destinations within this distance (0: the original model)"),
  TWOPT_DOUBLE("mult", mult, "event pool = mult * cars per rank + memory (+ --extramem)"),
  TWOPT_UINT("memory", optimistic_memory, "additional event buffers per rank"),
  TWOPT_END()
};

int main(int argc, char **argv) {
  tw_opt_add(app_opt);
  tw_init(&argc, &argv);

  ttl_lps = (tw_lpid)num_x * num_y;
  if (ttl_lps % tw_nnodes())
    tw_error(TW_LOC, "%llu intersections not divisible by %u ranks", (unsigned long long)ttl_lps, tw_nnodes());
  if (lookahead > mean_total) tw_error(TW_LOC, "Lookahead must be less than mean");
  if (ttl_lps > (1u << 20)) tw_error(TW_LOC, "more than 2^20 intersections: the tie-break priority would not be exact");
  if (trip_radius && (2 * trip_radius >= num_x || 2 * trip_radius >= num_y))
    tw_error(TW_LOC, "--trip-radius must be less than half the grid's width and height");
  if (trip_radius && percent_end > 0)
    tw_error(TW_LOC, "--trip-radius does not combine with --percent-end");
  /* Charades: Time mean (integer) = mean - lookahead */
  exp_mean = trunc(mean_total) - trunc(lookahead);
  g_tw_lookahead = lookahead;
  unsigned int nlp = (unsigned int)(ttl_lps / tw_nnodes());
  g_tw_events_per_pe = (unsigned int)(mult * (total_cars / tw_nnodes() + 1)) + optimistic_memory;

  tw_define_lps(nlp, sizeof(tr_msg));
  for (tw_lpid i = 0; i < g_tw_nlp; i++) tw_lp_settype(i, &tr_lps[0]);

  if (g_tw_mynode == 0) {
    printf("========================================\n");
    printf("Traffic (Charades-compatible)...........\n");
    printf("   Intersections..........%u x %u (%u per rank)\n", num_x, num_y, nlp);
    printf("   Cars...................%u\n", total_cars);
    printf("   Mean (total)...........%lf, lookahead %lf\n", mean_total, lookahead);
    if (trip_radius) printf("   Trips..................cars keep driving, destinations within %u\n", trip_radius);
    printf("========================================\n\n");
  }

  tw_run();
  if (g_tw_mynode == 0) printf("Trips finished on rank 0's LPs: %llu\n", finished_trips);
  tw_end();
  return 0;
}
