/*
 * phold-gauss: PHOLD with optional Gaussian destinations, written to match
 * Charades' models/phold/simple (phold.C / phold.h) draw for draw.
 *
 * Gaussian PHOLD is the benchmark of Bachan et al., PADS'24 (Devastator):
 * every event goes to an LP drawn from a normal distribution with standard
 * deviation --sigma (in LP ids, wrapping around) centred on the sender.
 *
 * Semantics (identical to Charades):
 *   --mean is the TOTAL mean delay: offset = lookahead + Exp(mean - lookahead),
 *     where tw_rand_exponential(rng, L) = -L * log(U), i.e. L is the mean
 *     (same definition in ROSS and Charades).
 *   Initial events: offset = lookahead + (mean - lookahead) * Exp(1).
 *   --sigma > 0: dest = (gid + llround(0 + N(0,1) * sigma)) mod total_LPs,
 *     wrapping; --remote is ignored. N(0,1) is Box-Muller with two uniforms,
 *     sqrt(-2 log U1) * sin(2 pi U2), exactly Charades' tw_rand_normal01.
 *     We do NOT use ROSS's tw_rand_normal_sd: it caches (U1, U2) in the
 *     stream and alternates sin/cos with a flip-flop that tw_rand_reverse_unif
 *     does not restore, so it consumes 0 or 2 uniforms per call and is not
 *     rolled back correctly by counting reverse calls.
 *   --sigma = 0: ROSS/Charades remote coin flip: U <= remote ? uniform dest
 *     in [0, total_LPs) : self.
 *   --int-ts=1 (default): offsets are truncated to integers, because Charades'
 *     Time is uint64_t and `Time offset = lookahead + exp(...)` truncates.
 *     --int-ts=0 keeps ROSS's continuous double timestamps.
 *   Events with ts >= end are dropped (both engines).
 *
 * RNG draws per event: sigma > 0: 2 (normal) + 1 (exponential) = 3;
 * sigma = 0: 1 (coin) + 1 (exponential) + 1 if remote (integer).
 * The reverse handler undoes exactly that many draws.
 *
 * LP layout: ROSS's default LINEAR mapping puts gids
 * [rank*nlp, (rank+1)*nlp) on each rank (contiguous blocks, same as
 * Charades' BlockMapper per PE). Each LP's single RNG stream is seeded with
 * stream id = gid by clcg4 with tw_rand_init(31, 41), as in Charades.
 */
#include <ross.h>
#include <math.h>

typedef struct { long dummy; } pg_state;
typedef struct { long dummy; } pg_msg;

static unsigned int nlp_per_pe = 8;
static unsigned long long num_lps = 0;   /* total; overrides --nlp if > 0 */
static double mean = 1000.0;             /* total mean delay (Charades --mean) */
static double lookahead = 5.0;           /* Charades default */
static double percent_remote = 0.1;      /* Charades default */
static double sigma = 0.0;
static unsigned int start_events = 1;
static unsigned int int_ts = 1;
static double mult = 1.4;
static unsigned int optimistic_memory = 100;

static double exp_mean;                  /* mean - lookahead */
static tw_lpid ttl_lps;

static tw_peid pg_map(tw_lpid gid) { return (tw_peid)(gid / g_tw_nlp); }

static inline double pg_offset(double raw) { return int_ts ? trunc(raw) : raw; }

static void pg_init(pg_state *s, tw_lp *lp) {
  (void)s;
  for (unsigned int i = 0; i < start_events; i++) {
    double off = pg_offset(lookahead + exp_mean * tw_rand_exponential(lp->rng, 1.0));
    tw_event_send(tw_event_new(lp->gid, off, lp));
  }
}

static void pg_event(pg_state *s, tw_bf *bf, pg_msg *m, tw_lp *lp) {
  (void)s; (void)m;
  tw_lpid dest = lp->gid;
  bf->c1 = 0;
  if (sigma > 0) {
    double u1 = tw_rand_unif(lp->rng);
    double u2 = tw_rand_unif(lp->rng);
    double z = sqrt(-2.0 * log(u1)) * sin(tw_opi * u2);
    long long n = (long long)ttl_lps;
    long long step = llround(0.0 + (z * sigma));
    dest = (tw_lpid)(((((long long)lp->gid + step) % n) + n) % n);
  } else if (tw_rand_unif(lp->rng) <= percent_remote) {
    bf->c1 = 1;
    dest = (tw_lpid)tw_rand_integer(lp->rng, 0, (long)ttl_lps - 1);
  }
  if (dest >= ttl_lps) tw_error(TW_LOC, "bad dest %llu", (unsigned long long)dest);

  double off = pg_offset(lookahead + tw_rand_exponential(lp->rng, exp_mean));
  tw_event_send(tw_event_new(dest, off, lp));
}

static void pg_event_rc(pg_state *s, tw_bf *bf, pg_msg *m, tw_lp *lp) {
  (void)s; (void)m;
  tw_rand_reverse_unif(lp->rng);       /* exponential */
  tw_rand_reverse_unif(lp->rng);       /* coin, or normal U2 */
  if (sigma > 0) tw_rand_reverse_unif(lp->rng);   /* normal U1 */
  if (bf->c1 == 1) tw_rand_reverse_unif(lp->rng); /* remote dest */
}

static void pg_commit(pg_state *s, tw_bf *bf, pg_msg *m, tw_lp *lp) {
  (void)s; (void)bf; (void)m; (void)lp;
}
static void pg_final(pg_state *s, tw_lp *lp) { (void)s; (void)lp; }

static tw_lptype pg_lps[] = {
  {(init_f)pg_init, (pre_run_f)NULL, (event_f)pg_event, (revent_f)pg_event_rc,
   (commit_f)pg_commit, (final_f)pg_final, (map_f)pg_map, sizeof(pg_state)},
  {0},
};

static const tw_optdef app_opt[] = {
  TWOPT_GROUP("Gaussian PHOLD (Charades-compatible)"),
  TWOPT_UINT("nlp", nlp_per_pe, "LPs per rank (ignored if --num-lps > 0)"),
  TWOPT_ULONGLONG("num-lps", num_lps, "total LPs (must divide evenly by ranks); 0 = use --nlp"),
  TWOPT_DOUBLE("mean", mean, "TOTAL mean delay incl. lookahead (Charades --mean)"),
  TWOPT_DOUBLE("lookahead", lookahead, "minimum delay; also sets g_tw_lookahead"),
  TWOPT_DOUBLE("remote", percent_remote, "remote probability (used only when sigma == 0)"),
  TWOPT_DOUBLE("sigma", sigma, "if > 0: Gaussian destinations, std. dev. in LP ids (ignores --remote)"),
  TWOPT_UINT("start-events", start_events, "initial events per LP"),
  TWOPT_UINT("int-ts", int_ts, "1: truncate offsets to integers like Charades' uint64 Time; 0: doubles"),
  TWOPT_DOUBLE("mult", mult, "event pool = mult * nlp * start-events + memory (+ --extramem)"),
  TWOPT_UINT("memory", optimistic_memory, "additional event buffers per rank"),
  TWOPT_END()
};

int main(int argc, char **argv) {
  tw_opt_add(app_opt);
  tw_init(&argc, &argv);

  if (num_lps > 0) {
    if (num_lps % tw_nnodes())
      tw_error(TW_LOC, "--num-lps=%llu not divisible by %u ranks", num_lps, tw_nnodes());
    nlp_per_pe = (unsigned int)(num_lps / tw_nnodes());
  }
  if (lookahead > mean) tw_error(TW_LOC, "Lookahead must be less than mean");
  exp_mean = mean - lookahead;
  ttl_lps = (tw_lpid)tw_nnodes() * nlp_per_pe;
  g_tw_lookahead = lookahead;
  g_tw_events_per_pe = (unsigned int)(mult * nlp_per_pe * start_events) + optimistic_memory;

  tw_define_lps(nlp_per_pe, sizeof(pg_msg));
  for (tw_lpid i = 0; i < g_tw_nlp; i++) tw_lp_settype(i, &pg_lps[0]);

  if (g_tw_mynode == 0) {
    printf("========================================\n");
    printf("Gaussian PHOLD (Charades-compatible)....\n");
    printf("   Total LPs..............%llu (%u per rank)\n",
           (unsigned long long)ttl_lps, nlp_per_pe);
    printf("   Start events...........%u\n", start_events);
    printf("   Lookahead..............%lf\n", lookahead);
    printf("   Mean (total)...........%lf\n", mean);
    if (sigma > 0) printf("   Sigma..................%lf (Gaussian destinations)\n", sigma);
    else           printf("   Remote.................%lf\n", percent_remote);
    printf("   Integer timestamps.....%u\n", int_ts);
    printf("========================================\n\n");
  }

  tw_run();
  tw_end();
  return 0;
}
