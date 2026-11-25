# EXP-067: Computational cost of matched fitness oracles and bounded-neighbourhood GNK baselines

## Status and scientific question

**Planned protocol; implementation pending. No RUN has been created or executed.**

2026-11-19 model-version note: the scientist subsequently selected a full-pairwise
default for candidate 0.6.0rc1. This protocol still explicitly targets the frozen
0.5.0rc1 ranked construction. Its A/B representation comparison is not a benchmark
of the new default. The guarded harness remains disabled; testing the new default
under HYP-010 requires an explicit protocol amendment rather than substituting
the current checkout into these arms. No design or frozen inputs are changed here.

This tests HYP-010 under AIM-001, motivated by RES-063 and the maintained
profile-conditioned model following EXP-065. Does this implementation offer a
reproducible advantage for online fitness queries after accounting for setup,
memory, workload and retained capabilities? This does not test biological
fidelity, learnability or searchability (HYP-011/012).

The prediction is a memory advantage over exact local tables, potentially at a
cost in warm query speed. A stateless GNK oracle may remove the memory advantage.
A null result or a time/memory trade-off is informative. No universal speedup or
new asymptotic class is predicted.

## What n means

**n = q^L is the number of possible complete sequences**, q=20, and L=log_q(n).
Q is the number actually queried; B is batch size. The smoothed profile changes
preferences, not the set of allowed sequences.

Standard NK with K partners per focal site has L local functions, each containing
q^(K+1) entries. Its local tables therefore contain q^(K+1) log_q(n) values,
not n complete-sequence scores. At fixed K and q this is O(log n) memory.
An online NK oracle does not require whole-landscape enumeration. The local-table
exponential is in K+1; if K=L-1, that representation grows to O(n log n).

[Complexity notes](complexity_notes.md) derive setup, normalization, storage,
complete queries and mutation queries for all arms, including retained Gram arrays,
support sorting/retries and public validation allocations. Fixed-complexity
O(log n) storage/query bounds are shared; constants and setup costs can differ.

## Four arms, two exact comparisons

| Arm | Representation | What the comparison establishes |
|---|---|---|
| A: factorized | Frozen Adapt-Env 0.5.0rc1, profile additive term, rank-eight centered contrasts, exact realized mutation normalization | Historical candidate oracle |
| B: exact local tables | Expand A's support functions into q^k tables; preserve profile, reference, coefficients and readout | A/B isolates representation cost for the same map |
| C: conventional GNK tables | Independent-entry GNK prescription implemented by a frozen keyed counter-to-Gaussian function, eagerly tabulated | Conventional bounded-support local-table baseline |
| D: counter GNK | Identical local value function to C, evaluated on demand without an entry cache | C/D isolates storage versus recomputation for the same map |

There is no Ehrlich arm. Landscapy's convenience NK constructor enumerates full
sequences; that is an API choice, not an NK requirement. It is not the primary
online baseline. C/D implement the bounded-support GNK prescription directly.
No arm enumerates the full landscape during timed work.

A/B and C/D are different ensembles. Arbitrary local tables contain lower-order
components; A's centered rank-limited products restrict the available tensors.
Matching topology, maximum order and expected mutation scale does not establish
matching expressivity or descriptors. A versus C/D is a qualified practical
comparison, not an exact acceleration of unrestricted GNK.

All arms share the profile, reference and support topology per seed/length,
float64 precision and canonical amino-acid ordering. Expand B with tensor outer
products and bounded per-support scratch, then release unnecessary factors/Grams.
Charge A's common construction plus conversion to B's cold total, reporting each.
C/D use the same topology algorithm without building unused factor fields. Give
B/C/D mixed-radix indices and incident-support updates; do not handicap them with
scalar cell loops or explicit arrays of every mutant sequence. A uses the frozen
public API. A future optimized A is a separately identified implementation.

## Normalization and counter contract

A/B retain exact realized interaction mutation RMS sigma_A. For C/D let h_S have
ideal independent mean-zero, variance-one local entries and normalized order
weights w_k. Set

    F(x) = sum_k c_k sum_{S: |S|=k} h_S(x_S)
    c_k = sqrt(L w_k / (2 k M_k)).

A uniform nonself mutation changes a support with probability k/L and has expected
squared entry difference two. Thus expected mutation energy in F is
sum_k (2 k M_k/L)c_k^2 = 1. Use G=sigma_A F and A/B's reference subtraction and
readout. This is ensemble-expected, not exact realized normalization. Overlapping
arbitrary tables may have nonzero realized cross-energies. Sampled RMS diagnostics
must not refit C/D. Exact normalization, if later added, needs a protocol amendment
and charged computation; it is not free for D.

Use documented Random123 Philox or Threefry, not an ad-hoc scalar hash. Before a
production RUN freeze rounds, word width, key, support/local-index counter packing,
domain separation, overflow semantics, open-interval uniform conversion, Gaussian
transform, backend and precision. Packing must be injective on the declared grid.
Check published known-answer vectors independently of C/D agreement. A finite
counter generator supplies a pseudorandom realization, not every possible iid
continuous-Gaussian table. C/D equality refers to the frozen prescription.
Unbounded entry memoization is prohibited; O(M) current-sequence index/value caches
for mutation updates are allowed and charged to memory.

## Design and inputs

The [config](config.yaml) is the prospective machine-readable design. Source:
Adapt-Env commit 745c00b3a369da33e7d7c7b9850ccbc212fac714. Preserve the existing
45-file release lock. Performance grid: L={32,128,512,1024}, orders {2,3,4}, target
degrees {6,11,7}, weights {0.35,0.45,0.20}, rank eight, r=0.5, steepness=0.5 and
reference occupancy=0.5. Readout values are a computational fixture, not a newly
calibrated long-protein regime.

Tile the four archived TrpB frequency rows in model-site order into a new L-row
array; use the unchanged profile converter with uniform mixture 0.01. Do not pass
repeated source-row indices, which the converter rejects. Preserve gap-labelled
source alphabet and non-gap conditioning; archive derived matrix/profile hashes.
This is a scaling fixture, not a full-length biological alignment. No MSA search
or re-estimation is required. Share the existing seed-selected computational
reference across arms.

Ten field seeds (0..9), paired across arms, have three fresh-process repetitions
per length/arm: 480 timed workers. Use cyclic rotations of A/B/C/D indexed by
(seed+repetition) modulo four; shuffle length blocks once using PCG64 seed 67001.
Workers run sequentially on one declared host, one compute thread, no GPU and
no competing benchmark jobs.

The L=4 locked `trpb_run178_r05_profile_shared_stop_intrinsic` preset is a
correctness anchor only. Its capacity-capped topology is excluded from scaling
fits and primary timing.

## Workloads and measurement boundaries

Primary cold workload: construct one oracle and answer Q=9,600 queries in 100
separate B=96 calls. Each batch interleaves 48 independent uniform sequences and
48 nonself single mutants of the shared reference. Choose mutant sites uniformly
and replacements uniformly among the other 19 alleles. Freeze one PCG64 query
trace per length/seed, independent of fitness; replay identical inputs in all arms.
This isolates oracle cost, without training an ML model or running an optimizer.

Measure setup and the first complete trace before warm-up. Cold oracle time
includes profile conversion, topology, generation/normalization, reference scoring,
table expansion or counter setup, and query scoring. Imports, environment startup
and trace generation are excluded but measured separately; report whole-worker
elapsed time too. Stream the trace in B-row chunks and retain only output checksums.
Record input I/O separately and report scoring-only and scoring-plus-I/O totals.
Internal batch_size alone does not bound validation memory for a Q-row public call.

Checkpoint the completed primary timing and memory records immediately. After
the cold trace, warm up with one identical trace. Measure scalar latency on the
first 64 trace sequences using B=1; measure full-trace throughput at B={96,1024},
allowing a final partial batch. Measure five warm repetitions; never choose the
fastest. The scalar subset avoids billions of unnecessary Python component visits.
Also evaluate the complete 19L nonself
neighborhood at each of 20 frozen uniform parents, using the public mutation API
or matched incident-update interface. Charge baseline scoring once per parent in
every arm. A cached delta alone is not a complete public call.

Report T(Q)=T_setup+T_queries(Q), using measured trace prefixes and warm throughput.
Label extrapolated amortization curves as projections. Separately report theoretical
array payload, unique owned ndarray bytes, steady RSS, setup peak RSS, primary-trace
peak RSS and whole-worker peak RSS. The primary memory endpoint is the maximum
RSS through setup and the first trace, measured before secondary work. All copying
diagnostics and correctness checks run in separate processes; being outside a
timing region would not reset a process high-water memory measurement. Capture the actual support-builder
path/retries/exchanges where instrumentation is available without modifying A;
otherwise mark them unavailable, not inferred from time.

## Resource envelope and failures

Per worker: 4 GiB absolute RSS, one compute thread, 1,800 seconds wall time.
Orchestration request: 8 GB, 24 hours. Set BLAS/OpenMP threads consistently. Before
production freeze hardware, OS, Python/NumPy, compiler/math libraries, power mode,
RSS sampler and external-interference rule. Do not pool macOS and cluster timings.
Account for child memory if a backend launches subprocesses.

Preflight allocations. Enforce RSS with a supported mechanism or a frequent parent
RSS monitor with recorded interval/overshoot on macOS; label the latter an
approximate cap. Virtual-address limits are not RSS limits. Stop and record cap/time
breaches, retaining the last measurement. Record swap/memory pressure; invalidate
blocks affected by swapping or competing jobs under the frozen rule. Retain their
logs and rerun all four arms of an invalidated seed/repetition block. Missing or
censored observations never become invented timings or ratios.

At L=1024 B/C local-table payload is 2,543,910,400 bytes (about 2.37 GiB), before
bookkeeping/scratch. The envelope does not automatically force tables to fail.
B has factors alive during expansion, so report that peak. Resource failure is an
outcome; do not reduce precision or shrink failed cells after seeing results.
A secondary-phase timeout/OOM does not erase completed, checkpointed primary metrics.
Mark each phase independently and keep whole-worker peak RSS descriptive, since
it may include larger secondary allocations.

## Correctness gates

Before timing, in separate correctness-only processes that exit before timed
workers launch (including the RMS diagnostic below):

1. Exhaust all 20^4 sequences and all local tables at the L=4 anchor for A/B and
   C/D. Compare latent additive/interaction/logit values, final fitness and local
   entries. Use independent counter known-answer tests.
2. For each large cell compare 1,024 independent queries and 128 mutations with
   direct full-score differences and reversed mutations. For first/median/last
   support in each order, test every pair state and boundary plus 256 seeded states
   at higher order. Validate full expansion exhaustively in the small anchor.
3. Test duplicate, reordered and differently chunked inputs. C/D local entries
   must be bit-identical on the pinned backend. Accumulated values/deltas require
   abs(error) <= 1e-10 + 1e-10*abs(expected); record maximum absolute/relative errors.
   Do not relax tolerance after timing. Saturated sigmoid outputs cannot substitute
   for latent checks.
4. Verify source/input hashes and all expected seeds/cells. Failed correctness
   gates halt performance interpretation and preserve diagnostics.

Sample 4,096 independent mutation differences per seed/length to report realized
RMS for all arms without refitting. This is a capability diagnostic, not evidence
of identical ensembles. Exhaustive full-state evaluation occurs only in the small
correctness fixture and is never charged as required online construction.

## Analysis, decisions and figures

Ten seeds are independent paired units. Take medians over three process repetitions
within seed; first take warm-repeat medians within process. Do not treat 30 workers
as 30 independent landscapes. Retain all raw observations/failures. An incomplete
paired set for a primary metric has no confirmatory ratio: show available data descriptively
and report completion fractions, without opportunistic seed replacement.

The confirmatory family has 12 A-versus-B/C/D contrasts: two lengths (512,1024),
two resources (cold setup-plus-9,600-query scoring time and setup-plus-primary-trace
peak RSS),
three comparators. Compute seed-level log(comparator/A) ratios and exponentiated
mean Student-t intervals, with Bonferroni two-sided family alpha 0.05 (individual
coverage 99.5833%). State the log-ratio normality assumption and show all ten
paired values. Bootstrap 95% intervals are descriptive sensitivity checks, not
replacement decision rules. Other lengths/workloads and scaling fits are descriptive.
Fit against L=log_20(n); display log10(n)=L*log10(20), avoiding numeric overflow.

A simultaneous lower bound above zero on log(comparator/A), equivalently above
one after exponentiation, supports a reproducible resource advantage
for that comparison/workload only. Report absolute saved seconds/bytes and the
other resource too; detectability alone is not practical importance. An acceptable
use case must complete its primary phase within the declared envelope. If resources disagree report
the Pareto trade-off, without a post-hoc scalar score or universal two-fold gate.
A/B establishes exact-map representation trade-offs. Broader practical claims must
include C and D and disclose different ensembles/normalization. If D is equally
practical, record that on-demand scalability is not distinctive to factorization.
No computational outcome establishes biological fidelity or HYP-011/012.

Planned figures: four panels for setup time, peak memory, warm throughput and total
cost versus Q, showing seeds, intervals, failures and theoretical curves distinctly;
a separate correctness/capability table. Preserve raw data, manifest, environment,
paired contrasts and scaling summary, with a reproducible figure notebook. These
are expected artifacts, not existing outputs. No ANA/RES is created before evidence.

## Execution readiness and amendments

run.sh intentionally refuses execution. Before a production RUN, implement/audit
adapters and measurement harness, freeze the counter contract and isolated dependency
lock, pass correctness-only checks, and record host/enforcement. Performance-driven
changes require a dated amendment retaining the original protocol. This planned
record does not claim a runnable benchmark exists. The inherited STACK-000 ID is
bookkeeping: its historical overlay is disabled and is not a validated benchmark
environment. No run has been submitted.

## Primary literature

- Brookes, Aghazadeh & Listgarten (2022), *On the sparsity of fitness functions and
  implications for learning*, [PNAS](https://doi.org/10.1073/pnas.2109649118): GNK
  local-neighborhood tables and interaction representations.
- Otwinowski & Plotkin (2014), *Inferring fitness landscapes by regression produces
  biased estimates of epistasis*, [PNAS](https://pmc.ncbi.nlm.nih.gov/articles/PMC4050575/):
  NK construction with 2^(K+1) entries per local binary function.
- Salmon, Moraes, Dror & Shaw (2011), *Parallel Random Numbers: As Easy as 1, 2, 3*,
  [SC11](https://thesalmons.org/john/random123/papers/random123sc11.pdf): established
  counter-based pseudorandom generation for the stateless comparator.

Big-O counts and normalization equations are our derivations from model definitions
and frozen source, not asymptotic theorems quoted from these papers. See the notes
for assumptions and source anchors.
