# TeNeS benchmark harness

A/B performance comparison harness for TeNeS. It runs the
`tenes_simple -> tenes_std -> tenes` pipeline over parameterized cases,
collects per-timer results from `output/timers.json`, and renders a
Markdown comparison report.

## Quick start (A/B comparison)

```bash
# 1. Build variant A (e.g. develop) and run the suite
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label baseline --tenes-dir build/src --tool-dir build/tool

# 2. Switch to variant B (e.g. cotengra branch), rebuild, and run again
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label cotengra --tenes-dir build/src --tool-dir build/tool

# 3. Compare
python3 benchmark/bench.py compare \
    benchmark/results/baseline benchmark/results/cotengra \
    --output report.md
```

Alternatively keep two build directories and switch `--tenes-dir`.

Multiple suite files can be given in one invocation; they run sequentially
under the same label (case names must be unique across the suites):

```bash
python3 benchmark/bench.py run \
    benchmark/suites/contraction.toml benchmark/suites/e2e.toml \
    --label baseline --tenes-dir build/src --tool-dir build/tool
```

## Timers

`tenes` always writes `output/timers.json`: cumulative wall times keyed by
hierarchical names.

- `total` — whole process (as `time.dat`'s "time all")
- `phase/*` — phase times (same values as `time.dat`); `count` is
  meaningless (always 1) for these entries
- `contract/<backend>/<N>x<M>` — contraction kernels per unit-cell shape;
  backends: `itps_ctm`, `itps_mf`, `density_ctm`
- `max_rank` / `min_rank` — per-rank extrema under MPI (load imbalance)

To add a new instrumentation point, place
`ScopedTimer scoped_timer("your/name");` (declared in `src/timer.hpp`)
inside the scope to measure. The harness and report pick it up
automatically; names are grouped by their first path component.

## Suites

Suite files (TOML) live in `benchmark/suites/`:

- `contraction.toml` — unit-cell x D x chi sweep (main suite)
- `e2e.toml` — representative end-to-end workloads
- `smoke.toml` — minimal case for CI smoke test
- `ci.toml` — reduced suite for GitHub Actions

A `[[case]]` either renders a `template` (a `simple.toml` with `${var}`
placeholders, parameters from `params` plus the direct product of `sweep`)
or copies an existing `input` file (a `simple.toml`). Derived parameters:
`Lsub = [L, W]` expands to `L`/`W`; `chi_ratio = r` yields
`chi = r * D * D`. Case names must be unique after expansion.

Note: the smallest supported square-lattice unit cell is 2x1 (`tenes_simple`
rejects 1x1 square lattices), so unit-cell sweeps in `contraction.toml`,
`ci.toml`, and `smoke.toml` all start at `[2, 1]`.

## Statistics and caveats

- Each case runs `repeat` times; the report shows median [min, max] and
  the ratio of medians (B/A). If the min-max intervals overlap, the row
  is marked "no sig. diff".
- If call counts differ between A and B ("count differs"), the algorithm
  changed and time ratios should be interpreted accordingly.
- Observables (`output/*.dat` except `time.dat` / `parameters.dat`) of
  `run_0` are compared with `np.isclose` (rtol 1e-3, atol 1e-4); any
  mismatch is reported as a WARNING at the top: faster-but-wrong must
  not go unnoticed.

## HPC (MPI) usage

Pass the launcher prefix and run both labels inside the *same* job to
avoid node-placement noise:

```bash
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label baseline --launcher "srun -n 64" ...
python3 benchmark/bench.py run benchmark/suites/contraction.toml \
    --label cotengra --launcher "srun -n 64" ...
```

`OMP_NUM_THREADS` and the launcher string are recorded in `meta.json`.

## GitHub Actions

The `Benchmark` workflow (manual trigger: workflow_dispatch) runs
`suites/ci.toml` on a shared runner and uploads the results directory as
an artifact. Shared-runner numbers are noisy; treat them as reference
values only and compare locally.
