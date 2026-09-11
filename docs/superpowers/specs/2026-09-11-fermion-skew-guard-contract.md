# Contract: fermion mode on skewed and one-row unit cells

Status: behaviour contract given to the test author (tests in test/fermion/skew_unfold.cpp,
test/input.cpp, test/python/test_tenes_std.py, test/python/test_fermion_models.py).
Branch `fermion`, base commit 03cd7cd6. Evidence for every claim below:
`docs/superpowers/notes/2026-09-11-fermion-skew-revisit.md` (read it first; it supersedes
`work/fermion/skew-validation/FINDINGS.md`, whose "skew breaks the fermionic
signs" conclusion is RETRACTED).

## 1. Background

Fermion mode currently refuses, in three layers, (a) any unit cell with
`skew != 0` and (b) any unit cell with `LX < 2` or `LY < 2`.

(a) was based on a 2026-08-20 measurement that is now known to be wrong: it
predates the CTM folding fix 3bef24a4, compared a one-row cell against a 2x2
cell, and read an ansatz restriction as a sign error: `[2,1]` skew 1 is exactly
the `[2,2]` skew-0 calculation restricted to T0 = T3, T1 = T2, and the simple
update stays in that symmetric subspace (the bosonic XY ferromagnet does the
same). (Corrected after review: the first version said the skewed cell cannot
represent the checkerboard state; that holds only for sublattice-swapping cells
such as `[2,2]` skew 1, not for `[2,1]` skew 1.) On HEAD, a skewed cell
and its unfolded skew-0 equivalent (definition in §4) give bit-identical simple
updates, mean-field measurements identical to the last bit, CTM measurements
equal to finite-chi precision (1e-13 at chi=32 on a converged state), and
single-bond full updates equal to 4e-15. The fermion layer never reads
coordinates; lattice geometry enters only through `SquareLattice::neighbor()` /
`other()`, which already implement skew for bosons.

(b) was added (8c9d8779) because a one-wide cell can make a site its OWN
nearest neighbour, and the simple update then writes the same site tensor twice
per bond (last write wins) — for fermions the parity ledger of that site can
also end up describing the discarded tensor. That hazard is real, but (b) is
broader than it: `L_sub = [2, 1]` with `skew = 1` (what `tenes_simple` builds
for a square lattice with `W = 1`) has no self-neighbour — every bond of site 0
goes to site 1 — and it is verified correct.

## 2. The new rule (all three layers)

A fermion-mode input is refused **if and only if some site of the unit cell is
its own nearest neighbour** in some direction (left, up, right, down), with the
neighbour resolved through the cell's periodic + skew boundary. With TeNeS'
conventions (`T(x, y) = T(x + skew, y + LY)` in `SquareLattice.hpp`) this is
exactly: `LX == 1` (horizontal self-neighbour, whatever the skew), or
`LY == 1` and `skew ≡ 0 (mod LX)` (vertical self-neighbour).

Examples (L_sub, skew) → verdict:

| accepted (newly) | still refused |
|---|---|
| [2,2] 1, [3,2] 1, [3,3] 2, [2,3] 1 | [1,1] any skew |
| [2,1] 1, [2,1] -1, [3,1] 1, [3,1] 2, [4,1] 2 | [1,2] any skew, [1,3] 0 |
| [2,2] 0 and every other cell accepted today | [2,1] 0, [3,1] 0, [2,1] 2, [3,1] 3, [2,1] -2 |

Negative skew is legal input (the bosonic ctest `Honeycomb_skew` uses
`skew = -1`). Note that the C++ `SquareLattice` member `skew` keeps the sign of
the input (it is `skew % LX` in C++ semantics), so do not assume it is
normalised to `[0, LX)`.

Everything else in the fermion guards is unchanged: parity metadata,
ground-state only, no RSVD / Simple_Gauge_Fix / correlation functions /
multisite / ops-form / beyond-nearest-neighbour, the parity-odd checks, etc.

### 2.1 Solver input layer — `validate_fermion_constraints` (src/iTPS/load_toml.cpp)

Throws `tenes::input_error` per the rule. The message (it goes through the
existing `throw_fermion_guard` wrapper, "fermion mode in this version does not
support <reason>; ...") must let the user see the cause: it contains the
unit-cell shape — both `L_sub` values and the `skew` value as the lattice holds
it — and states that a site would be its own nearest neighbour. Pin this with
flexible patterns: the words "own" and "neighbor" or "neighbour" (accept both
spellings), "L_sub", "skew", and the numbers. Accepted cells must pass the whole
`validate_fermion_constraints` call (with an otherwise valid fermion input:
parity metadata of the right length, etc.).

### 2.2 tenes_std — `Model._validate_fermion_mode_input` (tool/tenes_std.py)

Raises `RuntimeError` per the rule, evaluated on the `[tensor]` table
(`l_sub`, optional `skew` defaulting to 0). Message: contains "L_sub" and its
value, "skew" and its value, and the own-neighbour cause (both spellings
acceptable); it must not contain "M1" or "M2". For an accepted skewed fermion
cell the whole `Model(...)` construction succeeds and the emitted input.toml
carries `skew` unchanged (existing behaviour).

### 2.3 tenes_simple (tool/tenes_simple.py)

The fermionic models (`spinless fermion`, `hubbard`) on the square lattice
with `W = 1` are accepted: the returned lattice has `skew == 1`, the std.toml
text says `skew = 1` and has a one-row `L_sub` (`[L, 1]` in the tensor
section), and feeding that std.toml to tenes_std succeeds (fermion = true,
parity metadata present, `skew = 1` in the result). `tenes_simple` can never
produce a self-neighbour fermion cell (the square lattice asserts `L > 1`, and
`W = 1` always gets `skew = 1`), so it has no cell-shape rejection left at
all. Unchanged: non-square lattices are still refused for fermionic models;
`W >= 2` still gives `skew = 0`; bosonic models are untouched.

## 3. Existing tests whose premise is retracted

These assert the old behaviour and must be replaced (not merely deleted — keep
an equivalent net for what is still true):

- test/python/test_fermion_models.py: `TestFermionSkewGuard` (the three
  rejection tests and their "measured / wrong numbers" framing) and
  `TestHubbardSkewGuard.test_skewed_cell_is_rejected`. The acceptance and
  bosonic regression nets in those classes stay valid.
- test/python/test_tenes_std.py: `TestFermionSkewGuard` (all rejection tests)
  and `TestFermionUnitCellDimensionGuard` (its `[2,1]`/`[1,2]` rejections are
  still correct for skew 0, but the message requirements change per §2.2).
- Comments in test files that cite work/skew-validation/FINDINGS.md as
  evidence of a sign bug: correct them (cite docs/superpowers/notes/2026-09-11-fermion-skew-revisit.md).
- The C++ input tests have no skew/L_sub guard test today; add one (§2.1)
  where the other `validate_fermion_constraints` tests live (test/input.cpp
  uses `fermion_cell_toml()` and calls `validate_fermion_constraints`
  directly — follow that pattern).

## 4. Physics regression: covariance under unfolding (new C++ test)

This is the test that keeps skewed fermion cells correct in CI. It does not go
through `validate_fermion_constraints` (the `iTPS` constructor never calls it),
so it is GREEN on HEAD already; its teeth are established by mutation (§5).

**Unfolding.** A cell `[LX, LY]` with skew `s` (s ≢ 0 mod LX) is the same
infinite lattice as the cell `[LX, LY·m]` with skew 0, `m = lcm(LX, s)/s`
(with `s` taken modulo LX into [1, LX)). Unfolded site `(x, y)` holds the tensor
of skew-cell site `((x − s·⌊y/LY⌋) mod LX, y mod LY)`. For `[2,1]` skew 1 → `[2,2]`:
the unfolded sites 0,1,2,3 hold skew sites 0,1,1,0. Every skew-cell bond
(site, leg) has m images (the unfolded sites that hold `site`, same leg);
because no site is its own neighbour, the images touch pairwise disjoint site
tensors and bond weights.

**Required checks**, each with random, parity-even, site-distinct initial
tensors and site-distinct bond weights (copied onto the images), and the
physical parity `[0, 1]` (free spinless fermions; a gate from a hopping +
interaction Hamiltonian so the update is non-trivial):

1. *Simple update — exact.* Apply the skew cell's two-site gates bond by bond;
   in the unfolded cell apply each gate to all its images consecutively, in the
   same bond order. After ≥ 3 sweeps with an imaginary-time step large enough
   that the tensors change visibly, every unfolded site's `Tn`, bond weights
   (`lambda_tensor`) and virtual parity ledger equal those of the skew site it
   images, exactly (the same arithmetic runs on the same numbers; allow at
   most a few ulp if you must, and justify it).
2. *CTM measurement — finite-chi precision.* With the CTM environment of each
   cell computed, `measure_onesite()` and `measure_twosite()` of every unfolded
   site/bond equal those of the imaged skew site/bond. The difference is the
   finite-chi CTMRG residual, so ANCHOR the tolerance on what you observe
   (print it, set tol with a documented margin), and make sure the state is
   one whose CTM actually converges (precondition, not a hope). The two-site
   observables must include (i) one that is NOT symmetric under exchanging
   its two sites (e.g. n ⊗ (1 − n)), so that a wrongly oriented pair changes
   the value, and (ii) one with non-zero (odd, odd) matrix elements (hopping),
   so that a sign error changes the value; both horizontal and vertical bonds,
   including the bonds that wrap through the skewed boundary.
3. *Mean-field measurement — exact.* The same with the mean-field
   environment (`MeanField_Env = true`): agreement to rounding.
4. *Full update, one bond — CTM precision.* From a state with converged CTM
   in both cells, one `full_update` on skew-cell bond (site, leg) and on its
   FIRST image only in the unfolded cell: the two updated site tensors agree
   (observed 4e-15 on HEAD; anchor again), and the unfolded sites of the other
   image are untouched. One horizontal and one vertical bond.

**Cells.** At least `[2,1]` skew 1 → `[2,2]` (the tenes_simple `W = 1` cell)
and one cell with `LX = 3` and skew 1 (e.g. `[3,1]` skew 1 → `[3,3]`). The
LX = 3 cell matters: with LX = 2 a sign error in how skew is applied (x − s vs
x + s) is invisible because ±1 coincide mod 2, and no existing ctest (bosonic
ones included) has a skewed cell with LX > 2. If affordable, also `[2,2]`
skew 1 → `[2,4]` (skew with both sides ≥ 2). Keep the new C++ test under about
20 s total in the Debug preset (CI runners have two physical cores); D = 2,
small chi.

Helpers you will want: `iTPSTestAccessor` (test/test_fermion_common.hpp) gives
`Tn`, `lambda_tensor`, `finfo`, `update_reduced_density_environment`; `iTPS`
has public `simple_update(up)`, `full_update(up)`, `update_CTM()`,
`measure_onesite()`, `measure_twosite()`. The existing test
"fermion twosite measurement is translation invariant across wraps"
(test/test_fermion_layer.cpp) is the closest pattern. Register the file in
test/CMakeLists.txt (a new source of `test_fermion_layer`, or its own
executable if runtime says so).

## 5. Mutations the orchestrator will run against your tests

You cannot edit src/ or tool/, so the orchestrator applies these to a copy and
expects the named tests to go RED. Design so they do:

- M1: `SquareLattice::calc_neighbors` without its `skew != 0` branch → §4 tests.
- M2: `SquareLattice::index()` applying skew with the wrong sign
  (`x += skew * y_offset`) → §4 tests on the LX = 3 cell (CTM uses index()).
- M3: `Left_move_single` in src/iTPS/core/ctm_single.cpp iterating over
  `lattice.LY` instead of `lattice.LY_noskew` (the skewed column is longer
  than the cell) → §4 check 2 (the fermion CTM is the single-layer density
  CTM built from these moves).

Note that the §4 tests compare two cells that run the SAME code, so a
mutation that is not skew-specific (e.g. a sign error on every vertical bond)
corrupts both sides identically and is NOT expected to be caught by §4 — the
existing E2E and convention tests own those. §4 exists for the geometry that
only a skewed cell exercises.
- M4: the new C++ guard rejecting on `LY < 2` alone (i.e. old rule) → §2.1 tests.
- M5: the new tenes_std guard ignoring skew (i.e. old narrow-dimension rule) → §2.2 tests.
- M6: tenes_std guard checking only the vertical direction (missing LX == 1) → §2.2 tests.

## 6. Ground rules

- Modify ONLY test files: test/python/*.py, test/input.cpp,
  test/test_fermion_layer.cpp or new files under test/fermion/, and
  test/CMakeLists.txt for registration. Never touch src/, tool/, docs/.
- Do not run `git stash`, `git checkout`, `git reset`, `git restore`, or any
  command that changes files you did not write. Do not commit.
- Do not run black or clang-format; match the surrounding style by hand.
- Build with `cmake --build --preset gcc --target <target> -j 8` from the repo
  root (Debug, g++-16, MPI off). Run C++ test binaries with the current
  directory set to `work/fermion/skew-guard/run/` (they write output_* into the
  CWD); never from the repo root or the build tree. Python:
  `venv/bin/python3 -m pytest test/python/<file> -q` from the repo root.
- Deliverable: the tests AND a report at
  `work/fermion/skew-guard/test-author-report.md` — its absence is itself a
  defect. The report lists every test added/changed/removed, and for each one
  its expected state on HEAD (RED with the exact reason — the current guard's
  error text — or GREEN as a characterization test), the observed CTM/FU
  differences you anchored tolerances on, the runtime of the new C++ test, and
  anything in this contract you believe is wrong (say so; do not work around
  it silently).

## 7. Addendum after test authoring (2026-09-11)

What the test author found and how it was resolved:

- **Skew a non-zero multiple of LX never reaches the C++ guard.** The
  `SquareLattice` constructor computes `LY_noskew = LY * (lcm(LX, skew) / skew)`
  after reducing its `skew` parameter modulo LX, which is `0 / 0` for such a
  skew. On x86-64 this is SIGFPE inside `gen_lattice`; on arm64 it silently
  gives `LY_noskew = 0`. The bug predates this change and affects bosonic cells
  too (e.g. `[2,2]` skew 2). The §2 table rows "[2,1] 2, [3,1] 3, [2,1] -2"
  therefore hold for tenes_std (tested) but cannot be tested in the solver input
  layer, whose refusal test uses skew-0 cells. Fixing the constructor is outside
  this change.
- **§4 checks 2 and 4 use a different state from check 1.** On the state
  evolved by a few simple-update sweeps the chi = 8 CTM residual is 1e-5 to 3e-4
  and `[2,1]` skew 1 does not converge. Checks 2 and 4 use random,
  parity-even, site-distinct tensors with every odd virtual component scaled by
  0.3; there the CTM converges in 6-7 sweeps and the residual is at most 3e-11
  against observables of at least 1.7e-4. Tolerances: 1e-9 (CTM measurement)
  and 2e-9 (one full-update bond, with the ALS iteration count pinned so the two
  cells cannot stop one sweep apart). Checks 1 and 3 are exact on the evolved
  state.
- **"observed 4e-15" for the full update** was specific to the fold state at
  chi = 32; on the test state the worst case is 6.6e-11.
- Mutations M1-M6 and further variants were run against a reference
  implementation by the test author and against the final implementation by
  the orchestrator; every one turns the named tests RED.

## 8. Review round (2026-09-11)

### 8.1 SquareLattice with a skew that is a non-zero multiple of LX

The review showed the §7 constructor bug is reachable, not merely untestable:
tenes_std now accepts, for example, a fermion `[2,2]` skew 2 (no site is its own
neighbour), and on arm64 the solver then runs with `LY_noskew = 0` — the CTM
moves loop over nothing — and prints plausible wrong numbers without a warning
(E = -0.629 instead of -0.684 for the same lattice written as skew 0). Before
this change tenes_std refused every fermion skew != 0, which kept such inputs
out of the tool pipeline. The constructor is therefore fixed in this change.
Behaviour after the fix:

- `SquareLattice(X, Y, s)` with `s ≡ 0 (mod X)`, `s != 0` (e.g. `[2,2]` 2,
  `[2,1]` 2, `[2,1]` -2, `[3,1]` 3, `[3,2]` 6, `[1,2]` 1, `[1,1]` 1) is the same
  lattice as `SquareLattice(X, Y, 0)`: no crash, member `skew == 0`,
  `LX_noskew == X`, `LY_noskew == Y`, `N_UNIT_noskew == X*Y`, and every
  `neighbor()`, `other()` and `index()` equal to the skew-0 lattice's.
- Every other skew is unchanged: `LY_noskew = Y * lcm(X, r) / r` with
  `r = s mod X` in [1, X), and the member keeps `s % X` (C++ semantics, sign
  kept), which §2.1's tests already rely on.
- Consequences in fermion mode: `[2,2]` skew 2 is accepted and behaves as
  skew 0; `[2,1]` 2, `[2,1]` -2, `[3,1]` 3, `[1,1]` 1, `[1,2]` 5 now reach the
  §2.1 guard and are refused by it — so the §2 table's "still refused" rows
  become testable in the solver input layer.
- Bosonic cells get the same fix (they crashed on x86-64 and computed with
  `LY_noskew = 0` on arm64).

On the current code these constructions divide by zero: SIGFPE on x86-64, and
on arm64 (the development machine) no trap but `LY_noskew == 0`. The new tests
must be RED on the current code for that reason and must not depend on a trap
(a SIGFPE would kill the whole test binary on x86-64 CI before the fix; after
the fix there is nothing to trap).

### 8.2 Test-only corrections

- **Hollow ledger comparisons.** §4 checks 1 and 4 compare virtual parity
  ledgers, but on the current state every ledger stays even-first `[0, 1]`, so
  the comparison cannot fail. A reviewer's mutation that writes the SU ledger to
  a site found by in-cell coordinates (ignoring the skew) survived. Make the
  comparison bite: choose a state/gate/truncation for which at least one bond
  ledger after the sweeps differs from even-first, and assert that as a
  premise. The orchestrator will re-run that mutation.
- **Checkerboard wording.** Test comments that justify the retraction with "the
  skewed cell cannot represent the checkerboard state" are wrong for `[2,1]`
  skew 1 (see §1 as corrected).
- **Stale guard names.** Comments in test/input.cpp that name the removed
  "tensor.L_sub-dimensions guard" should name the self-neighbour guard.
