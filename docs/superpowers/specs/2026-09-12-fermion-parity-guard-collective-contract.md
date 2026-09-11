# Contract: collective parity guards (final review 2, major-1)

Branch `fermion`, after 6c5c883b. Prose only; the tests are yours.

## Problem

`tenes::fermion::enforce_even_parity` (the simple-update guard, now a public
function in src/fermion/fops.hpp; it used to sit in an anonymous namespace of
src/iTPS/core/simple_update.cpp and was moved without a behaviour change) decides
whether to throw from `parity_violation(a)`, which is the maximum over the
PROCESS-LOCAL slice, against `max_abs(a)`, which is reduced over all ranks. With
two or more MPI ranks, only the rank that owns the offending element throws; the
others carry on into the next collective call and the run hangs. The same pattern
sits in two debug-build checks: the layer checks at the top of
`doubled_pipeline_traced` (src/fermion/reduced.hpp) and `validate_block_diagonal`
(src/fermion/fops.hpp). The full-update and save/load guards were made collective
in 63b07273; these three were missed.

## Behaviour after the fix

1. A new collective predicate in src/fermion/fops.hpp, namespace `tenes::fermion`:

       template <class tensor>
       void require_even_parity(const ftensor<tensor>& a, const char* context);

   Every rank of the tensor's communicator computes the same two numbers: the
   largest magnitude in the parity-odd sector over ALL ranks, and `max_abs(a)`
   (all ranks). It throws `std::runtime_error` on every rank iff the first
   exceeds `1e-10 * max(1, second)`, and on no rank otherwise. The message starts
   with `context`, says the tensor is not parity even, and contains the magnitude
   and the threshold. It must be called by every rank (it is collective).
2. `enforce_even_parity(ftensor& a)` = `require_even_parity(a, <SU context>)`,
   then, when it did not throw, the parity-odd elements are set to zero on every
   rank (unchanged behaviour: round-off below the threshold is clipped). Its
   message keeps saying "fermion Simple_update_bond produced odd-parity elements"
   or equivalent.
3. The layer checks of `doubled_pipeline_traced` use `require_even_parity` with
   the contexts "doubled_pipeline_traced: bra layer" / "... ket layer" (debug
   builds only, as today).
4. `validate_block_diagonal(sorted, row_even, col_even, context)` (debug builds
   only, as today; `sorted` is a plain tensor) decides collectively too: the
   off-diagonal maximum and the scale are both reduced over the ranks before the
   comparison.
5. `parity_violation()` itself stays process-local (its Doxygen says so; other
   code relies on it).
6. With one rank, every decision is what it is today.

## What the tests should establish

Put them in test/fermion/mpi_reduction.cpp: that file is compiled into
`test_fermion_layer` (run at one rank everywhere) and, in MPI builds, is also
registered at two ranks as `test_fermion_layer_mpi2`
(`--source-file=*mpi_reduction.cpp`). Read its header comment and the existing
case first; note in particular that `get_value` on a non-owning rank returns
false without writing, so references must be built arithmetically, never read
back that way.

- `require_even_parity`: with an above-threshold element in the odd sector that
  is owned by exactly one rank (pick it at run time from `local_size()` /
  `global_index()` so the test works for any rank count), it throws on EVERY rank;
  with the same element below the threshold it throws on NO rank; a clean tensor
  never throws. At one rank these are trivially consistent; at two ranks the
  "owned by one rank" premise must actually hold — assert it (e.g. count owners
  with a reduction) so the test cannot become hollow if the layout changes.
- `enforce_even_parity`: the same throw-everywhere / throw-nowhere contract, and
  after a below-threshold call the odd-sector elements are zero on their owners
  while even-sector elements are untouched.
- `validate_block_diagonal`: the same collective decision, under `#ifndef NDEBUG`
  (it is a no-op in release builds; say so in the test).
- The layer checks of `doubled_pipeline_traced` are covered through
  `require_even_parity`; do NOT call `doubled_pipeline_traced` itself with a
  one-rank violation — before the fix the non-throwing rank would enter the
  collective contractions that follow and the test would hang instead of fail.

Design so that these mutations turn the two-rank run RED (not hang): each of the
three guards reverted to a rank-local violation; `require_even_parity` reducing
the violation but comparing it with a local scale; `enforce_even_parity` not
calling `require_even_parity`.

## Ground rules

- Modify ONLY test/fermion/mpi_reduction.cpp (and test/CMakeLists.txt only if a
  registration change is truly needed; the existing `--source-file` filter already
  picks up new cases in that file). Never touch src/, tool/, docs/.
- No git stash / checkout / reset / restore / clean / add / commit; no formatter
  on tree files; match the surrounding style by hand.
- Build: `cmake --build --preset gcc --target test_fermion_layer -j 8` (serial,
  Debug) and `cmake --build --preset gcc-mpi --target test_fermion_layer -j 8`
  (OpenMPI, Debug). Run from work/fermion/parity-guard/run/ (create it):
  `../../../../out-gcc/build/test/test_fermion_layer --source-file=*mpi_reduction.cpp`
  and `mpiexec -n 2 ../../../../out-gcc-mpi/build/test/test_fermion_layer --source-file=*mpi_reduction.cpp`.
- The functions of §1 and §4's collective form do not exist yet. Tests that call
  `require_even_parity` will not compile until the fix lands — that is the
  expected RED for them; say so. For `enforce_even_parity` (exists, rank-local
  today) and `validate_block_diagonal` (exists, rank-local today) the two-rank run
  must be RED now for the right reason (a rank that does not throw), and must not
  hang. If you want to see the compiled tests go GREEN, use a stand-in fix in a
  scratch copy (work/fermion/parity-guard/scratch/), never in src/.
- Deliverable: the tests and work/fermion/parity-guard/test-author-report.md (its
  absence is a defect): each test, its state now with the reason (compile error,
  or the exact failing check at two ranks), the SHA-256 of every file you changed,
  and anything in this contract you think is wrong.

## Addendum after test authoring

- The message format is pinned by an existing test: `test/fermion/fold_geometry.cpp`
  looks for "bra layer is not parity even" / "ket layer is not parity even", so the
  message reads `"<context> is not parity even; odd-parity elements up to max_abs=<v>
  (threshold <t>)"`. The simple-update context is "fermion Simple_update_bond: the
  updated tensor".
- A reverted-to-local scale is only visible when an element larger than 1 sits on a
  rank other than the deciding one (the threshold is floored at max(1, .)); the tests
  place a 4.0 there.
- The `require_even_parity` cases cannot compile before the fix, so the two-rank RED of
  the other cases was shown with a copy that drops them
  (work/fermion/parity-guard/scratch/mpi_reduction.stage1.cpp, not committed).
- Quote the doctest filter in zsh: `'--source-file=*mpi_reduction.cpp'`.
