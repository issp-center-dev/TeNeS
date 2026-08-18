# Spinless free fermions on the square lattice (experimental)

- how to run
    - `tenes input.toml`
    - takes about a minute; results are written to `output/`
- model
    - spinless free fermions on the square lattice (fermion mode, experimental)
    - H = -t sum_<ij> (c_i^dag c_j + h.c.) - mu sum_i n_i with t = 1, mu = 0 (half filling)
- input
    - fermion mode supports only the expert-mode input for now:
      `tenes_simple` / `tenes_std` cannot generate fermionic inputs yet,
      so `input.toml` (including the Trotter gates) is written by hand
    - the fermionic ingredients are `fermion = true` in `[parameter.general]`
      and `parity = [0, 1]` in `[[tensor.unitcell]]`; everything else is a
      standard expert-mode input and the exchange signs are generated
      automatically
- observables
    - number density n (exact: 0.5; this sample reproduces it to ~1e-9)
    - energy per site via the bond Hamiltonian
      (exact: -0.81056; this sample (D=2, chi=8) gives -0.7328, i.e. 7.8%;
      D=4, chi=32 gives -0.8034, i.e. 0.72%, but the two-site measurement
      is currently much slower at D=4)
- unitcell
    - 2x2
- notes
    - see the `fermion` entry of the `[parameter.general]` reference in the
      manual for the list of features not yet supported in fermion mode
    - away from half filling (mu != 0) the simple update needs considerably
      more imaginary time to converge (num_step ~ 4500 at tau = 0.01 for
      mu = -1)
    - the exact reference energy is the filled-Fermi-sea integral
      E = (2pi)^-2 int_{eps(k)<0} eps(k) d^2k with eps(k) = -2t(cos kx + cos ky) - mu
