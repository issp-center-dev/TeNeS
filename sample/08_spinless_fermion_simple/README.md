# Spinless free fermions on the square lattice (simple mode, experimental)

- how to run
    - `tenes_simple simple.toml`
    - `tenes_std std.toml`
    - `tenes input.toml`
- model
    - spinless free fermions on the square lattice (fermion mode, experimental)
    - H = -t sum_<ij> (c_i^dag c_j + h.c.) - mu sum_i n_i with t = 1, mu = 0 (half filling)
- expected results
    - D = 2 and chi = 8 give E ≈ -0.7328 and n = 0.5
    - the exact reference energy is -0.81056
    - this is numerically the same calculation as the expert-mode sample in
      `sample/07_spinless_fermion`
- notes
    - the current fermion support in `tenes_simple` / `tenes_std` is limited to
      nearest-neighbor bonds on the square lattice
    - skewed unit cells, including the one automatically selected for square
      lattices with W = 1, are not supported
    - correlation functions, correlation lengths, full update, and finite
      temperature calculations are not supported in fermion mode
