- how to run
    - simple.toml
        - `tenes_simple simple.toml`
        - `tenes_std std.toml`
        - `tenes input.toml`
    - run.py + basic.toml
        - `python3 run.py`
- model
    - Hardcore Bose-Hubbard model on triangular lattice
- Parameter
    - simple.toml
        - t/V = 0.4
        - mu/V = 3.0
        - superfluid phase
    - run.py + basic.toml
        - t/V = 0.1
        - mu/V = [-2, 3]
        - vacuum -> superfluid -> 1/3 solid -> supersolid
- Observables
    - number density
    - offdiag order parameter |b| as an indicator for the superfluid phase
    - Structure factor S(Q) with respect to sqrt(3) order as an indicator for the solid phase
        - calculated only by run.py
- Unitcell
    - 3x3
- Reference
    - S. Wessel and M. Troyer, PRL 95, 127205 (2005)
