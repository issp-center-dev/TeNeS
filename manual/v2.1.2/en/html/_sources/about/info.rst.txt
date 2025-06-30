Overview
=================
TeNeS (**Te** nsor **Ne** twork **S** olver) is an open-source program package for calculation of two-dimensional many-body quantum states based on the tensor network method.
This package calculates ground-state wavefunctions for user-defined Hamiltonian, and evaluates user-defined physical quantities such as magnetization and correlation functions.
TeNeS can calculate finite temperature quantity and real-time evolution as well as the ground-state quantity.
For predefined models and lattices, there is a tool that makes it easy for users to generate input files.
TeNeS uses an OpenMP/MPI hybrid parallelized tensor operation library and thus can deal with large-scale calculation by using massively parallel machines.

Developers
==================
TeNeS is developed by the following members.

- Tsuyoshi Okubo (Graduate School of Science, Univ. of Tokyo)
- Satoshi Morita (Faculty of Science and Technology, Keio University)
- Yuichi Motoyama (Institute for Solid State Physics, Univ. of Tokyo)
- Kazuyoshi Yoshimi (Institute for Solid State Physics, Univ. of Tokyo)
- Tatsumi Aoyama (Institute for Solid State Physics, Univ. of Tokyo)
- Takeo Kato (Institute for Solid State Physics, Univ. of Tokyo)
- Naoki Kawashima (Institute for Solid State Physics, Univ. of Tokyo)

Version information
======================

- ver. 2.1.2: released on 2025-06-30.
- ver. 2.1.0: released on 2024-02-28.
- ver. 2.0.0: released on 2023-11-17.
- ver. 2.0-beta: released on 2023-10-25.
- ver. 1.3.4: released on 2023-09-13.
- ver. 1.3.3: released on 2023-07-14.
- ver. 1.3.2: released on 2023-06-08.
- ver. 1.3.1: released on 2022-10-21.
- ver. 1.3.0: released on 2022-10-20.
- ver. 1.2.0: released on 2021-12-13.
- ver. 1.1.1: released on 2020-11-09.
- ver. 1.1.0: released on 2020-07-09.
- ver. 1.0.0: released on 2020-04-17.
- ver. 1.0-beta: released on 2020-03-30.
- ver. 0.1: released on 2019-12-04.

License
==================

This package is distributed under GNU General Public License version 3 (GPL v3) or later.

Papers
========

When you publish the results by using TeNeS, we would appreciate if you cite the following paper:

- `Y. Motoyama, T. Okubo, K. Yoshimi, S. Morita, T. Kato, and N. Kawashima, "TeNeS: Tensor Network Solver for Quantum Lattice Systems", Comput. Phys. Commun. 279, 108437 (2022) <https://www.sciencedirect.com/science/article/pii/S0010465522001564>`_
- `Y. Motoyama, T. Okubo, K. Yoshimi, S. Morita, T. Aoyama, T. Kato, and N. Kawashima, "TeNeS-v2: Enhancement for real-time and finite temperature simulations of quantum many-body systems", Comput. Phys. Commun. **315**, 109692 (2025) <https://www.sciencedirect.com/science/article/pii/S0010465525001948>`_

Copyright
==================

© *2019- The University of Tokyo. All rights reserved.*

This software was developed with the support of \"*Project for advancement of software usability in materials science*\" of The Institute for Solid State Physics, The University of Tokyo.
