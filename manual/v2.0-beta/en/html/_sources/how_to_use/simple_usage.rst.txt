.. highlight:: none

Usage of ``tenes_simple``
----------------------------

``tenes_simple`` is a tool that creates an input file of ``tenes_std`` for predefined models and lattices.


.. code:: bash

   $ tenes_simple simple.toml

- Takes a file as an argument
- Output an input file for ``tenes_std``
- Command line options are as follows

   - ``--help``

      - Show help message

   - ``--version``
     
      - Show version number

   - ``--output=filename``

      - Specify the output file name ``filename``
      - Default is ``std.toml``
      - File name cannot be the same as the input file name

   - ``--coordinatefile=coordfile``

      - Specify the output coordinate file name ``coordfile``
      - Default is ``coordinates.dat``
      - In a coordinate file, the first, second, and third columns denote site index, x coordinate, and y coordinate (in Cartesian), respectively.

   - ``--use-site-hamiltonian``

      - Onsite terms in Hamiltonian like Zeeman term will be output as site Hamiltonians
      - If not specified, these terms will be absorbed into the nearest neighbor bond Hamiltonians

The currently defined models and lattices are as follows:

- Model

   - Spin system

- Lattice

   - Square lattice
   - Triangular lattice
   - Honeycomb lattice
   - Kagome lattice

See :ref:`sec-simple-format` for details of the input file.
Below, a sample file for the S=1/2 Heisenberg model on the square lattice is shown.

::

   [lattice]
   type = "square lattice" # type of lattice
   L = 2                   # size of unitcell
   W = 2                   # size of unitcell
   virtual_dim = 3         # bond dimension
   initial = "antiferro"   # initial state

   [model]
   type = "spin" # type of model
   J = 1.0       # Heisenberg interaction

   [parameter]
   [parameter.general]
   is_real = true # use real tensor

   [parameter.simple_update]
   num_step = 1000  # number of steps
   tau = 0.01       # imaginary time step

   [parameter.full_update]
   num_step = 0    # number of steps
   tau = 0.01      # imaginary time step

   [parameter.ctm]
   dimension = 9       # bond dimension
