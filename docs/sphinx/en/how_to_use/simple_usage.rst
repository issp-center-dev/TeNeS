.. highlight:: none

Usage of ``tenes_simple``
----------------------------

``tenes_simple`` is a tool that creates an input file of ``tenes_std`` for predefined models, lattices.


.. code:: bash

   $ tenes_simple simple.toml

- Takes a file as an argument
- Output input file for ``tenes_std``
- Other command line options are as follows
   - ``--help``
      - Show help message
   - ``--output=filename``
      - Specify the output file name ``filename``
      - Default is ``std.toml``
      - File name cannot be the same as the input file name

The currently defined models and lattices are as follows:

- Model
   - Spin system
- Lattice
   - Square lattice
   - Triangular lattice
   - Honeycomb lattice

See :ref:`sec-simple-format` for deails of the input file.
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
