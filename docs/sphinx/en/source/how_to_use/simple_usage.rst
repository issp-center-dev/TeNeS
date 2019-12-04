.. highlight:: none

Usage of ``tenes_simple``
----------------------------

``tenes_simple`` is a tool that creates an input file of ``tenes`` for predefined models, lattices.

.. code:: bash

   $ ./tenes_simple --help
   usage: tenes_simple [-h] [-o OUTPUT] input

   Simple input generator for TeNeS

   positional arguments:
     input                 Input TOML file

   optional arguments:
     -h, --help            show this help message and exit
     -o OUTPUT, --output OUTPUT
                           Output TOML file

-  Take the input file name as an argument.
-  The command line options are:

   - ``help``
      - show help messages.
   - ``output``
      - output file name.
      - default name is ``input.toml``.
      - Cannot have the same file name as the input file name.
