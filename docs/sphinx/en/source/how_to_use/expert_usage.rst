.. highlight:: none

Usage of ``tenes``
------------------------------

If you prepare input files of ``tenes`` by yourself,
You can directry run ``tenes`` as follows:

.. code:: bash

 $ tenes --help
 TeNeS: TEnsor NEtwork Solver for 2D quantum lattice system
 
   Usage:
     tenes [--quiet] <input_toml>
     tenes --help
     tenes --version
 
   Options:
     -h --help       Show this help message.
     -v --version    Show the version.
     -q --quiet      Do not print any messages.

-  Take the input file name as an argument.
-  The command line options are:

   - ``help``
     - Show help messages.
   - ``version``
     - Show the version number.
   - ``quiet``
     - Do not print any messages to the standard output.

.. code:: bash

    tnsolve input.toml

