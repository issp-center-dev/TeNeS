.. highlight:: none

Expert usage
------------------------------

If you prepare input files of ``tenes`` by yourself,
You can directry run ``tenes`` as follows:

.. code:: bash

    $ tenes --help
    TeNeS: PEPS+CTM method for solving 2D quantum lattice system
    Usage: tenes [OPTIONS] input_toml

    Positionals:
      input_toml TEXT REQUIRED    Input TOML file

    Options:
      -h,--help                   Print this help message and exit
      -v,--version                Show version information

-  Take the input file name as an argument.
-  The command line options are:

   -  ``help``
   -  ``version``

.. code:: bash

    tnsolve input.toml

