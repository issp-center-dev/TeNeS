***************************
Usage
***************************

``tenes``, the main program of TeNeS, needs to create input file to define the model, order of operations, etc.
You can create the input file directly, but the following script is provided for ease of use:

- ``tenes_std`` : A tool that generates an input file to execute ``tenes``. An input file of ``tenes_std`` defines a lattice model etc. by yourself according to a predetermined format.

-  ``tenes_simple``: A tool that generates input files for ``tenes_std`` from another simpler input file which specifies lattice model predefined.
  
In order to simulate other models and/or lattices than predefined ones, you should create the input file of ``tenes_std`` and convert it.
Please see :doc:`../file_specification/index` for details on the input files of TeNeS.

The following sections describe how to use each script, and finally how to use ``tenes``.
  
.. toctree::
   :maxdepth: 2

   simple_usage
   standard_usage
   expert_usage
