***************************
Usage
***************************

``TeNeS`` needs to create several input files to define models, order of operations, etc. You can create and run the input file directly, but the following script is provided for ease of use:

-  ``tenes_simple.py``: A tool that generates input files for ``tenes`` by using a simple input file which specifies lattice model predefined.

- ``tenes_standard.py`` (from ver. 1.0):  A tool that generates input files for ``tenes`` by using input files for generating the lattice model.
  
The following sections describe how to use each script, and finally how to use ``tenes``.
If you want to work with other models or grids, you can do so by creating the input file for ``tenes`` directly.
See :doc:`../file_specification/index` for details on the input file of ``tenes``.
  
.. toctree::
   :maxdepth: 2

   simple_usage
   standard_usage
   expert_usage
