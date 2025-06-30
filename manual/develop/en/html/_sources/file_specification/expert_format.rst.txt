.. highlight:: none

.. _sec-expert-format:

Input file for ``tenes`` 
---------------------------------

-  File format is
   `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__
   format.
-  The input file has five sections: ``parameter``, ``tensor``, ``evolution``, ``observable``, ``correlation``.

``parameter`` section
========================

.. include:: ./parameter_section.rst


``tensor`` section
========================

.. include:: ./tensor_section.rst


``observable`` section
==========================

.. include:: ./observable_section.rst


``evolution`` section
========================

Specify the imaginary time evolution opetrators used in simple and full updates.
One-site and two-sites (nearest neighbor bond) operators can be defined.
This section has two subsections: ``simple`` and ``full``.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``group``,       "Group of the evolution operator", Integer (0-)
   ``site``,        "Index of site",                                                      Integer (0-)
   ``source_site``, "Index of source site",                                               Integer (0-)
   ``source_leg``,  "Direction from source site to target site",                          Integer (0-3)
   ``dimensions``,  "Dimension of a tensor of imaginary time evolution operator",         A list of integers
   ``elements``,    "Non-zero elements of a tensor of imaginary time evolution operator", String


``group`` specifies the group of the evolution operator.
It corresponds to the index of ``tau`` and ``num_steps`` in ``parameter.simple_update`` and ``parameter.full_update``.

``site`` is available for one-site operator, and ``source_site`` and ``source_leg`` are for two-site operator.

``source_leg`` is specified as an integer from 0 to 3.
Defined as ``0: -x, 1: + y, 2: + x, 3: -y`` in the clockwise order from the -x direction.

``dimensions`` is different from ``dim`` in ``observable`` section, so you need to specify the dimensions of all legs.
The order of the legs is ``source_initial, target_initial, source_final, target_final``, just like ``elements``.

Example :: 

    [evolution]
  
    # One site
    [[evolution.simple]]
    group = 0
    site = 0
    dimensions = [2, 2]
    elements = """
    0 0  1.0012507815756226 0.0
    1 1  0.9987507809245809 0.0
    """
  
    # Two site
    [[evolution.simple]]
    group = 0
    source_site = 0
    source_leg = 2
    dimensions = [2, 2, 2, 2]
    elements = """
    0 0 0 0  0.9975031223974601 0.0
    1 0 1 0  1.0025156589209967 0.0
    0 1 1 0  -0.005012536523536871 0.0
    1 0 0 1  -0.005012536523536871 0.0
    0 1 0 1  1.0025156589209967 0.0
    1 1 1 1  0.9975031223974601 0.0
    """


``correlation`` section
==========================

.. include:: ./correlation_section.rst

``correlation_length`` section
====================================

.. include:: ./correlation_length_section.rst
