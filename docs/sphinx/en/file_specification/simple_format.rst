.. highlight:: none

.. _sec-simple-format:

Input file for ``tense_simple`` 
---------------------------------

-  File format is
   `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__
   format.
-  The input file has four sections : ``model``, ``parameter``, ``lattice``, ``correlation`` .

   - The ``parameter`` section is copied to the standard mode input.

``model`` section
==========================

Specify the model to calculate.
Spin system (spin) is defined.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10

   ``type``, Model type, String, --


The parameter names such as interactions depend on the model type.

Spin system: spin
~~~~~~~~~~~~~~~~~~~~~

Spin system

.. math ::

 \mathcal{H} = \sum_{\langle ij \rangle}\left[\sum_\alpha^{x,y,z} J^\alpha_{ij} S^\alpha_i S^\alpha_j + B \left(\vec{S}_i\cdot\vec{S}_j\right)^2 \right] - \sum_i \left[ H S^z_i + \Gamma S^x_i - D \left(S^z_i\right)^2 \right]


The parameters of the one-body terms are defined as follows.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10

   ``S``, Magnituide of the local spin, Real (integer or half integer), 0.5
   ``H``, "longitudinal magnetic field :math:`H`", Real, 0.0
   ``G``, "Transverse magnetic field :math:`\Gamma` ", Real, 0.0
   ``D``, "On-site spin anisotropy :math:`D`", Real, 0.0

The exchange interaction :math:`J` can have a bond dependency.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10

   ``J0``, "Exchange interaction of 0th direction nearest neighbor bond", Real, 0.0
   ``J1``, "Exchange interaction of 1st direction nearest neighbor bond", Real, 0.0
   ``J2``, "Exchange interaction of 2nd direction nearest neighbor bond", Real, 0.0
   ``J0'``, "Exchange interaction of 0th direction next nearest neighbor bond", Real, 0.0
   ``J1'``, "Exchange interaction of 1st direction next nearest neighbor bond", Real, 0.0
   ``J2'``, "Exchange interaction of 2nd direction next nearest neighbor bond", Real, 0.0
   ``J0''``, "Exchange interaction of 0th direction third nearest neighbor bond", Real, 0.0
   ``J1''``, "Exchange interaction of 1st direction third nearest neighbor bond", Real, 0.0
   ``J2''``, "Exchange interaction of 2nd direction third nearest neighbor bond", Real, 0.0

The bond direction depends on the lattice defined in the ``lattice`` section.
For example, a square lattice can be defined for each of the two bond directions, x-direction (0) and y-direction (1).
By omitting the direction number, you can specify all directions at once.
You can also specify Ising-like interaction by adding one character of `xyz` at the end.
If the same bond or component is specified twice or more, an error will occur.

.. image:: ../../img/J.pdf
   :width: 400px
   :align: center

Biquadratic interaction  :math:`B` can also have a bond dependency like as :math:`J`.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10

   ``B0``, "Biquadratic interaction of 0th direction nearest neighbor bond", Real, 0.0
   ``B1``, "Biquadratic interaction of 1st direction nearest neighbor bond", Real, 0.0
   ``B2``, "Biquadratic interaction of 2nd direction nearest neighbor bond", Real, 0.0
   ``B0'``, "Biquadratic interaction of 0th direction next nearest neighbor bond", Real, 0.0
   ``B1'``, "Biquadratic interaction of 1st direction next nearest neighbor bond", Real, 0.0
   ``B2'``, "Biquadratic interaction of 2nd direction next nearest neighbor bond", Real, 0.0
   ``B0''``, "Biquadratic interaction of 0th direction third nearest neighbor bond", Real, 0.0
   ``B1''``, "Biquadratic interaction of 1st direction third nearest neighbor bond", Real, 0.0
   ``B2''``, "Biquadratic interaction of 2nd direction third nearest neighbor bond", Real, 0.0


:math:`S ^ z` and :math:`S ^ x` are automatically defined as one-site physical quantities used for physical quantity measurement.
If ``parameter.general.is_real = false``, :math:`S ^ y` is also defined. In addition, bond hamiltonian and :math:`S^z` ,  :math:`S^x`,  :math:`S^y` correlations with nearest neighobor bonds are automatically defined
as two-site physical quantities.

``lattice`` section
==========================

Specify the lattics to calculate.
Square, honeycomb, and triangular lattice are defined.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10

   ``type``, "lattice name (square, triangular or honeycomb lattice)", String, --
   ``L``, Unit cell size in x direction, Integer, --
   ``W``, Unit cell size in y direction, Integer, ``L``
   ``initial``, Inital tensor, String, "random"
   ``noise``, Noise for initial tensor, Real, 1e-2

The unit cell has a rectangular shape with the size of ``L`` times ``W``.
Sites in a unit cell are numbered sequentially from 0. They are arranged in order from the x direction.


``ex.) L = 2, W = 3``::

 y
 ^     4 5
 |     2 3
 .->x  0 1

``initial`` and ``noise`` are parameters that determine the initial state of the wave function.
In ``initial``, ``"ferro", "antiferro", "random"`` can be specified.
``noise`` is the amount of fluctuation given to the elements of the tensor.
If ``tensor_load`` is set in ``parameter.general``, ``initial`` is ignored.

square lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example of the ``type = square lattice``, we show the simple case of ``L=2``.
The definition of tensor alignment and nearest, next, and 3rd nearest neigbor bonds
is shown in :numref:`fig_square_1st`, :numref:`fig_square_2nd`, :numref:`fig_square_3rd`. 
The square represents the tensor, and the thin line represents the nearest neighbor bond.
Bold arrows and numbers represent two types of bonds.
The broken line represents one unit cell.

.. figure:: ../../img/Square_1st.pdf
   :width: 200px
   :align: center
   :name: fig_square_1st

   The nearest neighbor bond of the square lattice.

.. figure:: ../../img/Square_2nd.pdf
   :width: 200px
   :align: center
   :name: fig_square_2nd

   The 2nd nearest neighbor bond of the square lattice.

.. figure:: ../../img/Square_3rd.pdf
   :width: 200px
   :align: center
   :name: fig_square_3rd

   The 3rd nearest neighbor bond of the square lattice.

triangular lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example of the ``type = triangular lattice``, we show the simple case of ``L=2``.
The definition of tensor alignment and nearest, next, and 3rd nearest neigbor bonds
is shown in :numref:`fig_triangular_1st`, :numref:`fig_triangular_2nd`, :numref:`fig_triangular_3rd`.
The square represents the tensor, and the thin line represents the nearest neighbor bond.
Bold arrows and numbers represent two types of bonds.
The broken line represents one unit cell.

.. figure:: ../../img/Triangular_1st.pdf
   :width: 200px
   :align: center
   :name: fig_triangular_1st

   The nearest neighbor bond of the triangular lattice.

.. figure:: ../../img/Triangular_2nd.pdf
   :width: 200px
   :align: center
   :name: fig_triangular_2nd

   The 2nd nearest neighbor bond of the triangular lattice.

.. figure:: ../../img/Triangular_3rd.pdf
   :width: 200px
   :align: center
   :name: fig_triangular_3rd

   The 3rd neighbor bond of the triangular lattice.

honeycomb lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

As an example of the ``type = honeycomb lattice``, we show the simple case of ``L=2``.
The definition of tensor alignment and nearest, next, and 3rd nearest neigbor bonds
is shown in :numref:`fig_honeycomb_1st`, :numref:`fig_honeycomb_2nd`, :numref:`fig_honeycomb_3rd`.
In a honeycomb lattice, the unit cell size ``L`` in the x direction must be even number.
The square represents the tensor, and the thin line represents the nearest neighbor bond.
Bold arrows and numbers represent two types of bonds.
The broken line represents one unit cell.

.. figure:: ../../img/Honeycomb_1st.*
   :width: 200px
   :align: center
   :name: fig_honeycomb_1st
	  
   The nearest neighbor bond of the honeycomb lattice.

.. figure:: ../../img/Honeycomb_2nd.*
   :width: 200px
   :align: center
   :name: fig_honeycomb_2nd

   The 2nd nearest neighbor bond of the honeycomb lattice.

.. figure:: ../../img/Honeycomb_3rd.*
   :width: 200px
   :align: center
   :name: fig_honeycomb_3rd

   The 3rd nearest neighbor bond of the honeycomb lattice.


 
``parameter`` section
==========================

Parameters defined in this section is not used in ``tenes_simple`` but they are copied to the input file of ``tenes_std``.

.. include:: ./parameter_section.rst


``correlation`` section
==========================

For ``tenes_simple`` , correlation functions :math:`C = \langle A(0)B(r)\rangle` are not calculated by default.
For calculating correlation functions, they have to be specified in the same file format as the input file of ``tenes``.
For details, See ``correlation`` section :doc:`expert_format`.
