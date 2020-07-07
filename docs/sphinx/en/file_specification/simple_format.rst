.. highlight:: none

.. _sec-simple-format:

Input file for ``tenes_simple`` 
---------------------------------

-  File format is
   `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__
   format.
-  The input file has four sections : ``model``, ``parameter``, ``lattice``, ``correlation`` .

   - The ``parameter`` section is copied to the standard mode input.

``model`` section
==========================

Specify the model to calculate.
In this version, spin system (``"spin"``) is defined.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 5, 40, 10, 10

   ``type``, Model type, String, --


The parameter names such as interactions depend on the model type.

Spin system: ``"spin"``
~~~~~~~~~~~~~~~~~~~~~~~~~

Hamiltonian is described as

.. math ::

 \mathcal{H} = \sum_{\langle ij \rangle}\left[\sum_\alpha^{x,y,z} J^\alpha_{ij} S^\alpha_i S^\alpha_j + B \left(\vec{S}_i\cdot\vec{S}_j\right)^2 \right] - \sum_i \sum_\alpha^{x,y,z} h^\alpha S^\alpha_i - \sum_i D \left(S^z_i\right)^2

The parameters of the one-body terms are defined as follows.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 5, 40, 10, 10

   ``S``, Magnituide of the local spin, Real (integer or half integer), 0.5
   ``hx``, "Magnetic field along :math:`S^x`, :math:`h^x`", Real, 0.0
   ``hy``, "Magnetic field along :math:`S^y`, :math:`h^y`", Real, 0.0
   ``hz``, "Magnetic field along :math:`S^z`, :math:`h^z`", Real, 0.0
   ``D``, "On-site spin anisotropy :math:`D`", Real, 0.0

The exchange interaction :math:`J` can have a bond dependency.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 5, 40, 10, 10

   ``J0``, "Exchange interaction of **0th** direction **nearest neighbor** bond", Real, 0.0
   ``J1``, "Exchange interaction of **1st** direction **nearest neighbor** bond", Real, 0.0
   ``J2``, "Exchange interaction of **2nd** direction **nearest neighbor** bond", Real, 0.0
   ``J0'``, "Exchange interaction of **0th** direction **next nearest neighbor** bond", Real, 0.0
   ``J1'``, "Exchange interaction of **1st** direction **next nearest neighbor** bond", Real, 0.0
   ``J2'``, "Exchange interaction of **2nd** direction **next nearest neighbor** bond", Real, 0.0
   ``J0''``, "Exchange interaction of **0th** direction **third nearest neighbor** bond", Real, 0.0
   ``J1''``, "Exchange interaction of **1st** direction **third nearest neighbor** bond", Real, 0.0
   ``J2''``, "Exchange interaction of **2nd** direction **third nearest neighbor** bond", Real, 0.0

The bond direction depends on the lattice defined in the ``lattice`` section.
For a square lattice, for example, coupling constants along two bond directions can be defined, x-direction (0) and y-direction (1).
By omitting the direction number, you can specify all directions at once.
You can also specify Ising-like interaction by adding one character of `xyz` at the end.
If the same bond or component is specified twice or more, an error will occur.

To summarize,

.. image:: ../../img/J.*
   :width: 400px
   :align: center

The biquadratic interaction :math:`B` can also have a bond dependency like as :math:`J`.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 5, 40, 10, 10

   ``B0``, "Biquadratic interaction of **0th** direction **nearest neighbor** bond", Real, 0.0
   ``B1``, "Biquadratic interaction of **1st** direction **nearest neighbor** bond", Real, 0.0
   ``B2``, "Biquadratic interaction of **2nd** direction **nearest neighbor** bond", Real, 0.0
   ``B0'``, "Biquadratic interaction of **0th** direction **next nearest neighbor** bond", Real, 0.0
   ``B1'``, "Biquadratic interaction of **1st** direction **next nearest neighbor** bond", Real, 0.0
   ``B2'``, "Biquadratic interaction of **2nd** direction **next nearest neighbor** bond", Real, 0.0
   ``B0''``, "Biquadratic interaction of **0th** direction **third nearest neighbor** bond", Real, 0.0
   ``B1''``, "Biquadratic interaction of **1st** direction **third nearest neighbor** bond", Real, 0.0
   ``B2''``, "Biquadratic interaction of **2nd** direction **third nearest neighbor** bond", Real, 0.0


One-site operators :math:`S ^ z` and :math:`S ^ x` are automatically defined.
If ``parameter.general.is_real = false``, :math:`S ^ y` is also defined.
In addition, bond Hamiltonian

.. math ::

 \mathcal{H}_{ij} = \left[\sum_\alpha^{x,y,z} J^\alpha_{ij} S^\alpha_i S^\alpha_j + B \left(\vec{S}_i\cdot\vec{S}_j\right)^2 \right] 
 - \frac{1}{z} \left[ \sum_\alpha^{x,y,z} h^\alpha \left(S^\alpha_i + S^\alpha_j \right) + D \left(\left(S^z_i\right)^2 + \left(S^z_j\right)^2 \right) \right],

and spin correlations with nearest neighbor bonds :math:`S^\alpha_iS^\alpha_j` ( :math:`\alpha=x,y,z` ) are automatically defined as two-site operators.
In the bond Hamiltonian, one body terms (:math:`h^\alpha` and :math:`D` term) appear only in the nearest neighbor bonds, and :math:`z` is the number of the coodinate number.

``lattice`` section
==========================

Specify the lattices to calculate.
Square, triangular, honeycomb, and Kagome lattices are defined.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"

   ``type``,        "lattice name (square,                triangular or honeycomb lattice)", String, --
   ``L``,           Unit cell size in x direction,        Integer,                           --
   ``W``,           Unit cell size in y direction,        Integer,                           ``L``
   ``virtual_dim``, Bond dimension,                       Integer,                           --
   ``initial``,     Inital state,                         String,                            "random"
   ``noise``,       Noise for elements in initial tensor, Real,                              1e-2

``initial`` and ``noise`` are parameters that determine the initial state of the wave function.
If ``tensor_load`` is set in ``parameter.general``, ``initial`` is ignored.

- ``initial``

   - ``"ferro"`` 

      - Ferromagnetic state with :math:`S^z = S`

   - ``"antiferro"``

      - Antiferromagnetic state.
        For square lattice and honeycomb lattice, the Neel order state (:math:`S^z = S` for the A sublattice and :math:`S^z = -S` for the B sublattice.)
        For triangular lattice and kagome lattice, the 120 degree order state (spins on sites belonging to the A, B, and C sublattice are pointing to :math:`(\theta, \phi) = (0,0), (2\pi/3, 0)` and :math:`(2\pi/3, \pi)` direction, respectively.)

   - ``"random"``

      - Random state.

- ``noise``

  - The amount of fluctuation in the elements of the initial tensor

Square lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A square lattice ``type = "square lattice"`` consists of ``L`` sites in the :math:`(1,0)` direction and ``W`` sites in the :math:`(0,1)` direction.
As a concrete example, :numref:`fig_square_lattice` (a) shows the structure for ``L=3, W=3``.
In addition, the definitions of the first, second and third nearest neighbor bonds are shown in
:numref:`fig_square_lattice` (b), (c), and (d), respectively.
The blue line represents a bond of ``bondtype = 0`` and the red line represents a bond of ``bondtype = 1``.

.. figure:: ../../img/SquareLattice.*
   :width: 550px
   :align: center
   :name: fig_square_lattice

   Square lattice.
   (a) Site structure with ``L=3, W=3``
   (b) Nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 0 degree direction and ``bondtype=1`` (red) one in the 90 degree direction.
   (c) Second nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 45 degree direction and ``bondtype=1`` (red) one in the -45 degree direction.
   (d) Third nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 0 degree direction and ``bondtype=1`` (red) one in the 90 degree direction.


Triangular lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A triangular lattice ``type = "triangular lattice"`` consists of ``L`` sites in the :math:`(1,0)` direction and ``W`` sites in the :math:`(1/2, \sqrt{3}/2)` direction.
As a concrete example, :numref:`fig_triangular_lattice` (a) shows the structure for ``L=3, W=3``.
In addition, the definitions of the first, second and third nearest neighbor bonds are shown in
:numref:`fig_triangular_lattice` (b), (c), and (d), respectively.
The blue, red, and green lines represent bonds of ``bondtype = 0, 1``, and ``2``, respectively.
(e) shows the corresponding square TPS with ``L=3, W=3``.

.. figure:: ../../img/TriangularLattice.*
   :width: 550px
   :align: center
   :name: fig_triangular_lattice

   Triangular lattice.
   (a) Site structure with ``L=3, W=3``
   (b) Nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 0 degree direction, ``bondtype=1`` (red) one in the 60 degree direction, and ``bondtype=2`` (green) one in the 120 degree direction.
   (c) Second nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 90 degree direction, ``bondtype=1`` (red) one in the -30 degree direction, and ``bondtype=2`` (green) one in the 30 degree direction.
   (d) Third nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 0 degree direction, ``bondtype=1`` (red) one in the 60 degree direction, and ``bondtype=2`` (green) one in the 120 degree direction.
   (e) Corrensponding square TPS of the triangular lattice with ``L=3, W=3``.

Honeycomb lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a honeycomb lattice ``type = "honeycomb lattice"``, units consisting of two sites of coordinates :math:`(0, 0)` and :math:`(\sqrt{3}/2, 1/2)` are arranged with ``L`` units in the :math:`(\sqrt{3},0)` direction and ``W`` units in the :math:`(1/2, 3/2)` direction.
As a concrete example, :numref:`fig_honeycomb_lattice` (a) shows the structure for ``L=2, W=2``.
In addition, the definitions of the first, second and third nearest neighbor bonds are shown in
:numref:`fig_honeycomb_lattice` (b), (c), and (d), respectively.
The blue, red, and green lines represent bonds of ``bondtype = 0, 1``, and ``2``, respectively.
(e) shows the corresponding square TPS with ``L=2, W=2``.

.. figure:: ../../img/HoneycombLattice.*
   :width: 550px
   :align: center
   :name: fig_honeycomb_lattice

   Honeycomb lattice.
   (a) Site structure with ``L=2, W=2``. The dashed ellipse denotes one unit.
   (b) Nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 30 degree direction, ``bondtype=1`` (red) one in the 150 degree direction, and ``bondtype=2`` (green) one in the -90 degree direction.
   (c) Second nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the 120 degree direction, ``bondtype=1`` (red) one in the 60 degree direction, and ``bondtype=2`` (green) one in the 0 degree direction.
   (d) Third nearest neighbor bonds. ``bondtype=0`` (blue) bond extends in the -30 degree direction, ``bondtype=1`` (red) one in the -150 degree direction, and ``bondtype=2`` (green) one in the 90 degree direction.
   (e) Corrensponding square TPS of the honeycomb lattice with ``L=2, W=2``. Note that the most top-right red tensor in the honeycomb lattice moves to the most top-left position, and the boundary condition is skewed.


Kagome lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In a kagome lattice ``type = "kagome lattice"``, units consisting of three sites of coordinates :math:`(0, 0)`, :math:`(1, 0)`, and :math:`(1/2, \sqrt{3}/2)` are arranged with ``L`` units in the :math:`(2,0)` direction and ``W`` units in the :math:`(1,\sqrt{3})` direction.
As a concrete example, :numref:`fig_kagome_lattice` (a) shows the structure for ``L=2, W=2``.
In addition, the definitions of the first, second and third nearest neighbor bonds are shown in
:numref:`fig_kagome_lattice` (b), (c), and (d), respectively.
The blue and the red lines represent bonds of ``bondtype = 0``, and ``1``, respectively.
(e) shows the corresponding square TPS with ``L=2, W=2``.

.. figure:: ../../img/KagomeLattice.*
   :width: 550px
   :align: center
   :name: fig_kagome_lattice

   Kagome lattice.
   (a) Site structure with ``L=2, W=2``. The dashed circle denotes one unit.
   (b) Nearest neighbor bonds. ``bondtype=0`` (blue) bonds form upper triangle and ``bondtype=1`` (red) bonds form lowertriangle.
   (c) Second nearest neighbor bonds.
   (d) Third nearest neighbor bonds. ``bondtype=0`` (blue) bond passes over a site and ``bondtype=1`` (red) one does not.
   (e) Corresponding square TPS of the kagome lattice with ``L=2, W=2``. The white circles are the dummy tensors with bonds of dimension one.

 
``parameter`` section
==========================

Parameters defined in this section is not used in ``tenes_simple`` but they are copied to the input file of ``tenes_std``.

.. include:: ./parameter_section.rst


``correlation`` section
==========================

For ``tenes_simple`` , correlation functions :math:`C = \langle A(0)B(r)\rangle` are not calculated by default.
For calculating correlation functions, they have to be specified in the same file format as the input file of ``tenes``.
For details, See ``correlation`` section :doc:`expert_format`.
