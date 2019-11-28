.. highlight:: none

Input file for ``tense_simple`` 
---------------------------------

-  File format is
   `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__
   format.
-  The input file has five sections : ``model``, ``parameter``, ``lattice``, ``observable``, ``correlation`` .

-  In future, the file will be split by specifying a file name.

``parameter`` section
==========================

The contents of the ``parameter`` section are copied directly to the ``parameter`` section of the input file of ``tenes``.
You can also specify the imaginary time step size for the imaginary time evolution operator in simple and full updates in the subsections ``simple_update``, ``full_update``.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``tau``, imaginary time step in the imaginary time evolution operator, Real, 0.01

The following parameters are common to the ``tenes`` input file.

``parameter.tensor``
~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``D``,        The virtual bond dimensions of the central tensor,  Integer,   2
   ``CHI``,      The virtual bond dimensions of the angular transfer matrix,        Integer,   4
   ``save_dir``, Directory to write optimized tensors, Str, ""
   ``load_dir``, Directory to read initial tensor, Str, ""


- ``save_dir``
  - Store optimized tensors below this directory.
  - When it is empty, the tensors are not saved.
- ``load_dir``
  - Read various tensors from below this directory.
  - Must be same degree of parallelism as when saved.
  - Not read if it is empty.

``parameter.simple_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``num_step``,              Number of simple updates, Integer, 0
   ``inverse_lambda_cutoff``, cutoff of the mean field to be considered zero in the simple update, Real, 1e-12

``parameter.full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``num_step``,                 Number of full updates,  Integer,   0
   ``inverse_precision``,        Cutoff of singular values to be considered as zero when computing the pseudoinverse matrix with full update, Real,   1e-12
   ``convergence_epsilon``,      Convergence criteria for truncation optimization with full update, Real, 1e-12
   ``iteration_max``,            Maximum iteration number for truncation optimization on full updates,    Integer,   1000
   ``gauge_fix``,                Whether the tensor gauge is fixed, Boolean, true
   ``fastfullupdate``,           Whether the Fast full update is adopted, Boolean, true

``parameter.ctm``
~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``inverse_projector_cutoff``, Cutoff of singular values to be considered as zero when computing CTM projectors, Real,   1e-12
   ``convergence_epsilon``,      CTM convergence criteria,                                        Real,   1e-10
   ``iteration_max``,            Maximum iteration number of convergence for CTM,                           Integer,   100
   ``projector_corner``,         Whether to use only the 1/4 corner tensor in the CTM projector calculation,          Boolean, true
   ``use_rsvd``,                 Whether to replace SVD with Random SVD,                    Boolean, false
   ``rsvd_oversampling_factor``, ,                                                       Integer,   2


``parameter.random``
~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``seed``, Seed of the pseudo-random number generator used to initialize the tensor, Integer, 11

Example
~~~~~~~

::

    [parameter]
    [parameter.tensor]
    D  = 4     # tensor_dim
    CHI  = 16  # env_dim

    [parameter.simple_update]
    num_step = 1000

    [parameter.full_update]
    num_step = 1

    [parameter.ctm]
    iteration_max = 5


``lattice`` section
==========================

Specify the lattice information.
Square lattice and honeycomb lattice are defined as lattice types.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``type``, "Lattice name (square lattice or honeycomb lattice)", Str, --
   ``L_sub``, Unit cell size, An integer or a list of two integers, --


When a list of two integers is passed as ``L_sub``, the first element gives the value of ``Lx`` and the second one does ``Ly``.
If ``L_sub`` is an integer, Both ``Lx`` and ``Ly`` will have the same value.
A list of three or more elements causes an error.

Sites in a unit cell are indexed starting from 0.
These are arranged in order from the x direction.

Sites in a unit cell of ``L_sub = [2,3]`` are arranged as follows::

 y
 ^     4 5
 |     2 3
 .->x  0 1


Square lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

There are two types of bond, horizontal (0) and vertical (1) (corresponding to ``-`` and ``|`` in the below figure).

The unit cell for ``L_sub = 2`` is given as follows::

 0   1
 |   |
 2 - 3 - 2
 |   | 
 0 - 1 - 0


Honeycomb lattice
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Unit cell size (Each element of ``L_sub``) must be an even number.

There are 3 types of bonds: x, y, and z (corresponding to ``-``, ``~``, ``|``  in the below figure).
Each site with an even index has a rightward (x), a leftward (y), and an upward (z) bonds and
each site with an odd index has a leftward (x), a rightward (y), and a bottomward (z) bonds.

The unit cell for ``L_sub = 2`` is given as follows::

 0   1
     |
 2 ~ 3 - 2
 |   
 0 - 1 ~ 0


``model`` section
==========================

Specify the type of the model.
Spin system is only defined for ver. 0.1.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 15, 10

   ``type``, The type of the model, Str, --


Parameter names such as interaction depends on the model type.

Spin system
~~~~~~~~~~~~~~~~~~~~~

Spin system

.. math ::

 \mathcal{H} = \sum_{\langle ij \rangle}\left[\sum_\alpha^{x,y,z} J^\alpha_{ij} S^\alpha_i S^\alpha_j + B \left(\vec{S}_i\cdot\vec{S}_j\right)^2 \right] - \sum_i \left[ h S^z_i + \Gamma S^x_i - D \left(S^z_i\right)^2 \right]

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 10

   ``S``, "Magnituide of the local spin", Real, 0.5
   ``Jx``, "The x component of the exchange interaction J", Real or a list of Real, 1.0
   ``Jy``, "The y component of the exchange interaction J", Real or a list of Real, 1.0
   ``Jz``,"The z component of the exchange interaction J", Real or a list of Real, 1.0, Real or a list of Real, 1.0
   ``BQ``, "Biquadratic interaction B", Real or a list of Real, 0.0
   ``h``, "longitudinal magnetic field h", Real, 0.0
   ``G``, "Transverse magnetic field
 :math:`\Gamma` ", Real, 0.0
   ``D``, "On-site spin anisotropy D", Real, 0.0


By providing a list of exchange and biquadratic interactions, we can vary the magnitude of the interaction for each type of lattice bond.
If the number of elements in the list is less than the type of lattice bond, the remainder is filled in with the last element specified.


``observable`` section
==========================

By default, the local physical quantities used for physical quantities measurements: :math:`S^z`  and :math:`S^x` .
More detailed physical quantities measurements can be made by overwriting the format common to ``tenes``.
The following format is common to "tenes`.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``local_operator``,    Site opertor (ex. Sz),                      A list of Str
   ``hamiltonian``,       Bond hamiltonian,                           A list of Str
   ``hamiltonian_bonds``, Type of bond Hamiltonian and the set of bonds that act, Str

``local_operator``, ``hamiltonian``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Same as ``evolution.matrix`` .
The order you define corresponds exactly to the index of the operator Hamiltonian.

``hamiltonian_bonds``
~~~~~~~~~~~~~~~~~~~~~

Same as ``evolution.simple_update`` .

Example
~~~~~~~~

::

    [observable]
    local_operator = [
    """
      0.5  0.0
      0.0 -0.5
    """,
    """
      0.0 0.5
      0.5 0.0
    """,
    ]

    hamiltonian_bonds = """
    0 1 h 0
    3 2 h 0
    2 3 h 0
    1 0 h 0
    0 2 v 0
    3 1 v 0
    2 0 v 0
    1 3 v 0
    """

    hamiltonian = [
    """
      0.25   0.0    0.0     0.0
      0.0   -0.25   0.5     0.0  
      0.0    0.5   -0.25    0.0  
      0.0    0.0    0.0     0.25
    """,
    ]

``correlation`` section
==========================

For ``tenes_simple`` , correlation functions :math:`C = \langle A(0)B(r)\rangle` are not calculated by default.
For calculating correlation functions, it can be specified in the same file format as ``tenes``.
In the following, the parameters about correlation function are described.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``r_max``,    Maximum distance r of the correlation function, Integer
   ``operators``, "Numbers of operators A and B that measure correlation functions", A list for Integer


The operators defined in the ``observable`` section are used.

Example
~~~~~~~

::

    [correlation]
    r_max = 5
    operators = [[0,0], [0,1], [1,1]]
