.. highlight:: none

Input file for ``tenes`` 
---------------------------------


-  File format is
   `TOML <https://github.com/toml-lang/toml/blob/master/versions/ja/toml-v0.5.0.md>`__
   format.
-  The input file has five sections: ``parameter``, ``lattice``, ``evolution``, ``observable``, ``correlation``

   -   In the future we will be able to split it by specifying the file name.

``parameter`` section
========================

Various parameters such as the number of updates are specified in this section. 
Only the parameters in this section have default values.
This section has five sub sections, ``tensor``, ``simple_update``, ``full_update``, ``ctm``, ``random``.

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
========================

Specify the ''unit cell'' information.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``L_sub``, Unit cell size, A list of integer

``evolution`` section
========================

Define the imaginary time evolution opetrators used in simple and full updates.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``matrix``,        Matrix representation about the imaginary time evolution opetrators, A list of string
   ``simple_update``, The order of the bonds that act on the index of the imaginary time evolution operator in simple update, A list of string
   ``full_update``,   The order of the bonds that act on the index of the imaginary time evolution operator in full update,   A list of string

``matrix``
~~~~~~~~~~

- One matrix is defined by a list of string.
- Columns are separated by one or more blanks, and rows are separated by one or more newlines.
- The order defined corresponds exactly to the number of the matrix. This order numbers are used to specify ``*_update`` (0-origin).

``*_update``
~~~~~~~~~~~~

-  One row represents one operator action.
-  Each line consists of four fields: ``int int string int``.

   1. A site to which bond connects
   2. A site to which bond connects
   3. Horizontal (h) or Vertical (v)
   4. Operator number (0-origin)

Example
~~~~~~~

.. code:: 

    [evolution]
    simple_update = """
    0 1 h 0
    3 2 h 0
    2 3 h 0
    1 0 h 0
    0 2 v 0
    3 1 v 0
    2 0 v 0
    1 3 v 0
    """

    full_update = """
    0 1 h 0
    3 2 h 0
    2 3 h 0
    1 0 h 0
    0 2 v 0
    3 1 v 0
    2 0 v 0
    1 3 v 0
    """

    matrix = [
    """
    0.9975031223974601 0.0 0.0 0.0
    0.0 1.0025156589209967 -0.005012536523536887 0.0
    0.0 -0.005012536523536888 1.0025156589209967 0.0
    0.0 0.0 0.0 0.9975031223974601
    """
    ]

``observable`` section
==========================

In this section, the information about observing physical quantities is specified.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``local_operator``,    Site opertor (ex. Sz),                      A list of string
   ``hamiltonian``,       Bond hamiltonian,                           A list of string
   ``hamiltonian_bonds``, Type of bond Hamiltonian and the set of bonds that act, string

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

In the following, the parameters about the correlation function :math:`C = \langle A(0)B(r) \rangle` are described.

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
