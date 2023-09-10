.. highlight:: none

Set various parameters that appear in the calculation, such as the number of updates.
This section has five subsections: ``general``, ``simple_update``, ``full_update``,
``ctm``, ``random``.


``parameter.general``
~~~~~~~~~~~~~~~~~~~~~~~~~~

General parameters for ``tenes``.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 20, 30, 10, 10

   ``mode``,        "Calculation mode",                                        String, ``\"ground state\"``
   ``is_real``,     "Whether to limit all tensors to real valued ones",        Boolean, false
   ``iszero_tol``,  "Absolute cutoff value for reading operators",             Real,    0.0
   ``measure``,     "Whether to calculate and save observables",               Boolean, true
   ``measure_interval``, "Interval of measurement in time evolution process",  Integer or list of integers, 10
   ``output``,      "Directory for saving result such as physical quantities", String,  \"output\"
   ``tensor_save``, "Directory for saving optimized tensors",                  String,  \"\"
   ``tensor_load``, "Directory for loading initial tensors",                   String,  \"\"

- ``mode``

  - Specify the calculation mode
  - ``"ground state"``

    - Search for the ground state of the Hamiltonian
    - ``tenes_std`` calculates the imaginary time evolution operator :math:`U(\tau) = e^{-\tau H}` from the Hamiltonian :math:`H`

  - ``"time evolution"``

    - Calculate the time evolution of the observables from the initial state
    - ``tenes_std`` calculates the time evolution operator :math:`U(t) = e^{-it H}` from the Hamiltonian :math:`H`

  - ``"finite temperature"``

    - Calculate the finite temperature expectation values of the observables
    - ``tenes_std`` calculates the imaginary time evolution operator :math:`U(\tau) = e^{-\tau H}` from the Hamiltonian :math:`H`

- ``is_real``

  - When set to ``true``, the type of elements of the tensor becomes real. 
  - If one complex operator is defined at least,  calculation will end in errors before starting.

- ``iszero_tol``

  - When the absolute value of operator elements loaded is less than ``iszero_tol``, it is regarded as zero

- ``meaure``

  - When set to ``false``, the stages for measuring and saving observables will be skipped
  - Elapsed time ``time.dat`` is always saved

- ``measure_interval``

  - Specify the interval of measurement in time evolution process

- ``output``

  - Save numerical results such as physical quantities to files in this directory
  - Empty means ``"."`` (current directory)

- ``tensor_save``

  - Save optimized tensors to files in this directory
  - If empty no tensors will be saved

- ``tensor_load``

  - Read initial tensors from files in this directory
  - If empty no tensors will be loaded

``parameter.simple_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Parameters in the simple update procedure.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``tau``,           "(Imaginary) time step :math:`\tau` in (imaginary) time evolution operator", Real or list of real,    0.01
   ``num_step``,      "Number of simple updates",                                              Integer or list of integers, 0
   ``lambda_cutoff``, "cutoff of the mean field to be considered zero in the simple update",   Real,    1e-12
   ``gauge_fix``,     "Whether the tensor gauge is fixed",                                     Boolean, false
   ``gauge_maxiter``, "Maximum number of iterations for fixing gauge", Integer, 100
   ``gauge_converge_epsilon``, "Convergence criteria of iterations for fixing gauge", Real, 1e-2


- ``tau``

  - Specify the (imaginary) time step :math:`\tau` in (imaginary) time evolution operator

    - ``tenes_std`` uses it to calculate the imaginary time evolution operator :math:`e^{-\tau H}` from the Hamiltonian
    - ``tenes`` uses it to calculate the time of each measurement

      - For finite temperature calculation, note that the inverse temperature increase :math:`2\tau` at a step because :math:`\rho(\beta + 2\tau) = U(\tau)\rho(\beta)\bar{U}(\tau)`

  - When a list is specified, the time step can be changed for each group of time evolution operators

- ``num_step``

  - Specify the number of simple updates
  - When a list is specified, the number of simple updates can be changed for each group of time evolution operators

``parameter.full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~

Parameters in the full update procedure.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``tau``,                 "(Imaginary) time step :math:`\tau` in (imaginary) time evolution operator",                                       Real or list of reals,    0.01
   ``num_step``,            "Number of full updates",                                                                                      Integer or list of integers, 0
   ``env_cutoff``,          "Cutoff of singular values to be considered as zero when computing environment through full updates",          Real,    1e-12
   ``inverse_precision``,   "Cutoff of singular values to be considered as zero when computing the pseudoinverse matrix with full update", Real,    1e-12
   ``convergence_epsilon``, "Convergence criteria for truncation optimization with full update",                                           Real,    1e-6
   ``iteration_max``,       "Maximum iteration number for truncation optimization on full updates",                                        Integer, 100
   ``gauge_fix``,           "Whether the tensor gauge is fixed",                                                                           Boolean, true
   ``fastfullupdate``,      "Whether the fast full update is adopted",                                                                     Boolean, true

``parameter.ctm``
~~~~~~~~~~~~~~~~~

Parameters for corner transfer matrices, CTM.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``dimension``,                "Bond Dimension of CTM :math:`\chi`",                                                             Integer, 4
   ``projector_cutoff``,         "Cutoff of singular values to be considered as zero when computing CTM projectors",                          Real,    1e-12
   ``convergence_epsilon``,      "CTM convergence criteria",                                                                                  Real,    1e-6
   ``iteration_max``,            "Maximum iteration number of convergence for CTM",                                                           Integer, 100
   ``projector_corner``,         "Whether to use only the 1/4 corner tensor in the CTM projector calculation",                                Boolean, true
   ``use_rsvd``,                 "Whether to replace SVD with random SVD",                                                                    Boolean, false
   ``rsvd_oversampling_factor``, "Ratio of the number of the oversampled elements to that of the obtained elements in random SVD method", Real,    2.0
   ``meanfield_env``,            "Use mean field environment obtained through simple update instead of CTM", Boolean, false

For Tensor renomalization group approach using random SVD, please see the following reference, S. Morita, R. Igarashi, H.-H. Zhao, and N. Kawashima, `Phys. Rev. E 97, 033310 (2018) <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.97.033310>`_ .


``parameter.contraction``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Parameters for tensor contraction for calculating observables.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``optimize``, "Whether to optimize the contraction order of tensors", Boolean, true
   ``order_file``, "File of contraction order",                      String,  ""


When calculating physical quantities, it is necessary to contract the tensor network, and the computational cost depends on the order of contraction.
The contraction order is defined for each combination of the number of columns and rows of the bulk tensor in the network to be contracted, and dimensions of bonds.
The contraction order can be specified in three ways: dynamically during execution, using a predetermined order, or using the order provided by TeNeS.

``order_file`` specifies a file containing the contraction order already calculated.
The contraction order can be specified for each network shape, and for the specified network, it is used regardless of ``optimize``.

``optimize`` specifies whether to optimize dynamically.

- ``"always"``

   - always optimize the contraction order of tensors that are not specified in ``order_file``.

- ``"never"``

   - never optimize the contraction order of tensors that are not specified in ``order_file``.

- ``"automatic"`` (default value)

   - tensors that are not specified in ``order_file`` are optimized dynamically if the size of the virtual bond of the bulk tensor is the same, and otherwise optimized dynamically.

- ``"old"``

   - use the order used in versions before v1. This does not perform dynamic optimization.
   - ``order_file`` is not used.

The file specified by ``order_file`` can be generated by giving ``input.toml`` to the ``tenes_optimize_contraction_order`` command.
The filename is specified by ``order_file``.
When ``order_file`` is not specified, ``tenes_optimize_contraction_order`` stops with an error.

For details of the contraction order optimization algorithm, see
R. N. C. Pfeifer, J. Haegeman, and F. Verstraete: `Phys. Rev. E 90, 033315 (2014) <https://journals.aps.org/pre/abstract/10.1103/PhysRevE.90.033315>`_.

``parameter.random``
~~~~~~~~~~~~~~~~~~~~~

Parameters for random number generators.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``seed``, "Seed of the pseudo-random number generator used to initialize the tensor", Integer, 11

Each MPI process has the own seed as ``seed`` plus the process ID (MPI rank).

Example
~~~~~~~

::

  [parameter]
  [parameter.general]
  is_real = true
  [parameter.simple_update]
  num_step = 100
  tau = 0.01
  [parameter.full_update]
  num_step = 0  # No full update
  tau = 0.01
  [parameter.ctm]
  iteration_max = 10
  dimension = 9 # CHI
