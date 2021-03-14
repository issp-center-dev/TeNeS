.. highlight:: none

Set various parameters that appear in the calculation, such as the number of updates.
This section has five subsections: ``general``, ``simple_update``, ``full_update``,
``ctm``, ``random``.

Imaginary-time step :math:`\tau` for simple update ``parameter.simple_update.tau`` and that for full update ``parameter.full_update.tau`` are used only in standard mode ``tenes_std``, not used in ``tenes``.


``parameter.general``
~~~~~~~~~~~~~~~~~~~~~~~~~~

General parameters for ``tenes``.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 20, 30, 10, 10

   ``is_real``,     "Whether to limit all tensors to real valued ones",        Boolean, false
   ``iszero_tol``,  "Absolute cutoff value for reading operators",             Real,    0.0
   ``measure``,     "Whether to calculate and save observables",               Boolean, true
   ``output``,      "Directory for saving result such as physical quantities", String,  \"output\"
   ``tensor_save``, "Directory for saving optimized tensors",                  String,  \"\"
   ``tensor_load``, "Directory for loading initial tensors",                   String,  \"\"

- ``is_real``

  - When set to ``true``, the type of elements of the tensor becomes real. 
  - If one complex operator is defined at least,  calculation will end in errors before starting.

- ``iszero_tol``

  - When the absolute value of operator elements loaded is less than ``iszero_tol``, it is regarded as zero

- ``meaure``

  - When set to ``false``, the stages for measuring and saving observables will be skipped
  - Elapsed time ``time.dat`` is always saved

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

   ``tau``,           "Imaginary time step :math:`\tau` in imaginary time evolution operator", Real,    0.01
   ``num_step``,      "Number of simple updates",                                              Integer, 0
   ``lambda_cutoff``, "cutoff of the mean field to be considered zero in the simple update",   Real,    1e-12


..   ``gauge_fix``,     "Whether the tensor gauge is fixed",                                     Boolean, true

``parameter.full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~

Parameters in the full update procedure.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``tau``,                 "Imaginary time step :math:`\tau` in imaginary time evolution operator",                                       Real,    0.01
   ``num_step``,            "Number of full updates",                                                                                      Integer, 0
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
