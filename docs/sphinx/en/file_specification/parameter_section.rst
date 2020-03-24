.. highlight:: none

Set various parameters that appear in the calculation, such as the number of updates.
This section has four subsections: ``general``, ``simple_update``, ``full_update``,
``ctm``, ``random``.

Imaginary time increments for simple update ``parameter.simple_update.tau`` and that for full update ``parameter.full_update.tau`` are original parameters used in standard mode ``tenes_std``, not used in ``tenes``.


``parameter.general``
~~~~~~~~~~~~~~~~~~~~~~~~~~

Set general parametes for ``tenes``.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 20, 30, 10, 10

   ``is_real``, "Flag to limit all tensors to real", Boolean, false
   ``iszero_tol``, "Absolute cutoff value for loading tensors", Real, 0.0
   ``output``, "Directori name for outputting physical quantities etc.", String, \"output\"
   ``tensor_save``, "Directory name for saving optimized tensors", String, \"\"
   ``tensor_load``, "Directory name for loading initial tensors",       String, \"\"


- ``is_real``

  - When set to ``true``, the type of elements of the tensor becomes real. 
  - If one complex operator is defined at least,  calculation does not start. 

- ``iszero_tol``

  - When the absolute value of tensor elements loaded is less than ``iszero_tol``, it is regarded as zero

- ``output``

  - Save numerical results such as physical quantities to files in this directory
  - The default directory is the current directory.

- ``tensor_save``

  - Save optimized each tensor to files in this directory
  - Not save if empty 

- ``tensor_load``

  - Read each tensor from this directory
  - Must be set to the same degree of parallelism where the tensors were saved
  - Not load tensors if empty 

``parameter.simple_update``
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``tau``, "Imaginary time step in imaginary time evolution operator", Real, 0.01  
   ``num_step``,      "Number of simple updates",                                            Integer, 0
   ``lambda_cutoff``, "cutoff of the mean field to be considered zero in the simple update", Real,    1e-12

``parameter.full_update``
~~~~~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``tau``, "Imaginary time step in imaginary time evolution operator", Real, 0.01  	    
   ``num_step``,            "Number of full updates",                                                                                      Integer, 0
   ``env_cutoff``,          "Cutoff of singular values to be considered as zero when computing environment through full updates",          Real,    1e-12
   ``inverse_precision``,   "Cutoff of singular values to be considered as zero when computing the pseudoinverse matrix with full update", Real,    1e-12
   ``convergence_epsilon``, "Convergence criteria for truncation optimization with full update",                                           Real,    1e-12
   ``iteration_max``,       "Maximum iteration number for truncation optimization on full updates",                                        Integer, 1000
   ``gauge_fix``,           "Whether the tensor gauge is fixed",                                                                           Boolean, true
   ``fastfullupdate``,      "Whether the Fast full update is adopted",                                                                     Boolean, true

``parameter.ctm``
~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``dimension``, "Dimension of virtual bond in angular transfer matrix", Integer, 4
   ``projector_cutoff``,         "Cutoff of singular values to be considered as zero when computing CTM projectors",                          Real,    1e-12
   ``convergence_epsilon``,      "CTM convergence criteria",                                                                                  Real,    1e-10
   ``iteration_max``,            "Maximum iteration number of convergence for CTM",                                                           Integer, 100
   ``projector_corner``,         "Whether to use only the 1/4 corner tensor in the CTM projector calculation",                                Boolean, true
   ``use_rsvd``,                 "Whether to replace SVD with Random SVD",                                                                    Boolean, false
   ``rsvd_oversampling_factor``, "Ratio of the number of the oversampled elements to that of the obtained elements in the Random SVD method", Real,    2.0


``parameter.random``
~~~~~~~~~~~~~~~~~~~~~

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 30, 30, 10, 10 

   ``seed``, "Seed of the pseudo-random number generator used to initialize the tensor", Integer, 11

Each MPI process has the own seed as ``seed`` plus the process ID (MPI rank).

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
  
