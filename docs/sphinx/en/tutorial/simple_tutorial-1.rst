.. highlight:: none

Ising model with transverse magnetic field
---------------------------------------------

This section presents an example of calculation when the transverse magnetic field is applied to the Ising model.
By changing the variable  ``G`` in the input file,
the magnitude of the transverse magnetic field will be modified.
For example, when the transverse magnetic field is 0, the input file becomes

.. code::

   [parameter]
   [parameter.general]
   is_real = true

   [parameter.simple_update]
   num_step = 1000
   tau = 0.01

   [parameter.full_update]
   num_step = 0
   tau = 0.01

   [parameter.ctm]
   iteration_max = 10
   dimension = 10

   [lattice]
   type = "square lattice"
   L = 2
   W = 2
   virtual_dim = 2
   initial = "ferro"

   [model]
   type = "spin"
   Jz = -1.0
   Jx = 0.0
   Jy = 0.0
   G  = 0.0


In this case, since ``Jz = -1.0`` , the ferro magnetic state becomes the ground state at ``G=0``. 
When the input file name is ``simple.toml`` , type the following commands to execute ``tenes`` :
  
.. code:: bash

   $ tenes_simple simple.toml
   $ tenes_std std.toml
   $ tenes input.toml

Then, the following logs are output to the standard output.

.. code:: bash


	  Number of Processes: 1
	  Number of Threads / Process: 1
	  Tensor type: real
	  Start simple update
	  10% [100/1000] done
	  20% [200/1000] done
	  30% [300/1000] done
	  40% [400/1000] done
	  50% [500/1000] done
	  60% [600/1000] done
	  70% [700/1000] done
	  80% [800/1000] done
	  90% [900/1000] done
	  100% [1000/1000] done
	  Start calculating observables
	  Start updating environment
	  Start calculating onesite operators
	  Save onesite observables to output_0/onesite_obs.dat
	  Start calculating twosite operators
	  Save twosite observables to output_0/twosite_obs.dat
	  Save observable densities to output_0/density.dat
	  Save elapsed times to output_0/time.dat

	  Onesite observables per site:
	  Sz          = 0.5 0
	  Sx          = -1.28526262482e-13 0
	  Twosite observables per site:
	  hamiltonian = -0.5 0
	  SzSz        = 0.5 0
	  SxSx        = -1.7374919982e-18 0
	  SySy        = 1.73749202733e-18 0
	  Wall times [sec.]:
	  simple update = 3.545813509
	  full update   = 0
	  environmnent  = 0.123170523
	  observable    = 0.048149856

	  Done.


First, the information of parallelization and the tensors is displayed. Next, the execution status of the calculation process is displayed. After finishing the calculation, the expected values per site of the one-site operators ``Sz``, ``Sx`` and Hamiltonian, the nearest correlation ``SzSz``, ``SxSx``, ``SySy`` are output. Finally, the calculation time for each process is output (the unit is seconds). ``density.dat``, ``parameters.dat``, ``time.dat``, ``onesite_obs.dat``, and ``twosite_obs.dat`` are output to the output folder. For details on each output file, see File Format. For example, the value of ``<Sz>`` can be read from the second line of onesite_obs.dat. By changing ``G`` in increments of 0.1 from 0 to 3.0 and running ``tenes_simple`` and ``tenes``, the following result is obtained.
As an exapmle of the sample script, ``tutorial_example.py`` , ``tutorial_read.py`` are prepared in the ``sample/01_transverse_field_ising`` directory.
The calculation will be done by typing the following command:

.. code::

   $ python tutorial_example.py

For MacBook2017 (1.4 GHz Intel Core i7), the calculation was finished in  a few minutes. By typing the following command, G, energy, ``<Sz>`` and ``<Sx>`` are ouputted in the standard output:

.. code::

   $ python tutorial_read.py



.. image:: ../../img/tutorial_1_Sz_vs_G.*
   :width: 400px
   :align: center


As you can see from the figure, with increasing ``G``, the ``<Sz>`` decreases from ``0.5`` to ``0``, while the ``<Sx>`` increases from ``0`` to ```0.5``.
