.. highlight:: none

Ising model with transverse magnetic field
---------------------------------------------

This section presents a calculation of the transverse magnetic field Ising model as an example.
By changing the variable ``G`` in the input file,
the magnitude of the transverse magnetic field will be modified.
For example, when the transverse magnetic field is 0, the input file is

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


In this case, since ``Jz = -1.0`` , the ferro magnetic state manifests itself as the ground state at ``G=0``. 
When the input file name is ``simple.toml`` , type the following commands to execute ``tenes``
(before typing them, please install TeNeS and set PATH properly.):
  
.. code:: bash

   $ tenes_simple simple.toml
   $ tenes_std std.toml
   $ tenes input.toml

Then, the following logs are output:

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


First, the information of parallelization and the tensors (complex or not) is displayed.
Next, the execution status of the calculation process is displayed.
After finishing the calculation, the expected values per site of the one-site operators ``Sz``, ``Sx`` and two-site ones Hamiltonian, the nearest correlation ``SzSz``, ``SxSx``, ``SySy`` are output.
Finally, the calculation time for each process is output in units of seconds.
``density.dat``, ``parameters.dat``, ``time.dat``, ``onesite_obs.dat``, and ``twosite_obs.dat`` are saved to the output directory.
For details on each output file, see :ref:`sec-output-format`.
For example, the value of ``<Sz>`` can be read from ``onesite_obs.dat``.
By changing ``G`` in increments of 0.2 from 0 to 3.0 and running ``tenes_simple`` and ``tenes``, the following result is obtained.
As an example of the sample script, ``tutorial_example.py`` , ``tutorial_read.py`` are prepared in the ``sample/01_transverse_field_ising`` directory.
The calculation will be done by typing the following command:

.. code::

   $ python tutorial_example.py

For MacBook2017 (1.4 GHz Intel Core i7), the calculation was finished in a few minutes.
By typing the following command, G, energy, ``<Sz>`` and ``<Sx>`` are outputted in the standard output:

.. code::

   $ python tutorial_read.py


.. figure:: ../../img/tutorial_1_Sz_vs_G.*
   :name: fig_transverse
   :width: 400px
   :align: center
   
   ``G`` dependence of ``<Sz>`` and ``<Sx>``.

As seen from :numref:`fig_transverse` , with increasing ``G``, the ``<Sz>`` decreases from ``0.5`` to ``0``, while the ``<Sx>`` increases from ``0`` to ``0.5``.

Magnetization process of the Heisenberg model on triangular and square lattices
--------------------------------------------------------------------------------

Next, we introduce the calculation of the magnetization process of the
quantum Heisenberg model with spin :math:`S = 1/2` defined on a
triangular lattice. The Hamiltonian looks like this:

.. math::

   \begin{aligned}
   H = J \sum_{\langle i,j \rangle}\sum_{\alpha}^{x,y,z} {S}_i^{\alpha} {S}_j^{\alpha} - \sum_i h S_i^z\end{aligned}

Here, :math:`\langle i, j \rangle` represents the pair of adjacent lattices, and :math:`h` represents the magnitude of the external magnetic field applied in the :math:`z` direction. 
Let’s calculate the ground state of this model and find :math:`\langle S_z \rangle\equiv \frac{1}{N_u}\sum_i^{N_u} \langle S_i^z \rangle`, where :math:`N_u` is the total number of sites in the unit cell, as a function of the magnetic field :math:`h`. To do this, use the toml file ``basic.toml`` in the ``sample/05_magnetization`` directory and the python script ``tutorial_magnetization.py``. 
The ``basic.toml`` file contains model settings and parameters.

::

    [parameter]
    [parameter.general]
    is_real = true

    [parameter.simple_update]
    num_step = 200
    tau = 0.01

    [parameter.full_update]
    num_step = 0
    tau = 0.01

    [parameter.ctm]
    iteration_max = 10
    dimension = 10

    [lattice]
    type = "triangular lattice"
    L = 3
    W = 3
    virtual_dim = 4
    initial = "random"

    [model]
    type = "spin"
    J = 1.0

The lattice section specifies a triangular lattice, and the unit cell
size specifies :math:`3\times 3`. Here, in order to make the calculation
lighter, only ``simple update`` is performed, and the imaginary time
interval :math:`\tau` is assumed to be :math:`\tau = 0.01`. For
simplicity, :math:`J =1`. Using this basic setting file,
tutorial_magnetization.py calculates the magnetization when the magnetic
field is swept.

::

    import subprocess
    from os.path import join
    import numpy as np
    import toml

    num_h = 21
    min_h = 0.0
    max_h = 5.0
    num_step_table = [100, 200, 500, 1000, 2000]

    fout = open("magnetization.dat","w")
    for idx, h in enumerate(np.linspace(min_h, max_h, num=num_h)):
        print("Caclulation Process: {}/{}".format(idx+1, num_h))
        inum = 0
        num_pre = 0
        fout.write("{} ".format(h))
        for num_step in num_step_table:
            ns = num_step - num_pre
            print("Step numter: {}".format(num_step))
            with open("basic.toml") as f:
                dict_toml = toml.load(f)
            dict_toml["parameter"]["general"]["output"] = "output_{}_{}".format(idx,num_step)
            dict_toml["parameter"]["general"]["tensor_save"] = "tensor_save_{}_{}".format(idx,num_step)
            dict_toml["model"]["H"] = float(h)
            dict_toml["parameter"]["simple_update"]["num_step"] = ns
            if inum > 0:
                dict_toml["parameter"]["general"]["tensor_load"] = "tensor_save_{}_{}".format(idx,num_pre)
            with open("simple_{}_{}.toml".format(idx,num_step), 'w') as f:
                toml.dump(dict_toml, f)
            cmd = "tenes_simple simple_{}_{}.toml -o std_{}_{}.toml".format(idx,num_step,idx,num_step)
            subprocess.call(cmd.split())
            cmd = "tenes_std std_{}_{}.toml -o input_{}_{}.toml".format(idx,num_step,idx,num_step)
            subprocess.call(cmd.split())
            cmd = "tenes input_{}_{}.toml".format(idx,num_step)
            subprocess.call(cmd.split())
            with open(join("output_{}_{}".format(idx,num_step), "density.dat")) as f:
                lines = f.readlines()
                mag_sz = lines[0].split('=')[1].strip()
            fout.write("{} ".format(mag_sz))
            inum = inum + 1
            num_pre = num_step
        fout.write("\n")
    fout.close()

In this script, the magnetic field :math:`h` is changed in steps of
:math:`0.25` from :math:`0` to :math:`5`, and the ground state energy
and :math:`\langle S_z \rangle` are calculated and output to ``energy.dat``
and ``magnetization.dat``. In order to see what happens when the number of
time steps for simple update is changed, calculations are also performed
with :math:`100`, :math:`200`, :math:`500`, :math:`1000`, and
:math:`2000` steps for each magnetic field. In order to reduce the
amount of calculation, the information of the wave function obtained
with a small number of steps is stored in tensor_save, and this is used
as the initial state for the calculation of a larger number of steps.
For example, the python script first performs a calculation with the
number of time steps set to 100, and output the result. Then, it perform
a calculation with the number of time steps set to 200 using the wave
function at the end of the calculation of the number of steps 100. The
script consequently reduce the amount of the calculation by 100 steps
for the latter.

Let’s actually run it. After passing through a path to tenes in advance,
execute calculation by typing as follows.

::

    python tutorial_magnetization.py

The calculation will finish within a few hours if you use a notebook PC
using a single processor. After the calculation is completed, start up
gnuplot and type

::

    load 'plot.gp'

to obtain the magnetization curve as shown in the right panel of
:numref:`fig_triangular`. In a similar way,

::

    load 'plot_ene.gp'

we obtain the ground-state energy as shown in the left panel of
:numref:`fig_triangular` .

As can be seen from the result for a sufficiently large number of steps
(for example, 2000 steps), a plateau structure occurs in the
magnetization process at the magnetization of :math:`1/3` of the
saturation magnetization :math:`\langle S_z \rangle = 0.5`. On this
plateau, spins on the three lattices form a periodic magnetic structure
with :math:`\uparrow`, :math:`\uparrow`, :math:`\downarrow`, and a spin
gap is generated. This plateau structure is unique to triangular
lattices. This plateau structure is unique to the triangular lattice. To
see whether the accuracy of calculation is enough or not, it is helpful
to chekc the step dependence of energy. In principle, the ground-state
energy should decrease as the number of steps increases, but in some
magnetic fields, the calculated energy increases. This is a sign that
the calculation accuracy is not good. It is presumed that it is
necessary to increase the bond dimension.

.. figure:: ../../img/Fig_Triangular.pdf
   :name: fig_triangular
   :width: 800px

   Ground state energy (left figure) and magnetization (right figure) of the Heisenberg model on the triangular lattice.

Next, let’s perform the calculation for a model on a square lattice. Use the toml file ``basic_square.toml`` and the python script ``tutorial_magnetization_square.py`` in the ``sample/05_magnetization`` directory.
The content of ``basic_square.toml`` is the same as ``basic.toml`` except that the ``lattice`` section has been changed as follows.

::

    [lattice]
    type = "square lattice"
    L = 2
    W = 2
    \begin{lstlisting}

    To perform the calculation, type
    \begin{lstlisting}
    python tutorial_magnetization.py

After the calculation is completed, start up gnuplot and type

::

    load 'plot_square.gp'

Then, the magnetization curve shown in the right panel of
:numref:`fig_square` is obtained. In a similar way, if you type

::

    load 'plot_ene_square.gp'

we obtain the ground-state energy as shown in the left panel of
:numref:`fig_square`. The calculation is almost converged at 2000
steps, and it can be seen that the plateau structure does not appear
unlike the triangular lattice Heisenberg model. Since the energy
generally decreases as the number of steps is increased, it is assumed
that the calculation accuracy is sufficiently high.

.. figure:: ../../img/Fig_Square.*
   :name: fig_square
   :width: 800px

   Ground state energy (left figure) and magnetization (right figure) of the Heisenberg model on the square lattice.

[fig2]

