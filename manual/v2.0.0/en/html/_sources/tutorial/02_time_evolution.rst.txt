.. highlight:: none

Real-Time Evolution of the Transverse Field Ising Model
--------------------------------------------------------

Here, we introduce a calculation example for the real-time evolution of the Ising model on a square lattice when a transverse magnetic field, denoted by ``hx``, is applied.
The Hamiltonian is

.. math::

   \begin{aligned}
   H = J^z \sum_{\langle i,j \rangle} {S}_i^{z} {S}_j^{z} - h^x \sum_i S_i^x.
   \end{aligned}

Please note that the model is defined using spin operators of size 1/2, not Pauli operators.
The input and script files used in this tutorial can be found in ``sample/02_time_evolution``.

Initially, we compute the ground state (refer to the ``simple.toml`` file) which serves as our starting state. Specifically, it's set as:

.. literalinclude:: ../../../../sample/02_time_evolution/simple.toml

Given that ``Jz = -1.0``, the ground state becomes ferromagnetic. We use the ground state as the initial state, and save the state tensor with ``tensor_save = "save_tensor"``.

Next, we prepare the input file for the real-time evolution. This can be achieved by setting the ``mode`` to ``time``. Below is a sample input file (``simple_te_strong.toml``):

.. literalinclude:: ../../../../sample/02_time_evolution/simple_te_strong.toml

In this case, the transverse field is set to ``hx = 2.0``, and the time-step for evolution is ``tau = 0.01``.
Moreover, since we are utilizing the ground state as our initial condition, we load the state tensor with ``tensor_load = "save_tensor"``.

Once preparing the input file, we execute ``tenes_simple``, ``tenes_std``, and ``tenes`` in order. The results are saved in the ``output_te_strong`` directory.
Basically, the output is the same as the ground state, but with the addition of time in the first column.
For example, ``FT_density.dat`` records the expectation values of physical quantities over time::

    # The meaning of each column is the following: 
    # $1: time
    # $2: observable ID
    # $3: real
    # $4: imag
    # The meaning of observable IDs are the following: 
    # 0: Energy
    # 1: Sz              
    # 2: Sx              
    # 3: Sy              
    # 4: bond_hamiltonian
    # 5: SzSz            
    # 6: SxSx            
    # 7: SySy            

    0.00000000000000000e+00 0 -5.00184764052080899e-01  0.00000000000000000e+00
    0.00000000000000000e+00 1  4.99999945646528332e-01  0.00000000000000000e+00
    0.00000000000000000e+00 2  9.24306486797199186e-05  0.00000000000000000e+00
    0.00000000000000000e+00 3  2.34088935337348195e-06  0.00000000000000000e+00
    0.00000000000000000e+00 4 -5.00184764052080899e-01  3.47535331983321418e-21
    0.00000000000000000e+00 5  4.99999902788251294e-01 -8.46256269499545126e-22
    0.00000000000000000e+00 6  1.12653588020163689e-05  6.35907290717320676e-22
    0.00000000000000000e+00 7 -1.12840199341671039e-05 -2.06527532941704114e-21


The second column represents the type of physical quantity, and in this case, ``1`` represents the longitudinal magnetization :math:`m^z = \langle S^z \rangle`.
We can extract the time evolution of the magnetization by extracting rows with the second column equal to ``1``::

    awk '$2 == 1 {print $1, $3, $4}' output_te_strong/TE_density.dat > magnetization_strong.dat

For observing the time evolution with different transverse magnetic fields, we've also prepared sample input files named ``simple_te_middle.toml`` (``hx = 0.8``) and ``simple_te_weak.toml`` (``hx = 0.5``).
Additionally, there's a script named ``run.sh`` to execute these calculations in one go. Ensure that paths to tools like ``tenes`` are set correctly, and then execute the calculations with:

::

    sh run.sh

The computation will conclude in several seconds. Once done, launch gnuplot and enter:

::

    load 'plot.plt'

This will plot the temporal evolution of magnetization, :math:`S_z`. The result is displayed in :numref:`fig_tutorial7_timeevolution`.

.. figure:: ../../img/tutorial_07_timeevolution.*
    :name: fig_tutorial7_timeevolution
    :width: 600px

    Graph illustrating the real-time evolution of the Ising model. The vertical axis represents magnetization, and the horizontal axis represents time.


When the strength of the transverse-field exceeds the quantum phase transition point, the magnetization oscillates beyond 0 :ref:`[DQPT] <Ref-DQPT>`.

As time evolution progresses, the entanglement increases. At a certain point, the tensor network's capacity may be insufficient to express the wave function. In our case, the jump at ``t=4.25`` for ``hx=2.0`` indicates this issue. When applying this in practice, ensure no such discontinuities exist. If jumps are observed, steps like increasing the ``virtual_dimension`` might be necessary. For instance, adjusting it to ``virtual_dimension = 10`` and redoing the calculation as described above will eliminate the discontinuity, as can be seen in :numref:`fig_tutorial7_te_D10`.

.. figure:: ../../img/tutorial_07_timeevolution_D10.*
    :name: fig_tutorial7_te_D10
    :width: 600px

    Graph showcasing the real-time evolution of the Ising model. The vertical axis denotes magnetization, while the horizontal axis represents time. Results when ``virtual_dimension = 10`` are applied.


.. rubric:: Reference

.. _Ref-DQPT:
[DQPT]
M. Heyl, A. polkovnikov, and S. Kehrein, *Dynamical Quantum Phase Transitions in the Transverse-Field Ising Model*, Phys. Rev. Lett. **110**, 135704 (2013). `link <https://doi.org/10.1103/PhysRevLett.110.135704>`__
