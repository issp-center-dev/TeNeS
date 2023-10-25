.. highlight:: none

Real-Time Evolution of the Transverse Field Ising Model
--------------------------------------------------------

Here, we introduce a calculation example for the real-time evolution of the Ising model on a square lattice when a transverse magnetic field, denoted by ``hx``, is applied. The input and script files used in this tutorial can be found in ``sample/07_timeevolution``.

Initially, we compute the ground state (refer to the ``simple.toml`` file) which serves as our starting state. Specifically, it's set as:

.. literalinclude:: ../../../../sample/07_timeevolution/simple.toml

Given that ``Jz = -1.0``, the ground state becomes ferromagnetic. We use the ground state as the initial state, and save the state tensor with ``tensor_save = "save_tensor"``.

Next, we prepare the input file for the real-time evolution. This can be achieved by setting the ``mode`` to ``time``. Below is a sample input file (``simple_te_strong.toml``):

.. literalinclude:: ../../../../sample/07_timeevolution/simple_te_strong.toml

In this case, the transverse field is set to ``hx = 2.0``, and the time-step for evolution is ``tau = 0.01``. Moreover, since we are utilizing the ground state as our initial condition, we load the state tensor with ``tensor_load = "save_tensor"``. For observing the time evolution with different transverse magnetic fields, we've also prepared sample input files named ``simple_te_middle.toml`` and ``simple_te_weak.toml``. Additionally, there's a script named ``run.sh`` to execute these calculations in one go. Ensure that paths to tools like ``tenes`` are set correctly, and then execute the calculations with:

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

As time evolution progresses, the entanglement increases. At a certain point, the tensor network's capacity may be insufficient to express the wave function. In our case, the jump at ``t=4.25`` for ``hx=2.0`` indicates this issue. When applying this in practice, ensure no such discontinuities exist. If jumps are observed, steps like increasing the ``virtual_dimension`` might be necessary. For instance, adjusting it to ``virtual_dimension = 10`` and redoing the calculation as described above will eliminate the discontinuity, as can be seen in :numref:`fig_tutorial7_te_D10`.

.. figure:: ../../img/tutorial_07_timeevolution_D10.*
    :name: fig_tutorial7_te_D10
    :width: 600px

    Graph showcasing the real-time evolution of the Ising model. The vertical axis denotes magnetization, while the horizontal axis represents time. Results when ``virtual_dimension = 10`` are applied.