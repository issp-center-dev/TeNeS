.. highlight:: none

Finite Temperature Calculations for the Transverse Field Ising Model
----------------------------------------------------------------------

In this section, we present a calculation example of the ferromagnetic Ising model on a square lattice subjected to a transverse magnetic field, denoted by ``hx``, at finite temperatures. The input and script files used in this tutorial are located in the ``sample/08_finitetemperature`` directory. Below is a sample input file (`simple_ft_strong.toml`):

.. literalinclude:: ../../../../sample/08_finitetemperature/simple_ft_strong.toml

To perform finite temperature calculations, set the ``mode`` to ``finite``. Here, the transverse magnetic field is set to ``hx = 2.0`` with ``tau = 0.01`` (the inverse temperature step size is 2 times ``tau``). To observe the behavior at different transverse magnetic fields, we've provided additional sample input files: ``simple_ft_middle.toml``, ``simple_te_weak.toml``, and ``simple_ft_zero.toml``. Moreover, a script named ``run.sh`` has been set up to execute all these calculations simultaneously. Ensure you've added tools like ``tenes`` to your PATH, then initiate the calculations with:

::

    sh run.sh

The computation should complete in about a minute. To visualize the results, scripts have been prepared to plot energy, heat capacity, and magnetization (:math:`S_x`, :math:`S_z`): ``plot_e.plt``, ``plot_c.plt``, ``plot_mx.plt``, and ``plot_mz.plt``. Running the following:

::

    gnuplot -persist plot_e.plt
    gnuplot -persist plot_c.plt
    gnuplot -persist plot_mx.plt
    gnuplot -persist plot_mz.plt

will display plots for energy, heat capacity, and magnetizations (:math:`S_x` and :math:`S_z`). The resulting plots are illustrated in :numref:`fig_tutorial8_finitetemperature`. For comparison, results obtained using Quantum Monte Carlo calculations are also shown (using ``ALPS/looper``).

.. figure:: ../../img/tutorial_08_finitetemperature.*
    :name: fig_tutorial8_finitetemperature
    :width: 600px

    Graphs for the finite temperature calculations of the Ising model: (a) energy, (b) heat capacity, (c) :math:`S_x`, and (d) :math:`S_z`. The vertical axis represents the physical quantity, and the horizontal axis denotes temperature.