.. highlight:: none

Phase diagram of the hardcore boson model on a trianglar lattice 
-----------------------------------------------

Finally, let us consider a zero-temperature phase diagram of the hardcore boson model on a trianglar lattice.
The Hamiltonian of this model is given as

.. math::

   \begin{aligned}
   H = \sum_{\langle i,j \rangle} \Bigl[ -t (b_i^\dagger b_j + b_j^\dagger b_i) + V n_i n_j \Bigr] -\mu \sum_i n_i ,
   \end{aligned}

where \ :math:`\langle i, j\rangle`\  indicates a pair of the nearest-neighbor sites, \ :math:`\mu`\  is a chemical potential, \ :math:`t`\  is a hopping energy, \ :math:`V`\   is a strength of the nearest-neighbor interaction.
For a hardcore boson system, the maximum number of bosons at each site is restricted to 0 or 1.
It is known that several ordered phases characterized by two types of long-range order appear in this model :ref:`[Wessel] <Ref-Wessel>` .
One is a solid-like order which exists at a 1/3 filling, where one of three sites is filled in a \ :math:`\sqrt{3}\times\sqrt{3}`\  ordering with wave vector \ :math:`\boldsymbol{Q}=(4\pi/3,0)`\  (see the inset of :numref:`fig_tutorial6_hardcore_boson` ).
This long-range order is characterized by the structure factor \ :math:`S(\boldsymbol{Q})`\.
The other is a superfluid order which is characterized by the offdiagonal order parameter \ :math:`|\langle b \rangle|`\.

To perform calculation for this system, the user can use toml files named ``basic.toml`` , ``nn_obs.toml`` and a python script file ``run.py`` in the direction ``sample/06_hardcore_boson_triangular`` . Here, ``basic.toml`` specifies the model and its parameters. This file is almost the same as the triangular Heisenberg model described in the previous section and differs from it only in the section ``model`` in the last part, where the line ``type = "boson"`` specifies the hardcore boson model and ``t = 0.1``  , ``V = 1``  determines the strength of the hopping and nearest-neighbor interaction.

The ``nn_obs.toml`` file describes the structure factor to be measured.

.. literalinclude:: ../../../../sample/06_hardcore_boson_triangular/nn_obs.toml

Using this file, the structure factor \ :math:`S(\boldsymbol{Q})`\  can be calculated.

Let us execute calculation using the script ``run.py`` .
After setting of the paths, execute calculation by typing the following command:

::

    python run.py

The calculation will finish within several minutes or several tens of minutes.
After the calculation, start ``gnuplot`` and execute the following command:

::

    load 'plot.gp'

Then, we obtain a graph for the structure factor \ :math:`S(\boldsymbol{Q})`\  and the superfluid order parameter \ :math:`|\langle b \rangle|`\  as shown in :numref:`fig_tutorial6_hardcore_boson` (a).
We note that this calculation is not so accurate because the bond dimension used in the calculation is small.
By commenting out the four lines in the beginning of the script ``run.py`` , we can perform more accurate calculation at the expense of execution time.
The result is shown in :numref:`fig_tutorial6_hardcore_boson` (b).
From these figures, we find that there exists three phases for the ground state, i.e., (a) a superfluid phase (\ :math:`-0.5 \lesssim \mu/V \lesssim -0.2`\), (b) a solid phase (\ :math:`-0.2 \lesssim \mu/V \lesssim 2.4`\), and (c) a supersolid phase (\ :math:`2.4 \lesssim \mu/V`\).
This result is consistent with the previous work :ref:`[Wessel] <Ref-Wessel>` .

.. figure:: ../../img/tutorial_6_hardcore_boson.*
	:name: fig_tutorial6_hardcore_boson
	:width: 600px

	Phase diagram of the hardcore boson model on a triangular lattice. (a) The bond dimension is 2. (b) The bond dimension is 5. Inset: a particle pattern for a solid phase.

.. rubric:: Reference

.. _Ref-Wessel:

[Wessel] 
S. Wessel, M. Troyer, *Supersolid hard-core bosonson the triangular lattice*, Phys. Rev. Lett. **95**, 127205 (2005). `link <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.127205>`__.
