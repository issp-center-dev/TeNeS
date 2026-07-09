.. highlight:: none

Phase diagram of the hardcore boson model on a trianglar lattice 
--------------------------------------------------------------------

Finally, let us consider a zero-temperature phase diagram of the hardcore boson model on a trianglar lattice.
The Hamiltonian of this model is given as

.. math::

   \begin{aligned}
   H = \sum_{\langle i,j \rangle} \Bigl[ -t (b_i^\dagger b_j + b_j^\dagger b_i) + V n_i n_j \Bigr] -\mu \sum_i n_i ,
   \end{aligned}

where :math:`\langle i, j\rangle` indicates a pair of the nearest-neighbor sites, :math:`\mu` is a chemical potential, :math:`t` is a hopping energy, :math:`V` is a strength of the nearest-neighbor interaction.
For a hardcore boson system, the maximum number of bosons at each site is restricted to 0 or 1.
It is known that several ordered phases characterized by two types of long-range order appear in this model :ref:`[Wessel] <Ref-Wessel>`.
One is a superfluid order which is characterized by the offdiagonal order parameter :math:`|\langle b \rangle|`.
The other one is a solid-like order which exists at a 1/3 filling, where one of three sites is filled in a :math:`\sqrt{3}\times\sqrt{3}` ordering with wave vector :math:`\boldsymbol{Q}=(4\pi/3,0)` (see the inset of :numref:`fig_tutorial6_hardcore_boson`).
This long-range order is characterized by the structure factor :math:`S(\boldsymbol{Q}) = \sum_{ij}^{N_\text{sites}} \langle n_i n_j \rangle \exp[-i\boldsymbol{Q}\cdot(r_i-r_j)] /N_\text{sites}`.

To perform calculation for this system, the user can use toml files named ``basic.toml``, ``nn_obs.toml`` and a python script file ``run.py`` in the direction ``sample/05_hardcore_boson_triangular``. 
Here, ``basic.toml`` specifies the model and its parameters.
This file is almost the same as the triangular Heisenberg model described in the previous section and differs from it only in the section ``model`` in the last part, where the line ``type = "boson"`` specifies the hardcore boson model and ``t = 0.1``, ``V = 1`` determines the strength of the hopping and nearest-neighbor interaction.

To calculate the structure factor :math:`S(\boldsymbol{Q})`, the correlations of densities at all the pairs of sites (in the unitcell) :math:`\langle n_i n_j \rangle` are required.
These observables are not defined by ``tenes_simple``, users should define them by themselves.
``nn_obs.toml`` specifies them for :math:`3 \times 3` unitcell, and the content of the file is appended to ``std__XXX_YYY.toml`` in ``run.py``.

For larger unitcell, the computational cost increases, and hence the calculatation of the structure factor gets difficult.
On the other hand, TeNeS can calculate the correlation functions along the x and y directions (in the square lattice iTPS) at a low cost.
Structure factor can be calculated using these correlation functions.
Additionally, the Fourier transform of the density operator :math:`\langle n_i \rangle` is also available.
The ground state at 1/3 filling is three-fold degenerate, but one of them is selected in the finite bond dimension calculation, and calculated density has site dependence.

Let us execute calculation using the script ``run.py``.
After setting of the paths, execute calculation by typing the following command:

::

    python run.py

The calculation will finish within several minutes or several tens of minutes.
After the calculation, start ``gnuplot`` and execute the following command:

::

    load 'plot.gp'

Then, we obtain a graph as shown in :numref:`fig_tutorial6_hardcore_boson`.
:math:`S(\boldsymbol{Q})` is the structure factor calculated from the density correlations at all the pairs of sites in the unitcell,
:math:`S'(\boldsymbol{Q})` is the structure factor calculated from the density correlations along the x-axis,
:math:`n(\boldsymbol{Q})` is the Fourier transform of the densities,
and :math:`(|\langle b \rangle| + |\langle b^\dagger \rangle|)/2` is the superfluid order parameter.
We note that this calculation is not so accurate because the bond dimension used in the calculation is small.
By increasing the bond dimensions specified in the beginning of the script ``run.py``, we can perform more accurate calculation at the expense of execution time.
From these figures, we find that there exists three phases for the ground state, i.e., (a) a superfluid phase (:math:`-0.5 \lesssim \mu/V \lesssim -0.2`), (b) a solid phase (:math:`-0.2 \lesssim \mu/V \lesssim 2.4`), and (c) a supersolid phase (:math:`2.4 \lesssim \mu/V`).
This result is consistent with the previous work :ref:`[Wessel] <Ref-Wessel>` .

.. figure:: ../../img/tutorial_6_hardcore_boson.*
	:name: fig_tutorial6_hardcore_boson
	:width: 600px

	Phase diagram of the hardcore boson model on a triangular lattice. Inset shows a particle pattern in a solid phase.

.. rubric:: Reference

.. _Ref-Wessel:

[Wessel] 
S. Wessel, M. Troyer, *Supersolid hard-core bosonson the triangular lattice*, Phys. Rev. Lett. **95**, 127205 (2005). `link <https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.95.127205>`__.
