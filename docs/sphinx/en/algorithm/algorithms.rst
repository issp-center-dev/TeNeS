###########################
Algorithm
###########################

Tensor Network States
===========================

In TeNeS, we consider two-dimensional infinite tensor product states (iTNS), which are natural extension of iMPS to higher dimensions. We assume a square lattice tensor network with a translational symmetry, whose diagram is shown as

.. image:: ../../img/iTPS.*
   :align: center

and try to find an approximate ground state wavefunction of two-dimensional quantum many-body systems. Notice that square lattice tensor networks can represent lattices other than the square lattice, such as the honeycomb and the triangular lattices, by considering proper mapping.


Contraction of iTPS
===========================
In order to calculate expectation values over a TNS, :math:`\langle \Psi|O|\Psi\rangle/\langle \Psi|\Psi\rangle`, generally we need to contract tensor networks corresponding to :math:`\langle \Psi|O|\Psi\rangle` and :math:`\langle \Psi|\Psi\rangle`. For example, a tensor network corresponding to :math:`\langle \Psi|\Psi\rangle` is given by

.. image:: ../../img/iTPS_braket.*
   :align: center

which is often called as double layered tensor network. The contraction of a double layered tensor network often need huge computation cost. In the case of MPS (and iMPS), fortunately, we can contract it efficiently, *e.g*, by considering a transfer matrix consist of local tensors. However, in the case of TPS (and iTPS), exact contraction is impossible except for small finite size systems (or infinite cylinders) and we often use approximate contraction methods. Among several efficient methods for contracting iTPS in two-dimension, TeNeS supports corner transfer matrix renormalization group (CTMRG) method :ref:`[CTMRG] <Ref-CTMRG>`, which expresses an infinitely extended double layered tensor network by using *corner transfer matrices* and *edge tensors*.

When we simplify the double layered tensor network by using a locally contracted tensor,

.. image:: ../../img/double_tensor.*
   :align: center

a tensor network diagram for the corner transfer matrix representation is given as

.. image:: ../../img/CTM.*
   :align: center

A corner transfer matrix and an edge tensor are defined as

.. image:: ../../img/CandE.*
   :align: center

The accuracy of the corner transfer matrix representation is determined by the bond dimension :math:`\chi` of corner transfer matrices, which is indicated as thick lines in the diagrams.

In the CTMRG algorithm, we iteratively optimise corner transfer matrices and edge tensors by *absorbing* local tensors until they converges. For example, an absorbing procedure, so called *left move*, is described as follows:

.. image:: ../../img/LeftMove.*
   :align: center

The projectors in the above diagram is calculated in several ways :ref:`[CTMRG] <Ref-CTMRG>` and they reduces the degree of freedoms to :math:`\chi`.

When we consider iTPS with the bond dimension :math:`D` and CTMs with the bond dimension :math:`\chi`, the leading computation cost of CTMRG scales as :math:`O(\chi^2 D^6)` and :math:`O(\chi^3 D^4)`. Notice that the bond dimension of double layered tensor network becomes :math:`D^2` by using locally contracted tensors. Thus, typically we increase :math:`\chi` as :math:`\chi \propto O(D^2)`. In this setup, the leading computation cost of CTMRG algorithm is reduced to :math:`O(D^{10})`, while the memory usage scales :math:`O(D^{8})`. In order to achive the computation cost discussed abobe, we need to use a partial singular value decomposition (SVD)  (or the truncated SVD) technique. When we use the full SVD insted of the partial SVD, the computation cost becomes :math:`O(D^{12})`. 

Once we obtain the corner transfer matrices and edge tensors, we can also calculate :math:`\langle \Psi|O|\Psi\rangle` efficiently. For example, a local magnetization :math:`\langle \Psi|S^z_i|\Psi\rangle` is represented as

.. image:: ../../img/Sz.*
   :align: center


and similarly the nearest neighbor correlation :math:`\langle \Psi|S^z_iS^z_{i+1}|\Psi\rangle` is represented as

.. image:: ../../img/SzSz.*
   :align: center

Notice that by using the second representation, we can calculate expectation values of any two-site operators. Although we can generalize such a diagram for any operators, the computation cost to contract the tensor network becomes huge for larger clusters.

Optimization of iTPS
===========================
In order to use iTPS as variational wavefunctions for the ground state, we need to optimize it so that it give us the minimum energy expectation value,

.. math::
   E = \frac{\langle \Psi|\mathcal{H}|\Psi\rangle}{\langle \Psi|\Psi\rangle},

where :math:`\mathcal{H}` represents the Hamiltonian of the target system. Among two types of popular optimization algorithms, the imaginary evolution (ITE) and the variational optimization, we support the ITE in TeNeS. In TeNeS, we consider approximate ITE within the iTPS ansatz:

.. math::
   |\Psi^{\mathrm{iTPS}} \rangle  \simeq e^{-T \mathcal{H}} |\Psi_0\rangle,

where :math:`|\Psi_0 \rangle` is an arbitrary initial iTPS. If :math:`T` is sufficiently large, the left hand side, :math:`|\Psi^{\mathrm{iTPS}}\rangle`, is expected to be a good approximation of the ground state.

In TeNeS, we assume that the Hamiltonian can be represented as a sum of short range two-body interactions as

.. math::
   \mathcal{H} = \sum_{\{(i,j)\}}H_{ij},

and apply Suzuki-Trotter decomposition to the ITE operator with small time step :math:`\tau`:

.. math::
   e^{-\tau \mathcal{H}} = \prod_{\{(i,j)\}} e^{-\tau H_{ij}} + O(\tau^2).

We can also consider higher order Suzuki-Trotter decomposition. By using the Suzuki-Trotter decomposition form, the ITE is represented as 

.. math::
   e^{-T \mathcal{H}} |\Psi_0\rangle = \left( \prod_{\{(i,j)\}} e^{-\tau H_{ij}}\right)^{N_{\tau}} |\Psi_0\rangle + O(\tau),

where :math:`N_{\tau} = T/\tau` is the number of ITEs with sufficiently small :math:`\tau`. In order to simulate the right hand side of the equation, we divide :math:`\prod_{\{(i,j)\}}` into several subsets. In each subset, (local) ITE operators satisfy two properties: they commute with each other and they have the same translation symmetry with the iTPS ansatz. For example, in the case of two-site iMPS for the one-dimensional nearest-neighbor interaction Hamiltonian, we have two subsets:

.. image:: ../../img/iMPS_ITE.*
   :align: center

Then, we approximate the wavefunction after multiplication of each ITE-operator subset as an iTPS with the bond dimension :math:`D`:

.. math::
   |\Psi_{\tau}^{\mathrm{iTPS}} \rangle  \simeq \prod_{\{(i,j) \in \mathrm{subset}_n \}}e^{-\tau H_{ij}} |\Psi^{\mathrm{iTPS}}\rangle,

where :math:`\prod_{\{(i,j) \in \mathrm{subset}_n \}}` means the product of operators in the :math:`n\mathrm{th}` subset, and :math:`|\Psi_{\tau}^{\mathrm{iTPS}}\rangle` is a new iTPS. By using a diagram, it is represented as follows:

.. image:: ../../img/iMPS_ITE_iMPS.*
   :align: center

Notice that by applying :math:`e^{-\tau H_{ij}}` the bond dimension of the exact iTPS representation generally increases. In order to continue the simulation stably, we need to *truncate* the bond dimension to a constant :math:`D`.
	   
Naively, efficient truncation can be done by solving the minimization problem

.. math::
   \min \left \Vert |\Psi_{\tau}^{\mathrm{iTPS}} \rangle -\prod_{\{(i,j) \in \mathrm{subset}_n \}} e^{-\tau H_{ij}} |\Psi^{\mathrm{iTPS}}\rangle \right \Vert^2.

However, in practice, solving this minimization problem needs huge computation cost because it is a highly nonlinear problem due to the translational symmetry of iTPS. Thus, instead, we usually consider an alternative local problem where we apply only a local ITE operator and try to find optimal iTPS :math:`|\Psi_{\tau}^{\mathrm{iTPS}}\rangle` in which only a few local tensors are modified from the original :math:`|\Psi^{\mathrm{iTPS}}\rangle`. This minimization problem is written as 

.. math::
   \min \left \Vert |\Psi_{\tau}^{\mathrm{iTPS}} \rangle - e^{-\tau H_{ij}} |\Psi^{\mathrm{iTPS}}\rangle \right \Vert^2.

In the case of the nearest-neighbor interaction on the one-dimensional chain, the diagrams corresponding to this minimization problems are 

.. image:: ../../img/iMPS_ITE_local.*
   :align: center

The squared norm :math:`\left \Vert |\Psi_{\tau}^{\mathrm{iTPS}} \rangle - e^{-\tau H_{ij}} |\Psi^{\mathrm{iTPS}}\rangle \right \Vert^2` can be calculated by using, *e.g.*, CTMRG and we can solve the minimization problem easily :ref:`[ITE] <Ref-ITE>` Although this new iTPS breaks translational symmetry, we make translationally symmetric iTPS by *copying* updated local tensors to other parts so that the obtained iTPS can be considered as an approximated solution of the original minimization problem:

.. image:: ../../img/Copy.*
   :align: center

This ITE approach is often called as *full update*. The leading computation cost of the full update come from CTMRG and then it scales as :math:`O(D^{10})` or :math:`O(D^{12})` depending on SVD algorithms.

The *simple update* (or *simplified update*) is a cheaper version of ITE optimization. In order to avoid expensive environment calculation by CTMRG, we consider a part of the tensor network instead to treat the whole :ref:`[SimpleUpdate] <Ref-SimpleUpdate>` in the simple update. For example, in the case of the nearest-neighbor interaction, we consider the following local optimization problem:

.. image:: ../../img/Simple_opt.*
   :align: center

In this diagram, :math:`\lambda_i` represents a non-negative diagonal matrix considered to be a mean field  corresponding to the neglected environment beyond the bond :math:`i`. The definition of :math:`\lambda_i` will be given later. This optimization problem can be viewed as the low rank approximation of a matrix consisting of two tensors and a ITE operator, and then we can solve it by SVD. The procedure of the simple update is given in the following diagram:

.. image:: ../../img/Simple_update.*
   :align: center

The singular values obtained from the SVD of the matrix is used as the mean field :math:`\lambda` in the next step. The computation cost of the simple update is :math:`O(D^{5})`, if we use QR decomposition before we construct the matrix :ref:`[QR] <Ref-QR>`. Thus, it is much cheaper that that of the full update.

Although the computation cost of the simple update is cheaper than that of the full update, it is known that the simple update shows strong initial state dependence and it tends to overestimate the local magnetization. Thus, for complicated problems, we need to carefully check results obtained by the simple update. 



.. rubric:: References

.. _Ref-TNS:

[TNS]
R. Orús, *A practical introduction to tensor networks: Matrix product states and projected entangled pair states*, Annals. of Physics **349**, 117 (2014). `link <https://linkinghub.elsevier.com/retrieve/pii/S0003491614001596>`__; R. Orús, *Tensor networks for complex quantum systems*, Nature Review Physics **1**, 538 (2019). `link <https://doi.org/10.1038/s42254-019-0086-7>`__. 

.. _Ref-MPS:

[MPS]
U. Schollwcök, *The density-matrix renormalization group in the age of matrix product states*, Annals. of Physics **326**, 96 (2011). `link <https://linkinghub.elsevier.com/retrieve/pii/S0003491610001752>`__

.. _Ref-CTMRG:

[CTMRG]
T. Nishino and K. Okunishi, *Corner Transfer Matrix Renormalization Group Method*, J. Phys. Soc. Jpn. **65**, 891 (1996).; R. Orús and G. Vidal, *Simulation of two-dimensional quantum systems on an infinite lattice revisited: Corner transfer matrix for tensor contraction*, Phys. Rev. B **80**, 094403 (2009). `link <https://doi.org/10.1103/PhysRevB.80.094403>`__ ; P. Corboz *et al.*, *Competing States in the t-J Model: Uniform d-Wave State versus Stripe State*, Phys. Rev. Lett. **113**, 046402 (2014). `link <https://doi.org/10.1103/PhysRevLett.113.046402>`__

.. _Ref-ITE:

[ITE]
J. Jordan *et al.*, *Classical Simulation of Infinite-Size Quantum Lattice Systems in Two Spatial Dimensions*, Phys. Rev. Lett. **101**, 250602, (2008). `link <https://doi.org/10.1103/PhysRevLett.101.250602>`__; R. Orús and G. Vidal, *Simulation of two-dimensional quantum systems on an infinite lattice revisited: Corner transfer matrix for tensor contraction*, Phys. Rev. B **80**, 094403 (2009). `link <https://doi.org/10.1103/PhysRevB.80.094403>`__

.. _Ref-SimpleUpdate:

[SimpleUpdate]
H. G. Jiang *et al.*, *Accurate Determination of Tensor Network State of Quantum Lattice Models in Two Dimensions*, Phys. Rev. Lett. **101**, 090603 (2008). `link <https://doi.org/10.1103/PhysRevLett.101.090603>`__

.. _Ref-QR:

[QR]
L. Wang *et al.*, *Monte Carlo simulation with tensor network states*, Phys. Rev. B **83**, 134421 (2011). `link <https://doi.org/10.1103/PhysRevB.83.134421>`__
