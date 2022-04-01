.. highlight:: none

Definition of lattices, models, and operators using the standard mode
----------------------------------------------------------------------------

By using the standard mode, users can define own lattices, models, and operators.
In this section, we explain how to use the standard mode.

Definition of unit cell
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../../img/en_tutorial_1_Tensor.*
     :width: 400px
     :align: center

     ``[tensor]`` and ``[[tensor.unitcell]]``

Unit cells are defined using ``[tensor]`` and ``[[tensor.unitcell]]``::

   [tensor]
   L_sub = [2, 2]          # 2x2 unitcell
   skew = 0                # Displacement in x direction 
                           # when go beyond a y-direction boundary

   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # Bond dimensions (←, ↑, →, ↓)
   index = [0, 3]             # Indices of tensors in the unit cell
   physical_dim = 2           # Physical bond dimensions
   initial_state = [1.0, 0.0] # Initial state coefficients
   noise = 0.01               # Fluctuation of elements in initial tensor


The initial state :math:`\ket{\Psi}` is prepared as the direct product state of the per-site initial states :math:`\ket{\psi}_i` : :math:`\ket{\psi} = \otimes_i \ket{\psi_i}`.
:math:`\ket{\psi}_i` can be specified as follows, with the elements of the ``initial_state = [a0, a1, ..., a_{d-1}]``,

.. math::

   \Ket{\psi}_i \propto \sum_{k=0}^{d-1}a_k\Ket{k}



Definition of model and lattice (Bond Hamiltonian)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../../img/en_tutorial_1_Hamiltonian.*
     :width: 400px
     :align: center

     [[hamiltonian]]

TeNeS treats the Hamiltonian as the sum of bond Hamiltonians (two-site Hamiltonians).
Note that site Hamiltonians such as the Zeeman terms are also included in the bond Hamiltonian.

.. math::

   \begin{aligned}
   mathcal{H} = \sum_{i,j}\mathcal{H}_{i,j}\end{aligned}

A bond is a directed pair of two sites: source and target.
The Bond Hamiltonian is defined by its (non-zero) matrix elements and the bond it acts on.
To define matrix elements is to define a model, and to define bonds is to define a lattice.

In the input file of the standard mode, say ``std.toml``, each bond Hamiltonian is specified as ``[[hamiltonian]]``.
The bonds where the bond Hamiltonian acts are specified by ``bonds`` string::

   [[hamiltonian]]
   bonds = """ # Set of acting bonds (1 bond per line)
   0 1 0       # Row 1: Number of the source in the unit cell
   1 1 0       # Row 2: x coordinate(displacement) of target from the source
   2 1 0       # Row 3: y coordinate(displacement) of target from the source
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """

One line of three integers corresponds one bond.
The first integer is the index of the source site.
The other two integers are x and y displacement of the target site from the source site.
For example, ``0 1 0`` means the pair of the site 0 and the right neighbor (x+=1 and y+=0), the site 1,
and ``1 0 1`` means the pair of the site 1 and the top neighbor (x+=0 and y+=1), the site 3.

The dimension of the bond Hamiltonian, i.e., the number of states of the source and target sites, is specified by ``dim``, and
the non-zero elements of the bond Hamiltonian is defined by ``elements``::

   dim = [2, 2]      # Number of possible states of the acting bond [source, target]
   elements = """    # (nonzero) matrix elements of the Hamiltonian (one element per row)
   0 0 0 0 0.25 0.0  # Field 1: State of source before action
   1 0 1 0 -0.25 0.0 # Field 2: State of target before action
   0 1 1 0 0.5 0.0   # Field 3: State of source after action
   1 0 0 1 0.5 0.0   # Field 4: State of target after action
   0 1 0 1 -0.25 0.0 # Field 5: Real part of element
   1 1 1 1 0.25 0.0  # Field 6: Imaginary part of element
   """

One line of ``elements`` corresponds one element.
The first two integers are the states of the source and target sites **before** the Hamiltonian acts on,
and
the following two integers are the states of the source and target sites **after** the Hamiltonian acts on.
The remaining two numbers are the real and imaginary part of the element of the bond Hamiltonian.


Definition of operators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. figure:: ../../img/en_tutorial_1_Observable.*
     :width: 400px
     :align: center

     ``[[observable.onesite]]``



Operators whose expected values are finally computed are defined in ``[observable]``.
The current version of TeNeS can evaluate onesite and twosites operators.
Although the energy operator is just the Hamiltonian (sum of the bond Hamiltonians),
it should be defined in ``[observable]`` if users want to calculate.
For convenience, ``tenes_std`` copies ``[[hamiltonian]]`` as ``[[observable.twosites]]`` when twosites operator with ``group = 0`` is not defined.

For an example of onesite operator, the z-component of the spin operator ::

   S^z = \begin{pmatrix}
   0.5 & 0.0 \\
   0.0 & -0.5 \\
   \end{pmatrix}
 
is defined as follows::

   [[observable.onesite]] # onesite operator
   name = "Sz"            # Name
   group = 0              # 1-site operator identification number
   sites = []             # Indices of tensors on which the operator acts ([] means all)
   dim = 2                # Dimensions of operators
   elements = """         # Non-zero elements of operator matrix (one element per line)
   0 0 0.5 0.0            # Fields 1 and 2: before and after action
   1 1 -0.5 0.0           # Fields 3 and 4: Real and imaginary parts of the element
   """
   
Non-zero elements of the matrix can be specified in the similar way of bond Hamiltonians.

Twosites operators can be defined in the similar way how to define the bond Hamiltonian.
As an example of twosites operator, spin-spin correlation on nearest neighbor bonds :math:`S^z_i S^z_j` is defined as follows::

   [[observable].twosite]] # twosite operator
   name = "SzSz"           # Name
   group = 1               # Index of twosite operator (independent of indices of onesite)
   dim = [2, 2]            # Dimension
   bonds = """             # Bond that acts on
   0 1 0
   1 1 0
   2 1 0
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   ops = [0, 0] # When it can be written as a direct product of onesite operators,
                # their indices.
                # In this case, "Sz" is the onesite operator with index 0.
                # Matrix elements can also be written explicitly as elements
   

When the twosites operator is written as the direct product of the two onesite operators,
``ops`` can be used to specify them.

   
Example: Antiferromagnetic Heisenberg model in staggered field
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Let us consider the antiferromagnetic Heisenberg model in staggered field.
The Hamiltonian is as follows

.. math::

   \mathcal{H} = J \sum_{\braket{ij}} S_i \cdot S_j - h \sum_{i \in A} S_i^z + h \sum_{j \in B} S_j^z,

where :math:`\sum_{\braket{ij}}` is summation over the nearest neighbor bonds
and :math:`A` and :math:`B` are the sublattices of the square lattice.
When splitting the Hamiltonian into the summation of bond Hamiltonians, :math:`\mathcal{H} = \sum_{\braket{ij}} \mathcal{H}_{ij}`
the bond Hamiltonian can be written as follows

.. math::

   \begin{split}
   \mathcal{H}_{ij}
   &= J S_i \cdot S_j - \frac{h}{4} S_i^z \otimes I_j + \frac{h}{4} I_i \otimes S_j^z \\
   &= \begin{pmatrix}
   J/4 &&& \\
   & (-J+h)/4 & J/2 & \\
   & J/2 & (-J-h)/4 & \\
   &&& J/4 \\
   \end{pmatrix},
   \end{split}

where :math:`I` is an identity operator, and the order of basis is 
:math:`\ket{\uparrow\uparrow}, \ket{\downarrow\uparrow}, \ket{\uparrow\downarrow}, \ket{\downarrow\downarrow}`.
It should be note that we divide the strength of the field, :math:`h`, by the number of the bonds connected to each site, :math:`z=4` in order to aviod double-counting.

When :math:`J = 0, h = 1`, for example, an input file of ``tenes_std``, ``std.toml``, is as follows::

   [parameter]
   [parameter.general]
   is_real = true
   tensor_save = "tensor"
   [parameter.simple_update]
   num_step = 1000
   tau = 0.01
   [parameter.full_update]
   num_step = 0
   tau = 0.01
   [parameter.ctm]
   dimension = 4
   iteration_max = 100

   [tensor]
   type = "square lattice"
   L_sub = [2, 2]
   skew = 0

   [[tensor.unitcell]]
   virtual_dim = [2, 2, 2, 2]
   index = []
   physical_dim = 2
   noise = 0.01

   [[hamiltonian]]
   dim = [2, 2]
   bonds = """
   0 1 0
   3 1 0
   0 -1 0
   3 -1 0
   0 0 1
   3 0 1
   0 0 -1
   3 0 -1
   """
   elements = """
   1 0 1 0  1.0 0.0
   0 1 0 1 -1.0 0.0
   """

   [observable]
   [[observable.onesite]]
   name = "Sz"
   group = 0
   sites = []
   dim = 2
   elements = """
   0 0 0.5 0.0
   1 1 -0.5 0.0
   """

   [[observable.onesite]]
   name = "Sx"
   group = 1
   sites = []
   dim = 2
   elements = """
   1 0 0.5 0.0
   0 1 0.5 0.0
   """

   [[observable.twosite]]
   name = "SzSz"
   group = 1
   bonds = """
   0 1 0
   0 0 1
   1 1 0
   1 0 1
   2 1 0
   2 0 1
   3 1 0
   3 0 1
   """
   dim = [2,2]
   ops = [0,0]

Note that in ``bonds`` of ``[[hamiltonian]]``, source sites (the first column) is always sites belonging to the A sublattice, 0 and 3.
We can calculate this model and obtain results as ::
   
   $ tenes_std std.toml
   $ tenes input.toml

      ... skipped ...

   Onesite observables per site:
     Sz          = 0 0
     Sx          = -1.32597e-18 0
   Twosite observables per site:
     hamiltonian = -2 0
     SzSz        = -0.5 0

      ... skipped

Especially, the expectation values of onesite operators written in ``output/onesite_obs.dat`` are::

   0 0 5.00000000000000000e-01 0.00000000000000000e+00
   0 1 -5.00000000000000000e-01 0.00000000000000000e+00
   0 2 -5.00000000000000000e-01 0.00000000000000000e+00
   0 3 5.00000000000000000e-01 0.00000000000000000e+00
   1 0 -4.60377857579530558e-18 0.00000000000000000e+00
   1 1 -1.39327011854595808e-18 0.00000000000000000e+00
   1 2 -4.60726081547908400e-18 0.00000000000000000e+00
   1 3 5.30041788535222114e-18 0.00000000000000000e+00

This means that spins on the A sublattice are up and those on the B are down.
By imposing the staggered magnetic field (:math:`J = 0, h = 1`), a tensor product state representing the Neel state is obtained.
Tensors are saved into the ``tensor`` directory because we set ``tensor_save = "tensor"``,
and therefore we can use them as the initial states of another calculation by setting ``tensor_load = "tensor"``.
