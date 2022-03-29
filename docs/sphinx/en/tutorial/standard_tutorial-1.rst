.. highlight:: none

Definition of lattices, models, and operators using the standard mode
----------------------------

While the simple mode is a tool to generate Hamiltonian and unit cell information from parameters of predefined models and lattices, the standard mode enables the definition of lattices, models, and operators. In this section, we explain how to use the standard mode.

.. figure:: ../../img/en_tutorial_1_Network.*
     :name: fig_transverse
     :width: 400px
     :align: center

     Tensor network of 4-site unit cell

Definition of unit cell (lattice)
----------------------------

.. figure:: ../../img/en_tutorial_1_Tensor.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [tensor] and [[tensor.unitcell]]

Unit cells are defined using [tensor] and [[tensor.unitcell]]:

.. code::

   [tensor]/*What kind of lattice did you define? 
   type = "square lattice" # Ignored (remnant of simple.toml)
   L_sub = [2, 2] # \ :math:`2\times 2`\ unitcell
   skew = 0 # Displacement in \ :math:`x`\-direction 

   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # Bond dimensions (\ :math:`\leftarrow~\uparrow~rightarrow~downarrow`\ order)
   index = [0, 3] #　Number indicating which tensor is in the unit cell
   physical_dim = 2 # Physical bond dimensions
   initial_state = [1.0, 0.0] # Initial state coefficients
   noise = 0.01 # Initial tensor fluctuation


The total initial state\ :math:`\psi`\ can be written in the direct product state of the per-site initial state \ :math:`\psi_i`\:\ :math:`| \Psi \rangle = \otimes_i |\Psi_i\rangle`\.
\ :math:`\psi_i`\ can be written as follows, with the elements of the initial_state array from the front as\ :math:`a_0,a_1,\cdots,a_{d-1}`\

.. math::

   \begin{aligned}
   |\Psi_i\rangle \propto \sum_k^{d-1}a_k|k\rangle\end{aligned}



Definition of model (Bond Hamiltonian)
----------------------------

.. figure:: ../../img/en_tutorial_1_Hamiltonian.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [[hamiltonian]]


Hamiltonians handled by TeNeS are the sum of Bond Hamiltonians (2-site Hamiltonians)
(Site Hamiltonians such as magnetic fields are also incorporated as Bond Hamiltonians)

.. math::

   \begin{aligned}
   mathcal{H} = \sum_{i,j}\mathcal{H}_{i,j}\end{aligned}

Consider a bond to be a pair of source and target sites

The Bond Hamiltonian is defined by its matrix elements and the bond it acts on.
Defining a matrix element enables us to define a model.
Defining a bond enables us to define a lattice.
It is prohibited for source and target to be tensors of the same number.


Definition of Bond Hamiltonian in std.toml

Definition of Bond Hamiltonian acting bond

.. code::

   [[hamiltonian]]
   dim = [2, 2] # Number of possible states of the acting bond [source, target]
   bonds = """ # Set of acting bonds (1 bond per line)
   0 1 0 # Row 1: Number of the source in the unit cell
   1 1 0 # Row 2: \ :math:`x`\coordinate(displacement) of target from the source
   2 1 0 # Row 3: \ :math:`y`\coordinate(displacement) of target from the source
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """

0 1 0 is 0 and right neighbor (1) (\ :math:`x+=1, y+=0`\)
1 0 1 is 1 and upper neighbor (3) (\ :math:`x+=0, y+=1`\)
1 1 0 is 0 and right neighbor (0) 


Definition of the matrix elements of the Bond Hamiltonian operator

.. code::

   elements = """ # (nonzero) matrix elements of the Hamiltonian (one element per row)
   0 0 0 0 0.25 0.0 # Row 1: state of source before action
   1 0 1 0 -0.25 0.0 # Row 2: state of target before action
   0 1 1 0 0.5 0.0 # Row 3: state of source after action
   1 0 0 1 0.5 0.0 # Row 4: state of target after action
   0 1 0 1 -0.25 0.0 # Row 5: Real part of element
   1 1 1 1 0.25 0.0 # Row 6: Imaginary part of element
   """

0 0 0 0 0.25 0.0 is \ :math:`\langle 00|\mathcal{H}_b|00\rangle=0.25`\
0 1 1 0 0.25 0.0 is \ :math:`\langle 10|\mathcal{H}_b|01\rangle=0.5`\



Definition of operators
----------------------------

.. figure:: ../../img/en_tutorial_1_Observable.*
     :name: fig_transverse
     :width: 400px
     :align: center

     [[observable.onesite]]


Definition of the operator whose expected value is finally computed
Currently calculable for 1-site and 2-site operators

Energy operator = Bond Hamiltonian also needs to be specified 
(tenes_std automatically copy it as a 2-site operator with the number 0)

Expression of the 1-site operator is

.. math::

   \begin{aligned}
   S^z = \begin{pmatrix}
   0.5 & 0.0 \\ 0.0 & -0.5
   \end{pmatrix}\end{aligned}
 


.. code::

   [observable]
   [[observable.onesite]] # 1-site operator
   name = "Sz" # Name
   group = 0 # 1-site operator identification number
   sites = [] # Number of tensors on which the 1-site operator acts ([] means all)
   dim = 2 # Dimensions of 1-site operators
   elements = """ # Non-zero elements of 1-site operator matrix (1 element per row)
   0 0 0.5 0.0 # Rows 1 and 2: before and after action
   1 1 -0.5 0.0 # Rows 3 and 4: Real and imaginary parts of the element
   """
   


Definition of the operator whose expected value is finally computed
Currently calculable for 1-site and 2-site operators

Energy operator = Bond Hamiltonian also needs to be specified 
(tenes_std automatically copy it as a 2-site operator with the number 0)

Expression of the two-site operator is

.. math::

   \begin{aligned}
   S^z_i S^z_j
   \end{aligned}
 


.. code::

   [[observable].twosite]] # 2-site operator
   name = "SzSz" # Name
   group = 1 # 2-site operator identification number (Independent of 1-site)
   dim = [2, 2] # Dimension
   bonds = """ # Bond that acts (site pair)
   0 1 0
   1 1 0
   2 1 0
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   ops = [0, 0] # When it can be written as a direct product of 1-site operators, its identification number
                # This time, "Sz" is the 1-site operator with number 0
                # Matrix elements can also be written explicitly as elements
                # (Similar format to Bond Hamiltonian)
   


   
   
Hamiltonian of the two-dimensional Heisenberg model of antiferromagnet
----------------------------

.. figure:: ../../img/en_tutorial_1_2DHeisenberg.*
     :name: fig_transverse
     :width: 400px
     :align: center

     Two-dimensional Heisenberg model of antiferromagnet

.. code::

   [tensor]/*What kind of lattice did you define? 
   type = "square lattice" # Ignored (remnant of simple.toml)
   L_sub = [2, 2] # \ :math:`2\times 2`\ unitcell
   skew = 0 # Displacement in \ :math:`x`\-direction 

   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] # Bond dimensions (\ :math:`\leftarrow~\uparrow~rightarrow~downarrow`\　order)
   index = [0, 3] #　Number indicating which tensor is in the unit cell
   physical_dim = 2 # Physical bond dimensions
   initial_state = [1.0, 0.0] # Initial state coefficients
   noise = 0.01 # Initial tensor fluctuation
   
   [[tensor.unitcell]]
   virtual_dim = [4, 4, 4, 4] 
   index = [1, 2] 
   physical_dim = 2 
   initial_state = [0.0, 1.0] 
   noise = 0.01 
   
   [[hamiltonian]]
   dim = [2, 2] # Number of possible states of the acting bond [source, target]
   bonds = """ # Set of acting bonds (1 bond per line)
   0 1 0 # Row 1: Number of the source in the unit cell
   1 1 0 # Row 2: \ :math:`x`\coordinate(displacement) of target from the source
   2 1 0 # Row 3: \ :math:`y`\coordinate(displacement) of target from the source
   3 1 0
   0 0 1
   1 0 1
   2 0 1
   3 0 1
   """
   elements = """ # (nonzero) matrix elements of the Hamiltonian (one element per row), Assume J=-1, h=1
   0 0 0 0 0.0 0.0 # Row 1: state of source before action
   1 0 1 0 0.25 0.0 # Row 2: state of target before action
   0 1 1 0 0.5 0.0 # Row 3: state of source after action
   1 0 0 1 0.5 0.0 # Row 4: state of target after action
   0 1 0 1 -0.75 0.0 # Row 5: Real part of element
   1 1 1 1 0.0 0.0 # Row 6: Imaginary part of element
   """
