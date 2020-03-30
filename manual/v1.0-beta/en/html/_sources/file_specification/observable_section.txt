.. highlight:: none


Define various settings related to physical quantity measurement.
This section has two types of subsections, ``onesite`` and ``twosite``.

``observable.onesite``
~~~~~~~~~~~~~~~~~~~~~~~~~

Define one-body operators that indicate physical quantities defined at each site :math:`i`.

.. csv-table::
   :header: "Name", "Description", "type"
   :widths: 15, 30, 20

   ``name``,     "Operator name",       String
   ``group``,    "Identification number of operators",   Integer
   ``sites``,    "Site number",         Integer or a list of integer
   ``dim``,      "Dimension of an operator",       Integer
   ``elements``, "Non-zero elements of an operator", String

``name``  specifies an operator name.

``group`` specifies an identification number of one-site operators.

``sites`` specifies a site number where an operater acts on.
By using a list, multiple site the operators can be defined on site at the same time.
An empty list ``[]`` means all sites.

``dim`` specifies a dimension of an operator.

``elements`` is a string specifying the non-zero element of an operator.
One element is specified by one line consisting of two integers and two floating-point numbers separated by spaces.

- The first two integers are the state numbers before and after the act of the operator, respectively.
- The latter two floats indicate the real and imaginary parts of the elements of the operator, respectively.

Example
.......

As an example, the case of :math:`S^z` operator for S=1/2 

.. math::
   S^z = \left(\begin{array}{cc} 0.5 & 0.0 \\ 0.0 & -0.5 \end{array}\right)

is explained.


First, set the name to ``name = "Sz"`` and the identification number to ``group = 0``.

Next, if the same operator is used at all sites, set ``sites = []``.
Otherwise, for example, if there are sites with different spin length :math:`S`, specify a specific site number such as ``sites = [0,1]``.

The dimension of the operator is ``dim = 2``, because it is the size of the matrix shown above.

Finally, the operator element is defined.
When we label two basis on site as :math:`|\uparrow\rangle = |0\rangle` and :math:`|\downarrow\rangle = |1\rangle`,
non-zero elements of :math:`S^z` are represented as

::

  elements = """
  0 0   0.5 0.0
  1 1  -0.5 0.0
  """

As a result, :math:`S^z` operator for S=1/2 is defined as follows:
::

  [[observable.onesite]]
  name = "Sz"
  group = 0
  sites = []
  dim = 2
  elements = """
  0 0  0.5 0.0
  1 1  -0.5 0.0
  """


``observable.twosite``
~~~~~~~~~~~~~~~~~~~~~~~~~

Define two-body operators that indicate physical quantities defined on two sites.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``name``,     "Operator name",                      String
   ``group``,    "Identification number of operators", Integer
   ``bonds``,    "Bond",                               String
   ``dim``,      "Dimension of an operator",           Integer
   ``elements``, "Non-zero elements of an operator",   String
   ``ops``,      "Index of onesite operators",         A list of integer


``name``  specifies an operator name.

``group`` specifies an identification number of two sites operators.

``bonds`` specifies a string representing the set of site pairs on which the operator acts.
One line consisting of three integers means one site pair.

- The first integer is the number of the source site.
- The last two integers are the coordinates (dx, dy) of the other site (target site) from the source site.

    - Both dx and dy must be in the range :math:`-3 \ le dx \ le 3`.

``dim`` specifies a dimension of an operator. 
In other words, the number of possible states of the site where the operator acts on.
In the case of interaction between two :math:`S=1/2` spins, for example, ``dim = [2, 2]`` .

``elements`` is a string specifying the non-zero element of an operator.
One element consists of one line consisting of four integers and two floating-point numbers separated by spaces.

- The first two integers are the status numbers of the source site and target site **before** the operator acts on.
- The next two integers show the status numbers of the source site and target site **after** the operator acts on.
- The last two floats indicate the real and imaginary parts of the elements of the operator.

Using ``ops``, a two-body operator can be defined as a direct product of the one-body operators defined in ``observable.onesite``.
For example, if :math:`S^z` is defined as ``group = 0`` in ``observable.onesite``,  :math:`S ^ z_iS ^ z_j` can be expressed as ``ops = [0,0]``.

If both ``elements`` and ``ops`` are defined, the process will end in error.

Example
.......

As an example, the bond Hamiltonian for S=1/2 Heisenberg model on square lattice at ``Lsub=[2,2]`` 

.. math::
  \mathcal{H}_{ij} = S_i^z S_j^z + \frac{1}{2} \left[S_i^+ S_j^- + S_i^- S_j^+ \right]

is explained.

First, the name and identification number is set as ``name = "hamiltonian"`` and ``group = 0``.
``dim = [2,2]`` because the state of each site is a superposition of the two states 
:math:`|\uparrow\rangle` and :math:`|\downarrow\rangle`.

Next, let's define the bonds. In this case, the site numbers are assigned as  ::

  2 3
  0 1

The bond connecting 0 and 1 is represented as ``0 1 0`` because 1 is located at (1,0) from 0.
Similarly, The bond connecting 1 and 3 is represented as ``1 0 1`` because 3 is located at (0,1) from 1.

Finally, how to define the elements of the operator is explained.
First, the basis of the site is needed to be labeled. Here, we label :math:`|\uparrow\rangle` as 0 and :math:`|\downarrow\rangle` as 1.
Using this basis and label number, for example, one of diagonal elements :math:`\left\langle \uparrow_i \uparrow_j | \mathcal{H}_{ij} | \uparrow_i \uparrow_j \right\rangle = 1/4` is specified by ``0 0 0 0 0.25 0.0``.
Likewise, one of off-diagonal elements :math:`\left\langle \uparrow_i \downarrow_j | \mathcal{H}_{ij} | \downarrow_i \uparrow_j \right\rangle = 1/2` is specified by ``1 0 0 1 0.5 0.0``.

As a result, the Heisenberg Hamiltonian for S=1/2 is defined as follows:
::

  [[observable.twosite]]
  name = "hamiltonian"
  group = 0
  dim = [2, 2]
  bonds = """
  0 0 1
  0 1 0
  1 0 1
  1 1 0
  2 0 1
  2 1 0
  3 0 1
  3 1 0
  """
  elements = """
  0 0 0 0  0.25 0.0
  1 0 1 0  -0.25 0.0
  0 1 1 0  0.5 0.0
  1 0 0 1  0.5 0.0
  0 1 0 1  -0.25 0.0
  1 1 1 1  0.25 0.0
  """
