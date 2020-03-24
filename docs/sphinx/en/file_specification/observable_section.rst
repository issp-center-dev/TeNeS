.. highlight:: none


Define various settings related to physical quantity measurement.
This section has two types of subsections, ``onesite`` and ``twosite``.

``observable.onesite``
~~~~~~~~~~~~~~~~~~~~~~~~~

Define one-body operators that indicate physical quantities defined at each site.

.. csv-table::
   :header: "Name", "Description", "type"
   :widths: 15, 30, 20

   ``name``,     "Operator name",       String
   ``grouparrow``,    "Identification number of operators",   Integer
   ``sites``,    "Site number",         Integer or a list of integer
   ``dim``,      "Dimension of an operator",       Integer
   ``elements``, "Non-zero elements of an operator", String

``name``  specifies an operator name.

``grouparrow`` specifies an identification number of one-site operators.

``sites`` specifies a site number where an operater acts on.
By using a list, the operators can be defined multiple at the same time.
An empty list means all sites.

``dim`` specifies a dimension of an operator.

``elements`` is a string specifying the non-zero element of an operator.
One element consists of one line consisting of two integers and two floating-point numbers separated by spaces.
The first two are the state numbers before and after the operation of the operator, respectively.
The latter two indicate the real and imaginary parts of the elements of the operator, respectively.

Example
.......


As an example, the case of Sz operator for S=1/2 
.. math::
   S^z = \left(\begin{array}{cc} 0.5 & 0.0 \\ 0.0 & -0.5 \end{array}\right)

is explained.


First, set the name to ``name = "Sz"`` and the identification number to ``grouparrow = 0``.

Next, if the same operator is used at all sites, set ``sites = []``.
Otherwise, for example, if there are sites with different spin magnitudes, specify a specific site number such as ``sites = [0,1]``.

In this case, the dimension of the operator is ``dim = 2``, because it is the size of the matrix display shown above.

Finally, the operator element is defined.
For non-zero elements, the index (zero-based) and the elements can be ordered in order:

::

  elements = """
  0 0   0.5 0.0
  1 1  -0.5 0.0
  """

As a result, Sz operator for S=1/2 is defined as follows:
::

  [[observable.onesite]]
  name = "Sz"
  grouparrow = 0
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

   ``name``,     "Operator name",       String
   ``grouparrow``,    "Identification number of operators",   Integer
   ``bonds``,    "Bond",                   String
   ``dim``,      "Dimension of an operator",       Integer
   ``elements``, "Non-zero elements of an operator", String
   ``ops``,      "Identifcation number of onesite operator", A list of integer


``name``  specifies an operator name.

``grouparrow`` specifies an identification number of two sites operators.

``sites`` specifies a site number where an operater acts on.
By using a list, the operators can be defined multiple at the same time.
An empty list means all sites.


``bonds`` specifies a string representing the set of site pairs on which the operator acts.
One line consisting of three integers means one site pair.
The first integer is the number of the source site.
The last two integers are the coordinates (dx, dy) of the destination site (target) from the source site.
Both dx and dy must be in the range :math:`-3 \ le dx \ le 3`.

``dim`` specifies a dimension of an operator. In other words, the number of possible states of the site where the operator acts on.


``elements`` is a string specifying the non-zero element of an operator.
One element consists of one line consisting of four integers and two floating-point numbers separated by spaces.
The first two are the status numbers of the source site and target site before the operator acts on.
The next two show the status numbers of the source site and target site after the operator acts on.
The last two indicate the real and imaginary parts of the elements of the operator.

Using ``ops``, a two-body operator can be defined as a direct product of the one-body operators defined in ``observable.onesite``.
For example, if :math:`S ^ z` is defined as ``grouparrow = 0`` 
in ``observable.onesite``,  :math:`S ^ z_iS ^ z_j` can be expressed as ``ops = [0,0]``.

Example
.......

As an example, the bond Hamiltonian for S=1/2 Heisenberg model on square lattice at ``Lsub=[2,2]`` 

.. math::
  \mathcal{H}_{ij} = S_i^z S_j^z + \frac{1}{2} \left[S_i^+ S_j^- + S_i^- S_j^+ \right]

is explained.

First, the name and identification number is set as ``name = "hamiltonian"`` and ``grouparrow = 0``.
Since the state of each site is a superposition of the two states 
:math:`| \ uparrow \ rangle` and :math:`| \ downarrow \ rangle`, 
the dimension is 2, and ``dim = [2,2 ]``.


Next, let's see the bond. In this case, the site number is assigned as  ::

  2 3
  0 1

The bond connecting 0 and 1 is represented as ``0 1 0`` because 1 is located at (1,0) from 0.
Simiraly, The bond connecting 1 and 3 is represented as ``1 0 1`` because 3 is located at (0,1) from 1.

Finally, the elements of the operator is explained.
First, the basis of the site is needed to be labeled. Here, we set :math:`| \ uparrow \ rangle` to 0 and :math:` | \ downarrow \ rangle` to 1.
Using this basis and label number, for example, one of diagonal elements :math:`\left\langle \uparrow_i \uparrow_j | \mathcal{H}_{ij} | \uparrow_i \uparrow_j \right\rangle = 1/4` is specified by ``0 0 0 0 0.25 0.0``.
Likewise, one of off-diagonal elemensts :math:`\left\langle \uparrow_i \downarrow_j | \mathcal{H}_{ij} | \downarrow_i \uparrow_j \right\rangle = 1/2` is specified by ``1 0 0 1 0.5 0.0``.

As a result, the Heisenberg Hamiltonian for S=1/2 is defined as follows:
::

  [[observable.twosite]]
  name = "hamiltonian"
  grouparrow = 0
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

