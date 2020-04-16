.. highlight:: none

Specify the unit cell information (Information of bonds is given in the ``hamiltonian`` (``tenes_std``) and ``evolution`` (``tenes``)  sections.).
Unit cell has a shape of a rectangular with the size of ``Lx`` times ``Ly``.
``lattice`` section has an array of subsections ``unitcell`` .

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 20, 20, 15

   ``L_sub``, "Unit cell size", Integer or a list of integer, "--"
   ``skew``, "Shift value in skew boundary condition", Integer , 0

When a list of two integers is passed as ``L_sub``, the first element gives the value of ``Lx`` and the second one does ``Ly``.
A list of three or more elements causes an error.
If ``L_sub`` is an integer, Both ``Lx`` and ``Ly`` will have the same value.

Sites in a unit cell are indexed starting from 0.
These are arranged in order from the x direction.

.. figure:: ../../img/tensor_sec_fig1.*
   :width: 150px

   An example for ``L_sub = [2,3]``.

``skew`` is the shift value in the x direction when moving one unit cell in the y direction.

.. figure:: ../../img/tensor_sec_fig2.*
   :width: 400px

   An example for ``L_sub = [3,2], skew = 1`` (ruled line is a separator for unit cell).


``tensor.unitcell`` subsection
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The information of site tensors :math:`T_{ijkl\alpha}^{(n)}` is specified.
Here, :math:`i, j, k, l` indicates the index of the virtual bond, 
:math:`\alpha` indicates the index of the physical bond, 
and :math:`n` indicates the site number.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 20, 20

   ``index``,         "Site number",                                            Integer or a list of integer
   ``physical_dim``,  "Dimension of physical bond for a site tensor",           Integer
   ``virtual_dim``,   "Dimension of virtual bonds :math:`D` for a site tensor", Integer or a list of integer
   ``initial_state``, "Initial tensor",                                         a list of real
   ``noise``,         "Noise for initial tensor",                               Real


Multiple sites can be specified at once by setting a list to ``index``.
An empty list ``[]`` means all sites.

By setting a list to ``virtual_dim``, individual bond dimensions in four directions can be specified.
The order is left (-x), top (+y), right (+x), and bottom (-y).

An initial state of a system :math:`|\Psi\rangle` is represented as
the direct product state of the initial states at each site :math:`i` :math:`|\Psi_i\rangle`:

.. math::
   |\Psi\rangle = \otimes_i |\Psi_i\rangle,

where :math:`|\Psi_i\rangle = \sum_\alpha A_\alpha |\alpha\rangle_i` is the initial state at :math:`i` site.
Site tensors are initialized to realize this product state with some noise.
``initial_state`` specifies values of expansion coefficient :math:`A_\alpha`, which will be automatically normalized.
If an array consisting of only zeros is given, the initial state becomes random state.
The tensor itself is initialized such that all elements with a virtual bond index of 0 are :math:`T_{0000\alpha} = A_\alpha`.
The other elements are independently initialized by a uniform random number of ``[-noise, noise)``.
For example, in the case of :math:`S=1/2` , 
set ``initial_state = [1.0, 0.0]`` when you want to set the initial state as the state :math:`|\Psi_i\rangle = |\uparrow\rangle = |0\rangle`.
When you want to set the initial state as the state :math:`|\Psi_i\rangle = \left(|\uparrow\rangle + |\downarrow\rangle\right)/\sqrt{2}`, set ``initial_state = [1.0, 1.0]``.
