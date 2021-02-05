.. highlight:: none

In this section, the parameters about the site-site correlation function :math:`C = \left\langle A(\boldsymbol{r}_0)B(\boldsymbol{r}_0 + \boldsymbol{r}) \right\rangle` is specified.
If you omit this section, no correlation functions will be calculated.

Coordinates :math:`\boldsymbol{r}, \boldsymbol{r}_0` measured in the system of square lattice TNS.
For example, the coordinate of the right neighbor tensor is :math:`\boldsymbol{r} = (1,0)` and that of the top neighbor one is :math:`\boldsymbol{r} = (0,1)`.
TeNeS calculates the correlation functions along the positive direction of :math:`x` and :math:`y` axis, that is,

.. math::
   \boldsymbol{r} = (0,0), (1,0), (2,0), \dots, (r_\text{max}, 0), (0,1), (0,2), \dots, (0, r_\text{max})

The coordinate of each site of the unitcell is used as the center coordinate, :math:`\boldsymbol{r}_0`.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``r_max``,     "Maximum distance :math:`r` of the correlation function", Integer
   ``operators``, "Indices of operators A and B to be measured",            A list of integer

The operators defined in the ``observable.onesite`` section are used.

Example
~~~~~~~~

For example, if :math:`S^z` is defined as 0th operator and :math:`S^x` is defined as 1st one,
then :math:`S^z(0) S^z(r), S^z(0) S^x(r), S^x(0) S^x(r)` for :math:`0 \le r \le 5`
are measured by the following definition:

::

    [correlation]
    r_max = 5
    operators = [[0,0], [0,1], [1,1]]
