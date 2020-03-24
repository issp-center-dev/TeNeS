.. highlight:: none

In the following, the parameters about the correlation function :math:`C = \langle A(0)B(r) \rangle` is described.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``r_max``,     "Maximum distance :math:`r` of the correlation function",          Integer
   ``operators``, "Numbers of operators A and B that measure correlation functions", A list of integer

The operators defined in the ``observable.onesite`` section are used.

Example
~~~~~~~~

For example, if :math:`S ^ z` is defined as 0 and :math:`S ^ x` is defined as 1, 
then :math:`S ^ z (0) S ^ z (r), S ^ z (0) S ^ x (r), S ^ x (0) S ^ x (r)` 
are measured by the following definition:

::

    [correlation]
    r_max = 5
    operators = [[0,0], [0,1], [1,1]]
