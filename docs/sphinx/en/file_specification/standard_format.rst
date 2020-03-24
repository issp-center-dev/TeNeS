.. highlight:: none

.. _sec-std-format:

Input file for ``tenes_std``
---------------------------------

-  File format:
   `TOML <https://qiita.com/minoritea/items/c0de47b8beb813c655d4>`__ format
-  This file has 5 sections: ``parameter``, ``tensor``, ``hamiltonian``, ``observable``, ``correlation``

   - The four sections other than ``hamiltonian`` are identical to the ``tenes`` input file format, with the following exceptions, and are copied to the ``tenes`` input file.
  
   - By setting a real number for ``parameter.simple_update.tau`` and ``parameter.full_update.tau``, the imaginary time step for the imaginary time evolution operator can be specified.

``parameter`` section
===========================

.. include:: ./parameter_section.rst


``tensor`` section
===========================

.. include:: ./tensor_section.rst


``observable`` section
=============================

.. include:: ./observable_section.rst


``hamiltonian`` section
==============================

Let the whole Hamiltonian be the sum of the Bond Hamiltonian (two-site Hamiltonian). 
.. math::
  \mathcal{H} = \sum_{i,j} \mathcal{H}_{ij}

In ``hamiltonian`` section, each two-site Hamiltonian is defined.
The format is similar to that of the two-site operator made in ``observable.twosite``.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``bonds``,    "Bond",             String
   ``dim``,      "Dimension of an operator",       A list of integer
   ``elements``, "Non-zero elements of an operator", String

``bonds`` specifies a string representing the set of site pairs on which the operator acts.
One line consisting of three integers means one site pair.
The first integer is the number of the source site.
The last two integers are the coordinates (dx, dy) of the destination site (target) from the source site.

``dim`` specifies a dimension of an operator. In other words, the number of possible states of the site where the operator acts on.

``elements`` is a string specifying the non-zero element of an operator.
One element consists of one line consisting of four integers and two floating-point numbers separated by spaces.
The first two are the status numbers of the source site and target site before the operator acts on.
The next two show the status numbers of the source site and target site after the operator acts on.
The last two indicate the real and imaginary parts of the elements of the operator.


``correlation`` section
===========================

.. include:: ./correlation_section.rst
