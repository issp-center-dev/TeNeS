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

Let the whole Hamiltonian be the sum of the site Hamiltonian (one-site Hamiltonian) and bond Hamiltonian (two-site Hamiltonian). 

.. math::
  \mathcal{H} = \sum_i \mathcal{H}_i + \sum_{i,j} \mathcal{H}_{ij}

In ``hamiltonian`` section, each local Hamiltonian is defined.
The format is similar to that of the one-site and two-site operator specified in ``observable.onesite`` and ``observable.twosite``.

.. csv-table::
   :header: "Name", "Description", "Type"
   :widths: 15, 30, 20

   ``dim``,      "Dimension of an operator",         A list of integers
   ``sites``,    "Site",                             A list of integers
   ``bonds``,    "Bond",                             String
   ``elements``, "Non-zero elements of an operator", String

``dim`` specifies a dimension of an operator.
In other words, the number of possible states of the site where the operator acts on.
In the case of interaction between two :math:`S=1/2` spin, for example, ``dim = [2,2]`` .
``tenes_std`` judges whether a local Hamiltonian is site one or bond one from the number of integers in ``dims``; if one, a site Hamiltonian is defined and otherwise a bond one.

``sites``, a list of integers, specifies a set of sites where the site operator acts.
An empty list (``[]``) means all the sites.

``bonds`` specifies a string representing the set of site pairs on which the operator acts.
One line consisting of three integers means one site pair.

- The first integer is the number of the source site.
- The last two integers are the coordinates (dx, dy) of the destination site (target) from the source site.

``elements`` is a string specifying the non-zero element of an operator.
One element consists of one line consisting of two (site) or four (bond) integers and two floating-point numbers separated by spaces.

- For site Hamiltonian

    * The first integer is the index of the state of the site **before** the operator acts on.
    * The next one shows the index of the state of the site **after** the operator acts on.
    * The last two indicate the real and imaginary parts of the elements of the operator.

- For bond Hamiltonian

    * The first two integers are the indices of the states of the source site and target site **before** the operator acts on.
    * The next two show the indices of the states of the source site and target site **after** the operator acts on.
    * The last two indicate the real and imaginary parts of the elements of the operator.


``correlation`` section
===========================

.. include:: ./correlation_section.rst

``correlation_length`` section
==================================

.. include:: ./correlation_length_section.rst
