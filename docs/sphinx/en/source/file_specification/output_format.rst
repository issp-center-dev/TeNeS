.. highlight:: none

Output files
---------------------------------

Output files are generated in the ``output`` directry.


``parameters.dat``
=====================

Paramters in ``parameter`` and ``lattice`` sections are outputted.

``energy.dat``
==============

The energy of each site is output.

``site_obs.dat``
=================


-  The expected values of the site operator are outputted.
-  Each row consists of four columns.

   1. Index of the operator
   2. Index of the sites
   3. Real part of the expected value
   4. Imaginary part of the expected value

Example
~~~~~~~

::

    # $1: op_index
    # $2: site_index
    # $3: real
    # $4: imag

    0 0 1.92549465249573365e-02 0.00000000000000000e+00
    0 1 -1.92620814130195529e-02 0.00000000000000000e+00
    0 2 -1.95243093055922252e-02 0.00000000000000000e+00
    0 3 1.91619477632061150e-02 0.00000000000000000e+00
    1 0 4.07206063348768799e-01 0.00000000000000000e+00
    1 1 -4.07243511737157671e-01 0.00000000000000000e+00
    1 2 -4.07255967738734126e-01 0.00000000000000000e+00
    1 3 4.07308918791554009e-01 0.00000000000000000e+00

``neighbor_obs.dat``
======================

-  Nearest neighbor correlations for site operations are outputted.
-  Each row consists of four columns.

   1. Index of the operator
   2. Index of the sites
   3. Index of the sites
   4. Real part of the expected value
   5. Imaginary part of the expected value

::

    # $1: op_index
    # $2: source_site
    # $3: target_site
    # $4: real
    # $5: imag

    0 0 1 -7.05927615064968900e-02 0.00000000000000000e+00
    0 0 2 -7.27068456430051274e-02 0.00000000000000000e+00
    0 1 0 -7.13284385957392297e-02 0.00000000000000000e+00
    0 1 3 -7.19523349256113581e-02 0.00000000000000000e+00
    0 2 3 -7.12610364895483045e-02 0.00000000000000000e+00
    0 2 0 -7.19731507561011952e-02 0.00000000000000000e+00
    0 3 2 -7.05633558230210067e-02 0.00000000000000000e+00
    0 3 1 -7.26803750807340498e-02 0.00000000000000000e+00
    1 0 1 -1.85942869237103348e-01 0.00000000000000000e+00
    1 0 2 -1.87164731677545187e-01 0.00000000000000000e+00
    1 1 0 -1.86360382550076586e-01 0.00000000000000000e+00
    1 1 3 -1.86768451086366694e-01 0.00000000000000000e+00
    1 2 3 -1.86384181909805935e-01 0.00000000000000000e+00
    1 2 0 -1.86747576732693515e-01 0.00000000000000000e+00
    1 3 2 -1.85975089525013598e-01 0.00000000000000000e+00
    1 3 1 -1.87196522916879049e-01 0.00000000000000000e+00

``correlation.dat``
=====================

-  Correlation functions are outputted.
-  Each row consists of eight columns.

   1. Index of the left-hand operator
   2. Site index of the left operator
   3. Index of the right operator
   4. Site index of the right operator
   5. Unit cell offset of the right operator (x)
   6. Unit cell offset of the right operator (y)
   7. 

::

    # $1: op_index
    # $2: source_site
    # $3: target_site
    # $4: real
    # $5: imag

    0 0 1 -7.05927615064968900e-02 0.00000000000000000e+00
    0 0 2 -7.27068456430051274e-02 0.00000000000000000e+00
    0 1 0 -7.13284385957392297e-02 0.00000000000000000e+00
    0 1 3 -7.19523349256113581e-02 0.00000000000000000e+00
    0 2 3 -7.12610364895483045e-02 0.00000000000000000e+00
    0 2 0 -7.19731507561011952e-02 0.00000000000000000e+00
    0 3 2 -7.05633558230210067e-02 0.00000000000000000e+00
    0 3 1 -7.26803750807340498e-02 0.00000000000000000e+00
    1 0 1 -1.85942869237103348e-01 0.00000000000000000e+00
    1 0 2 -1.87164731677545187e-01 0.00000000000000000e+00
    1 1 0 -1.86360382550076586e-01 0.00000000000000000e+00
    1 1 3 -1.86768451086366694e-01 0.00000000000000000e+00
    1 2 3 -1.86384181909805935e-01 0.00000000000000000e+00
    1 2 0 -1.86747576732693515e-01 0.00000000000000000e+00
    1 3 2 -1.85975089525013598e-01 0.00000000000000000e+00
    1 3 1 -1.87196522916879049e-01 0.00000000000000000e+00

``correlation.dat``
=====================

-  Correlation functions are outputted.
-  Each row consists of eight columns.

   1. Index of the left-hand operator
   2. Site index of the left operator
   3. Index of the right operator
   4. Site index of the right operator
   5. Unit cell offset of the right operator (x)
   6. Unit cell offset of the right operator (y)
   7. Real part of the expected value
   8. Imaginary part of the expected value

Example
~~~~~~~

::

    # $1: left_op
    # $2: left_site
    # $3: right_op
    # $4: right_site
    # $5: offset_x
    # $6: offset_y
    # $7: real
    # $8: imag

    0 0 0 1 0 0 -7.05927615064967928e-02 0.00000000000000000e+00 
    0 0 0 0 1 0 1.19668843226761017e-02 0.00000000000000000e+00 
    0 0 0 1 1 0 -2.43086229320005863e-03 0.00000000000000000e+00 
    0 0 0 0 2 0 7.42729194528496308e-04 0.00000000000000000e+00 
    0 0 0 1 2 0 -4.38794819416885419e-04 0.00000000000000000e+00 
    0 0 0 2 0 0 -7.27068456430051135e-02 0.00000000000000000e+00 
    0 0 0 0 0 1 1.23339845746621279e-02 0.00000000000000000e+00 
    0 0 0 2 0 1 -2.50111186244407349e-03 0.00000000000000000e+00 
    0 0 0 0 0 2 7.54607806587391516e-04 0.00000000000000000e+00 
    0 0 0 2 0 2 -4.47734559969679546e-04 0.00000000000000000e+00 
    1 0 1 1 0 0 -1.85942869237103237e-01 0.00000000000000000e+00 
    ...
    1 3 1 1 0 3 -1.65874245891461547e-01 0.00000000000000000e+00

``time.dat``
=====================

The calculation time is outputted.
