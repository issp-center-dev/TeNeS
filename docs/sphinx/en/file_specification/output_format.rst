.. highlight:: none

.. _sec-output-format:

Output files
---------------------------------

Output files are generated in the ``output`` directry.


``parameters.dat``
=====================

Paramters in the ``parameter`` and ``lattice`` sections defined in the input file are outputted.

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
-  Each row consists of five columns.

   1. Index of the operator
   2. Index of the sites
   3. Index of the sites
   4. Real part of the expected value
   5. Imaginary part of the expected value

- In addition, norm of the wave function :math:`\langle \Psi | \Psi \rangle` is outputted as an operator with index of -1.

   - If the imaginary part is finite, something is wrong. A typical cause is that the bond dimension of the CTM is too small.

::

    # $1: op_index
    # $2: source_site
    # $3: target_site
    # $4: real
    # $5: imag

   0 0 4.08513219339284250e-01 -1.30821773764633317e-18
   0 1 -4.08513219339288636e-01 -6.09055344460395354e-20
   0 2 -4.08513219339289801e-01 -6.46773648122808868e-19
   0 3 4.08513219339286304e-01 -1.46983151052043956e-19
   1 0 -5.63321380968339106e-05 -1.06329584429579662e-18
   1 1 5.63321015724943072e-05 8.78344073011021613e-18
   1 2 5.63321015892866812e-05 4.36590453150828335e-18
   1 3 -5.63321381136249429e-05 -1.14509908259615756e-18
   2 0 7.95013947121512355e-06 -5.82669256660366468e-18
   2 1 -7.95019526521416614e-06 9.61281832688426906e-18
   2 2 -7.95019527785552810e-06 -9.34047696321987429e-18
   2 3 7.95013948581251836e-06 -6.49339357741872464e-19
   -1 0 1.04968851690174758e+00 7.97972798949331263e-17
   -1 1 1.04968851690174780e+00 7.65446733774766130e-17
   -1 2 1.04968851690174647e+00 8.32667268468867405e-17
   -1 3 1.04968851690174803e+00 7.41594285980085033e-17

``correlation.dat``
=====================

-  Correlation functions are outputted.
-  Each row consists of eight columns.

   1. Index of the left operator
   2. Site index of the left operator
   3. Index of the right operator
   4. Site index of the right operator
   5. Unit cell offset of the right operator (x)
   6. Unit cell offset of the right operator (y)
   7. Real part of the expected value
   8. Imaginary part of the expected value

- In addition, norm of the wave function :math:`\langle \Psi | \Psi \rangle` is outputted as an operator with index of -1.

   - If the imaginary part is finite, something is wrong. A typical cause is that the bond dimension of the CTM is too small.

Example
~~~~~~~

::

   # $1: op_group
   # $2: source_site
   # $3: dx
   # $4: dy
   # $5: real
   # $6: imag

   0 0 0 1 -3.30408104727482554e-01 -3.63538996091175880e-19
   0 0 1 0 -3.26902334621655521e-01 -1.28557778331473411e-19
   0 1 0 1 -3.30408104727482110e-01 6.13195629489298286e-18
   0 1 1 0 -3.28820570518176758e-01 5.98724951760379135e-18
   0 2 0 1 -3.32375821733345012e-01 -5.42272048973129865e-18
   0 2 1 0 -3.26902334621652579e-01 9.69166076872613868e-20
   0 3 0 1 -3.32375821733344956e-01 5.07748884268378299e-18
   0 3 1 0 -3.28820570518176702e-01 4.86902738935337153e-18
   1 0 0 1 -1.87348767102901825e-01 4.90760305979372382e-19
      ... skipped ...
   -1 3 1 0 1.07465536687797147e+00 7.74120351154650166e-17


``correlation_length.dat``
===========================


1. Direction (``0: x, 1: y``)
2. When direction is ``0`` it is :math:`y` coodinate, and otherwise :math:`x` coordinate
3. Correlation length :math:`\xi = 1/e_1`

The 4th and the subsequent columns show the logarithm of the absolute value of the eigenvalues of the transfer matrix, :math:`e_i = -\log\left|\lambda_i/\lambda_0\right|` (:math:`i>0`).
This information may be used to estimate the bond dimension dependence of the correlation length.
See PRX **8**, 041033 (2018) and PRX **8**, 031030 (2018) for more information.

Example
~~~~~~~~~~~~

::

   # $1: direction
   # $2: col or row index
   # $3: correlation length
   # $4-: eigenvalues e_i = -log|t_i/t_0|
   #      where i > 0 and t_i is i-th largest eigenvalue of T

   0 0 7.19213553469021272e-01 1.39040761283856007e+00 1.44013584036962405e+00 1.53522220522654251e+00
   0 1 7.19303527237354912e-01 1.39023369430805133e+00 1.39042786247674610e+00 1.53457094348925005e+00
   1 0 7.26232546918431754e-01 1.37696940772377285e+00 1.39968879441491767e+00 1.51923157420858113e+00
   1 1 7.26095712518373015e-01 1.37722890076244076e+00 1.38699264750702023e+00 1.52016493301531241e+00

``time.dat``
=====================

The calculation time is outputted.
