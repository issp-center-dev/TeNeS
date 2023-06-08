.. highlight:: none

.. _sec-output-format:

Output files
---------------------------------

Output files are generated in the ``output`` directry.


``parameters.dat``
=====================

Paramters in the ``parameter`` and ``lattice`` sections defined in the input file are outputted.

``density.dat``
================

The expectation value per site of each observable is outputted.
When the name of the operator (``name``) is an empty, the index of the operator is written.

Example
~~~~~~~~~~

::

   Sz          =  6.11647102532908438e-03  0.00000000000000000e+00
   Sx          = -1.18125085038094907e-01  0.00000000000000000e+00
   hamiltonian = -5.43684776153081639e-01  0.00000000000000000e+00
   SzSz        = -3.16323622995942133e-01  0.00000000000000000e+00
   SxSx        = -8.55704529153783616e-02  0.00000000000000000e+00
   SySy        = -1.41790700241760936e-01  0.00000000000000000e+00

``onesite_obs.dat``
======================


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

``twosite_obs.dat``
======================

-  Expectation values for two-site operations are outputted.
-  Each row consists of six columns.

   1. Index of the two-site operator
   2. Index of the source site
   3. x coordinate of the target site from the source site
   4. y coordinate of the target site from the source site
   5. Real part of the expected value
   6. Imaginary part of the expected value

- In addition, norm of the wave function :math:`\langle \Psi | \Psi \rangle` is outputted as an operator with index of -1.

   - If the imaginary part is finite, something is wrong. A typical cause is that the bond dimension of the CTM is too small.

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

``correlation.dat``
=====================

-  Correlation functions :math:`C^{\alpha \beta}_i(x,y) \equiv \langle \hat{A}^\alpha(x_i,y_i) \hat{A}^\beta(x_i+x,y_i+y) \rangle` are outputted.
-  Each row consists of eight columns.

   1. Site index of the left operator
   2. Index of the left operator
   3. Index of the right operator
   4. Site index of the right operator
   5. Unit cell offset of the right operator (x)
   6. Unit cell offset of the right operator (y)
   7. Real part of the expected value
   8. Imaginary part of the expected value


   1. Index of the left operator :math:`\alpha`
   2. Index of the left site :math:`i`
   3. Index of the right operator :math:`\beta`
   4. x coordinate of the right site :math:`x`
   5. y coordinate of the right site :math:`y`
   6. Real part :math:`\mathrm{Re}C`
   7. Imaginary part :math:`\mathrm{Im}C`

Example
~~~~~~~~~~~~

::

   # $1: left_op
   # $2: left_site
   # $3: right_op
   # $4: right_dx
   # $5: right_dy
   # $6: real
   # $7: imag

   0 0 0 1 0 -1.71759992763061836e-01 1.36428299157186382e-14 
   0 0 0 2 0 1.43751794649139675e-01 -1.14110668277268192e-14 
   0 0 0 3 0 -1.42375391377041444e-01 1.14103263451826963e-14 
   0 0 0 4 0 1.41835919840103741e-01 -1.11365361507372103e-14 
   0 0 0 5 0 -1.41783912096811515e-01 1.12856813523671142e-14 
   0 0 0 0 1 -1.72711348845767942e-01 1.40873628493918905e-14 
   0 0 0 0 2 1.43814797743900907e-01 -1.17958665742991377e-14 
   0 0 0 0 3 -1.42415176172922653e-01 1.22109610917000360e-14 
   0 0 0 0 4 1.41838862178711583e-01 -1.19321507524565005e-14 
   0 0 0 0 5 -1.41792935491960648e-01 1.23094733264734764e-14 
   1 0 1 1 0 -7.95389427681298805e-02 6.15901595234210079e-15 
   1 0 1 2 0 2.01916094009441903e-02 -1.27162373457160362e-15 
   ... Skipped ...
   2 3 2 0 5 -1.41888376278899312e-03 -2.38672137694415560e-16 

``correlation_length.dat``
===========================

The correlation length :math:`\xi` is outputted.

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

Example
~~~~~~~~~~~

::

   time simple update = 1.64429
   time full update   = 0
   time environmnent  = 0.741858
   time observable    = 0.104487
