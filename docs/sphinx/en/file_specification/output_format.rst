.. highlight:: none

.. _sec-output-format:

Output files
---------------------------------

Output files are generated in the ``output`` directry.


For all modes
==============


``parameters.dat``
~~~~~~~~~~~~~~~~~~~~~

Paramters in the ``parameter`` and ``lattice`` sections defined in the input file are outputted.

Example::

   simple_num_step = [10]
   simple_tau = [0.01]
   simple_inverse_lambda_cutoff = 1e-12
   simple_gauge_fix = 0
   simple_gauge_maxiter = 100
   simple_gauge_convergence_epsilon = 0.01

   full_num_step = [0]
   full_inverse_projector_cutoff = 1e-12
   full_inverse_precision = 1e-12
   full_convergence_epsilon = 1e-06
   full_iteration_max = 100
   full_gauge_fix = true
   full_fastfullupdate = true

   ctm_dimension = 10
   ctm_inverse_projector_cutoff = 1e-12
   ctm_convergence_epsilon = 1e-06
   ctm_iteration_max = 10
   ctm_projector_corner = true
   use_rsvd = false
   rsvd_oversampling_factor = 2
   meanfield_env = true

   mode = ground state
   simple
   Lcor = 0
   seed = 11
   is_real = 0
   iszero_tol = 0
   measure = 1
   tensor_load_dir = 
   tensor_save_dir = save_tensor
   outdir = output

   Lsub = [ 2 , 2 ]
   skew = 0

   start_datetime =  2023-06-08T16:41:50+09:00

``time.dat``
~~~~~~~~~~~~~~~~~~~~~

The calculation time is outputted.

Example::

   time simple update = 1.64429
   time full update   = 0
   time environmnent  = 0.741858
   time observable    = 0.104487


For ground state calculation mode
====================================

``density.dat``
~~~~~~~~~~~~~~~~

The expectation value per site of each observable is outputted.
When the name of the operator (``name``) is an empty, the index of the operator is written.
``Energy`` means the summation of ``site hamiltonian`` and ``bond hamiltonian``.

Example::

   Energy           = -5.00499902760266346e-01  0.00000000000000000e+00
   site hamiltonian = -4.99999945662006270e-04  0.00000000000000000e+00
   Sz               =  4.99999945662006284e-01  0.00000000000000000e+00
   Sx               =  9.24214061616647275e-05  0.00000000000000000e+00
   Sy               = -2.34065881671767322e-06  0.00000000000000000e+00
   bond hamiltonian = -4.99999902814604325e-01  2.22346094146706503e-21
   SzSz             =  4.99999902814604380e-01 -1.80051315353166456e-21
   SxSx             =  1.12631053560300631e-05  6.08792260271591701e-21
   SySy             = -1.12817627661272438e-05  4.76468712680822333e-21

``onesite_obs.dat``
~~~~~~~~~~~~~~~~~~~~~~

-  The expected values of the site operator :math:`\langle\hat{A}^\alpha_i\rangle = \langle\Psi | \hat{A}^\alpha_i | \Psi \rangle / \langle\Psi | \Psi \rangle` are outputted.
-  Each row consists of four columns.

   1. Index of the operator :math:`\alpha`
   2. Index of the sites :math:`i`
   3. Real part of the expected value :math:`\mathrm{Re}\langle\hat{A}^\alpha_i\rangle`
   4. Imag part of the expected value :math:`\mathrm{Im}\langle\hat{A}^\alpha_i\rangle`

- In addition, norm of the wave function :math:`\langle \Psi | \Psi \rangle` is outputted as an operator with index of -1.

   - If the imaginary part is finite, something is wrong. A typical cause is that the bond dimension of the CTM is too small.

Example::

   # The meaning of each column is the following: 
   # $1: op_group
   # $2: site_index
   # $3: real
   # $4: imag
   # The names of op_group are the following: 
   # 0: site hamiltonian
   # 1: Sz              
   # 2: Sx              
   # 3: Sy              
   # -1: norm

   0 0 -4.99999945520001373e-04 0.00000000000000000e+00
   0 1 -4.99999967900088089e-04 0.00000000000000000e+00
   0 2 -4.99999894622883147e-04 0.00000000000000000e+00
   0 3 -4.99999974605052581e-04 0.00000000000000000e+00
   1 0 4.99999945520001376e-01 0.00000000000000000e+00
   1 1 4.99999967900088049e-01 0.00000000000000000e+00
   1 2 4.99999894622883134e-01 0.00000000000000000e+00
   1 3 4.99999974605052522e-01 0.00000000000000000e+00
      ... Skipped ...
   -1 3 1.00000000000000044e+00 0.00000000000000000e+00

``twosite_obs.dat``
~~~~~~~~~~~~~~~~~~~~~~

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

Example::

   # The meaning of each column is the following: 
   # $1: op_group
   # $2: source_site
   # $3: dx
   # $4: dy
   # $5: real
   # $6: imag
   # The names of op_group are the following: 
   # 0: bond hamiltonian
   # 1: SzSz            
   # 2: SxSx            
   # 3: SySy            
   # -1: norm

   0 0 0 1 -2.49999925774909121e-01 3.38316768671362694e-21
   0 0 1 0 -2.49999967989907063e-01 4.24343236807659553e-22
   0 1 0 1 -2.49999972903562101e-01 -2.06825262200104597e-25
   0 1 1 0 -2.49999957625646446e-01 2.06789370628128221e-24
   0 2 0 1 -2.49999931343147630e-01 3.11801499860976615e-28
   0 2 1 0 -2.49999939447834718e-01 1.65429596395607220e-24
      ... Skipped ...
   -1 3 1 0 1.00000000000000067e+00 0.00000000000000000e+00

``correlation.dat``
~~~~~~~~~~~~~~~~~~~~~

-  Correlation functions :math:`C^{\alpha \beta}_i(x,y) \equiv \langle \hat{A}^\alpha(x_i,y_i) \hat{A}^\beta(x_i+x,y_i+y) \rangle` are outputted.
-  Each row consists of seven columns.

   1. Index of the left operator :math:`\alpha`
   2. Index of the left site :math:`i`
   3. Index of the right operator :math:`\beta`
   4. x coordinate of the right site :math:`x`
   5. y coordinate of the right site :math:`y`
   6. Real part :math:`\mathrm{Re}C`
   7. Imaginary part :math:`\mathrm{Im}C`

Example::

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
~~~~~~~~~~~~~~~~~~~~~~~~~~~

The correlation length :math:`\xi` is outputted. Each row consists of 3+n columns.

1. Direction (``0: x, 1: y``)
2. When direction is ``0`` it is :math:`y` coodinate, and otherwise :math:`x` coordinate
3. Correlation length :math:`\xi = 1/e_1`

The 4th and the subsequent columns show the logarithm of the absolute value of the eigenvalues of the transfer matrix, :math:`e_i = -\log\left|\lambda_i/\lambda_0\right|` (:math:`i>0`).
This information may be used to estimate the bond dimension dependence of the correlation length.
See PRX **8**, 041033 (2018) and PRX **8**, 031030 (2018) for more information.

Example::

   # The meaning of each column is the following: 
   # $1: direction 0: +x, 1: +y
   # $2: y (dir=0) or x (dir=1) coorinates
   # $3: correlation length xi = 1/e_1 
   # $4-: eigenvalues e_i = -log|t_i/t_0|
   #      where i > 0 and t_i is i-th largest eigenvalue of T

   0 0 2.18785686529154477e-01 4.57068291744370647e+00 4.57068291744370647e+00 4.88102462824739991e+00
   0 1 2.20658864940629751e-01 4.53188228022952533e+00 4.53188228022952533e+00 4.56359469233104953e+00
   1 0 2.23312072254469030e-01 4.47803824443704013e+00 4.47803824443704013e+00 6.03413555039678595e+00
   1 1 2.00830966658579996e-01 4.97931178960083720e+00 4.97931178960083720e+00 5.08813099309339911e+00


For time evolution mode
=========================

``TE_density.dat``
~~~~~~~~~~~~~~~~~~~

The expectation value per site of each obesrvable is outputted.
Each row consists of four columns.

1. Time :math:`t`
2. Operator ID :math:`\alpha`
3. Real part of the expected value :math:`\mathrm{Re}\langle\hat{A}^\alpha_i\rangle`
4. Imag part of the expected value :math:`\mathrm{Im}\langle\hat{A}^\alpha_i\rangle`

Example::

   # The meaning of each column is the following: 
   # $1: time
   # $2: observable ID
   # $3: real
   # $4: imag
   # The meaning of observable IDs are the following: 
   # 0: Energy
   # 1: site hamiltonian
   # 2: Sz              
   # 3: Sx              
   # 4: Sy              
   # 5: bond hamiltonian
   # 6: SzSz            
   # 7: SxSx            
   # 8: SySy            

   0.00000000000000000e+00 0 -5.00684745572451129e-01  0.00000000000000000e+00
   0.00000000000000000e+00 1 -6.84842757985213292e-04  0.00000000000000000e+00
   0.00000000000000000e+00 2  4.99999945661913914e-01  0.00000000000000000e+00
   0.00000000000000000e+00 3  9.24214061616496842e-05  0.00000000000000000e+00
      ... Skipped ...
   4.99999999999993783e+00 8  2.54571641402435656e-01  3.25677610112348483e-17


``TE_onesite_obs.dat``
~~~~~~~~~~~~~~~~~~~~~~~~

The expected values of the site operators :math:`\langle\hat{A}^\alpha_i\rangle = \langle\Psi | \hat{A}^\alpha_i | \Psi \rangle / \langle\Psi | \Psi \rangle` are outputted.
Each row consists of five columns.

1. Time :math:`t`
2. Index of the operator :math:`\alpha`
3. Index of the sites :math:`i`
4. Real part of the expected value :math:`\mathrm{Re}\langle\hat{A}^\alpha_i\rangle`
5. Imag part of the expected value :math:`\mathrm{Im}\langle\hat{A}^\alpha_i\rangle`

- In addition, norm of the wave function :math:`\langle \Psi | \Psi \rangle` is outputted as an operator with index of -1.

   - If the imaginary part is finite, something is wrong. A typical cause is that the bond dimension of the CTM is too small.


Example::

   # The meaning of each column is the following: 
   # $1: time
   # $2: op_group
   # $3: site_index
   # $4: real
   # $5: imag
   # The names of op_group are the following: 
   # 0: site hamiltonian
   # 1: Sz              
   # 2: Sx              
   # 3: Sy              
   # -1: norm

   0.00000000000000000e+00 0 0 -6.43318936197596913e-04 0.00000000000000000e+00
   0.00000000000000000e+00 0 1 -6.73418200262321655e-04 0.00000000000000000e+00
   0.00000000000000000e+00 0 2 -9.89240026254938282e-04 0.00000000000000000e+00
   0.00000000000000000e+00 0 3 -4.33393869225996210e-04 0.00000000000000000e+00
   0.00000000000000000e+00 1 0 4.99999945519898625e-01 0.00000000000000000e+00
   0.00000000000000000e+00 1 1 4.99999967900020936e-01 0.00000000000000000e+00
   0.00000000000000000e+00 1 2 4.99999894622765451e-01 0.00000000000000000e+00
      ... Skipped ...
   4.99999999999993783e+00 -1 3 9.99999999999999667e-01 0.00000000000000000e+00

``TE_twosite_obs.dat``
~~~~~~~~~~~~~~~~~~~~~~~~

-  Expectation values for two-site operations are outputted.
-  Each row consists of six columns.

   1. Time :math:`t`
   2. Index of the two-site operator
   3. Index of the source site
   4. x coordinate of the target site from the source site
   5. y coordinate of the target site from the source site
   6. Real part of the expected value
   7. Imaginary part of the expected value

- In addition, norm of the wave function :math:`\langle \Psi | \Psi \rangle` is outputted as an operator with index of -1.

   - If the imaginary part is finite, something is wrong. A typical cause is that the bond dimension of the CTM is too small.


Example::

   # The meaning of each column is the following: 
   # $1: time
   # $2: op_group
   # $3: source_site
   # $4: dx
   # $5: dy
   # $6: real
   # $7: imag
   # The names of op_group are the following: 
   # 0: bond hamiltonian
   # 1: SzSz            
   # 2: SxSx            
   # 3: SySy            
   # -1: norm

   0.00000000000000000e+00 0 0 0 1 -2.49999925774803150e-01 -1.01660465821037727e-20
   0.00000000000000000e+00 0 0 1 0 -2.49999967989888300e-01 4.23516895582898471e-22
   0.00000000000000000e+00 0 1 0 1 -2.49999972903488521e-01 -6.20403358955599675e-25
   0.00000000000000000e+00 0 1 1 0 -2.49999957625561042e-01 4.13590865617858526e-25
   0.00000000000000000e+00 0 2 0 1 -2.49999931343070220e-01 8.27316466562544801e-25
      ... Skipped ...
   4.99999999999993783e+00 -1 3 1 0 9.99999999999999445e-01 1.38777878078144568e-17

``TE_correlation.dat``
~~~~~~~~~~~~~~~~~~~~~~~

-  Correlation functions :math:`C^{\alpha \beta}_i(x,y) \equiv \langle \hat{A}^\alpha(x_i,y_i) \hat{A}^\beta(x_i+x,y_i+y) \rangle` are outputted.
-  Each row consists of eight columns.

   1. Time :math:`t`
   2. Index of the left operator :math:`\alpha`
   3. Index of the left site :math:`i`
   4. Index of the right operator :math:`\beta`
   5. x coordinate of the right site :math:`x`
   6. y coordinate of the right site :math:`y`
   7. Real part :math:`\mathrm{Re}C`
   8. Imaginary part :math:`\mathrm{Im}C`

Example::

   # The meaning of each column is the following: 
   # $1: time
   # $2: left_op
   # $3: left_site
   # $4: right_op
   # $5: right_dx
   # $6: right_dy
   # $7: real
   # $8: imag
   # The names of operators are the following: 
   # 0: site hamiltonian
   # 1: Sz              
   # 2: Sx              
   # 3: Sy              

   0.00000000000000000e+00 0 0 0 1 0 1.83422488349707711e-04 1.90382762094233524e-20 
   0.00000000000000000e+00 0 0 0 2 0 8.30943360551218668e-07 -4.19695835411528090e-23 
   0.00000000000000000e+00 0 0 0 3 0 4.12158436385765748e-07 -1.04903226091485958e-23 
   0.00000000000000000e+00 0 0 0 4 0 4.13819451426396547e-07 1.74438421668770658e-23 
   0.00000000000000000e+00 0 0 0 5 0 4.33224506806043380e-07 -8.71850465073480394e-24 
      ... Skipped ...
   4.99999999999993783e+00 2 3 2 0 5 3.96301355731331212e-02 -1.37659660157453792e-18 


``TE_correlation_length.dat``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The correlation length :math:`\xi` is outputted. Each row consists of 4+n columns.

1. Time :math:`t`
2. Direction (``0: x, 1: y``)
3. When direction is ``0`` it is :math:`y` coodinate, and otherwise :math:`x` coordinate
4. Correlation length :math:`\xi = 1/e_1`

The 5th and the subsequent columns show the logarithm of the absolute value of the eigenvalues of the transfer matrix, :math:`e_i = -\log\left|\lambda_i/\lambda_0\right|` (:math:`i>0`).
This information may be used to estimate the bond dimension dependence of the correlation length.
See PRX **8**, 041033 (2018) and PRX **8**, 031030 (2018) for more information.

Example::

   # The meaning of each column is the following: 
   # $1: time
   # $2: direction 0: +x, 1: +y
   # $3: y (dir=0) or x (dir=1) coorinates
   # $4: correlation length xi = 1/e_1 
   # $5-: eigenvalues e_i = -log|t_i/t_0|
   #      where i > 0 and t_i is i-th largest eigenvalue of T

   0.00000000000000000e+00 0 0 2.18785686529220424e-01 4.57068291744232891e+00 4.57068291744232891e+00 4.88102462824919758e+00
   0.00000000000000000e+00 0 1 2.20658864940612931e-01 4.53188228022987083e+00 4.53188228022987083e+00 4.56359469232955917e+00
   0.00000000000000000e+00 1 0 2.23312072254560540e-01 4.47803824443520515e+00 4.47803824443520515e+00 6.03413555040836602e+00
   0.00000000000000000e+00 1 1 2.00830966658709920e-01 4.97931178959761578e+00 4.97931178959761667e+00 5.08813099310449513e+00
   9.99999999999999917e-02 0 0 2.02379048126702904e-01 4.94122296382149528e+00 4.94122296382149617e+00 6.74309974506451315e+00
   9.99999999999999917e-02 0 1 2.20416567580991346e-01 4.53686404327366777e+00 4.53686404327366777e+00 6.18101616573088020e+00
   9.99999999999999917e-02 1 0 2.12137154053103655e-01 4.71393143960851368e+00 4.71393143960851368e+00 7.17220113786375002e+00
   9.99999999999999917e-02 1 1 1.90367314703518503e-01 5.25300260476656966e+00 5.25300260476656966e+00 7.61893825410630487e+00
   2.00000000000000039e-01 0 0 1.96835348300227503e-01 5.08038829730281805e+00 5.08038829730281805e+00 7.35176717846311778e+00
   2.00000000000000039e-01 0 1 2.02355022722768896e-01 4.94180963014702801e+00 4.94180963014702801e+00 6.57691315725687975e+00
   2.00000000000000039e-01 1 0 2.05314677188187883e-01 4.87057239986509760e+00 4.87057239986509760e+00 7.90951918842309798e+00
   2.00000000000000039e-01 1 1 1.63323696507474692e-01 6.12281023136305169e+00 6.12281023136305169e+00 7.83104916294462416e+00
      ... Skipped ...
   4.99999999999993783e+00 1 1 4.61585992965019176e-01 2.16644355600232430e+00 2.16644355600232430e+00 2.29497956495965427e+00
