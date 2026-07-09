FAQ
-------------------

Q1. How can I set the number of full updates?

A1. The full update improves the accuracy of the calculation, but its computational time becomes longer. The number of full updates should be determined by balancing computer performance, bond dimensions, and required computational accuracy. One way to see what happens is to just do a simple update first. Also, if you do a simple update before the full update and approach the ground state, the full update can be done efficiently. However, if you have a complex quantum state (such as a spin liquid) and the simple update does not approach the correct ground state, you must do a full update from the beginning.

Q2. How can I set the number of the simple updates?

A2. Increasing the number of simple updates over time should bring the simulated quantum state closer to the ground state, but this is not necessarily the case; if the bond dimension is small, the calculation accuracy may deteriorate during the update. To see if the calculation works, plot the ground state energy against the number of simple updates. It is a good idea to increase the number of updates and adopt the minimum energy as the calculation result. Another strategy is to increase the number of the simple update so that the energy is almost the same, but that is not necessarily the minimum energy.

Q3. How do I get the bond dimension?

A3. Increasing the bond dimension improves calculation accuracy but increases calculation time. It is necessary to determine the bond dimension by considering the balance of the computer resources and the accuracy required for the physical quantity desired. It is also important to change parameter.ctm.dimension together when changing the bond dimension assigned by lattice.virtual_dim. Typically, the latter takes a value greater than or equal to the square of the former.

Q4. How should I test the correctness of the calculated ground state?

A4. It is difficult to guarantee that the calculated ground state is correct, but the easiest way is to check whether lower energy states are obtained or not by using multiple seed numbers. It is also useful to calculate the initial configuration guessed from several candidate ground states and compare their energies, though its computational cost is expensive. It is also important to change the shape of the unit cell to check for other low-energy states. Although there is a size issue, it is recommended to compare with other methods such as exact diagonalization. (The exact diagonalization can be easily performed using H :math:`\Phi` .)

