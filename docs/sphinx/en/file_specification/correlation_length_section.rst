.. highlight:: none

This section describes how to calculate the correlation length :math:`\xi`.

.. csv-table::
   :header: "Name", "Description", "Type", "Default"
   :widths: 15, 30, 20, 20

   ``measure``,                  "Whether to calculate :math:`xi` or not",  Bool,    true
   ``num_eigvals``,              "The number of eigenvalues of the transfer matrix to be calculated", Integer, 4
   ``maxdim_dense_eigensolver``, "Maximum dimension of the transfer matrix where the diagonalization method for dense matrices is used", Integer, 200
   ``arnoldi_maxdim``,           "Dimension of the Hessenberg matrix generated by the Arnoldi method",                                     Integer, 50
   ``arnoldi_restartdim``,       "The number of the initial vectors generated by the restart process of the IRA method",                                     Integer, 20
   ``arnoldi_maxiterations``,    "Maximum number of iterations in the IRA method",                                  Integer, 1
   ``arnoldi_rtol``,             "Relative tolerance used in the Arnoldi method",   Float,   1e-10

The correlation length :math:`\xi` will be calculated from the dominant eigenvalues of the transfer matrices.
If the dimension of the transfer matrix is less than or equal to ``maxdim_dense_eigensolver``, an eigensolver for dense matrices (LAPACK's ``*geev`` routines) will be used.
If not, an iterative method, the implicit restart Arnoldi method (IRA method), will be used.

In the IRA method, a Hessenberg matrix with the size of ``arnoldi_maxdim`` is generated by the Arnoldi process.
Its eigenvalues are approximants of the first ``arnoldi_maxdim`` eigenvalues of the original matrix.
If not converged, the IRA method restarts the Arnoldi process with the newly generated ``arnoldi_restartdim`` initial vectors.
In the many cases of the transfer matrices, such a process is not necessary (``arnoldi_maxiterations = 1``).
