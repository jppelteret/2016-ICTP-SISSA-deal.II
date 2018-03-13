#  Lab 09 - matrix-free computations with MPI
## Numerical Solution of PDEs Using the Finite Element Method
### MHPC P2.13_seed

**Martin Kronbichler** <kronbichler@lnm.mw.tum.de>
and
**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

1. Run the included step-37 example. Compile the code in release mode (``make
   release'') and run with 1 and 2 MPI processes. Record the solver time and
   the number of iterations.

2. This example uses a geometric multigrid preconditioner that gives constant
   iteration numbers as the mesh is refined, as opposed to the
   ``PreconditionIdentity`` that was used in the early labs. Check the number
   of iterations with the conjugate gradient method.

3. Enable the code ``compare_performance_matrix_vector()`` in the ``run()``
   function. Study the steps in this function. Which of the matrix-based and
   the matrix-free product is faster on your machine?

4. Change the polynomial degree from ``2`` (set at the top of the file through
   the variable ``degree_finite_element``) to ``1`` and ``3``. Observe the
   difference in the performance of the matrix-vector product. Be careful with
   the large problem sizes and degree ``3`` as the sparse matrix could run
   out of RAM memory.

5. Compare the ratio between the sparse matrix-vector with the Trilinos sparse
   matrix and the matrix-free ``LaplaceOperator`` when run on ``1`` processor
   core with the ratio when using all processors on your system. What could be
   an explanation of the difference in the ratio?
