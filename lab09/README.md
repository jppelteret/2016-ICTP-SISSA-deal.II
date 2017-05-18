#  Lab 08 - MPI
## Numerical Solution of PDEs Using the Finite Element Method
### MHPC P2.13_seed

**Martin Kronbichler** <kronbichler@lnm.mw.tum.de>
and
**Luca Heltai** <luca.heltai@sissa.it>

* * * * *

1. Repeat the basic MPI commands from the file mpihello/main.cc

2. Run the included step-40 using ``mpirun -n 4 ./step-40`` and look at the
   graphical output.

3. Similar to shown in lecture, visualize the view of the mesh from each
   individual processor using ``GridOut::write_svg`` and the "global"
   mesh. Use 3 MPI tasks for this.

4. Now create a simple mesh (hyper_cube refined twice globally), run with two
   MPI tasks and print locally owned, locally active, and locally relevant
   ``IndexSet`` for each task.

5. Switch to release mode (``make release``), decide on a global refinement
   level that takes in the order of 30-60 seconds to solve, and study assembly
   and solve time with 1,2,4,8,12,16 MPI tasks. Which is the fastest, do the
   timings make sense based on how many cores your machine has?

6. Play with the test problem by switching to 3d and changing the geometry to
   something interesting. Your choice!
