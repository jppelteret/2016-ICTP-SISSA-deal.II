#  Lab 06 - Higher Order Mappings
## Numerical Solution of PDEs Using the Finite Element Method
### MHPC P2.13_seed

**Martin Kronbichler** <kronbichler@lnm.mw.tum.de>
and
**Luca Heltai** <luca.heltai@sissa.it>

* * * * *


1.  The topic of this lab session is a modified version of step-4 made
    available for you
    <https://www.dealii.org/8.5.0/doxygen/deal.II/step_4.html>

2.  For more information in higher order mappings see step-10\
    <https://www.dealii.org/8.5.0/doxygen/deal.II/step_10.html>

3.  Learn how deal.II passes the mapping to the functions that involve
    integrals and evaluation of positions, namely `FEValues` and
    `VectorTools::interpolate_boundary_values`.

4.  Run the program and check the graphical and text output.

5.  Adjust the right-hand side and solution to get the manufactured solution
    $$u(x,y) = r^2 \sin(3\theta)\cos(\frac{1}{2} \pi r)$$ and apply zero
    boundary conditions. You can use wolframalpha.com to compute $- \triangle
    u$. Make sure the L2 errors are converging. Also look at the error field
    (in addition to the numerical solution called `solution` and the
    analytical solution `analytic_solution`) in ParaView.

6.  What mapping degree is required to get optimal convergence of the L2
    error based on the polynomial degree of the finite element space?

7.  Try getting H1 convergence to work correctly too.

8.  Change the right hand side and the class `SolutionValues` according to the
    analytic solution $$u(x) = \sin(\pi x )\cdot\cos(\pi y)$$ from lab05, but
    now solved on the circle. Record the convergence rates as you increase the
    polynomial degree. Look at the error field in the ParaView output. Where
    is the error largest?

9.  In the asymptotic regime, the highest possible convergence rate on the
    cirlce appears to be 3.5 irrespective of the degree of the
    polynomial. Change the mesh to `GridGenerator::hyper_shell` and set the
    manifold of all cells to the spherical manifold, not just the faces on the
    boundary. Ensure that you get optimal convergence rates $$h^{p+1}$$ in the
    L2 norm for pth-degree polynomials.