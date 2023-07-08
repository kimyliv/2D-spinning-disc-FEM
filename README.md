# 2D-spinning-disc-FEM

The program calculates displacement, stresses and strains for a 2D spinning disc (neglible thickness).
Axisymmetric_disc_FEM.m is the main matlab file.

The program finds the solution by using Finite Element Method with 2 - node elements. The integration 
is done using numerical integration (Gaussian quadrature).

The system of equations needed for the FEM formulation is derived using the principle of virtual work,
which states that the sum of internal and external work is zero at equilibrium.
