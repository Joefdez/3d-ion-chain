# Heat transport in a three-dimensional trapped ion chain


Fortran 90 code which simulates the dynamics of a three-dimensional trapped ion chain (basically the three-dimensional version of the code available in another repository on this page).

The edges of the chain are in contact with heat baths, inducing heat transport between them. 

As the dynamics of the laser-ion interaction is stochastic, the equations of motion must be solved a large number of times to achieve convergence. The code is parallelized using MPI to accelerate convergence. 

The code has been optimized for speed. Memory optimizations can still be implemented, which would significantly improve performance and parallelizability. 

This code was used to generate science results presented in my Master's dissertation, and which contributed to a paper called "Delocalization and heat transport in multidimensional trapped ion systems", published in Physical Review E (pre-print available at https://arxiv.org/pdf/1901.08870.pdf). 
