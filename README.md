# Circuit-Simulation-Program-for-SPICE-netlists

This project is a C program for circuit simulation of SPICE netlists.

In particular, it supports the following operations: 
* DC Analysis
* DC Sweep (using a .DC command)
* Transient Analysis (using a .TRAN command)
* AC Analysis (using a .AC command)

For the solution of the occuring linear systems, the following methods are available:
* LU decomposition
* Cholesky factorization (including SPD in the .OPTIONS command)
* Bi-conjugate gradient (including ITER in the .OPTIONS command)
* Conjugate gradient (including ITER and SPD in the .OPTIONS command)

For the last two methods which are iterative, it is also possible to define a tolerance threshold by including ITOL=[value] in the .OPTIONS command.

There are implementations for a couple of additional methods for transient analysis, namely:
* Trapezoidal (activated with METHOD=TR in the .OPTIONS command)
* Backward-Euler (activated with METHOD=BE in the .OPTIONS command)

All the above methods are implemented for sparse matrices as well and can be used by including SPARSE in the .OPTIONS command. First you need to unzip CXSparse.zip in the libs folder in order to try this alternative.

The tests folder contains an indicative example of a netlist used for AC analysis. You can run this example using the given Makefile.
