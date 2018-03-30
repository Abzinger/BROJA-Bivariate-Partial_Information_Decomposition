# BROJA-Brivariate-Partial_Information_Decomposition
This repository provides code of the experiments performed in the paper A. Makkeh, D.O. Theis, R. Vicente, [*Bivariate Partial Information Decomposition: The Optimization Perspective*](http://www.mdpi.com/1099-4300/19/10/530) Entropy 2017, 19, 530-570;[doi:10.3390/e19100530](dx.doi.org/10.3390/e19100530).

The repository contains different approcahes to compute the Bertschinger-Rauh-Olbrich-Jost-Ay bivariate Partial Information Decomposition, **BROJA PID** (N. Bertschinger, J. Rauh, E. Olbrich, J. Jost, N. Ay, *Quantifying Unique Information*. Entropy 2014, 16, 2161-2183; [doi:10.3390/e16042161](dx.doi.org/10.3390/e16042161).):

* implementations in Julia that work together with several Julia packages, see below for details.

* implementations in Python that work together with the Python package [cvxopt](https://github.com/cvxopt/cvxopt) 

## Julia
The Julia implementation is build using the [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) Julia package which is part of [JuliaOpt](https://github.com/JuliaOpt/MathProgBase.jl). 

#### Dependencies 
The Julia implementation needs several packages and the user can run the file `Julia/Pkg_installer.jl` which installs all the required depedencies

#### Computing BROJA PID
 The Julia impelementation used different approaches to compute *BROJA PID*.  The following is the list of the implemented approaches:

* **Zero Order Optimization:**
  
  - [Artelys Knitro:](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2.) We refer by "Knitro_As" to their Active Set Method.

* **First Order Optimization:** 
  
  - [SCS](https://github.com/JuliaOpt/SCS.jl)(version 1.2.7.) stands for Splitting Conic Solver. It is a numerical optimization package for solving large-scale convex cone problems.
  
  - Projected Gradient descent: is our implementation of projected gradient decent in `Julia/graddesc.jl`  
  
  - [Artelys Knitro:](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2.) We refer by "Knitro_IpCG" to their IPM with conjugate gradient decesnt.

* **Second Order Optimization:** 
  - [*Mosek*](https://github.com/JuliaOpt/Mosek.jl)(version 8.0): is an optimization suite which offers algorithms for a vast range of convex optimization problems.
  
  - [Ipopt](https://github.com/JuliaOpt/Ipopt.jl)(version 3.0): is a software for nonlinear optimization which can also deal with non-convex objectives.At its heart is an Interior Point Method (as the name indicates), which is enhanced by approaches to ensure convergence even in the non-convex case.
  
  - [ECOS](https://github.com/JuliaOpt/ECOS.jl)(version Nov 8, 2016) is a lightweigt numerical software for solving convex cone programs, using an Interior Point approach.
  
  - [Artelys Knitro:](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2.) is an optimization suite which offers four algorithms for general Convex Programming. We refer by "Knitro_Ip" to their standard Interior Point Method.  
  
  -[Artelys Knitro:](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2.) "Knitro_SQP" designates their Sequential Quadratic Programming Method  

## Python 
* **Second Order Optimization**
- [CVXOPT](https://github.com/cvxopt/cvxopt)(version 1.1.9) is written and maintained by Andersen, Dahl, and Vandenberghe. It transforms the general Convex Problems with nonlinear objective function into an epigraph form, before it deploys an Interior Point method.

#### Performing experiments

