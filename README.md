# BROJA-Brivariate-Partial_Information_Decomposition
This repository provides code of the experiments performed in the paper A. Makkeh, D.O. Theis, R. Vicente, [*Bivariate Partial Information Decomposition: The Optimization Perspective*](http://www.mdpi.com/1099-4300/19/10/530) Entropy 2017, 19, 530-570; [doi:10.3390/e19100530](dx.doi.org/10.3390/e19100530).

The repository contains different approcahes to compute the Bertschinger-Rauh-Olbrich-Jost-Ay bivariate Partial Information Decomposition, **BROJA PID** (N. Bertschinger, J. Rauh, E. Olbrich, J. Jost, N. Ay, *Quantifying Unique Information*. Entropy 2014, 16, 2161-2183; [doi:10.3390/e16042161](dx.doi.org/10.3390/e16042161).):

* implementations in Julia that work together with several Julia packages, see below for details.

* implementations in Python that work together with the Python package [cvxopt](https://github.com/cvxopt/cvxopt).

## Julia
The Julia implementation is build using the [MathProgBase](https://github.com/JuliaOpt/MathProgBase.jl) Julia package which is part of [JuliaOpt](https://github.com/JuliaOpt/MathProgBase.jl). 

### Dependencies 
The Julia implementation needs several packages and the user can run the file `Julia/Pkg_installer.jl` which installs all the required depedencies

### Computing BROJA PID
 The Julia impelementation used different approaches to compute *BROJA PID*.  The following is the list of the implemented approaches:

* **Zero Order Optimization:**
  
  - [Artelys Knitro](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2): We refer by "Knitro_As" to their Active Set Method.

* **First Order Optimization:** 
  
  - [SCS](https://github.com/JuliaOpt/SCS.jl)(version 1.2.7): stands for Splitting Conic Solver. It is a numerical optimization package for solving large-scale convex cone problems.
  
  - Projected Gradient descent: is our implementation of projected gradient decent in `Julia/graddesc.jl`  
  
  - [Artelys Knitro](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2): We refer by "Knitro_IpCG" to their IPM with conjugate gradient decesnt.

* **Second Order Optimization:** 
  - [*Mosek*](https://github.com/JuliaOpt/Mosek.jl)(version 8.0): is an optimization suite which offers algorithms for a vast range of convex optimization problems.
  
  - [Ipopt](https://github.com/JuliaOpt/Ipopt.jl)(version 3.0): is a software for nonlinear optimization which can also deal with non-convex objectives.At its heart is an Interior Point Method (as the name indicates), which is enhanced by approaches to ensure convergence even in the non-convex case.
  
  - [ECOS](https://github.com/JuliaOpt/ECOS.jl)(version Nov 8, 2016): is a lightweigt numerical software for solving convex cone programs, using an Interior Point approach.
  
  - [Artelys Knitro](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2): is an optimization suite which offers four algorithms for general Convex Programming. We refer by "Knitro_Ip" to their standard Interior Point Method.  
  
  - [Artelys Knitro](https://github.com/JuliaOpt/KNITRO.jl)(version 10.2): "Knitro_SQP" designates their Sequential Quadratic Programming Method  

## Python 

### Dependencies 
The python implementation requires the following external packages to be insatlled:
* `numpy`
* `cvxopt`
* `gurobipy`

### Computing BROJA PID

* **Second Order Optimization**
  - [CVXOPT](https://github.com/cvxopt/cvxopt)(version 1.1.9): is written and maintained by Andersen, Dahl, and Vandenberghe. It transforms the general Convex Problems with nonlinear objective function into an epigraph form, before it deploys an Interior Point method. We use two ways to compute PID with the aid of **cvxopt**
    - Compute the solution of BROJA PID and use `cvxopt` as an *interior point solver*. It is implemented in `Python/cvxopt.py`
    - Transform the problem into a *Geometric Program* and use `cvxopt` as a *geometric programming solver*. It is implemented in `Python/cvxopt_geo_solve.py`
  
 * **Active set Optimization**

   - [Gurobi](https://www.gurobi.com/documentation/7.5/quickstart_windows/py_python_interface)(version 7.5), [CVXOPT](https://github.com/cvxopt/cvxopt)(version 1.1.9): Enhancement of the `Python/cvxopt.py` which have bails out when solution vanishes at any point, i.e., *q(x,y,z) = 0* for some *x,y,z*. The ad-hoc technique fixes *q(x,y,z) =0* and then checks whether the solution is optimal by checking *KKT* first order optimality condotionswhich is a *Linear Program* using gurobi.
  

## Experiments

* **Stored Instances**

  - Gates (Pure): These are paradagmatic gates in which the solution can be computed theoretically and they are well structured and considered to be small instances. The gates used are `XOR`, `AND`, `UNQ`, `RDN`, `XORAND`, `RDNUNQ`, `RDNUNQXOR`. For more details check [*Bivariate Partial Information Decomposition: The Optimization Perspective*](http://www.mdpi.com/1099-4300/19/10/530). The isntances are stored as follows: `PDFs/[gate_name]-0.0-0.dens`.
  - Gates (Noisy): These are paradagmatic gates distorted with noise. There are five different levels of noise. For more details check [*Bivariate Partial Information Decomposition: The Optimization Perspective*](http://www.mdpi.com/1099-4300/19/10/530). The isntances are stored as follows: `PDFs/Data/[gate_name]-0.0-10000.dens`, `PDFs/Data/[gate_name]-0.001-1000000.dens`,`PDFs/Data/[gate_name]-0.01-100000.dens`,`PDFs/Data/[gate_name]-0.05-100000.dens`, and `PDFs/Data/[gate_name]-0.1-10000.dens`. 
  - Uniformly at random: *X,Y,Z* are uniformly drawn at random in *{0,1}* with *|X|={2,...,5}* and *|Y|=|Z|={|X|,..., 5|X|}* 
    - The instances which are idependently drawn are stored as follows:`PDFs/broad_indep--[|X|]-[|Y|]-[|Z|].dens`
    - The instances which are dependently drawn are stored as follows:`PDFs/broad_dep--[|X|]-[|Y|]-[|Z|].dens`
  - Discretized Gaussian: These are generated using 20 *3*x*3* (randomly) generated covariance matrices of the gaussians random variables *X,Y,Z*. Then, the resulting distribution is discretized over different boxes *[0,b]* where *b\in{0,0.25,0.5,0.75,1}*.  For more details check [*Bivariate Partial Information Decomposition: The Optimization Perspective*](http://www.mdpi.com/1099-4300/19/10/530). The isntances are stored as follows: `PDFs/Data/gauss-[matrix_number]-[b].dens`.
 
 * **Generate new Instances**
 
 To generate new instances you can use the following:
   
   - `Python/time_series` generates *noisy gate* instances  and store them into files
   ``` 
   # Example: The following command generates `XOR` gate instances 
   # with noise `0.1` 
   # probability is estimated by sampling `10000` outputs of the gate 
   
   python3 time_series.py -w XOR 0.1 10000
   ```
   - `Julia/gaussians.jl` is a julia file which generates discretized gaussian instances using the following: 
         - `random_Σ()` to generates random 3x3 covarient matrix 
         - `make_apxgaussian_pdf(Σ,s)` discritize via the desired box persicion where *Σ* is the covariance matrix and *s* is the box size 
         - Note that the *pdf* should be stored approperiatly for later compuations see the julia script `generate_guassians.jl` for as an example to how the user should store the *pdfs*
 
 
   * **Preform the experiments**
   
