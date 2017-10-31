# toy_code

### 1. Gaussian Mixture Model

#### Files
* ```GMM_EM.m```: expectation maximization for Gaussian mixture model.
* ```GMM_VI.m```: variational inference for Gaussian mixture model.
* ```GMM_MCMC.m```: Gibbs sampling for Gaussian mixture model.
* ```data.mat```: input sample.

### 2. Monte Carlo Solver

#### Files

* ```MonteCarlo_Solver.cpp```: Monte Carlo solution of integral equations in the form of:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
![equation](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%5CLARGE%20I%20%3D%20%5Cint_%7Bx_1%7D%5E%7Bx_2%7D%5Cint_%7By_1%7D%5E%7By_2%7D%5Cint_%7Bz_1%7D%5E%7Bz_2%7Dx%5E%7Bx_0%7Dy%5E%7By_0%7D&plus;z_0e%5E%7B-z%7Ddxdydz)

#### Run

For example:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
![equation](https://latex.codecogs.com/gif.latex?%5Cinline%20%5Cdpi%7B100%7D%20%5CLARGE%20I%20%3D%20%5Cint_%7B0%7D%5E%7B2%7D%5Cint_%7B-1%7D%5E%7B1%7D%5Cint_%7B1%7D%5E%7B1.5%7Dx%5E%7B2%7Dy%5E%7B4%7D&plus;3e%5E%7B-z%7Ddxdydz)

```
$ g++ -o mc.o MonteCarlo_Solver.cpp 
$ ./mc.o 1000000 2 0 2 4 -1 1 3 1 1.5
```

### 3. Strassen Solver Solver

#### Files

* ```Strassen_Solver.cpp```: Strassen solver of matrix multiplication C = A * B.
* ```a.txt```: contains matrix A.
* ```b.txt```: contains matrix B.

#### Run

```
$ g++ -o strassen.o Strassen_Solver.cpp 
$ ./strassen.o 
```


