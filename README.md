# toy_code

### 1. Monte Carlo Solver

#### Files

* ```MonteCarlo_Solver.cpp```: Monte Carlo solution of integral equations in the form of:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
![equation](https://latex.codecogs.com/gif.latex?%5Cdpi%7B120%7D%20%5Chuge%20I%20%3D%20%5Cint_%7Bx_1%7D%5E%7Bx_2%7D%5Cint_%7By_1%7D%5E%7By_2%7D%5Cint_%7Bz_1%7D%5E%7Bz_2%7Dx%5E%7Bx_0%7Dy%5E%7By_0%7D&plus;z_0e%5E%7B-z%7Ddxdydz)

#### Run

For example:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
![equation](http://latex.codecogs.com/gif.latex?%24%24I%3D%5Cint_%7B0%7D%5E%7B2%7D%5Cint_%7B-1%7D%5E%7B1%7D%5Cint_%7B1%7D%5E%7B1.5%7Dx%5E%7B2%7Dy%5E%7B4%7D&plus;%7B3%7De%5E%7B-z%7Ddxdydz%24%24)

```
$ g++ -o mc.o MonteCarlo_Solver.cpp 
$ ./mc.o 1000000 2 0 2 4 -1 1 3 1 1.5
```

$-b \pm \sqrt{b^2 - 4ac} \over 2a$
$x = a_0 + \frac{1}{a_1 + \frac{1}{a_2 + \frac{1}{a_3 + a_4}}}$
$\forall x \in X, \quad \exists y \leq \epsilon$

