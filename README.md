# toy_code

### 1. ```MonteCarlo_Solver.cpp```

#### Description 

Monte Carlo solution of integral equations in the form of:

![equation](http://latex.codecogs.com/gif.latex?I%3D%5Cint_%7Bx_1%7D%5E%7Bx_2%7D%5Cint_%7By_1%7D%5E%7By_2%7D%5Cint_%7Bz_1%7D%5E%7Bz_2%7Dx%5E%7Bx_0%7Dy%5E%7By_0%7D&plus;%7Bz_0%7De%5E%7B-z%7Ddxdydz) 

#### Run

```
$ g++ -o mc.o MonteCarlo_Solver.cpp 
$ ./mc.o 1000000 2 0 2 1 -1 1 1 1 1.5
```

