# toy_code

### 1. ```MonteCarlo_Solver.cpp```

#### Description 

Monte Carlo solution of integral equations in the form of:

\begin{equation}
I = sum_k Xk 
\end{equation}

**Run**

```
$ g++ -o mc.o MonteCarlo_Solver.cpp 
$ ./mc.o 1000000 2 0 2 1 -1 1 1 1 1.5
```

