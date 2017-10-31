# toy_code

### 1. ```MonteCarlo_Solver.cpp```

#### Description 

Monte Carlo solution of integral equations in the form of:

![equation](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

\begin{equation}
I = sum_k Xk 
\end{equation}
$ \sum_{\forall i}{x_i^{2}} $

#### Run

```
$ g++ -o mc.o MonteCarlo_Solver.cpp 
$ ./mc.o 1000000 2 0 2 1 -1 1 1 1 1.5
```

