/***********************************************************
Function:
Monte Carlo Solution of Integral Equations

Jie Yang 2014

inputs: 
niter = number of iterations
x0, x1, x2, y0, y1, y2, z0, z1, z2 = 
        integral parameters
        - I = triple integral of (x^x0*y^y0+z0*e^(-z))
        - x ranges from x1 to x2
        - y from y1 to y2
        - z from z1 to z2
************************************************************/

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;

/*--------------------------------
 Define Monte Carlo Solver
--------------------------------*/

double MonteCarloSolver(int niter, double x0, double x1, double x2, double y0, double y1, double y2, double z0, double z1, double z2)

{
    double I = 0;      // I is the integral value using Monte Carlo
    
    for (int i = 0; i < niter; i++)
    {
        /*----------------------------------------------------------
         - The integral is convert to unit cube
         - In unit cube, x y z are all range from 0 to 1
         - rand() generates a PRN ranging from 0 to RAND_MAX=32767
         ---------------------------------------------------------*/
        
        double randx = double(rand())/double(RAND_MAX);
        double randy = double(rand())/double(RAND_MAX);
        double randz = double(rand())/double(RAND_MAX);
        
        double scalar = (x2-x1)*(y2-y1)*(z2-z1);
        
        I = I + scalar*(pow((x2-x1)*randx+x1,x0)*pow((y2-y1)*randy+y1,y0)+z0*exp(-1*((z2-z1)*randz+z1)));
    }
    I = I/niter;
    
    return I;
    
}

/*--------------------------------
 Main Function
 --------------------------------*/

int main(int argc, char **argv)
{
    
    int     niter   = atoi(argv[1]);
    double  x0      = stold(argv[2]);
    double  x1      = stold(argv[3]);
    double  x2      = stold(argv[4]);
    double  y0      = stold(argv[5]);
    double  y1      = stold(argv[6]);
    double  y2      = stold(argv[7]);
    double  z0      = stold(argv[8]);
    double  z1      = stold(argv[9]);
    double  z2      = stold(argv[10]);
    
    double I = MonteCarloSolver(niter, x0, x1, x2, y0, y1, y2, z0, z1, z2);
    
    cout << "Monte Carlo Solution of Integral Equations: " << endl;
    cout << "iteration = " << niter << "  I = " << I << endl;
    
    return 0;
}
