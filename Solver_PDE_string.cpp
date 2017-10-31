/***********************************************************
 Function:
 Solve partial differential equation for a vibrating string.
 Jie Yang 2014
************************************************************/

# include <iostream>
# include <iomanip>
# include <fstream>
# include <math.h>
using namespace std;
int main()
{
    // set parameters
    double h=0.1;  // step length for x
    double k=0.01;  // step length for t
    double T=2;     // maximum time t to solve until
    double p=10;    // graphing rate
    
    int tmax = ceil(T/k)+1;
    int xmax = ceil(1.0/h)+1;
    
    int m = ceil(0.5/h); // medium of x
    
    double rho = pow(k,2.0)/pow(h,2.0);
    
    // t corresponds to row, x corresponds to column
    double **u = new double *[tmax]();
    for (int i=0;i<tmax;i++)
        u[i] = new double[xmax]();
    
    // boundary condition
    // when x=0, u=0; when x=1, u=0
    for (int i=0;i<tmax;i++)
    {
        u[i][0]=0;
        u[i][xmax-1]=0;
    }

    // I take function f give in the assignment as input
    double * f = new double [xmax];
    for (int i=0;i<m;i++)
    {
        f[i]=i*h;
    }
    for (int i=m;i<xmax;i++)
    {
        f[i]=1-i*h;
    }
    
    // initial condition
    for (int i=0;i<xmax;i++)
    {
        u[0][i]=f[i];
        u[1][i]=f[i];
    }
    
    // solve the PDE
    for (int i=1;i<tmax-1;i++)
    {
        for (int j=1;j<xmax-1;j++)
        {
            u[i+1][j] = rho*(u[i][j+1]+u[i][j-1])+2*(1-rho)*(u[i][j])-u[i-1][j];
        }
    }
    
    ofstream file_PDE_string_result;
    file_PDE_string_result.open("./PDE_string_result.txt");
    for (int i=0;i<tmax;i=i+p)
	   {
           
           for (int j=0;j<xmax;j++)
           {
              file_PDE_string_result<<u[i][j]<<'\t';
           }
           file_PDE_string_result<<'\n';
       }
    cout<<"saved to PDE_string_result.txt"<<endl;
    
    return 0;
    
}

