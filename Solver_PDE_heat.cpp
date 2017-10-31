/***********************************************************
 Function:
 Solve partial differential equation for a heat.
 Jie Yang 2014
 ************************************************************/

# include <iostream>
# include <fstream>
# include <math.h>
# define PI 3.1415926
using namespace std;

double * Initial_Condition(int xMax,double h)
{
    double * f = new double [xMax];
    for (int i=0;i<xMax;i++)
    {
        f[i] = sin(PI*double(i)*h);
    }
    return f;
}


int main()
{
    // set parameters
    double h_Stable=0.1;
    double h_Unstable=0.1;

    double k_Stable=0.005;
    double k_Unstable=0.1;

    double T=2;
    double p=0.1;
    
    double sigma_Stable = k_Stable/pow(h_Stable,2.0);
    double sigma_Unstable = k_Unstable/pow(h_Unstable,2.0);

    
    int tMax_Stable = ceil(T/k_Stable)+1;
    int xMax_Stable = ceil(1.0/h_Stable)+1;
    
    int tMax_Unstable = ceil(T/k_Unstable)+1;
    int xMax_Unstable = ceil(1.0/h_Unstable)+1;
    
    //***************************************************************************
    // Stable parameter
    double ** u_Stable_Analytical = new double *[tMax_Stable]();
    for (int i=0;i<tMax_Stable;i++)
        u_Stable_Analytical[i] = new double[xMax_Stable]();
    
    for (int i=0;i<tMax_Stable;i++)
    {
        for (int j=0;j<xMax_Stable-1;j++)
            u_Stable_Analytical[i][j]=exp(-pow(PI,2.0)*i*k_Stable)*sin(PI*j*h_Stable);
        u_Stable_Analytical[i][xMax_Stable-1]=0;
    }
   double **u_Stable_Program = new double *[tMax_Stable]();
    for (int i=0;i<tMax_Stable;i++)
        u_Stable_Program[i] = new double[xMax_Stable]();
 
    // boundary condition
    for (int i=0;i<tMax_Stable;i++)
    {
        u_Stable_Program[i][0]=0;
        u_Stable_Program[i][xMax_Stable-1]=0;
    }
    
    // initial condition
    double * f_Stable = Initial_Condition(xMax_Stable,h_Stable);
    
    for (int i=0;i<xMax_Stable;i++)
    {
        u_Stable_Program[0][i]=f_Stable[i];
    }
    u_Stable_Program[0][xMax_Stable-1]=0;
    
    // recurrence
    for (int i=0;i<tMax_Stable-1;i++)
    {
        for (int j=1;j<xMax_Stable-1;j++)
        {
            u_Stable_Program[i+1][j] = sigma_Stable*(u_Stable_Program[i][j+1]+u_Stable_Program[i][j-1])+(1-2*sigma_Stable)*(u_Stable_Program[i][j]);
        }
    }
        
    //***************************************************************************
    // Unstable parameter
    double ** u_Unstable_Analytical = new double *[tMax_Unstable]();
    for (int i=0;i<tMax_Unstable;i++)
        u_Unstable_Analytical[i] = new double[xMax_Unstable]();
    
    for (int i=0;i<tMax_Unstable;i++)
    {
        for (int j=0;j<xMax_Unstable-1;j++)
            u_Unstable_Analytical[i][j]=exp(-pow(PI,2.0)*i*k_Unstable)*sin(PI*j*h_Unstable);
        u_Unstable_Analytical[i][xMax_Unstable-1]=0;
    }
    
    double **u_Unstable_Program = new double *[tMax_Unstable]();
    for (int i=0;i<tMax_Unstable;i++)
        u_Unstable_Program[i] = new double[xMax_Unstable]();
    
    // boundary condition
    for (int i=0;i<tMax_Unstable;i++)
    {
        u_Unstable_Program[i][0]=0;
        u_Unstable_Program[i][xMax_Unstable-1]=0;
    }
    
    // initial condition
    double * f_Unstable = Initial_Condition(xMax_Unstable,h_Unstable);
    
    for (int i=0;i<xMax_Unstable;i++)
    {
        u_Unstable_Program[0][i]=f_Unstable[i];
    }
    u_Unstable_Program[0][xMax_Unstable-1]=0;
    
    // recurrence
    for (int i=0;i<tMax_Unstable-1;i++)
    {
        for (int j=1;j<xMax_Unstable-1;j++)
        {
            u_Unstable_Program[i+1][j] = sigma_Unstable*(u_Unstable_Program[i][j+1]+u_Unstable_Program[i][j-1])+(1-2*sigma_Unstable)*(u_Unstable_Program[i][j]);
        }
    }
    
    // save value to file
    ofstream file_Stable_Analytical;
    ofstream file_Unstable_Analytical;
    file_Stable_Analytical.open("PDE_heat_Stable_Analytical.txt");
    file_Unstable_Analytical.open("PDE_heat_Unstable_Analytical.txt");
    
    for (int i=0;i<tMax_Stable;i=i+p/k_Stable)
    {
        for (int j=0;j<xMax_Stable;j++)
        {
            file_Stable_Analytical<<u_Stable_Analytical[i][j]<<'\t';
        };
        file_Stable_Analytical<<'\n';
    }
    cout<<"saved to PDE_heat_Stable_Analytical.txt"<<endl;
    
    for (int i=0;i<tMax_Unstable;i=i+p/k_Unstable)
    {
        for (int j=0;j<xMax_Unstable;j++)
        {
            file_Unstable_Analytical<<u_Unstable_Analytical[i][j]<<'\t';
        };
        file_Unstable_Analytical<<'\n';
    }
    cout<<"saved to PDE_heat_Unstable_Analytical.txt"<<endl;
    
    ofstream file_Stable_Program;
    ofstream file_Unstable_Program;
    file_Stable_Program.open("PDE_heat_Stable_Program.txt");
    file_Unstable_Program.open("PDE_heat_Unstable_Program.txt");
    
    for (int i=0;i<tMax_Stable;i=i+p/k_Stable)
    {
        for (int j=0;j<xMax_Stable;j++)
        {
            file_Stable_Program<<u_Stable_Program[i][j]<<'\t';
        };
        file_Stable_Program<<'\n';
    }
    cout<<"saved to PDE_heat_Stable_Program.txt"<<endl;
    
    for (int i=0;i<tMax_Unstable;i=i+p/k_Unstable)
    {
        for (int j=0;j<xMax_Unstable;j++)
        {
            file_Unstable_Program<<u_Unstable_Program[i][j]<<'\t';
        };
        file_Unstable_Program<<'\n';
    }
    cout<<"saved to PDE_heat_Unstable_Program.txt"<<endl;

    // error file
    ofstream file_Stable_Error;
    ofstream file_Unstable_Error;
    file_Stable_Error.open("PDE_heat_Stable_Error.txt");
    file_Unstable_Error.open("PDE_heat_Unstable_Error.txt");
    
    for (int i=0;i<tMax_Stable;i=i+p/k_Stable)
    {
        for (int j=0;j<xMax_Stable;j++)
        {
            file_Stable_Error<<u_Stable_Analytical[i][j]-u_Stable_Program[i][j]<<'\t';
        };
        file_Stable_Error<<'\n';
    }
    cout<<"saved to PDE_heat_Stable_Error.txt"<<endl;
    
    for (int i=0;i<tMax_Unstable;i=i+p/k_Unstable)
    {
        for (int j=0;j<xMax_Unstable;j++)
        {
            file_Unstable_Error<<u_Unstable_Analytical[i][j]-u_Unstable_Program[i][j]<<'\t';
        };
        file_Unstable_Error<<'\n';
    }
    cout<<"saved to PDE_heat_Untable_Error.txt"<<endl;
    
    // average error
    double * Stable_Error_AVG = new double [tMax_Stable];
    double * Unstable_Error_AVG = new double [tMax_Unstable];

    
    ofstream file_Stable_Error_AVG;
    ofstream file_Unstable_Error_AVG;

    file_Stable_Error_AVG.open("./PDE_heat_Stable_Error_AVG.txt");
    file_Unstable_Error_AVG.open("./PDE_heat_Unstable_Error_AVG.txt");
    
    for (int i=0;i<tMax_Stable;i=i+p/k_Stable)
    {
        Stable_Error_AVG[i]=0;
        for (int j=0;j<xMax_Stable;j++)
        {
            Stable_Error_AVG[i]=Stable_Error_AVG[i]+(u_Stable_Analytical[i][j]-u_Stable_Program[i][j])/double(xMax_Stable);
        };
        file_Stable_Error_AVG<<Stable_Error_AVG[i]<<endl;
    }
    cout<<"saved to PDE_heat_Stable_Error_AVG.txt"<<endl;
    
    for (int i=0;i<tMax_Unstable;i=i+p/k_Unstable)
    {
        Unstable_Error_AVG[i]=0;
        for (int j=0;j<xMax_Unstable;j++)
        {
            Unstable_Error_AVG[i]=Unstable_Error_AVG[i]+(u_Unstable_Analytical[i][j]-u_Unstable_Program[i][j])/double(xMax_Unstable);
        };
        file_Unstable_Error_AVG<<Unstable_Error_AVG[i]<<endl;
    }
    cout<<"saved to PDE_heat_Unstable_Error_AVG.txt"<<endl;
    
    
    return 0;
    
}

