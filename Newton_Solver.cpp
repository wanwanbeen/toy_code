/***********************************************************
 Function:
 Newton solver of a function f
 
 Jie Yang 2014
 
 inputs:
 val, w1, w2, w3 = function parameters
 - f = w1*exp(-w2t)cos(w3t)
 - for f = val, solve t
************************************************************/

#include<iostream>
#include<math.h>
#include<stdlib.h>
#include<time.h>
#define PI 3.1415926

using namespace std;

double InputFunction(double t, double w1, double w2, double w3, double val)
{
    double f;
    f=w1*exp(-w2*t)*cos(w3*t)-val;
    return f;
}

double DiffFunction(double t,double w1, double w2, double w3, double val, double step)
{
    double f1;
    f1 = (InputFunction(t+step, w1, w2, w3, val)-InputFunction(t, w1, w2, w3, val))/step;
    return f1;
}

double NewtonSolver(double t,double w1, double w2, double w3, double val, double step, double e)
{
    int LoopTime = 0;
    int LoopMax = 100;
    double f=InputFunction(t, w1, w2, w3, val);
    double f1=DiffFunction(t, w1, w2, w3, val, step);
    cout<<"t0="<<t<<" "<<"f0="<<f+val<<endl;
    while ((((f>=e)||(f<=-e))&&LoopTime<LoopMax))
    {
        t = t - f/f1;
        f = InputFunction(t, w1, w2, w3, val);
        f1 = DiffFunction(t, w1, w2, w3, val, step);
        cout<<"t="<<t<<" "<<"f="<<f+val<<endl;

        LoopTime++;
    };
    return t;
}

int main(int argc, char **argv)
{
    double step     = 0.00001;
    double e        = 0.00001;
    srand (time(NULL));
    double t0       = double(rand()%100)/100;
    double w1       = stold(argv[1]);
    double w2       = stold(argv[2]);
    double w3       = stold(argv[3]);
    double val      = stold(argv[4]);
    
    double t = NewtonSolver(t0, w1, w2, w3, val, step, e);
    double f = InputFunction(t, w1, w2, w3, val);
    cout<<"Solution: "<<"t_result="<<t<<" "<<"f_result="<<f+val<<endl;
    return 0;
}
