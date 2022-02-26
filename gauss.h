#include<iostream>
#include<math.h>
using namespace std;

struct gauss
{
    double f1(double x);
    double f2(double x, double y);
    double Gauss1d( double (*function)(double), int pc );
    double Gauss2d( double (*function)(double, double), int pc);

};

double gauss::f1(double x)
{
    return 5*(x*x) + 3*x + 6;
}

double gauss::f2(double x, double y)
{
    return 5*(x*x)*(y*y) + 3*x*y + 6;
}

double gauss::Gauss1d( double (*function)(double), int pc )
{
    double result = 0.0;
    if(pc == 2)
    {
        double x = 1.0/sqrt(3.0), w1 = 1;
        result = w1 * function(x);
        result += w1 * function(-x);
    }
    else if (pc == 3)
    {
        
        double x = sqrt(3.0/5.0), x1 = 0, w1 = 5.0/9.0, w2 = 8.0/9.0;
        result = w2 * function(x1);
        result += w1 * function(x);
        result += w1 * function(-x);
    }
    else
    {
        result = NULL;
    }
    return result;
}

double gauss::Gauss2d( double (*function)(double, double), int pc)
{
    double result = 0.0;
    if(pc == 2)
    {
        double x1 = 1.0/sqrt(3), w1 = 1;
        result = w1*w1*function(-x1, -x1);
        result += w1*w1*function(x1, -x1);
        result += w1*w1*function(x1, x1);
        result += w1*w1*function(-x1, -x1);
    }
    else if(pc == 3)
    {
        double x1 =sqrt(3.0/5.0), x2 = 0, w1 = 5.0/9.0, w2 = 8.0/9.0;
        result = w1*w1*function(-x1, -x1);
        result += w1*w2*function(x2, -x1);
        result += w1*w1*function(x1, -x1);
        result += w1*w2*function(-x1, x2);
        result += w2*w2*function(x2, x2);
        result += w1*w2*function(x1, x2);
        result += w1*w1*function(-x1, x1);
        result += w1*w2*function(x2, x1);
        result += w1*w1*function(x1, x1);
    }
    else
    {
        result = NULL;
    }
    return result;
}