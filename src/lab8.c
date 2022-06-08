#include <stdio.h>
#include <math.h>

double f0( double x, double l, double t )
{
    return sqrt(1.0+l*l*t*t*cos(t*x)*cos(t*x));
}

double Integral(double a, double b, double (*f)(double x, double y, double z), 
double eps, double l, double t);

int main()
{
    double a=0.0, b, eps=0.005, l, t;

    scanf("%lf %lf %lf", &l, &b, &t);
    printf("%.2f\n", Integral(a, b, f0, eps, l, t));

    return 0;
}

/* Your function will be put here */
#define SIMPSON(interval, f0, f1, f2) (interval * (f0 + 4 * f1 + f2) / 6)

double Integral(double a, double b, double (*f)(double x, double y, double z), 
double eps, double l, double t)
{
    double h = a + (b - a) / 4, fs[5], res1, res2;
    for (int i = 0; i < 5; i++) {
        fs[i] = f(a + i * h, l, t) / 100;
        // printf("f(%lf) = %lf\n", a + i * h, fs[i]);
    }
    res1 = SIMPSON(2 * h, fs[0], fs[2], fs[4]);
    res2 = SIMPSON(h, fs[0], fs[1], fs[2]) + SIMPSON(h, fs[2], fs[3], fs[4]);
    // printf("res1: %lf, res2: %lf\n", res1, res2);
    if (fabs((res1 - res2) / 16) < eps) {
        return res2;
    }
    else {
        return Integral(a, a + 2 * h, f, eps / 2, l, t) + Integral(a + 2 * h, b, f, eps / 2, l, t);
    }
}