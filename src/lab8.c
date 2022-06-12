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
#define N 1000000

double Simpson(double fs[], int n, double h)
{
    double res = fs[0] + fs[n];
    for (int i = 1; i < n; i += 2) {
        res += 4 * fs[i];
    }
    for (int i = 2; i < n; i += 2) {
        res += 2 * fs[i];
    }
    res = res * h / 3;
    return res;
}

double Integral(double a, double b, double (*f)(double x, double y, double z), 
double eps, double l, double t)
{
    double T, h, fs[N + 1], integral_T;
    T = M_PI / t;
    h = T / N;
    for (int i = 0; i <= N; i++) {
        fs[i] = f(a + i * h, l, t);
    }
    integral_T = Simpson(fs, N, h);
    int count = floor((b - a) / T), r = ceil((b - a - count * T) / h);
    if (r % 2) {
        r++;
    }
    return (integral_T * count + Simpson(fs, r, h)) / 100;
}