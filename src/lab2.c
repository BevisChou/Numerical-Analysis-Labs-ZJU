#include <stdio.h>
#include <math.h>

#define ZERO 1e-13 /* X is considered to be 0 if |X|<ZERO */
#define MAXN 11    /* Max Polynomial Degree + 1 */

double Polynomial_Root(int n, double c[], double a, double b, double EPS);

int main()
{
    int n;
    double c[MAXN], a, b;
    double EPS = 0.00005;
    int i;

    scanf("%d", &n);
    for (i=n; i>=0; i--) 
        scanf("%lf", &c[i]);
    scanf("%lf %lf", &a, &b);
    printf("%.4f\n", Polynomial_Root(n, c, a, b, EPS));

    return 0;
}

/* Your function will be put here */
#include <math.h>
#define MAXI 1000000
#define N 1000

double Calc(int n, double c[], double x)
{
    double res = 0;
    for (int i = n; i > 0; i--) {
        res = (res + c[i]) * x;
    }
    res += c[0];
    return res;
}

void Derivative(int n, double c[], double res[])
{
    for (int i = 1; i <= n; i++) {
        res[i - 1] = i * c[i];
    }
}

int Newton(int n, double c[], double *x, double EPS)
{
    // printf("---\ninitial point: %lf\n", *x);
    double dc[MAXN - 1], last, temp;
    Derivative(n, c, dc);
    for (int i = 0; i < MAXI; i++) {
        last = *x;
        temp = Calc(n - 1, dc, last);
        if (temp == 0) { 
            return 0;
        }
        *x = last - Calc(n, c, last) / temp;
        // printf("iter %d: x: %.10f\n", i, *x);
        if (fabs(*x - last) < EPS) {
            return 1;
        }
    }
    return 0;
}

double Polynomial_Root(int n, double c[], double a, double b, double EPS)
{
    double res, h = (b - a) / N;
    for (int i = 0; i <= N; i++) {
        res = a + i * h;
        if (Newton(n, c, &res, EPS)) {
            // printf("result: %e\n", res);
            return res;
        }
    }
    // printf("failed\n");
    return 0;
}

/*
case 1:
2 1 1.5 -1
0 2

case 2:
2 1 0 0
-1 1

case 3:
4 1 0 0 0 0
-1000000 1000000

*/