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
#include <stdlib.h>
#include <math.h>

#define MAXI 1000
#define N 1000

double fs[N + 1];

double calc(int n, double c[], double x);

int cmp(const void *lhs, const void *rhs) 
{
    int *l = (int *) lhs;
    int *r = (int *) rhs;
    if (fs[*l] < fs[*r]) return -1;
    else if (fs[*l] > fs[*r]) return 1;
    return 0;
}

void sample(int n, double c[], double a, double h, int idx[])
{
    for (int i = 0; i <= N; i++) {
        idx[i] = i;
        fs[i] = fabs(calc(n, c, a + i * h));
    }
    qsort(idx, N + 1, sizeof(int), cmp);
}

double calc(int n, double c[], double x)
{
    double res = 0;
    for (int i = n; i > 0; i--) {
        res = (res + c[i]) * x;
    }
    res += c[0];
    return res;
}

void derivative(int n, double c[], double res[])
{
    for (int i = 1; i <= n; i++) {
        res[i - 1] = i * c[i];
    }
}

int newton(int n, double c[], double d1c[], double d2c[], double x, double *res, double EPS)
{
    double d0 = calc(n, c, x), d1 = calc(n - 1, d1c, x), temp = d1 * d1 - d0 * calc(n - 2, d2c, x);
    if (fabs(temp) < ZERO) {
        return 1;
    }
    *res = x - d0 * d1 / temp;
    if (fabs(d0 / d1) < EPS) {
        return -1;
    }
    else {
        return 0;
    }
}

int trial(int n, double c[], double d1c[], double d2c[], double *x, double EPS)
{
    double ps[3], d0, d1;
    ps[0] = *x;
    // printf("---\nx0: %e\n", ps[0]);
    for (int i = 0; i < MAXI; i++) {
        if (fabs(calc(n, c, ps[0])) < ZERO) {
            *x = ps[0];
            // printf("ps[0] - f(%lf) = %lf\n", ps[0], calc(n, c, ps[0]));
            return 0;
        }
        // printf("before newton 1\n");
        switch (newton(n, c, d1c, d2c, ps[0], &ps[1], EPS)) {
            case 1:
                // printf("newton 1 failed\n");
                return 1;
            case -1:
                *x = ps[1];
                return 0;
            default:
                if (fabs(calc(n, c, ps[1])) < ZERO) {
                    *x = ps[1];
                    return 0;
                }
        }
        // printf("ps[1]: %lf\n", ps[1]);
        switch (newton(n, c, d1c, d2c, ps[1], &ps[2], EPS)) {
            case 1:
                // printf("newton 2 failed\n");
                return 1;
            case -1:
                *x = ps[2];
                return 0;
            default:
                if (fabs(calc(n, c, ps[2])) < ZERO) {
                    *x = ps[2];
                    return 0;
                }
        }
        d0 = ps[1] - ps[0];
        d1 = ps[2] - ps[1];
        if (fabs(d1 - d0) < ZERO) {
            // printf("ps[2]: %lf - to fail\n", ps[2]);
            return 1;
        }
        ps[0] = ps[0] - d0 * d0 / (d1 - d0);
        // printf("iter %d: %e\n", i, ps[0]);
    }
    return 1;
}

double Polynomial_Root(int n, double c[], double a, double b, double EPS)
{
    int idx[N + 1];
    double h, d1c[MAXN - 1], d2c[MAXN - 2], res;
    h = (b - a) / N;
    sample(n, c, a, h, idx);
    derivative(n, c, d1c);
    derivative(n - 1, d1c, d2c);
    for (int i = 0; i <= N; i++) {
        res = a + idx[i] * h;
        if (trial(n, c, d1c, d2c, &res, EPS) == 0) {
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

case 4:
3 1 4 0 -10
1 2

case 5:
10 1 -10 45 -120 210 -252 210 -120 45 -10 1
-1000000 1000000

case 6:
9 1 -9 36 -84 126 -126 84 -36 9 -1
-1000000 1000000

case 7:
8 1 -12 62 -180 321 -360 248 -96 16
-1000000 1.5

*/