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

int newton(int n, double c[], double dc[], double x, double *res)
{
    double d0 = calc(n, c, x), d1 = calc(n - 1, dc, x);
    if (fabs(d1) < ZERO) {
        return 1;
    }
    *res = x - d0 / d1;
    return 0;
}

int trial(int n, double c[], double dc[], double *x, double EPS)
{
    double ps[3], d0, d1;
    ps[0] = *x;
    for (int i = 0; i < MAXI; i++) {
        if (fabs(calc(n, c, ps[0])) < ZERO) {
            *x = ps[0];
            return 0;
        }
        if (newton(n, c, dc, ps[0], &ps[1]) != 0 || newton(n, c, dc, ps[1], &ps[2]) != 0) {
            return 1;
        }
        d0 = ps[1] - ps[0];
        d1 = ps[2] - ps[1];
        if (fabs(d1) < EPS) {
            *x = ps[2];
            return 0;
        }
        ps[0] = ps[0] - d0 * d0 / (d1 - d0);
    }
    return 1;
}

void swap(double *lhs, double *rhs) {
    double temp = *lhs;
    *lhs = *rhs;
    *rhs = temp;
}

double Polynomial_Root(int n, double c[], double a, double b, double EPS)
{
    if (a > b) {
        swap(&a, &b);
    }
    int idx[N + 1];
    double h, dc[MAXN - 1], res;
    h = (b - a) / N;
    sample(n, c, a, h, idx);
    derivative(n, c, dc);
    for (int i = 0; i <= N; i++) {
        res = a + idx[i] * h;
        if (trial(n, c, dc, &res, EPS) == 0 && res >= a && res <= b) {
            return res;
        }
    }
    return 0;
}
