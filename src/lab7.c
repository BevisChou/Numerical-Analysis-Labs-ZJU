#include <stdio.h>
#include <math.h>

#define MAX_m 200
#define MAX_n 5

double f1(double x)
{
    return sin(x);
}

double f2(double x)
{
    return exp(x);
}

int OPA( double (*f)(double t), int m, double x[], double w[], double c[], double *eps );

void print_results( int n, double c[], double eps)
{    
    int i;

    printf("%d\n", n);
    for (i=0; i<=n; i++)
        printf("%12.4e ", c[i]);
    printf("\n");
    printf("error = %9.2e\n", eps);
    printf("\n");
}

int main()
{
    int m, i, n;
    double x[MAX_m], w[MAX_m], c[MAX_n+1], eps;

    m = 90;
    for (i=0; i<m; i++) {
        x[i] = 3.1415926535897932 * (double)(i+1) / 180.0;
        w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f1, m, x, w, c, &eps);
    print_results(n, c, eps);

    m = 200;
    for (i=0; i<m; i++) {
        x[i] = 0.01*(double)i;
        w[i] = 1.0;
    }
    eps = 0.001;
    n = OPA(f2, m, x, w, c, &eps);
    print_results(n, c, eps);

    return 0;
}

/* Your function will be put here */
#include <string.h>

double buf_m[MAX_m], buf_n[MAX_n + 1], cache[MAX_n + 1], cache_phi[MAX_n + 1][MAX_m];

double sample(double c[], int degree, double x)
{
    double res = 0;
    for (int i = degree; i >= 0; i--) {
        res = res * x + c[i];
    }
    return res;
}

void poly_sum(double lhs[], double rhs[], double dst[])
{
    // lhs might be the same as dst
    for (int i = 0; i <= MAX_n; i++) {
        dst[i] = lhs[i] + rhs[i];
    }
}

void poly_scale(double p[], double s, double dst[])
{
    // p might be the same as dst
    for (int i = 0; i <= MAX_n; i++) {
        dst[i] = s * p[i];
    }
}

void poly_elevate(double p[], double dst[])
{
    // p might be the same as dst
    for (int i = MAX_n; i > 0; i--) {
        dst[i] = p[i - 1];
    }
    dst[0] = 0;
}

double calc_inner_prod(int m, double w[], double x[], double y[])
{
    double res = 0, buf[MAX_m];
    for (int i = 0; i < m; i++) {
        buf[i] = w[i] * x[i] * y[i];
    }
    // TODO: sort buf?
    for (int i = 0; i < m; i++) {
        res += buf[i];
    }
    return res;
}

double calc_err(double (*f)(double t), double c[], int degree, int m, double w[], double x[])
{
    for (int i = 0; i < m; i++) {
        // we already have buf_m[i] = f(x[i])
        buf_m[i] -= sample(c, degree, x[i]);
    }
    return calc_inner_prod(m, w, buf_m, buf_m);
}

double calc_phi_coef(double (*f)(double t), double phi[], int degree, int m, double w[], double x[])
{
    for (int i = 0; i < m; i++) {
        cache_phi[degree][i] = sample(phi, degree, x[i]);
        buf_m[i] = f(x[i]);
    }
    double res = calc_inner_prod(m, w, cache_phi[degree], buf_m);
    cache[degree] = calc_inner_prod(m, w, cache_phi[degree], cache_phi[degree]);
    res /= cache[degree];
    return res;
}

void calc_phi(double phis[][MAX_n + 1], int degree, int m, double w[], double x[])
{
    memset(phis[degree], 0, sizeof phis[degree]);
    if (degree == 0) {
        phis[0][0] = 1;
        return;
    }
    // calculate bk
    for (int i = 0; i < m; i++) {
        buf_m[i] = x[i] * cache_phi[degree - 1][i];
    }
    double bk = calc_inner_prod(m, w, buf_m, cache_phi[degree - 1]) / cache[degree - 1];
    poly_elevate(phis[degree - 1], phis[degree]);
    poly_scale(phis[degree - 1], -bk, buf_n);
    poly_sum(phis[degree], buf_n, phis[degree]);
    if (degree == 1) {
        return;
    }
    // calculate ck
    double ck = calc_inner_prod(m, w, buf_m, cache_phi[degree - 2]) / cache[degree - 2];
    poly_scale(phis[degree - 2], -ck, buf_n);
    poly_sum(phis[degree], buf_n, phis[degree]);
}

int OPA( double (*f)(double t), int m, double x[], double w[], double c[], double *eps )
{
    memset(c, 0, MAX_n * sizeof(double));
    double err, phis[MAX_n + 1][MAX_n + 1];
    int d = -1;
    do {
        d++;
        calc_phi(phis, d, m, w, x);
        poly_scale(phis[d], calc_phi_coef(f, phis[d], d, m, w, x), buf_n);
        poly_sum(c, buf_n, c);
        err = calc_err(f, c, d, m, w, x);
    } while (err > *eps && d < MAX_n);
    *eps = err;
    return d;
}