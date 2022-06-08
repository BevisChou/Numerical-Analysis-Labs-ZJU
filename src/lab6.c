#include <stdio.h>

#define MAX_N 10

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[]);

double S( double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[] );

int main()
{
    int n, Type, m, i;
    double x[MAX_N], f[MAX_N], a[MAX_N], b[MAX_N], c[MAX_N], d[MAX_N];
    double s0, sn, Fmax, t0, tm, h, t;

    scanf("%d", &n);
    for (i=0; i<=n; i++) 
        scanf("%lf", &x[i]);
    for (i=0; i<=n; i++) 
        scanf("%lf", &f[i]);
    scanf("%d %lf %lf %lf", &Type, &s0, &sn, &Fmax);

    Cubic_Spline(n, x, f, Type, s0, sn, a, b, c, d);
    for (i=1; i<=n; i++)
        printf("%12.8e %12.8e %12.8e %12.8e \n", a[i], b[i], c[i], d[i]);

    scanf("%lf %lf %d", &t0, &tm, &m);
    h = (tm-t0)/(double)m;
    for (i=0; i<=m; i++) {
        t = t0+h*(double)i;
        printf("f(%12.8e) = %12.8e\n", t, S(t, Fmax, n, x, a, b, c, d));
    }

    return 0;
}

/* Your functions will be put here */
void Solve_Tridiagonal(int n, double A[][MAX_N], double y[], double x[])
{
    double L[MAX_N][MAX_N], z[MAX_N], U[MAX_N][MAX_N];
    L[0][0] = A[0][0];
    U[0][1] = A[0][1] / L[0][0];
    z[0] = y[0] / L[0][0];
    for (int i = 1; i < n; i++) {
        L[i][i - 1] = A[i][i - 1];
        L[i][i] = A[i][i] - L[i][i - 1] * U[i - 1][i];
        U[i][i + 1] = A[i][i + 1] / L[i][i];
        z[i] = (y[i] - L[i][i - 1] * z[i - 1]) / L[i][i];
    }
    L[n][n - 1] = A[n][n - 1];
    L[n][n] = A[n][n] - L[n][n - 1] * U[n - 1][n];
    z[n] = (y[n] - L[n][n - 1] * z[n - 1]) / L[n][n];
    x[n] = z[n];
    for (int i = n - 1; i >= 0; i--) {
        x[i] = z[i] - U[i][i + 1] * x[i + 1];
    }
}

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[])
{
    double h[MAX_N], A[MAX_N][MAX_N], y[MAX_N];
    for (int i = 0; i <= n; i++) {
        a[i] = f[i];
    }
    for (int i = 0; i < n; i++) {
        h[i] = x[i + 1] - x[i];
    }
    switch (Type) {
        case 1:
            A[0][0] = 2 * h[0];
            A[0][1] = h[0];
            y[0] = 3 * ((a[1] - a[0]) / h[0] - s0);
            A[n][n - 1] = h[n - 1];
            A[n][n] = 2 * h[n - 1];
            y[n] = 3 * (sn - (a[n] - a[n - 1]) / h[n - 1]);
            break;
        case 2:
            A[0][0] = 2;
            A[0][1] = 0;
            y[0] = s0;
            A[n][n - 1] = 0;
            A[n][n] = 2;
            y[n] = sn;
            break;
    }
    for (int i = 1; i < n; i++) {
        A[i][i - 1] = h[i - 1];
        A[i][i] = 2 * (h[i - 1] + h[i]);
        A[i][i + 1] = h[i];
        y[i] = 3 * ((a[i + 1] - a[i]) / h[i] - (a[i] - a[i - 1]) / h[i - 1]);
    }
    Solve_Tridiagonal(n, A, y, c);
    for (int i = 0; i < n; i++) {
        b[i] = (a[i + 1] - a[i]) / h[i] - (c[i + 1] + 2 * c[i]) * h[i] / 3;
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);  
    }
    for (int i = n; i > 0; i--) {
        a[i] = a[i - 1];
        b[i] = b[i - 1];
        c[i] = c[i - 1];
        d[i] = d[i - 1];
    } 
}

double S( double t, double Fmax, int n, double x[], double a[], double b[], double c[], double d[] )
{
    if (t < x[0] || t > x[n]) {
        return Fmax;
    }
    int i;
    for (i = 1; t > x[i]; i++);
    double h = t - x[i - 1];
    return ((d[i] * h + c[i]) * h + b[i]) * h + a[i];
}