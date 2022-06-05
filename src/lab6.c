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
// #include <string.h>

#define TIMES2(num) (num + num)
#define TIMES3(num) (num + num + num)

void Print_Matrix(int size, double M[][MAX_N]);
void Print_Vector(int size, double V[]);

void Solve_Tridiagonal(int size, double A[][MAX_N], double y[], double x[])
{
    double L[MAX_N][MAX_N], z[MAX_N], U[MAX_N][MAX_N];
    // memset(L, 0, sizeof L); 
    // memset(U, 0, sizeof U); 
    L[0][0] = A[0][0];
    // printf("L[0][0]: %12.8e\n", L[0][0]);
    U[0][1] = A[0][1] / L[0][0];
    // printf("U[0][1]: %12.8e\n", U[0][1]);
    z[0] = y[0] / L[0][0];
    // printf("z[0]: %12.8e\n", z[0]);
    for (int i = 1; i < size - 1; i++) {
        L[i][i - 1] = A[i][i - 1];
        // printf("L[%d][%d]: %12.8e\n", i, i - 1, L[i][i - 1]);
        L[i][i] = A[i][i] - L[i][i - 1] * U[i - 1][i];
        // printf("L[%d][%d]: %12.8e\n", i, i, L[i][i]);
        U[i][i + 1] = A[i][i + 1] / L[i][i];
        // printf("U[%d][%d]: %12.8e\n", i, i + 1, U[i][i + 1]);
        z[i] = (y[i] - L[i][i - 1] * z[i - 1]) / L[i][i];
        // printf("z[%d]: %12.8e\n", i, z[i]);
    }
    L[size - 1][size - 2] = A[size - 1][size - 2];
    // printf("L[%d][%d]: %12.8e\n", size - 1, size - 2, L[size - 1][size - 2]);
    L[size - 1][size - 1] = A[size - 1][size - 1] - L[size - 1][size - 2] * U[size - 2][size - 1];
    // printf("L[%d][%d]: %12.8e\n", size - 1, size - 1, L[size - 1][size - 1]);
    z[size - 1] = (y[size - 1] - L[size - 1][size - 2] * z[size - 2]) / L[size - 1][size - 1];
    // printf("z[%d]: %12.8e\n", size - 1, z[size - 1]);
    x[size - 1] = z[size - 1];
    // printf("x[%d]: %12.8e\n", size - 1, x[size - 1]);
    for (int i = size - 2; i >= 0; i--) {
        x[i] = z[i] - U[i][i + 1] * x[i + 1];
        // printf("x[%d]: %12.8e\n", i, x[i]);
    }
}

void Cubic_Spline(int n, double x[], double f[], int Type, double s0, double sn, double a[], double b[], double c[], double d[])
{
    double h[MAX_N], da[MAX_N], dc[MAX_N], A[MAX_N][MAX_N], hy[MAX_N], y[MAX_N];
    // memset(A, 0, sizeof A); 
    for (int i = 0; i <= n; i++) {
        a[i] = f[i];
    }
    for (int i = 0; i < n; i++) {
        h[i] = x[i + 1] - x[i];
        da[i] = a[i + 1] - a[i];
        hy[i] = da[i] / h[i];
    }
    switch (Type) {
        case 1:
            A[0][0] = h[0] + h[0];
            A[0][1] = h[0];
            A[n][n - 1] = h[n - 1];
            // A[n][n] = 2 * h[n - 1];
            A[n][n] = TIMES2(h[n - 1]);
            y[0] = TIMES3(hy[0] - s0);
            y[n] = TIMES3(sn - hy[n - 1]);
            break;
        case 2:
            A[0][0] = 2;
            A[n][n] = 2;
            y[0] = s0;
            y[n] = sn;
            break;
    }
    for (int i = 1; i < n; i++) {
        A[i][i - 1] = h[i - 1];
        // A[i][i] = 2 * (h[i - 1] + h[i]);
        A[i][i] = TIMES2(h[i - 1] + h[i]);
        A[i][i + 1] = h[i];
        y[i] = TIMES3(hy[i] - hy[i - 1]);
    }
    // Print_Matrix(n + 1, A);
    // Print_Vector(n + 1, y);
    Solve_Tridiagonal(n + 1, A, y, c);
    // Print_Vector(n + 1, c);
    for (int i = 0; i < n; i++) {
        b[i] = hy[i] - (c[i + 1] + TIMES2(c[i])) * h[i] / 3;
        // printf("b[%d]: %12.8e\n", i, b[i]);
        d[i] = (c[i + 1] - c[i]) / (3 * h[i]);  
        // printf("d[%d]: %12.8e\n", i, d[i]);
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
    for (i = 1; t < x[i - 1] || t > x[i]; i++);
    double h = t - x[i - 1];
    return ((d[i] * h + c[i]) * h + b[i]) * h + a[i];
}

void Print_Matrix(int size, double M[][MAX_N]) {
    printf("---MATRIX %d * %d---\n", size, size);
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            printf("%12.8e ", M[i][j]);
        }
        printf("\n");
    }
}

void Print_Vector(int size, double V[]) {
    printf("---VECTOR %d---\n", size);
    for(int i = 0; i < size; i++) {
        printf("%12.8e ", V[i]);
    }
    printf("\n");
}

/*
case 1:
2
0.0 1.0 2.0
0.0 1.0 2.0
1 1.0 1.0 0.0
0.0 3.0 2

case 2:
2
1.0 2.0 3.0
2.0 3.0 5.0
1 2.0 1.0 0.0
0.0 3.0 2

case 3:
2
1.0 2.0 3.0
1.0 1.0 0.0
2 0.0 0.0 0.0
0.0 3.0 2

*/