#include <stdio.h>

#define Max_size 10000 /* max number of dishes */

void Price( int n, double p[] );

int main()
{
    int n, i;
    double p[Max_size];

    scanf("%d", &n);
    for (i=0; i<n; i++) 
        scanf("%lf", &p[i]);
    Price(n, p);
    for (i=0; i<n; i++)
        printf("%.2f ", p[i]);
    printf("\n");

    return 0;
}

/* Your function will be put here */
#include <string.h>

#define Max_iteration 10

void Solve_Tridiagonal(double A[][3], int n, double y[], double x[])
{
    double alpha[Max_size], gamma[Max_size], beta[Max_size], z[Max_size];
    // L_{i,i} = alpha[i], L_{i+1,i} = gamma[i+1] (i = 0, ..., n-1)
    // U_{i,i+1} = beta[i] (i = 0, ..., n-2)

    alpha[0] = A[0][1];
    beta[0] = A[0][2] / alpha[0];
    for (int i = 1; i < n - 1; i++) {
        gamma[i] = A[i][0];
        alpha[i] = A[i][1] - gamma[i] * beta[i - 1];
        beta[i] = A[i][2] / alpha[i];
    }
    gamma[n - 1] = A[n - 1][0];
    alpha[n - 1] = A[n - 1][1] - gamma[n - 1] * beta[n - 2];

    z[0] = y[0] / alpha[0];
    for (int i = 1; i < n; i++) {
        z[i] = (y[i] - gamma[i] * z[i - 1]) / alpha[i];
    }
    x[n - 1] = z[n - 1];
    for (int i = n - 2; i >= 0; i--) {
        x[i] = z[i] - beta[i] * x[i + 1];
    }
}

void Sum(double lhs[], double rhs[], int n, double dst[])
{
    for (int i = 0; i < n; i++) {
        dst[i] = lhs[i] + rhs[i];
    }
}

void Price( int n, double p[] )
{
    double y[Max_size], x[Max_size];
    memcpy(y, p, n * sizeof(double));
    memset(p, 0, n * sizeof(double));
    if (n == 2) {
        p[0] = (2 * y[0] - y[1]) / 3;
        p[1] = (2 * y[1] - y[0]) / 3;
        return;
    }

    double A[Max_size][3];
    for (int i = 1; i < n - 1; i++) {
        A[i][0] = 0.5;
        A[i][1] = 2;
        A[i][2] = 0.5;
    }
    A[0][0] = A[n - 1][2] = 0;
    A[0][1] = A[n - 1][1] = 2;
    A[0][2] = A[n - 1][0] = 0.5;

    for (int i = 0; i < Max_iteration; i++) {
        Solve_Tridiagonal(A, n, y, x);
        Sum(p, x, n, p);
        memset(y, 0, sizeof y);
        y[0] = -0.5 * x[n - 1];
        y[n - 1] = -0.5 * x[0];
    }
}