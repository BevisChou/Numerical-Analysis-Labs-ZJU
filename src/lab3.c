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

void Price( int n, double p[] )
{
    double price[Max_size];
    memcpy(price, p, n * sizeof(double));
    if (n == 2) {
        p[0] = (2 * price[0] - price[1]) / 3;
        p[1] = (2 * price[1] - price[0]) / 3;
        return;
    }

    double temp[Max_size], A[Max_size][3], coef;
    // temp could be of size 3, but would induce extra complexity
    for (int i = 1; i < n - 1; i++) {
        A[i][0] = 0.5;
        A[i][1] = 2;
        A[i][2] = 0.5;
    }

    memset(temp, 0, sizeof temp);
    A[0][0] = 0;
    temp[0] = 2;
    temp[1] = temp[n - 1]= 0.5;
    for (int i = n - 2; i > 0; i--) {
        temp[i] -= A[i][1];
        temp[i - 1] -= A[i][0];
        price[0] -= price[i];
        coef = 2 * temp[i];
        temp[i] = 0.5;
        temp[i - 1] /= coef;
        if (i > 2) {
            temp[1] /= coef;
        }
        if (i > 1) {
            temp[0] /= coef;
        }
        price[0] /= coef;
    }
    A[0][1] = temp[0];
    A[0][2] = temp[1];

    memset(temp, 0, sizeof temp);
    temp[0] = temp[n - 2] = 0.5;
    temp[n - 1] = 2;
    A[n - 1][2] = 0;
    for (int i = 1; i < n - 1; i++) {
        temp[i] -= A[i][1];
        temp[i + 1] -= A[i][2];
        price[n - 1] -= price[i];
        coef = 2 * temp[i];
        temp[i] = 0.5;
        temp[i + 1] /= coef;
        if (i < n - 3) {
            temp[n - 2] /= coef;
        }
        if (i < n - 2) {
            temp[n - 1] /= coef;
        }
        price[n - 1] /= coef;
    }
    A[n - 1][0] = temp[n - 2];
    A[n - 1][1] = temp[n - 1];

    Solve_Tridiagonal(A, n, price, p);
}