#include <stdio.h>

#define MAX_SIZE 10

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN);

int main()
{
    int n, MAXN, m, i, j, k;
    double a[MAX_SIZE][MAX_SIZE], v[MAX_SIZE];
    double lambda, TOL;

    scanf("%d", &n);
    for (i=0; i<n; i++) 
        for (j=0; j<n; j++) 
            scanf("%lf", &a[i][j]);
    scanf("%lf %d", &TOL, &MAXN);
    scanf("%d", &m);
    for (i=0; i<m; i++) {
        scanf("%lf", &lambda);
        for (j=0; j<n; j++)
            scanf("%lf", &v[j]);
        switch (EigenV(n, a, &lambda, v, TOL, MAXN)) {
            case -1: 
                printf("%12.8f is an eigenvalue.\n", lambda );
                break;
            case 0:
                printf("Maximum number of iterations exceeded.\n");
                break;
            case 1:
                printf("%12.8f\n", lambda );
                for (k=0; k<n; k++)
                    printf("%12.8f ", v[k]);
                printf("\n");
                break;
        }
    }

    return 0;
}

/* Your function will be put here */
#include <math.h>
#include <string.h>

void PrintMatrix(int n, double a[][MAX_SIZE])
{
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
}

void PrintVector(int n, double v[])
{
    for (int i = 0; i < n; i++) {
        printf("%lf ", v[i]);
    }
    printf("\n");
}

void SwapDouble(double *a, double *b)
{
    double temp = *a;
    *a = *b;
    *b = temp;
}

void SwapVec(int n, double a[], double b[])
{
    double temp[MAX_SIZE];
    memcpy(temp, a, n * sizeof(double));
    memcpy(a, b, n * sizeof(double));
    memcpy(b, temp, n * sizeof(double));
}

int Solve(int n, double a[][MAX_SIZE], double y[], double res[])
{
    int i, j, k;
    double copy[MAX_SIZE][MAX_SIZE], coef;
    for (i = 0; i < n; i++) {
        memcpy(copy[i], a[i], n * sizeof(double));
    }
    for (i = 0; i < n; i++) {
        if (copy[i][i] == 0) {
            for (j = i + 1; j < n && copy[j][i] == 0; j++);
            if (j == n) {
                return 1;
            }
            SwapVec(n, copy[i], copy[j]);
            SwapDouble(&y[i], &y[j]);
        }
        for (j = i + 1; j < n; j++) {
            coef = copy[j][i] / copy[i][i];
            for (k = i; k < n; k++) {
                copy[j][k] -= coef * copy[i][k];
            }
            y[j] -= coef * y[i];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        for (j = i + 1; j < n; j++) {
            y[i] -= copy[i][j] * res[j];
        }
        res[i] = y[i] / copy[i][i];
    }
    return 0;
}

void Transpose(int n, double a[][MAX_SIZE], double res[][MAX_SIZE])
{
    // variable res might point to the same location as variable a does
    int i, j;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            SwapDouble(&a[i][j], &a[j][i]);
        }
    }
}

int Inverse(int n, double a[][MAX_SIZE], double res[][MAX_SIZE])
{
    double y[MAX_SIZE];
    for (int i = 0; i < n; i++) {
        memset(y, 0, sizeof y);
        y[i] = 1;
        if(Solve(n, a, y, res[i]) != 0) {
            return 1;
        }
    }
    Transpose(n, res, res);
    // PrintMatrix(n, res);
    return 0;
}

void Multiply(int n, double a[][MAX_SIZE], double v[], double res[])
{
    int i, j;
    for (i = 0; i < n; i++) {
        res[i] = 0;
        for (j = 0; j < n; j++) {
            res[i] += a[i][j] * v[j];
        }
    }
}

int Normalize(int n, double v[])
{
    int i, res = -1;
    double max = 0, val;
    for (i = 0; i < n; i++) {
        if ((val = fabs(v[i])) > max) {
            max = val;
            res = i;
        }
    }
    val = v[res];
    for (i = 0; i < n; i++) {
        v[i] /= val;
    }
    return res;
}

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN)
{
    int i, idx;
    double p = *lambda, l, b[MAX_SIZE][MAX_SIZE], last[MAX_SIZE];
    for (i = 0; i < n; i++) {
        a[i][i] -= p;
    }
    if (Inverse(n, a, b) != 0) {
        return -1;
    }
    idx = Normalize(n, v);
    for (int i = 0; i < MAXN; i++) {
        memcpy(last, v, n * sizeof(double));
        Multiply(n, b, last, v);
        l = *lambda;
        *lambda = 1 / v[idx] + p;
        idx = Normalize(n, v);
        if (fabs(*lambda - l) < TOL) {
            return 1;
        }
    }
    return 0;
}