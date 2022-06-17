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

void MatCopy(int n, double src[][MAX_SIZE], double dst[][MAX_SIZE])
{
    for (int i = 0; i < n; i++) {
        memcpy(dst[i], src[i], n * sizeof(double));
    }
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

int MaxIdx(int n, double v[])
{
    int res = -1;
    double max = 0, val;
    for (int i = 0; i < n; i++) {
        if ((val = fabs(v[i])) > max) {
            max = val;
            res = i;
        }
    }
    return res;
}

int Normalize(int n, double v[])
{
    int idx = MaxIdx(n, v);
    double val = v[idx];
    for (int i = 0; i < n; i++) {
        v[i] /= val;
    }
    return idx;
}


int Solve(int n, double m[][MAX_SIZE], double b[], double res[])
{
    int i, j, k;
    double a[MAX_SIZE][MAX_SIZE], y[MAX_SIZE], coef;
    MatCopy(n, m, a);
    memcpy(y, b, n * sizeof(double));
    for (i = 0; i < n; i++) {
        if (a[i][i] == 0) {
            for (j = i + 1; j < n && a[j][i] == 0; j++);
            if (j == n) {
                return 1;
            }
            SwapVec(n, a[i], a[j]);
            SwapDouble(&y[i], &y[j]);
        }
        for (j = i + 1; j < n; j++) {
            coef = a[j][i] / a[i][i];
            for (k = i; k < n; k++) {
                a[j][k] -= coef * a[i][k];
            }
            y[j] -= coef * y[i];
        }
    }
    for (i = n - 1; i >= 0; i--) {
        for (j = i + 1; j < n; j++) {
            y[i] -= a[i][j] * res[j];
        }
        res[i] = y[i] / a[i][i];
    }
    return 0;
}

void Transpose(int n, double a[][MAX_SIZE])
{
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
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
    Transpose(n, res);
    return 0;
}

void Multiply(int n, double a[][MAX_SIZE], double v[], double res[])
{
    for (int i = 0; i < n; i++) {
        res[i] = 0;
        for (int j = 0; j < n; j++) {
            res[i] += a[i][j] * v[j];
        }
    }
}

void Difference(int n, double lhs[], double rhs[], double res[])
{
    for (int i = 0; i < n; i++) {
        res[i] = lhs[i] - rhs[i];
    }
}

int EigenV(int n, double a[][MAX_SIZE], double *lambda, double v[], double TOL, int MAXN)
{
    int i, idx;
    double p = *lambda, m[MAX_SIZE][MAX_SIZE], inv[MAX_SIZE][MAX_SIZE], temp[MAX_SIZE];
    MatCopy(n, a, m);
    for (i = 0; i < n; i++) {
        m[i][i] -= p;
    }
    if (Inverse(n, m, inv) != 0) {
        return -1;
    }
    idx = Normalize(n, v);
    for (i = 0; i < MAXN; i++) {
        memcpy(temp, v, n * sizeof(double));
        Multiply(n, inv, temp, v);
        *lambda = 1 / v[idx] + p;
        idx = Normalize(n, v);
        Difference(n, v, temp, temp);
        if (fabs(temp[MaxIdx(n, temp)]) < TOL) {
            return 1;
        }
    }
    return 0;
}