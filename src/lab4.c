#include <stdio.h>
#include <math.h>

#define MAX_SIZE 10
#define bound pow(2, 127)
#define ZERO 1e-9 /* X is considered to be 0 if |X|<ZERO */

enum bool { false = 0, true = 1 };
#define bool enum bool

int Jacobi( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN );

int Gauss_Seidel( int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN );

int main()
{
    int n, MAXN, i, j, k;
    double a[MAX_SIZE][MAX_SIZE], b[MAX_SIZE], x[MAX_SIZE];
    double TOL;

    scanf("%d", &n);
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++)
            scanf("%lf", &a[i][j]);
        scanf("%lf", &b[i]);
    }
    scanf("%lf %d", &TOL, &MAXN);

    printf("Result of Jacobi method:\n");
    for ( i=0; i<n; i++ )
        x[i] = 0.0;
    k = Jacobi( n, a, b, x, TOL, MAXN );
    switch ( k ) {
        case -2:
            printf("No convergence.\n");
            break;
        case -1: 
            printf("Matrix has a zero column.  No unique solution exists.\n");
            break;
        case 0:
            printf("Maximum number of iterations exceeded.\n");
            break;
        default:
            printf("no_iteration = %d\n", k);
            for ( j=0; j<n; j++ )
                printf("%.8f\n", x[j]);
            break;
    }
    printf("Result of Gauss-Seidel method:\n");
    for ( i=0; i<n; i++ )
        x[i] = 0.0;
    k = Gauss_Seidel( n, a, b, x, TOL, MAXN );
    switch ( k ) {
        case -2:
            printf("No convergence.\n");
            break;
        case -1: 
            printf("Matrix has a zero column.  No unique solution exists.\n");
            break;
        case 0:
            printf("Maximum number of iterations exceeded.\n");
            break;
        default:
            printf("no_iteration = %d\n", k);
            for ( j=0; j<n; j++ )
                printf("%.8f\n", x[j]);
            break;
    }

    return 0;
}

/* Your function will be put here */
#include <string.h>

enum Type { J = 0, GS = 1 };
#define Type enum Type

double temp[MAX_SIZE];

bool Elementary_Transform(int n, double a[][MAX_SIZE], double b[])
{
    int i, j, idx;
    double val, max;
    for (i = 0; i < n; i++) {
        idx = i;
        max = fabs(a[i][i]);
        for (j = i + 1; j < n; j++) {
            if ((val = fabs(a[j][i])) > max) {
                max = val;
                idx = j;
            }
        }
        if (idx != i) {
            // this indicates that variable max is non-zero
            memcpy(temp, a[i], n * sizeof(double));
            memcpy(a[i], a[idx], n * sizeof(double));
            memcpy(a[idx], temp, n * sizeof(double));
            val = b[i];
            b[i] = b[idx];
            b[idx] = val;
        }
        else if (max == 0) {
            idx = -1;
            max = 0;
            for (j = 0; j < i; j++) {
                if ((val = fabs(a[j][i])) > max) {
                    max = val;
                    idx = j;
                }
            }
            if (idx != -1) {
                for (j = 0; j < n; j++) {
                    a[i][j] += a[idx][j];
                }
                b[i] += b[idx];
            }
            else {
                return false;
            }
        }
    }
    return true;
}

double Infinity_Norm(int n, double v[])
{
    double ans = -1, val;
    for (int i = 0; i < n; i++) {
        if ((val = fabs(v[i])) > ans) {
            ans = val;
        }
    }
    return ans;
}

void Difference(int n, double lhs[], double rhs[], double dst[])
{
    for (int i = 0; i < n; i++) {
        dst[i] = lhs[i] - rhs[i];
    }
}

int Procedure(int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN, Type t)
{
    if (!Elementary_Transform(n, a, b)) {
        return -1;
    }
    int i, j, ans = 0;
    double norm, sum, *p = (t == J ? temp : x);
    do {
        ans++;
        memcpy(temp, x, n * sizeof(double));
        for (i = 0; i < n; i++) {
            sum = 0;
            for (j = 0; j < i; j++) {
                sum += a[i][j] * p[j];
            }
            for (j = i + 1; j < n; j++) {
                sum += a[i][j] * p[j];
            }
            x[i] = (b[i] - sum) / a[i][i];
            if (fabs(x[i]) > bound) {
                return -2;
            }
        }
        Difference(n, x, temp, temp);
        if (Infinity_Norm(n, temp) < TOL) {
            return ans;
        }
    } while (ans < MAXN);
    return 0;
}

int Jacobi(int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN)
{
    return Procedure(n, a, b, x, TOL, MAXN, J);
}

int Gauss_Seidel(int n, double a[][MAX_SIZE], double b[], double x[], double TOL, int MAXN)
{
    return Procedure(n, a, b, x, TOL, MAXN, GS);
}