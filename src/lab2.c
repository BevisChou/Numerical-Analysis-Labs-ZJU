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
#define SAMPLING_DENSITY 10
#define MAX_TRIAL 1000

enum Sign { P = 1, N = -1, Z = 0 };

static inline enum Sign Get_Sign(long double num)
{
    if(num > ZERO)
        return P;
    else if(num < -ZERO)
        return N;
    else
        return Z;
}

static inline long double Calculate_Polynomial(double x, double c[], int n)
{
    long double ans = 0;
    for(int i = n; i >= 0; i--)
        ans = ans * x + c[i];
    return ans;
}

static inline void Calculate_Derivative_Coefficients(double c[], int n, double res[])
{
    for(int i = n; i > 0; i--)
        res[i - 1] = c[i] * i;
}

long double Polynomial_Root_Bisection(double c[], int n, double a, double b, double EPS);
long double Polynomial_Root_Newton(double c[], int n, double a, double b, double EPS, int recursive_level);

double Polynomial_Root(int n, double c[], double a, double b, double EPS) 
{
    enum Sign left_sign = Get_Sign(Calculate_Polynomial(a, c, n)), right_sign = Get_Sign(Calculate_Polynomial(b, c, n));
    if(left_sign == Z)
        return a;
    if(right_sign == Z)
        return b;
    if(left_sign * right_sign == N)
        return Polynomial_Root_Bisection(c, n, a, b, EPS);
    else
        return Polynomial_Root_Newton(c, n, a, b, EPS, 2);
}

long double Polynomial_Root_Bisection(double c[], int n, double a, double b, double EPS)
{
    long double left = a, right = b, ans, poly;
    enum Sign left_sign = Get_Sign(Calculate_Polynomial(a, c, n)), ans_sign;
    do {
        ans_sign = Get_Sign(poly = Calculate_Polynomial(ans = (left + right) / 2, c, n));
        if(ans_sign == left_sign)
            left = ans;
        else
            right = ans;
    } while(fabsl(poly) > ZERO);
    return Polynomial_Root_Newton(c, n, left, right, EPS, 1);
}

typedef struct Sample {
    long double x;
    long double poly;
} Sample;

int compare (const void * lhs, const void * rhs)
{
    const Sample *sample_l = (const Sample *) lhs, *sample_r = (const Sample *) rhs;
    return fabsl(sample_l->poly) > fabsl(sample_r->poly) ? 1 : -1;
}

long double Polynomial_Root_Newton(double c[], int n, double a, double b, double EPS, int recursive_level)
{
    if(recursive_level == 0) 
        return 0;

    Sample samples[SAMPLING_DENSITY];
    long double step = (b - a) / SAMPLING_DENSITY, x, poly, d_poly;
    double derivative_c[n - 1];
    Calculate_Derivative_Coefficients(c, n, derivative_c);

    // Sampling and sorting.
    for(int i = 0; i < SAMPLING_DENSITY; i++)
        samples[i].poly = Calculate_Polynomial(samples[i].x = a + i * step, c, n);
    qsort(samples, SAMPLING_DENSITY, sizeof(Sample), compare);
    
    // Try each sample point.
    for(int i = 0; i < SAMPLING_DENSITY; i++)
    {
        int count = 0;
        x = samples[i].x;
        poly = samples[i].poly;
        d_poly = Calculate_Polynomial(x, derivative_c, n - 1);
        while(fabsl(poly) > EPS && count < MAX_TRIAL)
        {
            x = x - poly / d_poly;
            poly = Calculate_Polynomial(x, c, n);
            d_poly = Calculate_Polynomial(x, derivative_c, n - 1);
            count++;
        }
        if(count < MAX_TRIAL)
            return x;
    }

    // Try each interval.
    Sample candidates[SAMPLING_DENSITY];
    for(int i = 0; i < SAMPLING_DENSITY; i++)
        candidates[i].poly = Calculate_Polynomial(
            candidates[i].x = Polynomial_Root_Newton(
                c,
                n, 
                samples[i].x, 
                samples[i].x + step,
                EPS,
                recursive_level - 1), 
            c, 
            n);
    qsort(candidates, SAMPLING_DENSITY, sizeof(Sample), compare);
    return candidates[0].x;
}
