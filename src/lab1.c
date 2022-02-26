#include <stdio.h>

void Series_Sum( double sum[] );

int main()
{
    int i;
    double x, sum[3001];
    
    Series_Sum( sum );

    x = 0.0;
    for (i=0; i<3001; i++)
        printf("%6.2f %16.12f\n", x + (double)i * 0.10, sum[i]);

    return 0;
}

/* Your function will be put here */

#include <string.h>

void Series_Sum( double sum[] )
{
    long double res[3001];
    memset(res, 0, sizeof(res));
    int i, j, k = 6950;
    unsigned long long factor = (unsigned long long)k * (k + 1) * (k + 2) * (k + 3) * (k + 4);
    for (i = k; i > 0; i--)
    {
        for(j = 0; j <= 3000; j++)
            res[j] += (long double)1 / (i + 0.1 * j) / factor;
        factor = factor / (i + 4) * (i - 1);
    }
    for(int j = 0; j <= 3000; j++)
        sum[j] = (((res[j] * (4 - 0.1 * j) + 1.0 / 96) * (3 - 0.1 * j) + 1.0 / 18) * (2 - 0.1 * j) + 1.0 / 4) * (1 - 0.1 * j) + 1;
}