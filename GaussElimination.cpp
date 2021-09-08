// GaussElimination.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <iostream>
#include <vector>
#include <time.h>
using namespace std;
#define N 10000



// Random number generator
float Random(int low, int high)
{
    
    float num = ((float)rand() / RAND_MAX) * (float)(18.0);      
    num = num - 9;
    return num;
}

// Swapping row of matrix
void SwapRow(float** mat, int i, int j)
{   
    for (int k = 0; k <= N; k++)
    {
        float temp = mat[i][k];
        mat[i][k] = mat[j][k];
        mat[j][k] = temp;
    }
}

// Perform forward elimination
int ForwardElimination(float** mat)
{
    for (int k = 0; k < N; k++)
    {
        // Initialize maximum value and index for pivot
        int max = k;
        int max_value = mat[max][k];

        /* Find greater pivot value */
        for (int i = k + 1; i < N; i++)
        {  if (abs(mat[i][k]) > max_value)
			{   
				max_value = mat[i][k];
				max = i;
			}
		}
       
        // Main diagonal element is zero, it is singula. Nothing to do return. 
        if (!mat[k][max])
			return k; 

         // Swap row with max value row.    
        if (max != k)
            SwapRow(mat, k, max);

        for (int i = k + 1; i < N; i++)
        {
            /* factor f to set current row kth element to 0,
            * and subsequently remaining kth column to 0 */
            float f = mat[i][k] / mat[k][k];

            /* subtract fth multiple of corresponding kth
            row element*/
            for (int j = k + 1; j <= N; j++)
                mat[i][j] -= mat[k][j] * f;

            // filling lower triangular matrix with zeros
            mat[i][k] = 0;
        }
    }
    return -1;
}

void BackSubsitution(float** mat)
{
    float x[N]; // An array to store solution

    /* Start calculating from last equation up to the
    first */
    for (int i = N - 1; i >= 0; i--)
    {
        /* start with the RHS of the equation */
        x[i] = mat[i][N];

        /* Initialize j to i+1 since matrix is upper
        triangular*/
        for (int j = i + 1; j < N; j++)
        {
            /* subtract all the lhs values
            * except the coefficient of the variable
            * whose value is being calculated */
            x[i] -= mat[i][j] * x[j];
        }

        /* divide the RHS by the coefficient of the
        unknown being calculated */
        x[i] = x[i] / mat[i][i];
    }
  
}

void GaussianElimination(float** arr)
{
    // Create RREF Form
    int Singular = ForwardElimination(arr);

    // Matrix is singular
    if ( Singular != -1)
    {
        return;
    }
	BackSubsitution(arr);

}

int main()
{
	int i=0;
    float** arr = (float**)malloc(N * sizeof(float*));
    for (i = 0; i < N; i++)
        arr[i] = (float*)malloc((N+1) * sizeof(float));
  
    clock_t start = clock();    
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j <= N; j++)
        {
           arr[i][j] =(float)Random(0, 18);
        }  
    }
    clock_t end = clock();

	// Time taken to create random array
    double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    cout << cpu_time_used;
    start = clock();
    GaussianElimination(arr);
    end = clock();
    cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
    cout << "\n";
    cout << cpu_time_used;
    for (int i = 0; i < N; ++i)
        delete[] arr[i];
    delete[] arr;
    return 0;
}

