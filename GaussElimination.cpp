// GaussElimination.cpp : This file contains the 'main' function. Program execution begins and ends there.

#include <stdlib.h>
#include <vector>
#include <time.h>
#include <chrono>
#include <sstream>
#include <stdio.h>
#include<iostream>


using namespace std;
//#define N 10000


int N = 1000;
// Random number generator

float Random(int low, int high)
{

    float num = ((float)rand() / RAND_MAX) * (float)(18.0);
    num = num - 9;
    return num;
}

// Swapping row of matrix
inline void SwapRow(float** mat, int i, int j)
{
   /* for (int k = i; k <= N; k++)
    {
        float temp = mat[i][k];
        mat[i][k] = mat[j][k];
        mat[j][k] = temp;
    }*/
  
      
}

// Perform forward elimination
int ForwardElimination(float** mat)
{
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point end;
    string sComputationTime;
    double computationTime = 0.0;
    FILE* complete = fopen("C:\\temp\\PracticalForwardElimination.txt", "a");
    start = chrono::high_resolution_clock::now();
    for (int k = 0; k < N; k++)
    {
        // Initialize maximum value and index for pivot
        int max = k;
        int max_value = mat[max][k];

        /* Find greater pivot value */
        for (int i = k + 1; i < N; i++)
        {
            if (abs(mat[i][k]) > max_value)
            {
                max_value = mat[i][k];
                max = i;
            }
        }

        // Main diagonal element is zero, it is singula. Nothing to do return. 
        if (!mat[k][k])
            return k;

        // Swap row with max value row.    
        if (max != k)
        {
            float* temp = mat[k];
            mat[k] = mat[max];
            mat[max] = temp;

            //SwapRow(mat, k, max);
        }
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
    end = chrono::high_resolution_clock::now();
    computationTime = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    sComputationTime = to_string(computationTime);

    fprintf(complete, "%ld, %s\n", N, sComputationTime.c_str());
    fclose(complete);
    return -1;
}

void BackSubsitution(float** mat)
{
    float x[10000]; // An array to store solution
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point end;
    string sComputationTime;
    double computationTime = 0.0;
    FILE* complete = fopen("C:\\temp\\PracticalBackWardElimination.txt", "a");
    /* Start calculating from last equation up to the
    first */
    start = chrono::high_resolution_clock::now();
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
    end = chrono::high_resolution_clock::now();
    computationTime = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    sComputationTime = to_string(computationTime);

    fprintf(complete, "%ld, %s\n", N, sComputationTime.c_str());
    fclose(complete);
}

void GaussianElimination(float** arr)
{
    //this is trial
    // Create RREF Form
    int Singular = ForwardElimination(arr);

    // Matrix is singular
    if (Singular != -1)
    {
        return;
    }
    BackSubsitution(arr);

}

int main()
{



    long long unsigned int LIMIT = 100000000;
    long long unsigned int j = 0;
    long long unsigned int temp = 0;
    double loopTime = 0.0;
    double divisionTime = 0.0;
    double multiplicationTime = 0.0;
    double additionTime = 0.0;
    chrono::high_resolution_clock::time_point endC;
    
    // Measure Loop Time
    chrono::high_resolution_clock::time_point startC = chrono::high_resolution_clock::now();
    while (j < LIMIT) 
    {
        j++;
    }
    endC = chrono::high_resolution_clock::now();
    loopTime = chrono::duration_cast<chrono::nanoseconds>(endC - startC).count();

    j = 0;
    temp = 1;
    // Measure Multiplication Time
    startC = chrono::high_resolution_clock::now();
    while (j < LIMIT) 
    {
        j++;
        temp = temp * 2;
    }
    endC = chrono::high_resolution_clock::now();
    multiplicationTime = chrono::duration_cast<chrono::nanoseconds>(endC - startC).count() - loopTime;
    multiplicationTime = multiplicationTime / (double)LIMIT;

    j = 0;
    temp = LIMIT;
    // Measure Division Time
    startC = chrono::high_resolution_clock::now();
    while (j < LIMIT) 
    {
        j++;
        temp = temp / 2;
    }
    endC = chrono::high_resolution_clock::now();
    divisionTime = chrono::duration_cast<chrono::nanoseconds>(endC - startC).count() - loopTime;
    divisionTime = divisionTime / (double)LIMIT;

    temp = 1;
    j = 0;
    // Addition Time 
    startC = chrono::high_resolution_clock::now();
    while (j < LIMIT) {
        j++;
        temp = temp + 2;
    }
    endC = chrono::high_resolution_clock::now();
    additionTime = chrono::duration_cast<chrono::nanoseconds>(endC - startC).count() - loopTime;
    additionTime = additionTime / (double)LIMIT;

    FILE* compute = fopen("C:\\temp\\ComputationTime.txt", "a");
    fprintf(compute, "\nLoopTime = %f, MultiplicationTime = %f, DivisionTime = %f, AdditionTime = %f", loopTime, multiplicationTime, divisionTime, additionTime);
    fclose(compute);
    cout << multiplicationTime;
    FILE* operationCount = fopen("C:\\temp\\totalOperationTime.txt", "a");
    double totalDivisionStime;
    double totalMultiplicationTime;
    double totalAdditionTime;
    double totalTime;
    for (int i = 1000; i <= 10000; i = i + 1000)
    {
        totalDivisionStime = ((((double)i * (double)(i - 1)) / 2) * divisionTime)/ 1000000;
        totalMultiplicationTime = ((((double)i * (double)(i - 1) * (double)(2 * i - 1)) / 6) * multiplicationTime)/ 1000000;
        totalAdditionTime = ((((double)i * (double)(i - 1) * (double)(2 * i - 1)) / 6) * additionTime)/ 1000000;
        totalTime = totalDivisionStime + totalMultiplicationTime + totalAdditionTime;
        fprintf(operationCount, "\nN = %d, totalTime = %lf,  MultiplicationTime = %lf, DivisionTime = %lf, AdditionTime = %lf", i,totalTime, totalMultiplicationTime, totalDivisionStime, totalAdditionTime);
        fflush(operationCount);
    }
    fclose(operationCount);
    double computationTime = 0.0;
    double computationTimeS = 0.0;
    string forwardComputation;
    string backSubstitution;
    chrono::high_resolution_clock::time_point start;
    chrono::high_resolution_clock::time_point end;
    std::chrono::milliseconds gaussComputationTime;
    string sComputationTime;
    j = 0;
    temp = 1;
    startC = chrono::high_resolution_clock::now();
    while (j < LIMIT) {
        j++;
        temp = temp * 2;
    }
    endC = chrono::high_resolution_clock::now();
    FILE* complete = fopen("C:\\temp\\Gauss.txt", "a");
    FILE* theoriticalForward = fopen("C:\\temp\\TheoryForward.txt", "a");
    FILE* theoriticalBackward = fopen("C:\\temp\\TheoryBackward.txt", "a");
    computationTimeS = chrono::duration_cast<chrono::nanoseconds>(endC - startC).count();
    computationTimeS = computationTimeS - loopTime;
    computationTimeS = computationTimeS / LIMIT;
    //computationTimeS /= 3;
    for (N = 1000; N <= 10000; N = N + 1000)
    {
        int i = 0;
        forwardComputation = to_string((((double)2 / 3) * pow(N, 3) * computationTimeS) / 1000000);
        backSubstitution = to_string(((double)pow(N, 2) * computationTimeS) / 1000000);
        fprintf(theoriticalForward, " %d, %s\n", N, forwardComputation.c_str());
        fprintf(theoriticalBackward, "%d, %s\n", N, backSubstitution.c_str());

        // Allocate memory for random array
        float** arr = (float**)malloc(N * sizeof(float*));
        for (i = 0; i < N; i++)
            arr[i] = (float*)malloc((N + 1) * sizeof(float));


        for (int i = 0; i < N; i++)
        {
            for (int j = 0; j <= N; j++)
            {
                arr[i][j] = (float)Random(0, 18);
            }
        }

        start = chrono::high_resolution_clock::now();
        GaussianElimination(arr);
        //this is trial
  
        end = chrono::high_resolution_clock::now();
        computationTime = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        sComputationTime = to_string(computationTime);

        fprintf(complete, "%ld, %s\n", N, sComputationTime.c_str());


        for (int i = 0; i < N; ++i)
        {
            free(arr[i]);
        }
        free(arr);
        fflush(theoriticalBackward);
        fflush(theoriticalForward);
        fflush(complete);
    }
    fclose(theoriticalBackward);
    fclose(theoriticalForward);
    fclose(complete);
    return 0;
}

