#include "spkmeans.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#define FLT_MAX 3.402823466e+38F
#define FLT_MIN 1.1754943508e-38F
#define BUFFER_SIZE 10000
#define INITIAL_ARRAY_SIZE 1000
#define EPSILON 0.001

#define ASSERT_ARGS(expr) \
    if (!(expr)){             \
         printf("%s", "Invalid Input!"); \
         exit(0);\
    }

#define ASSERT_ERROR(expr) \
    if (!(expr)){             \
         printf("%s", "An Error Has Occured"); \
         exit(0);\
    }

/**
 * main, a shell function for the spectral clustering algorithm implementation
 */
int main(int argc, char *argv[]) {
//    checkArgs(argc, argv, 0);
    int i;
    double **arr, **Atag, **V, **jacobiMatrix, *jacobiArray;
    /*****/
//    jacobiMatrix = (double**) calloc(4,sizeof (double*));


    arr = (double**) calloc(3,sizeof (double*));
   // Atag = (double**) calloc(3,sizeof (double*));
   // V = (double**) calloc(3,sizeof (double*));
    for(i = 0; i < 3; i++) {
        arr[i] =  calloc(3, sizeof(double));
     //   Atag[i] =  calloc(3, sizeof(double));
       // V[i] =  calloc(3, sizeof(double));
        //V[i][i] = 1;
    }
 //   for(i = 0; i < 4; i++) {
 //       jacobiMatrix[i] = calloc(3, sizeof(double));
  //  }
    arr[0][0] = 3;
    arr[0][1] = 2;
    arr[0][2] = 4;
    arr[1][0] = 2;
    arr[1][1] = 0;
    arr[1][2] = 2;
    arr[2][0] = 4;
    arr[2][1] = 2;
    arr[2][2] = 3;

    printTest(arr, 3, 3);


    printf("\n");
    jacobiMatrix = calcJacobi(3, &arr);
    printf("\n jacobi Matrix: \n");
    printTest(jacobiMatrix, 4, 3);
 //   printf("\n Values: \n");
  //  printTest(Atag);
   // printf("\n Vectors: \n");
   // printTest(V);
    //printf("end");

    return 0;
}

/**
 * Checks the validity of arguments passed, reads and parses info in files
 * and starts the implementation of spectral k-means clustering
 * @param argc - # of arguments passed to main
 * @param origArgv - the arguments passed to main
 * @param isCAPI - used in order to distinguish between
 *                 C-based and C-API-based call to this function
 */
double **checkArgs(int argc, char **origArgv, int isCAPI, int *returnRowCount,
                   int *returnColCount) {
    int k, d, arraySize, *dArraySizeInfo;
    double *datap_array, **datapoint, **ret_matrix;
    enum goalEnum goal;
    char *ptr;

    ASSERT_ARGS(argc > 3)

    k = strtol(origArgv[1], &ptr, 10);

    /** checking whether k>=0,
     * strcmp with "0" used for the case which strtol fails and returns 0*/
    ASSERT_ARGS(k > 0 || !strcmp(origArgv[1], "0"))

    goal = checkGoal(origArgv[2]);
    ASSERT_ARGS(goal != INVALID)

    dArraySizeInfo = calloc(2, sizeof(int));

    /** process the datapoint file, and return info about their amount (arraySize)
     * and number of properties (d) */
    processDatapoints(origArgv[3], &datap_array, &datapoint, &dArraySizeInfo);

    d = dArraySizeInfo[0];
    arraySize = dArraySizeInfo[1];
    free(dArraySizeInfo);

    ASSERT_ARGS(k < arraySize)

    /** initialize the process to achieve the provided goal */
    ret_matrix = initProcess(&k, d, arraySize, goal, &datapoint, isCAPI);

    /** Determining row and col size for return matrix */
    switch (goal) {
        case jacobi:{
            *returnRowCount = arraySize + 1;
            *returnColCount = arraySize;
            break;
        }
        case spk:{
            *returnRowCount = arraySize;
            *returnColCount = k;
            break;
        }
        default:{
            *returnRowCount = arraySize;
            *returnColCount = arraySize;
            break;
        }
    }

    /** Free allocated memory and terminate*/
    free(datapoint);
    free(datap_array);

    return ret_matrix;
}

/**
 * Initializes datapoint-related arrays on a row-by-row basis
 * Mostly copied from HW1
 * @param filename - contains path to datapoints file
 * @param datap_array - pointer to 1D-support array containing datapoint values
 * @param datapoint - pointer to 2D array built upon datap_array
 * @param dArraySizeInfo - array used to return d and arraySize info back to
 *                         the checkArgs function
 */
void processDatapoints(char *filename, double **datap_array, double ***datapoint,
                       int **dArraySizeInfo) {
    char *token, *ptr, currRow[BUFFER_SIZE];
    FILE *fptr;
    int i, d, arraySize, tail = 0;

    /*opens file and reads first line to determine d*/
    fptr = fopen(filename, "r");
    ASSERT_ERROR(fptr != NULL)

    fgets(currRow, BUFFER_SIZE, fptr);
    ASSERT_ERROR(currRow != NULL)

    for (i = 0, d = 1; currRow[i]; i++) {
        d += (currRow[i] == ',');
    }

    arraySize = INITIAL_ARRAY_SIZE;

    *datap_array = calloc(d * arraySize, sizeof(double));
    ASSERT_ERROR(*datap_array != NULL)

    /** Populates datap_array with info from the file */
    do {/* as long as there are still lines to be read */
        /*loop to read line */
        token = strtok(currRow, ","); /* read the current line */
        while (token != NULL) {
            (*datap_array)[tail] = strtod(token, &ptr);
            token = strtok(NULL, ",");
            tail++;

            /*in-case we reached the edge of the current datap_array,
             * increase the arraySize by twice */
            if (tail == d * arraySize) {
                arraySize *= 2;
                *datap_array = realloc(*datap_array,
                                       d * arraySize * sizeof(double));
                ASSERT_ERROR(*datap_array != NULL)
            }
        }
    } while (fgets(currRow, BUFFER_SIZE, fptr));

    fclose(fptr);

    /* cut the datap_array to its intended arraySize */
    if (tail < d * arraySize - 1) {
        *datap_array = realloc(*datap_array, tail * sizeof(double));
        ASSERT_ERROR(*datap_array != NULL)
        arraySize = tail / d;
    }

    free(token);

    *datapoint = arrayToTwoDimArray(datap_array, arraySize, d);

    *dArraySizeInfo[0] = d;
    *dArraySizeInfo[1] = arraySize;

}

/**
 * Initializes datapoints and related arrays. Mostly copied from HW1.
 * Main function in C-API implementation.
 * @param k - pointer to amount of clusters the data needs to divide into.
 *            if 0, need to apply Eigengap Heuristic
 * @param d - dimension of datapoints
 * @param goal - goal of call to main (spk / wam / ddg / lnorm / jacobi)
 * @param datapoint - pointer to 2D-array containing all datapoints.
 *                    In case of goal = 'jacobi', it contains a symmetric
 *                    nXn matrix
 * @param init_centroids - pointer to 1D array containing the indices of
 *                         initial centroids provided by K-Means++ algorithm,
 *                         relevant only in C-API implementation
 * @param isCAPI - used in order to distinguish between
 *                 C-based and C-API-based call to this function
 */
double **
initProcess(int *k, int d, int arraySize, enum goalEnum goal,
            double ***datapoint, int isCAPI) {

    double **ret_matrix;

    ret_matrix = goalBasedProcess(k, d, arraySize, datapoint, goal, isCAPI);

    /** Freeing memory and returning/printing the result */

    if (!isCAPI) {
        printResult(&ret_matrix, goal);
        free(ret_matrix);
        return NULL;
    } else {
        /*Will be freed by C-API*/
        return ret_matrix;
    }
}

/**
 * Takes the datapoint info, and returns relevant 1D support array according to
 * the goal, later will be converted to matrix form according to needs.
 * @param k - number of clusters to assign. if 0, we need to use EigGap Heuristic
 * @param d - dimension of datapoints
 * @param arraySize - amount of datapoints
 * @param datapoint - pointer to 2D-array containing all datapoints
 * @param goal - the goal of the computation
 * @param init_centroids - pointer to 1D array containing the indices of
 *                         initial centroids provided by K-Means++ algorithm,
 *                         relevant only in C-API implementation
 * @return relevant matrix according to goal
 */
double **
goalBasedProcess(int *k, int d, int arraySize, double ***datapoint,
                 enum goalEnum goal, int isCAPI) {

    double **weightedAdjMatrix, **diagDegMatrix,
            **normLaplacianMatrix, **spkMatrix, **jacobiMatrix;

    /**
     * in case of goal == jacobi, we get a symmetric matrix that we need to
     * apply the Jacobi algorithm on
     */
    if (goal == jacobi){
        jacobiMatrix = calcJacobi(arraySize, datapoint);
        return jacobiMatrix;
    }

    weightedAdjMatrix = calcWeightedAdjMatrix(d, arraySize, datapoint);
    if (goal == wam) {
        return weightedAdjMatrix;
    }

    diagDegMatrix = calcDiagDegMatrix(arraySize, &weightedAdjMatrix);
    if (goal == ddg) {
        freeMatrix(&weightedAdjMatrix);
        return diagDegMatrix;
    }

    normLaplacianMatrix = calcNormLaplacian(arraySize, &weightedAdjMatrix,
                                            &diagDegMatrix);
    if (goal == lnorm) {
        freeMatrix(&weightedAdjMatrix);
        freeMatrix(&diagDegMatrix);
        return normLaplacianMatrix;
    }

    jacobiMatrix = calcJacobi(arraySize, &normLaplacianMatrix);

    /** TODO: if Python, we should return to implement K-Means++*/
    spkMatrix = calcSpectralClusters(k, arraySize, 0, &jacobiMatrix);

    freeMatrix(&weightedAdjMatrix);
    freeMatrix(&diagDegMatrix);
    freeMatrix(&normLaplacianMatrix);
    freeMatrix(&jacobiMatrix);
    return spkMatrix;
}

/**
 * Prints the result of the algorithm according to goal
 * wam / ddg / lnorm - print a nXn matrix
 * jacobi - first line is eigenvalues,
 *          next are nXn eigenvecor matrix
 * spk - first line is k-means++ init_centroid indices
 *       next are the k centroids after k-means algorithm
 * @param retArray - pointer to matrix to be printed to user
 * @param goal - distinguishes between spk-printing and matrix-printing
 */
void printResult(double ***retArray, enum goalEnum goal) {
    /* TODO: implement function */
}

/**
 * Checks goal string and returns enum accordingly
 * @param string - represents goal of the call to main
 * @return appropriate enum for goal, or INVALID if invalid
 */
static enum goalEnum checkGoal(char *string) {
    if (!strcmp(string, "spk")) {
        return spk;
    }
    if (!strcmp(string, "wam")) {
        return wam;
    }
    if (!strcmp(string, "ddg")) {
        return ddg;
    }
    if (!strcmp(string, "lnorm")) {
        return lnorm;
    }
    if (!strcmp(string, "jacobi")) {
        return jacobi;
    }
    return INVALID;
}

/**
 * Takes a n*m-sized array, and coverts it into "matrix" form
 * if allocation fails, returns NULL
 * @param array - pointer to a n*m-sized array
 * @param n - # of rows
 * @param m - # of columns
 * @return a 2D array referencing this array
 */
static double **arrayToTwoDimArray(double **array, int n, int m) {
    double **twoDimArray;
    int i;

    twoDimArray = calloc(n, sizeof(double *));
    ASSERT_ERROR(twoDimArray != NULL)

    for (i = 0; i < n; i++) {
        twoDimArray[i] = (*array) + i * m;
    }

    return twoDimArray;
}
/**
 * Returns a 2D matrix of size n X m
 * @param n - # of rows
 * @param m - # of cols
 */
static double ** createMatrix(int n, int m) {
    double *array;
    array = calloc(n*m, sizeof (double ));
    ASSERT_ERROR(array != NULL)
    return arrayToTwoDimArray(&array, n, m);
}

/**
 * Gets an empty matrix and puts 1 in the diagonal, thus making it
 * an identity matrix (I)
 * @param emptyMatrix - a pointer to a new, empty matrix
 * @param matrixSize - size of given nXn matrix
 */
void makeIntoIdentityMatrix(double ***emptyMatrix, int matrixSize) {
    int i;
    for (i = 0; i < matrixSize; i++){
        (*emptyMatrix)[i][i] = 1;
    }
}

/**
 * Frees 1D support matrix and 2D structure built on top of that
 * (In accordance to createMatrix matrix initialization)
 * @param matrix - pointer to matrix we want to free
 */
void freeMatrix(double ***matrix) {
    free(*matrix[0]);
    free(*matrix);
}

/**
 * Takes the datapoint info, and returns the weighted adj matrix.
 * @param d - dimension of datapoints
 * @param arraySize - amount of datapoints
 * @param datapoint - pointer to 2D-array containing all datapoints
 */
double ** calcWeightedAdjMatrix(int d, int arraySize, double ***datapoint) {

    int i, j;
    double value, **weightedAdjMatrix;

    /** creating 2-D array */
    weightedAdjMatrix = createMatrix(arraySize, arraySize);

    /** calculating the weights */
    for (i = 0; i < arraySize; i++) {
        for (j = i + 1; j < arraySize; j++) {
            value = calcWeight(d, datapoint, i, j);
            weightedAdjMatrix[i][j] = value;
            weightedAdjMatrix[j][i] = value;
        }
    }

    return weightedAdjMatrix;
}

/**
 * Takes the datapoint info, and returns the weighted adj (i,j).
 * @param d - dimension of datapoints
 * @param datapoint - 2D-array containing all datapoints
 * @param i - the first datapoint
 * @param j - the second datapoint
 * @return weighted adj (i,j))
 */
double calcWeight(int d, double ***datapoint, int i, int j) {

    int t;
    double value = 0;

    /** calculating the sum of the distances */
    for (t = 0; t < d; t++) {
        value += pow(((*datapoint)[i][t] - (*datapoint)[j][t]), 2);
    }

    value = exp(sqrt(value) / (-2));
    return value;
}

/**
 * Takes the weightedAdjMatrix, and returns the diagonal degree matrix
 * corresponding to it
 * @param arraySize - amount of datapoints
 * @param weightedAdjMatrix - 2-D array containing all the weights

 */
double ** calcDiagDegMatrix(int arraySize, double ***weightedAdjMatrix) {

    int i, j;
    double sumofweights, **diagDegMatrix;

    diagDegMatrix = createMatrix(arraySize, arraySize);

    for (i = 0; i < arraySize; i++) {
        sumofweights = 0;
        for (j = 0; j < arraySize; j++) {
            sumofweights += (*weightedAdjMatrix)[i][j];
        }
        (diagDegMatrix)[i][i] = sumofweights;
    }

    return diagDegMatrix;
}

/**
 * Takes the weighted and the diagonal matrices info,
 * and returns the normLaplacian matrix.
 * @param arraySize - amount of datapoints
 * @param weightedAdjMatrix - pointer to 2-D array containing all the weights
 * @param diagDegMatrix - pointer to 2-D array containing the weights sum
 *                        for each datapoint
 */
double ** calcNormLaplacian(int arraySize, double ***weightedAdjMatrix,
                         double ***diagDegMatrix) {

    int i, j;
    double value, *sqrtDegMatrix, **normLaplacianMatrix;

    normLaplacianMatrix = createMatrix(arraySize, arraySize);

    /** 1-D array represents sqrt of diagDegMatrix */
    sqrtDegMatrix = calloc(arraySize, sizeof(double));
    ASSERT_ERROR(sqrtDegMatrix != NULL)

    for (i = 0; i < arraySize; i++) {
        sqrtDegMatrix[i] = 1 / sqrt((*diagDegMatrix)[i][i]);
    }

    /**  sqrtDegMatrix[i][j] == sqrtDeg[i] * weight[i][j] * sqrtDeg[j]
     *                       == sqrtDegMatrix[j][i]
     * convince yourself it's true :) */

    for (i = 0; i < arraySize; i++) {
        (normLaplacianMatrix)[i][i] = 1;
        for (j = i + 1; j < arraySize; j++) {
            value = -(sqrtDegMatrix[i] * (*weightedAdjMatrix)[i][j] * sqrtDegMatrix[j]);
            (normLaplacianMatrix)[i][j] = value;
            (normLaplacianMatrix)[j][i] = value;
        }
    }

    free(sqrtDegMatrix);

    return normLaplacianMatrix;
}

/**
 * Takes a matrix and returns the result of running the Jacobi algorithm on it.
 * @param arraySize - size of the nXn matrix
 * @param inputMatrix - pointer to 2D array representation of input matrix
 */
double ** calcJacobi(int arraySize, double ***inputMatrix) {

    int row, col, converge;
    double c, s;
    double **A, **Atag, **V, **jacobiMatrix;

    /** initialization */
    Atag = createMatrix(arraySize, arraySize);
    V = createMatrix(arraySize, arraySize);

    makeIntoIdentityMatrix(&V, arraySize);

    /** run until convergence -  */

    A = (*inputMatrix);

    do {
        findmatrixP(arraySize, &A, &c, &s, &row, &col);         // P
        updateAtag(arraySize, &Atag, &A, c, s, row, col);         // A' = P^T * A * P
        updateV(arraySize, &V, c, s, row, col);               // V *= P
        converge = convergenceCheck(arraySize, &A, &Atag);

        if(!converge) {
            copymatrix(arraySize, &Atag, &A);                // A = Atag
        }
    } while (!converge);              // as long as delta > epsilon

    //

    /** Atag has the A"A
 *  V has the V"A */

    jacobiMatrix = copyJacoby(arraySize, &Atag, &V);

    freeMatrix(&Atag);
    freeMatrix(&V);

    return jacobiMatrix;
}


/*
double ** calcJacobi(int arraySize, double ***inputMatrix) {

    int row, col, c, s, converge;

    double **A, **P, **Atag, **V, **jacobiMatrix;

    ** initialization
    P = createMatrix(arraySize, arraySize);
    Atag = createMatrix(arraySize, arraySize);
    V = createMatrix(arraySize, arraySize);

    makeIntoIdentityMatrix(&P, arraySize);
    makeIntoIdentityMatrix(&V, arraySize);

    ** run until convergence -

    A = (*inputMatrix);

    do {
        findmatrixP(arraySize, &A, &P, &row, &col);         // P
        updateAtag(arraySize, &Atag, &A, &P, row, col);         // A' = P^T * A * P
        updateV(arraySize, &V, &P, row, col);               // V *= P
        converge = convergenceCheck(arraySize, &A, &Atag);

        if(!converge) {
            copymatrix(arraySize, &Atag, &A);                // A = Atag
        }
    } while (!converge);              // as long as delta > epsilon

 //

    ** Atag has the A"A
 *  V has the V"A

    jacobiMatrix = copyJacoby(arraySize, &Atag, &V);
    return jacobiMatrix;
}
*/


double **copyJacoby(int arraySize, double ***Atag, double ***V) {
    int i, j;
    double **jacobiMatrix;

    jacobiMatrix = createMatrix(arraySize+1, arraySize);

    for(j = 0; j < arraySize; j++) {// copy the EigenValues
        jacobiMatrix[0][j] = (*Atag)[j][j];
    }

    for(i = 1; i < arraySize + 1; i++) { // copy the EigenVectors
        for(j = 0; j < arraySize; j++){
            jacobiMatrix[i][j] = (*V)[i-1][j];
        }
    }

    return jacobiMatrix;
}


/**
 * Check weather the the index of the maximum odd-diag value in A.
 * @param arraySize - amount of datapoints
 * @param A - 2-D symmetric array
 * @param Atag - pointer to the row index of the max value in A (off diag), to be calculated
 * @return 1 if converges, else 0
 */
int convergenceCheck(int arraySize, double ***A, double ***Atag) {
    double offA, offAtag, epsilon = EPSILON;
    offA = calcOff(arraySize, A);
    offAtag = calcOff(arraySize, Atag);

    return (offA - offAtag) <= epsilon;
}

/**
 * Find the sum of squares of all off-diagonal elements of matrix
 * @param arraySize - amount of datapoints
 * @param matrix - 2-D symmetric array
 * @return the calculated sum
 */
double calcOff(int arraySize, double ***matrix){
    double off = 0;
    int i, j;
    for(i = 0; i < arraySize; i++){
        for(j = i+1; j < arraySize; j++){
            off += (2 * pow((*matrix)[i][j],2));
        }
    }
    return off;
}

/**
 * Gets 2 matrices, and cupy the 1'st to the 2'nd a matrix of size (N * N).
 * @param arraySize - amount of datapoints
 * @param matrix1 - source matrix
 * @param matrix2 - matrix to be updated
 * TODO if ill copy " symmetric matrix, will that be better ?
 */
void copymatrix(int arraySize, double ***matrix1, double ***matrix2) {

    int i,j;

    for (i = 0; i < arraySize; i++) {
        (*matrix2)[i][i] = (*matrix1)[i][i];
        for (j = i+1; j < arraySize; j++) {
            (*matrix2)[i][j] = (*matrix1)[i][j];
            (*matrix2)[j][i] = (*matrix1)[j][i];
        }
    }
}

/**
 * Calculate the rotation matrix P using the matrix A
 * @param arraySize - amount of datapoints
 * @param A - 2-D symmetric array
 * @param P - 2-D array to be calculated
 * @param row - pointer to the row index of the max value in A (off diag), to be calculated
 * @param col - pointer to the col index of the max value in A (off diag), to be calculated
 */
void findmatrixP(int arraySize, double ***A, double *c, double *s, int *row, int *col) {

    double teta, t, sign = -1;
    findMaxOffDiag(arraySize, A, row, col);

    /** update P values */
    teta = ((*A)[*col][*col] - (*A)[*row][*row]) / (2 * (*A)[*row][*col]);

    if (teta >= 0) {
        sign = 1;
    }
    t = sign / (fabs(teta) + sqrt(pow(teta, 2) + 1));
    (*c) = 1 / sqrt(pow(t, 2) + 1);
    (*s) = t * (*c);
}

/**
 * Find the index of the maximum odd-diag value in A.
 * @param arraySize - amount of datapoints
 * @param A - 2-D symmetric array
 * @param row - pointer to the row index of the max value in A (off diag), to be calculated
 * @param col - pointer to the col index of the max value in A (off diag), to be calculated
 */
void findMaxOffDiag(int arraySize, double ***A, int *row, int *col) {

    int i, j;
    double curmax = FLT_MIN;

    for (i = 0; i < arraySize; i++) {
        for (j = i + 1; j < arraySize; j++) {
            if ((*A)[i][j] > curmax) {
                curmax = (*A)[i][j];
                (*row) = i;
                (*col) = j;
            }
        }
    }
}

/**
 *Calculate the matrix Atag using the current A and P, A' = P^T * A * P.
 * @param arraySize - amount of datapoints
 * @param Atag - 2-D array to be calculated
 * @patam p -  2-D symmetric array
 * @param row - row index of the max value in Atag (off diag)
 * @param col - col index of the max value in Atag (off diag)
 */
void updateAtag(int arraySize, double ***Atag,double ***A, double c, double s, int row, int col) {

    int r;
    double ri, rj, ij, ii, jj;

    /** update each Atag[r,i]/Atag[r,j], for : r != i,j
     * can cut computation in half because the matrix is symmetric */
    for (r = 0; r < arraySize; r++) {
        ri = (*A)[r][row];
        rj = (*A)[r][col];
        if ((r != row) && (r != col)) {

            (*Atag)[r][row] = (c * ri) - (s * rj);
            (*Atag)[row][r] = (*Atag)[r][row];
            (*Atag)[r][col] = (c * rj) + (s * ri);
            (*Atag)[col][r] = (*Atag)[r][col];
        }
    }
    /**  Atag[i,i] / Atag[j,j]/ Atag[i,j] / Atag[j,i] */
    ij = (*A)[row][col];
    ii = (*A)[row][row];
    jj = (*A)[col][col];

    (*Atag)[row][row] = (c * c * ii) + (s * s * jj) - (2 * s * c * ij);
    (*Atag)[col][col] = (s * s * ii) + (c * c * jj) + (2 * s * c * ij);
    (*Atag)[row][col] = 0;
    (*Atag)[col][row] = 0;
}

/**
 *Calculate the current matrix V, which is the Eigenvectors matrix. V *= P
 * @param arraySize - amount of datapoints
 * @param V - 2-D array to be calculated
 * @patam P -  2-D symmetric array
 * @param row - row index of the max value in Atag (off diag)
 * @param col - col index of the max value in Atag (off diag)
 */
void updateV(int arraySize, double ***V, double c, double s, int row, int col) {

    int r;
    double ri, rj;

    for (r = 0; r < arraySize; r++) {
        ri = (*V)[r][row];
        rj = (*V)[r][col];
        /** V'[r][row] = (V[r][row] * P[row][row]) + (V[r][col] * P[col][row])
         * V'(r,i) = (V(r,i) * c) + ((V(r,j) * (-s)) */
        (*V)[r][row] = (ri * c) + (rj * -s);

        /** V'[r][col] = (V[r][row] * P[row][col]) + (V[r][col] * P[col][col])
         * V'(r,j) = (V(r,i) * s) + ((V(r,j) * c) */
        (*V)[r][col] = (ri * s) + (rj * c);
    }
}

double **
calcSpectralClusters(int *k, int arraySize, int isCAPI,
                     double ***jacobiMatrix) {

    double **U, **spkMatrix;

    //TODO: sort jacobi

    //TODO: determine k (if k==0) AND change k variable accordingly, important for CAPI

    //TODO: get k first eigenvectors from jacobi into U, normalize U

    //TODO: K-Means algorithm - if CAPI - need to return to python for K-Means++
    //TODO:                             - need to get back result from python
    //TODO: Need to cut off function when returning to python, then start new
    //TODO: method with KMeans implementation and finito

    //TODO: free memory and finish

    return NULL;
}

double ** KMeansAlgorithm(int k, int *d, int arraySize, double ***datapoints,
                          int isCAPI, int **init_centroids){
    return NULL;
}


void printTest(double **matrix, int n, int m){
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}