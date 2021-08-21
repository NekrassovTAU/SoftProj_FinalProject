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
#define MAX_ITER 300

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

#define BEEP(num) \
    printf("%s %f%s", "SPKMEANS Check #" , num, "\n");

/**
 * main, a shell function for the spectral clustering algorithm implementation
 */
int main(int argc, char *argv[]) {
    int rowCount, colCount;
    checkArgs(argc, argv, 0, &rowCount, &colCount);
/*
    int i;
    double **arr, **Atag, **V, **jacobiMatrix, *jacobiArray;

    arr = (double**) calloc(3,sizeof (double*));

    for(i = 0; i < 3; i++) {
        arr[i] =  calloc(3, sizeof(double));
    }

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
      jacobiMatrix = calcJacobi(3, &arr);
      printf("\n\n");
      printTest(jacobiMatrix, 4, 3);

      sortJacobi(3, &jacobiMatrix);
      printf("\n\n");
      printTest(jacobiMatrix, 4, 3);
*/
/*
    double **weightdAdjMatrix, **diagDegMatrix, **normLaplacianMatrix, **jacobiMatrix;

    // check the calc functions directly
    weightdAdjMatrix = TesterToWeight();

    TesterToSortJacobi();
*/
    return 0;
}

/**
 * Checks the validity of arguments passed, reads and parses info in files
 * and starts the implementation of spectral k-means clustering
 * @param argc - # of arguments passed to main
 * @param origArgv - the arguments passed to main
 * @param isCAPI - used in order to distinguish between
 *                 C-based and C-API-based call to this function
 * Rest of @params are used to convey data regarding output size and datapoints
 * dimensions
 */
double **checkArgs(int argc, char **origArgv, int isCAPI, int *returnRowCount,
                   int *returnColCount) {
    int k, d, arraySize, *dArraySizeInfo;
    double *datap_array, **ret_matrix, **datapoints;
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
    ASSERT_ERROR(dArraySizeInfo != NULL)


    /** process the datapoint file, and return info about their amount (arraySize)
     * and number of properties (d) */
    processDatapoints(origArgv[3], &datap_array, &datapoints, &dArraySizeInfo);

    d = dArraySizeInfo[0];
    arraySize = dArraySizeInfo[1];
    free(dArraySizeInfo);

    if (goal == spk){
        ASSERT_ARGS(k < arraySize)
    }

    /** initialize the process to achieve the provided goal */
    ret_matrix = initProcess(&k, d, arraySize, goal, &datapoints, isCAPI);

    /** Determining row and col size for return matrix */
    determineRowAndCol(goal, k, arraySize, returnRowCount, returnColCount);

    /** Free allocated memory and terminate*/

    return ret_matrix;
}

/**
 *
 * @param goal - goal of call to algorithm
 * @param k - amount of clusters (used in case of goal == spk)
 * @param arraySize - amount of datapoints
 * @param rowCount - the amount of rows in the return matrix (updated in this
 *                   method).
 * @param colCount - the amount of columns in the return matrix (updated in this
 *                   method).
 */
void determineRowAndCol(enum goalEnum goal, int k, int arraySize, int *rowCount,
                        int *colCount) {
    switch (goal) {
        case jacobi:{
            *rowCount = arraySize + 1;
            *colCount = arraySize;
            break;
        }
        case spk:{
            *rowCount = arraySize;
            *colCount = k;
            break;
        }
        default:{
            *rowCount = arraySize;
            *colCount = arraySize;
            break;
        }
    }
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

    copyDatapoint(datapoint, datap_array, arraySize, d);

    free(*datap_array);

    (*dArraySizeInfo)[0] = d;
    (*dArraySizeInfo)[1] = arraySize;
}

/**
 *
 * @param datapoint - 2D array to which we will input the datapoints info
 * @param datap_array - 1D array in which we inputed info in the first-place
 * @param arraySize - amount of datapoints
 * @param d - dimension of datapoints
 */
static void
copyDatapoint(double ***datapoint, double **datap_array, int arraySize, int d) {
    int i, j, count;
    count = 0;

    *datapoint = createMatrix(arraySize, d);
    for (i = 0 ; i < arraySize ; i++){
        for (j = 0 ; j < d ; j++){
            (*datapoint)[i][j] = (*datap_array)[count];
            count++;
        }
    }
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
        printResult(&ret_matrix, goal, arraySize, *k, d);
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

    spkMatrix = calcSpectralClusters(k, arraySize, isCAPI, &jacobiMatrix, datapoint);

    freeMatrix(&weightedAdjMatrix);
    freeMatrix(&diagDegMatrix);
    //freeMatrix(&normLaplacianMatrix);
    //freeMatrix(&jacobiMatrix);

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
void printResult(double ***retArray, enum goalEnum goal, int arraySize, int k, int d) {

    if((goal == wam) || (goal == lnorm)){
        printRegular(retArray, arraySize, arraySize);
    }
    if(goal == jacobi){
        printJacobi(retArray, arraySize);
    }
    if(goal == ddg) {
        printDiagonal(retArray, arraySize);
    }
    if(goal == spk){
        printRegular(retArray, k, k);
    }
}

/**
 * Prints matrix of size rows * cols
 * @param retArray - pointer to matrix to be printed to user
 * @param rows - number of rows in the matrix
 * @param cols - number of columns in the matrix
 */
void printRegular(double ***retArray, int rows, int cols){

    int i, j;

    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            printf("%.4f", (*retArray)[i][j]);
            if (j != cols - 1){   /* not last in the line */
                printf("%s", ",");
            }
            else{
                printf("\n");
            }
        }
    }
}
/**
 * Prints the Jacobi matrix with eigenValues & eigenVecrtors as lines
 * @param retArray - pointer to matrix to be printed to user
 * @param arraySize - number of rows in the matrix
 */
void printJacobi(double ***retArray, int arraySize) {

    int i, j;

    /** print the eigenValues - first line */
    for (j = 0; j < arraySize; j++) {
        printf("%.4f", (**retArray)[j]);
        if (j != arraySize - 1) {   /* not last component of the cluster*/
            printf("%s", ",");
        } else {
            printf("\n");
        }
    }

    /** print the eigenVectors - each verctor as a row */
    for (j = 0; j < arraySize; j++) {
        for (i = 1; i < arraySize + 1; i++) {
            printf("%.4f", (*retArray)[i][j]);
            if (i != arraySize) {
                printf("%s", ",");
            } else {
                printf("\n");
            }
        }
    }
}

/**
 * Prints the diagonal matrix, filled with 0's
 * @param retArray - pointer to matrix to be printed to user
 * @param arraySize - size of the given matrix
 */
void printDiagonal(double ***retArray, int arraySize) {

    int i, j;

    for(i = 0; i < arraySize; i++){
        for(j = 0; j < arraySize; j++){
            if(i == j){
                printf("%.4f", (**retArray)[i]);
            }
            else{
                printf("0.0000");
            }
            if (j != arraySize - 1){   /* not last component of the cluster*/
                printf("%s", ",");
            }
            else{
                printf("\n");
            }
        }
    }
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
double ** createMatrix(int n, int m) {
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
 * @return the weighted adjacency matrix
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
 * @return the weight between the i,j data points
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
 * Takes the weighted matrix, and returns representation of the diagonal degree matrix
 * @param arraySize - amount of datapoints
 * @param weightedAdjMatrix - pointer to 2D-array containing all the weights
 * @return 2-D array of size (1 * arraySize) represents the diag of the diagonal degree matrix
 */
double **calcDiagDegMatrix(int arraySize, double ***weightedAdjMatrix) {

    int i, j;
    double sumOfWeights, **diagDegMatrix;

    diagDegMatrix = createMatrix(1, arraySize);

    for (i = 0; i < arraySize; i++) {
        sumOfWeights = 0;
        for (j = 0; j < arraySize; j++) {
            sumOfWeights += (*weightedAdjMatrix)[i][j];
        }
        diagDegMatrix[0][i] = sumOfWeights;
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
 * @return The norm laplacian matrix
 */
double ** calcNormLaplacian(int arraySize, double ***weightedAdjMatrix,
                         double ***diagDegMatrix) {

    int i, j;
    double value, **normLaplacianMatrix;

    normLaplacianMatrix = createMatrix(arraySize, arraySize);


    for (i = 0; i < arraySize; i++) {
        (*diagDegMatrix)[0][i] = 1 / sqrt((*diagDegMatrix)[0][i]);
    }

    /**  sqrtDegMatrix[i][j] == sqrtDeg[i] * weight[i][j] * sqrtDeg[j]
     *                       == sqrtDegMatrix[j][i]
     * convince yourself it's true :) */

    for (i = 0; i < arraySize; i++) {
        (normLaplacianMatrix)[i][i] = 1;
        for (j = i + 1; j < arraySize; j++) {

            value = -(*diagDegMatrix)[0][i] * (*weightedAdjMatrix)[i][j] * (*diagDegMatrix)[0][j];
            (normLaplacianMatrix)[i][j] = value;
            (normLaplacianMatrix)[j][i] = value;
        }
    }

    return normLaplacianMatrix;
}

/**
 * Takes a symmetrical matrix, and returns the result of running the Jacobi algorithm on it.
 * @param arraySize - size of the nXn matrix
 * @param inputMatrix - pointer to 2D array, which is symmetrical matrix
 * @return The jacobi matrix
 */
double ** calcJacobi(int arraySize, double ***inputMatrix) {

    int row, col, converge, iterations = 0;
    double c, s, offA, offAtag, epsilon = EPSILON;
    double **A, **V, **jacobiMatrix;

    /** initialization */
    V = createMatrix(arraySize, arraySize);
    makeIntoIdentityMatrix(&V, arraySize);

    /** run until convergence -  */
    A = (*inputMatrix);
    offAtag = calcOff(arraySize, &A);

    do {
        offA = offAtag;

        findmatrixP(arraySize, &A, &c, &s, &row, &col); // P
        updateAtag(arraySize, &A, c, s, row, col, &offAtag); // A' = P^T * A * P
        updateV(arraySize, &V, c, s, row, col);  // V *= P

        converge = (offA - offAtag) < epsilon ? 1 : 0;
        iterations++;

     } while (!converge && (iterations < 100)); // as long as (delta > epsilon) or number of iterations is under 100

    /** Atag has the eigenValues
      *  V has the eifenVectors */

    jacobiMatrix = copyJacoby(arraySize, &A, &V);

    freeMatrix(&A);
    /**  seems fine now TODO: FIX NEXT LINE. PROBLEM IN FREE V - when i run threw terminal it doesn't get to the return */
    freeMatrix(&V);

    return jacobiMatrix;
}


/**
 * Copies the eigenValues and eigenVectors into the final jacobiMatrix and returns it
 * @param arraySize - amount of datapoints
 * @param Atag - pointer to 2-D array containing the eigenValues on the diag
 * @param V - pointer to 2-D array containing the eigenVectors as columns
 * @return The jacobi matrix
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
 /** TODO: delete if not necessary */
int convergenceCheck(int arraySize, double ***A, double ***Atag) {
    double offA, offAtag, epsilon = EPSILON;
    offA = calcOff(arraySize, A);
    offAtag = calcOff(arraySize, Atag);

    return (offA - offAtag) <= epsilon;
}

/**
 * Find the sum of squares of all off-diagonal elements of the given matrix
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
 */
/** TODO: delete if not necessary */
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
 * Calculate the paramaters c and s in the matrix P using the matrix A
 * @param arraySize - amount of datapoints
 * @param A - 2-D symmetric array
 * @param c - pointer to the value in the index (row,row) in matrix P
 * @param s - pointer to the value in the index (row,col) in matrix P
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
            if (fabs((*A)[i][j]) > curmax) {
                curmax = fabs((*A)[i][j]);
                (*row) = i;
                (*col) = j;
            }
        }
    }
}

/**
 *Calculate the matrix Atag using the current A and P (based on c&s) : A' = P^T * A * P
 * & update the offAtag accordingly to the changes made
 * @param arraySize - amount of datapoints
 * @param A - 2-D array to be calculated based on the given one
 * @param c - the value of P[row][row]
 * @param s - the value of P[row][col]
 * @param row - row index of the max value in Atag (off diag)
 * @param col - col index of the max value in Atag (off diag)
 * @param offAtag - pointer to the current off(A) value
 */
void updateAtag(int arraySize, double ***A, double c, double s, int row, int col, double *offAtag) {

    int r;
    double ri, rj, ij, ii, jj;

    /** update each Atag[r,i]/Atag[r,j], for : r != i,j
     * can cut computation in half because the matrix is symmetric */
    for (r = 0; r < arraySize; r++) {

        if ((r != row) && (r != col)) {
            ri = (*A)[r][row];
            rj = (*A)[r][col];
            (*A)[r][row] = (c * ri) - (s * rj);
            (*A)[row][r] = (*A)[r][row];
            (*A)[r][col] = (c * rj) + (s * ri);
            (*A)[col][r] = (*A)[r][col];

            *offAtag += 2 * ((*A)[r][row] * (*A)[r][row] + (*A)[r][col] * (*A)[r][col]) - 2 * (ri * ri + rj * rj);
        }
    }
    /**  Atag[i,i] / Atag[j,j]/ Atag[i,j] / Atag[j,i] */
    ij = (*A)[row][col];
    ii = (*A)[row][row];
    jj = (*A)[col][col];

    (*A)[row][row] = (c * c * ii) + (s * s * jj) - (2 * s * c * ij);
    (*A)[col][col] = (s * s * ii) + (c * c * jj) + (2 * s * c * ij);
    (*A)[row][col] = 0;
    (*A)[col][row] = 0;

    *offAtag -= 2 * (ij * ij);
}

/**
 *Calculate the current matrix V, which is the Eigenvectors matrix. V *= P (based on c&s)
 * @param arraySize - amount of datapoints
 * @param V - 2-D array to be calculated
 * @param c - the value of P[row][row]
 * @param s - the value of P[row][col]
 * @param row - row index of the max value in Atag (off diag)
 * @param col - col index of the max value in Atag (off diag)
 */
void updateV(int arraySize, double ***V, double c, double s, int row, int col){

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
calcSpectralClusters(int *k, int arraySize, int isCAPI, double ***jacobiMatrix,
                     double ***datapoints) {

    double **U, **T, **spkMatrix, **combined;
    int i, *init_centroids;

    /** create (2 * n) matrix, which includes the eigenValues in the 1'st line,
      * and the original index of each eigenValue in the line below  */
    combined = createMatrix(2, arraySize);

    /** sort jacobi's eigenValues*/
    sortEigenValues(arraySize, jacobiMatrix, &combined);

    /** determine k (if k==0) AND change k variable accordingly, important for CAPI */
    if(*k == 0){
        *k = getK(arraySize, &combined);
    }

    /** create U with the first k sorted jacobi's eigenVectors & normalize it*/
    U = getMatrixU(*k, arraySize, jacobiMatrix, combined);
    freeMatrix(&combined);

    normalizeU(*k, arraySize, &U);

    T = U;

    //TODO: K-Means algorithm - if CAPI - need to return to python for K-Means++
    if (isCAPI){
        return T;
    }

    init_centroids = calloc((*k), sizeof(int));
    for (i=0 ; i < (*k); i++){
        init_centroids[i] = i;
    }

    spkMatrix = KMeansAlgorithm(*k, arraySize, &T, &init_centroids);

    //TODO: free memory and finish
    free(init_centroids);
    freeMatrix(&T);
    return spkMatrix;
}

double **
KMeansAlgorithm(int k, int arraySize, double ***T, int **init_centroids) {
    int max_iter, update, *datap_cluster_assignment, *countArray;
    double **t_centroids, **sumArrayHead;

    max_iter = MAX_ITER;
    update = 1;

    t_centroids = centroidsFromList(T, init_centroids, k);

    datap_cluster_assignment = calloc( arraySize , sizeof (int));

    sumArrayHead = createMatrix(k, k);

    countArray = calloc(k, sizeof(int));
    ASSERT_ERROR(countArray != NULL)

    while(max_iter > 0 && update){
        update = updateCentroidsPerDatap(T, &t_centroids,
                                         &datap_cluster_assignment, k, k,
                                         arraySize, &sumArrayHead, &countArray);

        max_iter--;
    }

    free(datap_cluster_assignment);
    freeMatrix(&sumArrayHead);
    free(countArray);

    return t_centroids;

}

/**
 * Takes init_centroids indices and translates them to actual centroids list,
 * according to pointList provided
 * @param pointList - list of points from which we take the centroids
 * @param init_centroids - list of indices of initial centroids
 * @param k - amount of centroids (and their dimension in this case)
 * @return 2D array with the chosen centroids copied into them
 */
double **centroidsFromList(double ***pointList, int **init_centroids, int k) {

    int i, j;
    double **centroids;

    centroids = createMatrix(k , k);

    for (i = 0 ; i < k ; i++){
        for (j = 0 ; j < k ; j++){
            centroids[i][j] = (*pointList)[(*init_centroids)[i]][j];
        }
    }

    return centroids;
}

int updateCentroidsPerDatap(double ***datapoint, double ***centroid,
                            int **datap_cluster_assignment, int d, int k,
                            int arraySize, double ***sumArrayHead,
                            int **countArray) {
    int i, j, v, min_cluster, update, currCluster;
    double dist, min_dist, new_value;

    update = 0;

    for (i = 0 ; i < arraySize ; i++){
        min_dist = FLT_MAX, min_cluster = -1;

        for (j = 0; j < k; j++ ){
            dist = 0;
            for (v = 0; v < d; v++){
                dist += ((*datapoint)[i][v] - (*centroid)[j][v])*((*datapoint)[i][v] - (*centroid)[j][v]);
            }

            if (min_dist > dist){
                min_dist = dist;
                min_cluster = j;
            }
        }
        if((*datap_cluster_assignment)[i] != min_cluster){ /* there is a change in one or more of the data points cluster assignment */
            update = 1;
        }
        (*datap_cluster_assignment)[i] = min_cluster;
    }

    /* loop to initialize sum/counter*/

    for(i = 0; i < arraySize; i++){ /*count and sum up all the sizes*/
        currCluster = (*datap_cluster_assignment)[i];
        (*countArray)[currCluster]++;
        for(v = 0; v < d; v++){
            (*sumArrayHead)[currCluster][v] += (*datapoint)[i][v];
        }
    }

    /*update the new clusters and initialize to 0*/
    for(j = 0; j < k; j++) { /* each loop for different cluster*/
  //      printf(" %d ", (*countArray)[j]);
        for (v = 0; v < d; v++) { /* each loop for opponent of the current cluster*/

            new_value = (*sumArrayHead)[j][v] / (*countArray)[j];
            (*centroid)[j][v] = new_value;
            (*sumArrayHead)[j][v] = 0;
        }
        (*countArray)[j] = 0;
    }

    return update;
}


/**
 * Takes the jacobi matrix, and update the combined matrix to possess the eigenValues in sorted order.
 * @param arraySize - size of the nXn matrix
 * @param jacobiMatrix - pointer to the jacobi matrix
 * @param combined - pointer to 2D array, represents pairs of eigenValue & original index
 */
void sortEigenValues(int arraySize, double ***jacobiMatrix, double ***combined) {
    int i;
    double **tmp;

    tmp = createMatrix(2, arraySize);
    /** initialization of combined with the initial values */
    for (i = 0; i < arraySize; i++){
        (*combined)[0][i] = (*jacobiMatrix)[0][i];
        (*combined)[1][i] = i;
    }

    mergeSort1(combined, &tmp, 0, arraySize - 1); // sort the eigenValues

    freeMatrix(&tmp);
}


/**
 * Takes the jacobi matrix, sort it's eigenValues, and update the combined matrix.
 * @param combined - pointer to 2D array, includes pairs of eigenValue & originalValue
 * @param tmp - pointer to 2D Auxiliary array to combined
 * @param low - the lower index of the subArray to be sorted
 * @param high - the higher index of the subArray to be sorted
 */
void mergeSort1(double ***combined, double ***tmp, int low, int high){
    int mid;

    if(low < high){
        mid = (low + high) / 2;
        mergeSort1(combined, tmp, low, mid);
        mergeSort1(combined, tmp, mid + 1, high);
        merge1(combined, tmp, low, mid, high);
    }
    else{
        return;
    }
}


/**
 * Takes the combined 2D array, and merge 2 sorted subArrays of it, based on the given indices
 * @param combined - pointer to 2D array, includes pairs of eigenValue & originalValue
 * @param tmp - pointer to 2D Auxiliary array to combined
 * @param low - the lower index of the 1'st subArray to be merged
 * @param mid - the higher index of the 1'st subArray to be merged
 * @param high - the higher index of the 2'nd subArray to be merged
 */
void merge1(double ***combined, double ***tmp, int low, int mid, int high){
    int l1, l2, i;

    for(l1 = low, l2 = mid+1, i = low; l1 <= mid && l2 <= high; i++) {

        if ((*combined)[0][l1] <= (*combined)[0][l2]) {

            (*tmp)[0][i] = (*combined)[0][l1];
            (*tmp)[1][i] = (*combined)[1][l1];

            l1++;
        }
        else {

            (*tmp)[0][i] = (*combined)[0][l2];
            (*tmp)[1][i] = (*combined)[1][l2];
            l2++;
        }
    }

    while(l1 <= mid){  // copy the rest of the "1'st array"
        (*tmp)[0][i] = (*combined)[0][l1];
        (*tmp)[1][i] = (*combined)[1][l1];

        i++;
        l1++;
    }

    while(l2 <= high){ // copy the rest of the "2'nd array"
        (*tmp)[0][i] = (*combined)[0][l2];
        (*tmp)[1][i] = (*combined)[1][l2];
        i++;
        l2++;
    }

    for(i = low; i <= high; i++){    // copy the eigenValues from the temp array, back to the original jacobiMatrix
        (*combined)[0][i] = (*tmp)[0][i];
        (*combined)[1][i] = (*tmp)[1][i];
    }
}


/** TODO: delete if not necessary */
void sortEigenVectors(int arraySize, double ***jacobiMatrix, double **combined) {
    int i, j, column;
    double **pointer, **tmpMatrix = createMatrix(arraySize + 1, arraySize);

    for(j = 0; j < arraySize; j++){                 // use temp t
        column = (int) combined[1][j];
        for(i = 0; i < arraySize + 1; i++){
                tmpMatrix[i][j] = (*jacobiMatrix)[i][column];
            }
        }

    for(i = 0; i < arraySize + 1; i++){                 // use temp t
        for(j = 0; j < arraySize; j++) {
            (*jacobiMatrix)[i][j] = tmpMatrix[i][j];
        }
    }

    /** switch pointers **/
    /**
     * pointer = (*jacobiMatrix);
     * jacobiMatrix = &tmpMatrix;
     * freeMatrix(&pointer);
                                */
}


/**
 * Returns K based on the biggest difference between two adjacent eigenValues
 * @param arraySize - amount of datapoints
 * @param combined - jacobi matrix which includes eigenValues and EigenVectors
 */
int getK(int arraySize, double ***combined) {

    int j, k = 0;
    double cur, max = (-1);
    for(j = 0; j < arraySize / 2; j++) {
        cur = fabs((*combined)[0][j] - (*combined)[0][j + 1]);
        if(cur > max){
            max = cur;
            k = j + 1;
        }
    }
    return k;
}


/**
 * Returns matrix U (n*k), which is the first K columns (eigenVectors) in jacobiMatrix, based on combined
 * @param k - number of initial clusters
 * @param arraySize - amount of datapoints
 * @param jacobiMatrix - jacobi matrix which includes eigenValues and EigenVectors
 * @param combined - 2D array contatins pairs of sorted eigenValues & original index
 * @return the matrix U
 */
double **getMatrixU(int k, int arraySize, double ***jacobiMatrix, double **combined){

    int i, j, col;
    double **U = createMatrix(arraySize, k);

    /** take the first K eigenVectors if the jacobi was sorted by the eigenValue */
    for(j = 0; j < k; j++){
        col = (int) combined[1][j]; // The original index of the j eigenValue is sorted order
        for(i = 0; i < arraySize; i++){
            U[i][j] = (*jacobiMatrix)[i+1][col];
        }
    }

    return U;
}


/**
 * Normalize the given matrix U by rows
 * @param k - number of initial clusters
 * @param arraySize - amount of datapoints
 * @param U - matrix to be normalized by rows
 */
void normalizeU(int k, int arraySize, double ***U){
    int i, j;
    double valueOfRow;

    for(i = 0; i < arraySize; i++){
        valueOfRow = 0;
        for(j = 0; j < k; j++){
            valueOfRow += pow((*U)[i][j], 2);
        }

        valueOfRow = sqrt(valueOfRow);
        /*if(valueOfRow == 0){ // in case we divide in 0
            continue;
        }
         */
        if(valueOfRow != 0){
            for (j = 0; j < k; j++) {
                (*U)[i][j] /= valueOfRow;   // T(i,j) = U(i,j) / valueOfRow
            }
        }
    }
}


/********************** Testers **********************/


/**
 * Tester to print n*m Matrix
 */
void printTest(double **matrix, int n, int m){
    int i, j;
    for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
}

/** sortJacobi Tester */
void TesterToSortJacobi(){
    double **matrix;

    /** initialization */
    matrix = randomMatrix(9,8);

    /** Print the random matrix */
    printf("\nThe Random Matrix is:\n");
    printTest(matrix, 9, 8);

    /** sort the matrix by the first line */
//    sortJacobi(8, &matrix);

    /** print the sorted Matrix */
    printf("\nThe Sorted Matrix is:\n");
    printTest(matrix, 9, 8);

}

/** the method calculate weighted Adj Matrix on RANDOM Matrix */
double **TesterToWeight(){
    double **matrix, **weightedAdjMatrix;
    int d = 3, arraySize = 6;

    matrix = randomMatrix(arraySize, d);
    printf("\nRandom Matrix:\n");
    printTest(matrix, arraySize, d);

    weightedAdjMatrix = calcWeightedAdjMatrix(d, arraySize, &matrix);
    printf("\nWeighted Adj Matrix:\n");
    printTest(weightedAdjMatrix, arraySize, arraySize);

    return weightedAdjMatrix;
}



/** creates random n*m matrix */
double **randomMatrix(int n, int m) {

    int i, j;
    double **matrix = createMatrix(n, m);

    for(i = 0; i < n; i++){
        for(j = 0; j < m; j++){
            matrix[i][j] = rand();
        }
    }
    return matrix;
}

void print1Darray(int *array, int arraySize){

    int i;

    for(i = 0; i < arraySize; i++){
        printf("%d,", array[i]);
    }
}

