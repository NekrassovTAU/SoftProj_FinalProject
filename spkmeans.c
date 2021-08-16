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

#define BEEP(num) \
    printf("%s %d%s", "SPKMEANS Check #" , num, "\n");

/**
 * main, a shell function for the spectral clustering algorithm implementation
 */
int main(int argc, char *argv[]) {
    checkArgs(argc, argv, 0, NULL, NULL);
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

void determineRowAndCol(enum goalEnum goal, int k, int arraySize, int *rowCount,
                        int *colCount);

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
    ASSERT_ERROR(dArraySizeInfo != NULL)


    /** process the datapoint file, and return info about their amount (arraySize)
     * and number of properties (d) */
    processDatapoints(origArgv[3], &datap_array, &datapoint, &dArraySizeInfo);

    d = dArraySizeInfo[0];
    arraySize = dArraySizeInfo[1];
    free(dArraySizeInfo);

    ASSERT_ARGS(k < arraySize)

    BEEP(1)

    /** initialize the process to achieve the provided goal */
    ret_matrix = initProcess(&k, d, arraySize, goal, &datapoint, isCAPI);

    BEEP(2)

    /** Determining row and col size for return matrix */
    determineRowAndCol(goal, k, arraySize, returnRowCount, returnColCount);

    /** Free allocated memory and terminate*/
    free(datapoint);
    free(datap_array);

    return ret_matrix;
}

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

    *datapoint = arrayToTwoDimArray(datap_array, arraySize, d);

    (*dArraySizeInfo)[0] = d;
    (*dArraySizeInfo)[1] = arraySize;
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
        printf("\n*\n*\n*array size is : %d\n*\n*\n*\n*\n***", arraySize);
        printResult(&ret_matrix, goal, arraySize);
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
    BEEP(3)

    spkMatrix = calcSpectralClusters(k, arraySize, 0, &jacobiMatrix);

    BEEP(4)

    freeMatrix(&weightedAdjMatrix);
//    free(diagDegMatrix);
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
void printResult(double ***retArray, enum goalEnum goal, int arraySize) {
    /* TODO: implement function */
    printf("check2");

    if((goal == wam) || (goal == lnorm)){
        printRegular(retArray, arraySize);
    }
    if(goal == jacobi){
        printJacobi(retArray, arraySize);
    }
    if(goal == ddg) {
        printDiagonal(retArray, arraySize);
    }
    if(goal == spk){
        // TODO
    }
}

/** regular (rows * columns) Matrix */
void printRegular(double ***retArray, int arraySize){

    int i, j;

    for(i = 0; i < arraySize; i++){
        for(j = 0; j < arraySize; j++){
            printf("%.4f", (*retArray)[i][j]);
            if (j != arraySize - 1){   /* not last in the line */
                printf("%s", ",");
            }
            else{
                printf("\n");
            }
        }
    }
}

/** Transporse the matrix (beside the first line) while printing */
void printJacobi(double ***retArray, int arraySize) {
    printf("check3");
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


/** print the diagonal matrix - fill with 0's */
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

    value = exp(sqrt(value) / (-2)); // currently exp sets value to 0 - work it out

    return value;
}

/**
 * Takes the weightedAdjMatrix, and returns the diagonal degree matrix
 * corresponding to it
 * @param arraySize - amount of datapoints
 * @param weightedAdjMatrix - 2-D array containing all the weights

 */
 /** 2-D array */
/* double ** calcDiagDegMatrix(int arraySize, double ***weightedAdjMatrix) {

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
*/

double **calcDiagDegMatrix(int arraySize, double ***weightedAdjMatrix) {

    int i, j;
    double sumofweights, **diagDegMatrix;

    diagDegMatrix = createMatrix(1, arraySize);

    for (i = 0; i < arraySize; i++) {
        sumofweights = 0;
        for (j = 0; j < arraySize; j++) {
            sumofweights += (*weightedAdjMatrix)[i][j];
        }
        diagDegMatrix[0][i] = sumofweights;
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
 * Takes a matrix and returns the result of running the Jacobi algorithm on it.
 * @param arraySize - size of the nXn matrix
 * @param inputMatrix - pointer to 2D array representation of input matrix
 * @return The jacobi matrix
 */
double ** calcJacobi(int arraySize, double ***inputMatrix) {

    int row, col, converge, iterations = 0;
    double c, s, epsilon = EPSILON;
    double offA, offAtag;
    double **A, **Atag, **V, **jacobiMatrix;

    /** initialization */
    Atag = createMatrix(arraySize, arraySize);
    V = createMatrix(arraySize, arraySize);

    makeIntoIdentityMatrix(&V, arraySize);

    /** run until convergence -  */
    A = (*inputMatrix);
    copymatrix(arraySize, &A, &Atag); // Atag = A
    offAtag = calcOff(arraySize, &Atag);

    printTest(*inputMatrix, arraySize, arraySize);
    do {
        offA = offAtag;

        findmatrixP(arraySize, &A, &c, &s, &row, &col); // P
        printf("\n**** c = %f, s = %f, row = %d, col = %d\n", c,s,row,col);
        updateAtag(arraySize, &Atag, &A, c, s, row, col); // A' = P^T * A * P
        printf("Atag is :\n");
        printTest(Atag, arraySize,arraySize);
        updateV(arraySize, &V, c, s, row, col);  // V *= P
        printf("\nV is :\n");
        printTest(V, arraySize,arraySize);

        /** instead of calculating both just calculate one, and update the sec at the first line */
        offAtag = calcOff(arraySize, &Atag);
        converge = (offA - offAtag) < epsilon ? 1 : 0;
       // converge = convergenceCheck(arraySize, &A, &Atag);
        iterations++;
        if(!converge && (iterations < 100)) {

            copymatrix(arraySize, &Atag, &A); // A = Atag
        }
     } while (!converge && (iterations < 100)); // as long as (delta > epsilon) or number of iterations is under 100

    /** Atag has the A"A
      *  V has the V"A */
    printf("\n\n eigenValues are:\n\n");
    printTest(Atag,arraySize,arraySize);
    printf("\n\n eigenVectors are:\n\n");
    printTest(V,arraySize,arraySize);
    jacobiMatrix = copyJacoby(arraySize, &Atag, &V);

    freeMatrix(&Atag);
    /** TODO: FIX NEXT LINE. PROBLEM IN FREE V - when i run threw terminal it doesn't get to the return */
    // freeMatrix(&V);

    return jacobiMatrix;
}


/**
 * Copies the eigenValues and eigenVectors into the final jacobiMatrix and returns it
 * @param arraySize - amount of datapoints
 * @param Atag - pointer to 2-D array containing the eigenValues on the diag
 * @param diagDegMatrix - pointer to 2-D array containing the eigenVectors as columns
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
            if ((*A)[i][j] > curmax) {
                curmax = (*A)[i][j];
                (*row) = i;
                (*col) = j;
            }
        }
    }
}

/**
 *Calculate the matrix Atag using the current A and P (based on c&s) : A' = P^T * A * P.
 * @param arraySize - amount of datapoints
 * @param Atag - 2-D array to be calculated
 * @param c - the value of P[row][row]
 * @param s - the value of P[row][col]
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
calcSpectralClusters(int *k, int arraySize, int isCAPI,
                     double ***jacobiMatrix) {

    double **U, **T, **spkMatrix, **combined = createMatrix(2, arraySize);;

    /** sort jacobi's eigenValues*/
    sortJacobi(arraySize, jacobiMatrix, &combined);

    /** determine k (if k==0) AND change k variable accordingly, important for CAPI */
    if(*k == 0){
        *k = getK(arraySize, jacobiMatrix);
    }
    /** sort jacobi's first K eigenVectors into U */
    /** get k first eigenvectors from jacobi into U, normalize U */
    U = getMatrixU(*k, arraySize, jacobiMatrix, combined);
    normalizeU(*k, arraySize, &U);
    T = U;

    //TODO: K-Means algorithm - if CAPI - need to return to python for K-Means++
    //TODO:                             - need to get back result from python
    //TODO: Need to cut off function when returning to python, then start new
    //TODO: method with KMeans implementation and finito

    //TODO: free memory and finish

    return T;
}

double ** KMeansAlgorithm(int k, int *d, int arraySize, double ***datapoints,
                          int isCAPI, int **init_centroids){
    return NULL;
}


/**
 * Sort the jacobiMatrix eigenValues by ascending order, and the eigenVectors respectively
 * @param arraySize - amount of datapoints
 * @param jacobiMatrix - jacobi matrix which includes eigenValues and EigenVectors
 */
 /*
void sortJacobi(int arraySize, double ***jacobiMatrix) {
    int i;
    int *newIndex = (int *) calloc(arraySize, sizeof(int));        // keeps the new index of each eigenValue
    double *eigenValuesSorted = calloc(arraySize, sizeof(double)); // temp array to keep the eigenvalues sorted

    for (i = 0; i < arraySize; i++){
        newIndex[i] = -1;
    }

    mergeSort(*jacobiMatrix, &eigenValuesSorted, &newIndex, 0, arraySize - 1); // sort the eigen values

    printf("\n The newIndex 1-D array: \n ");
    print1Darray(newIndex, arraySize);
    printf("\neigenValues sorted: \n");
    printTest(*jacobiMatrix, arraySize+1, arraySize);

    sortEigenVectors(arraySize, jacobiMatrix, newIndex);                                // sort the eigenVectors respectively

    free(newIndex);
    free(eigenValuesSorted);
}
*/

/**
 * sort the eigenValues (first row) of the jacobiMatrix, and keeps the have been changed
 * @param jacobiMatrix - jacobi matrix which includes eigenValues and EigenVectors
 * @param eigenValuesSorted - 2-D array to be calculated
 * @param newIndex - pointer to 1-D array for keeping the new indices of the eigenValues that changed place
 * @param low - the index of the leftMost eigenValue to be sorted
 * @param high - the index of the rightMost eigenValue to be sorted
 */
 /*
void mergeSort(double **jacobiMatrix, double **eigenValuesSorted, int **newIndex, int low, int high){
    int mid;

    if(low < high){
        mid = (low + high) / 2;
        mergeSort(jacobiMatrix, eigenValuesSorted, newIndex, low, mid);
        mergeSort(jacobiMatrix, eigenValuesSorted, newIndex, mid + 1, high);
        merge(jacobiMatrix, eigenValuesSorted, newIndex, low, mid, high);
    }
    else{
        return;
    }
}
*/
/**
 *  Merge between two soretd slices - (low,mid) & (mid+1, high)
 * @param jacobiMatrix - jacobi matrix which includes eigenValues and EigenVectors
 * @param eigenValuesSorted - pointer to 1-D array for keeping the eigenValues sorted
 * @param newIndex - pointer to 1-D array for keeping the new indices of the eigenValues that changed place
 * @param low - the index of the leftMost eigenValue in the "1'st array" to be merged
 * @param mid - the index of the right most eigenValue in the "1'st array" to be merged
 * @param high - the index of the right most eigenValue in the "2'nd array" to be merged
 */
 /*
void merge(double **jacobiMatrix, double **eigenValuesSorted, int **newIndex, int low, int mid, int high){
    int l1, l2, i;

    for(l1 = low, l2 = mid+1, i = low; l1 <= mid && l2 <= high; i++) {

        if ((*jacobiMatrix)[l1] <= (*jacobiMatrix)[l2]) {

            (*eigenValuesSorted)[i] = (*jacobiMatrix)[l1];
            if ((i == l1 && (*newIndex)[l1] != -1) || i != l1) { // i=l1 but we update the index not for the first time, or i!=l1 so update it anyway
                (*newIndex)[l1] = i;
            }
            l1++;
        }
        else {

            (*eigenValuesSorted)[i] = (*jacobiMatrix)[l2];
            if ((i == l2 && (*newIndex)[l2] != -1) || i != l2) { // i=l2 but we update the index not for the same time, or i!=l2 so update it anyway
                (*newIndex)[l2] = i;
            }
            l2++;
        }
    }

    while(l1 <= mid){  // copy the rest of the "1'st array"
        (*eigenValuesSorted)[i] = (*jacobiMatrix)[l1];
        if ((i == l1 && (*newIndex)[l1] != -1) || i != l1) { // i=l1 but we update the index not for the same time, or i!=l1 so update it anyway
            (*newIndex)[l1] = i;
        }
        i++;
        l1++;
    }

    while(l2 <= high){ // copy the rest of the "2'nd array"
        (*eigenValuesSorted)[i] = (*jacobiMatrix)[l2];
        if ((i == l2 && (*newIndex)[l2] != -1) || i != l2) {  // i=l2 but we update the index not for the same time, or i!=l2 so update it anyway
            (*newIndex)[l2] = i;
        }
        i++;
        l2++;
    }

    for(i = low; i <= high; i++){    // copy the eigenValues from the temp array, back to the original jacobiMatrix
        (*jacobiMatrix)[i] = (*eigenValuesSorted)[i];
    }
}
*/
/*********************************************************************************************************/


void sortJacobi(int arraySize, double ***jacobiMatrix, double ***combined) {
    int i;

//    double **combined = createMatrix(2, arraySize);
    double **tmp = createMatrix(2, arraySize);
    for (i = 0; i < arraySize; i++){
        (*combined)[0][i] = (*jacobiMatrix)[0][i];
        (*combined)[1][i] = i;
    }

    mergeSort(*jacobiMatrix, combined, &tmp, 0, arraySize - 1); // sort the eigen values

    /** Alternative way : sort only the first K eigenVectors (so find K first)
    sortEigenVectors(arraySize, jacobiMatrix, combined);                                // sort the eigenVectors respectively
    */
}



void mergeSort(double **jacobiMatrix, double ***combined, double ***tmp, int low, int high){
    int mid;

    if(low < high){
        mid = (low + high) / 2;
        mergeSort(jacobiMatrix, combined, tmp, low, mid);
        mergeSort(jacobiMatrix, combined, tmp, mid + 1, high);
        merge(jacobiMatrix, combined, tmp, low, mid, high);
    }
    else{
        return;
    }
}



void merge(double **jacobiMatrix, double ***combined, double ***tmp, int low, int mid, int high){
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

/*********************************************************************************************************/

/**
 *  Sort the eigenVectors (columns) of the jacobiMatrix respectively to the eigenValues
 * @param arraySize - amount of datapoints
 * @param jacobiMatrix - jacobi matrix which includes eigenValues and EigenVectors
 * @param newIndex - 1-D array with che changes in eigenValues indices
 */
 /*
void sortEigenVectors(int arraySize, double ***jacobiMatrix, int *newIndex){
    int i, j;
    double **tmpMatrix = createMatrix(arraySize, arraySize);

    for(j = 0; j < arraySize; j++){                 // use temp matrix to save the eigenVectors that their eigenValues changed place earlier
        if(newIndex[j] != -1){     // switch
           for(i = 0; i < arraySize; i++){
               tmpMatrix[i][newIndex[j]] = (*jacobiMatrix)[i+1][j];
           }
        }
    }

    // could change the pointer instead
    for(j = 0; j < arraySize; j++){   // update the places of the eigenvectors in the jacobiMatrix.
        if(newIndex[j] != -1){ // the current column changed place during the sort
            for(i = 0; i < arraySize; i++){
                (*jacobiMatrix)[i+1][j] = tmpMatrix[i][j];
            }
        }
    }

    freeMatrix(&tmpMatrix);
}
*/

/**
 * Returns K based on the biggest difference between two adjacent eigenValues
 * @param arraySize - amount of datapoints
 * @param jacobiMatrix - jacobi matrix which includes eigenValues and EigenVectors
 */
int getK(int arraySize, double ***jacobiMatrix) {

    int j, k = 0;
    double cur, max = (-1);
    for(j = 0; j < arraySize / 2; j++) {
        cur = fabs((*jacobiMatrix)[0][j] - (*jacobiMatrix)[0][j+1]);
        if(cur > max){
            max = cur;
            k = j + 1;
        }
    }
    return k;
}


/**
 * Returns matrix U (n*k), which is the first K columns (eigenVectors) in jacobiMatrix
 * @param k - number of initial clusters
 * @param arraySize - amount of datapoints
 * @param jacobiMatrix - jacobi matrix which includes eigenValues and EigenVectors
 */
double **getMatrixU(int k, int arraySize, double ***jacobiMatrix, double **combined){

    int i, j, col;
    double **U = createMatrix(arraySize, k);

    /** take the first K eigenVectors if the jacobi was sorted by the eigenValue */
    for(j = 0; j < k; j++){
        col = combined[1][j]; // The original index of the j eigenValue is sorted order
        for(i = 0; i < arraySize; i++){
            U[i][j] = (*jacobiMatrix)[i+1][col];
        }
    }
    /** JacobiMatrix is sorted
    for(i = 0; i < arraySize; i++){
        for(j = 0; j < k; j++){
            U[i][j] = (*jacobiMatrix)[i+1][j];
        }
    }
     */
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
        for(j = 0; j < arraySize; j++) {
            (*U)[i][j] /= valueOfRow;   // T(i,j) = U(i,j) / valueOfRow
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

