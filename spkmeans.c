/**
 * Authors: Daniel Nekrassov, Lior Grinberg
 *
 * Note: Throughout comments, we use 2D-array and matrix interchangeably
 * */



#include "spkmeans.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>


#define FLT_MAX 3.402823466e+38F
#define FLT_MIN 1.1754943508e-38F

#define BUFFER_SIZE 1000
#define INITIAL_ARRAY_SIZE 100
#define EPSILON 1.0e-15F
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


/**
 * main, a shell function for the spectral clustering algorithm implementation
 * Input:
 * k - Number of required clusters. If equal 0, we use the heuristic
 * goal - spk / wam / ddg / lnorm / jacobi
 * file_name - path to file containing N observations (.txt/.csv format)
 */
int main(int argc, char *argv[]) {
    int rowCount, colCount; /* dummy variables, relevant only for C-API */
    startCSPKMeansExec(argc, argv, 0, &rowCount, &colCount);
    return 0;
}

/**
 * Checks the validity of arguments passed, reads and parses info in files
 * and starts the implementation of spectral k-means clustering
 * @param argc - # of arguments passed to main
 * @param origArgv - the arguments passed to main
 * @param isCAPI - used in order to distinguish between
 *                 C-based and C-API-based call to this function
 * @param returnRowCount - returns to C-API the number of rows in the
 *                         result matrix
 * @param returnColCount - returns to C-API the number of columns in the
 *                         result matrix
 * @return returns the result matrix to C-API for printing / K-Means++
 */
double **startCSPKMeansExec(int argc, char **origArgv, int isCAPI,
                            int *returnRowCount, int *returnColCount) {

    /** Initialization & Argument Checking*/
    int k, d, arraySize;
    double **ret_matrix, **datapoints;
    enum goalEnum goal;
    char *ptr;

    /*Check if we received right amount of arguments*/
    ASSERT_ARGS(argc > 3)

    /*Check if goal is valid*/
    goal = checkGoal(origArgv[2]);
    ASSERT_ARGS(goal != INVALID)

    /*check k only in case of goal = spk, since it is only needed then*/
    if (goal == spk){
        k = strtol(origArgv[1], &ptr, 10);

        /** checking whether k>=0,
         * strcmp with "0" used for the case which strtol fails and returns 0*/
        ASSERT_ARGS(k > 0 || !strcmp(origArgv[1], "0"))
    }
    /*if goal != spk, we do not care about the value of k
     * not need to check argument due to 2.7.12 in document*/
    else{
        k = 0;
    }

    /** process the datapoint file, and return info about their amount
     * (arraySize) and number of properties (d) */
    processDatapoints(origArgv[3], &datapoints, &d,
                      &arraySize);

    if (goal == spk){
        ASSERT_ARGS(k < arraySize)
    }

    /** initialize the process to achieve the provided goal */
    ret_matrix = initProcess(&k, d, arraySize, goal, &datapoints, isCAPI);

    /** Determining row and col size for return matrix */
    determineRowAndCol(goal, k, arraySize, returnRowCount,
                       returnColCount, isCAPI);

    /** Free memory and terminate*/
    freeMatrix(&datapoints);


    return ret_matrix;
}

/**
 * Checks goal string and returns enum accordingly
 * @param string - represents goal of the call to main
 * @return appropriate enum for goal, or INVALID if invalid
 */
enum goalEnum checkGoal(char *string) {
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
 * Determines the size of the return matrix according to goal and other metrics
 * @param goal - goal of call to algorithm
 * @param k - amount of clusters (used in case of goal == spk)
 * @param arraySize - amount of datapoints
 * @param rowCount - the amount of rows in the return matrix (updated in this
 *                   method).
 * @param colCount - the amount of columns in the return matrix (updated in this
 *                   method).
 * @param isCAPI - used in order to distinguish between
 *                 C-based and C-API-based call to this function
 */
void determineRowAndCol(enum goalEnum goal, int k, int arraySize, int *rowCount,
                        int *colCount, int isCAPI) {
    switch (goal) {
        case jacobi:{
            *rowCount = arraySize + 1;
            *colCount = arraySize;
            break;
        }
        case spk:{
            if (isCAPI){
                *rowCount = arraySize;
            }
            else{
                *rowCount = k;
            }
            *colCount = k;
            break;
        }
        case ddg:{
            *rowCount = 1;
            *colCount = arraySize;
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
 * @param d - pointer to int-type parameter used to save the dimension of
 *            the datapoints
 * @param arraySize - pointer to int-type parameter used to save the amount of
 *            datapoints in the file
 */
void processDatapoints(char *filename, double ***datapoint, int *d,
                       int *arraySize) {

    double *datap_array;
    char *token, *ptr, currRow[BUFFER_SIZE];
    FILE *fptr;
    int i, tail = 0;

    /*opens file and reads first line to determine d*/
    fptr = fopen(filename, "r");
    ASSERT_ERROR(fptr != NULL)

    fgets(currRow, BUFFER_SIZE, fptr);
    ASSERT_ERROR(currRow != NULL)

    for (i = 0, *d = 1; currRow[i]; i++) {
        *d += (currRow[i] == ',');
    }

    *arraySize = INITIAL_ARRAY_SIZE;

    datap_array = calloc((*d) * (*arraySize), sizeof(double));
    ASSERT_ERROR(datap_array != NULL)

    /** Populates datap_array with info from the file */
    do {/* as long as there are still lines to be read */
        /*loop to read line */
        token = strtok(currRow, ","); /* read the current line */
        while (token != NULL) {
            datap_array[tail] = strtod(token, &ptr);
            token = strtok(NULL, ",");
            tail++;

            /*in-case we reached the edge of the current datap_array,
             * increase the arraySize by twice */
            if (tail == (*d) * (*arraySize)) {
                *arraySize *= 2;
                datap_array = realloc(datap_array,
                                      (*d) * (*arraySize) * sizeof(double));
                ASSERT_ERROR(datap_array != NULL)
            }
        }
    } while (fgets(currRow, BUFFER_SIZE, fptr));

    fclose(fptr);

    /* cut the datap_array to its intended arraySize */
    if (tail < (*d) * (*arraySize) - 1) {
        datap_array = realloc(datap_array, tail * sizeof(double));
        ASSERT_ERROR(datap_array != NULL)
        *arraySize = tail / (*d);
    }

    free(token);

    /* Copying datap_array to a matrix-format for ease of use */
    copyDatapoint(datapoint, &datap_array, *arraySize, *d);

    free(datap_array);
}

/**
 * Copies datapoint data saved in a 1D-array to a 2D-array ("matrix").
 * @param datapoint - 2D array to which we will input the datapoints info
 * @param datap_array - 1D array in which we inputed info in the first-place
 * @param arraySize - amount of datapoints
 * @param d - dimension of datapoints
 */
void
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
 * Shell function for goalBasedProcess, either prints its results (if C)
 *                                      or returns it to C-API
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
 * @return returns the result matrix to C-API for printing / K-Means++
 */
double **
initProcess(int *k, int d, int arraySize, enum goalEnum goal,
            double ***datapoint, int isCAPI) {

    double **ret_matrix;

    /** Calculates needed results according to goal*/

    ret_matrix = goalBasedProcess(k, d, arraySize, datapoint, goal, isCAPI);

    /** Freeing memory and returning/printing the result */

    /*if C, we need to print the result, else we return in to CAPI*/
    if (!isCAPI) {
        printResult(&ret_matrix, goal, arraySize, *k);
        freeMatrix(&ret_matrix);
        return NULL;
    } else {
        /*Will be freed by C-API*/
        return ret_matrix;
    }
}

/**
 * Prints the result of the algorithm according to goal
 * wam / ddg / lnorm - print a nXn matrix
 * jacobi - first line is eigenvalues,
 *          next are nXn eigenvecor matrix
 * spk - k centroids after k-means algorithm
 * @param ret_matrix - pointer to matrix to be printed to user
 * @param goal - distinguishes between spk-printing and matrix-printing
 */
void printResult(double ***ret_matrix, enum goalEnum goal, int arraySize,
                 int k) {

    int rows, cols;

    determineRowAndCol(goal, k, arraySize, &rows, &cols, 0);
    fixZeros(ret_matrix, rows, cols);

    if(goal == jacobi){
        printJacobi(ret_matrix, arraySize);
    }
    else if(goal == ddg) {
        printDiagonal(ret_matrix, arraySize);
    }
    else {
        printRegular(ret_matrix, rows, cols);
    }

}

/**
 * Replace all the values in the matrix which are in the range of
 * (-0.00005 - 0) with the value 0.0
 * @param matix - pointer to 2-D array to be updated
 * @param rows - number of rows of the given matrix
 * @param cols - number of cols of the given matrix
 */

void fixZeros(double ***matrix, int rows, int cols) {

    int i, j;

    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            if((*matrix)[i][j] < 0 && (*matrix)[i][j] > -0.00005){
                (*matrix)[i][j] = 0.0f;
            }
        }
    }
}

/**
 * Prints matrix of size rows * cols
 * @param ret_matrix - pointer to matrix to be printed to user
 * @param rows - number of rows in the matrix
 * @param cols - number of columns in the matrix
 */
void printRegular(double ***ret_matrix, int rows, int cols){

    int i, j;

    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            printf("%.4f", (*ret_matrix)[i][j]);
            if (j != cols - 1){   /* not last in the line */
                printf("%s", ",");
            }
            else if (i != rows - 1){
                printf("\n");
            }
        }
    }
}


/**
 * Prints the Jacobi matrix with eigenValues & eigenVecrtors as lines
 * @param ret_matrix - pointer to matrix to be printed to user
 * @param arraySize - number of rows in the matrix
 */
void printJacobi(double ***ret_matrix, int arraySize) {

    int i, j;

    /** print the eigenValues - first line */
    for (j = 0; j < arraySize; j++) {
        printf("%.4f", (**ret_matrix)[j]);
        if (j != arraySize - 1) {   /* not last component of the cluster*/
            printf("%s", ",");
        } else {
            printf("\n");
        }
    }

    /** print the eigenVectors - each vector as a row */
    for (j = 0; j < arraySize; j++) {
        for (i = 1; i < arraySize + 1; i++) {
            printf("%.4f", (*ret_matrix)[i][j]);
            if (i != arraySize) {
                printf("%s", ",");
            } else if (i != arraySize - 1){
                printf("\n");
            }
        }
    }
}

/**
 * Prints the diagonal matrix, filled with 0's
 * @param ret_matrix - pointer to matrix to be printed to user
 * @param arraySize - size of the given matrix
 */
void printDiagonal(double ***ret_matrix, int arraySize) {

    int i, j;

    for(i = 0; i < arraySize; i++){
        for(j = 0; j < arraySize; j++){
            if(i == j){ /* on the diagonal */
                printf("%.4f", (**ret_matrix)[i]);
            }
            else{
                printf("0.0000");
            }
            if (j != arraySize - 1){   /* not last component of the cluster */
                printf("%s", ",");
            }
            else if (i != arraySize - 1){
                printf("\n");
            }
        }
    }
}

/**
 * Takes the datapoint info, and returns relevant matrix according to
 * the goal, later will be converted to matrix form according to needs.
 * @param k - number of clusters to assign.if 0, we need to use EigGap Heuristic
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

    spkMatrix = calcSpectralClusters(k, arraySize, isCAPI, &jacobiMatrix);

    freeMatrix(&weightedAdjMatrix);
    freeMatrix(&diagDegMatrix);
    freeMatrix(&normLaplacianMatrix);
    freeMatrix(&jacobiMatrix);

    return spkMatrix;
}

/**
 * Returns a 2D matrix of size n X m
 * @param n - # of rows
 * @param m - # of cols
 */
double ** createMatrix(int n, int m) {
    double *array, **twoDimArray;
    int i;

    array = calloc(n*m, sizeof (double ));
    ASSERT_ERROR(array != NULL)

    twoDimArray = calloc(n, sizeof(double *));
    ASSERT_ERROR(twoDimArray != NULL)

    for (i = 0; i < n; i++) {
        twoDimArray[i] = array + i * m;
    }

    return twoDimArray;
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
    free((*matrix)[0]);
    free(*matrix);
}

/**
 * Takes the datapoint info, and returns the weighted adjacency matrix.
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
            /* calculate the weight between the i,j vectors */
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
 * @param datapoint - pointer to 2D-array containing all datapoints
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
 * Takes the weighted matrix, and returns representation of the diagonal
 * degree matrix
 * @param arraySize - amount of datapoints
 * @param weightedAdjMatrix - pointer to 2D-array containing all the weights
 * @return 2-D array of size (1 * arraySize) represents the diag of the diagonal
 *         degree matrix
 */
double **calcDiagDegMatrix(int arraySize, double ***weightedAdjMatrix) {

    int i, j;
    double sumOfWeights, **diagDegMatrix;

    diagDegMatrix = createMatrix(1, arraySize);

    /** calculation of the values on the diagonal */
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

    /** calculate the D^(-0.5) matrix */
    for (i = 0; i < arraySize; i++) {
        (*diagDegMatrix)[0][i] = 1 / sqrt((*diagDegMatrix)[0][i]);
    }

    /** calculating the normLaplacianMatrix using the next formula :
     *  sqrtDegMatrix[i][j] = sqrtDeg[i] * weight[i][j] * sqrtDeg[j]
     *                       == sqrtDegMatrix[j][i]                    */
    for (i = 0; i < arraySize; i++) {
        (normLaplacianMatrix)[i][i] = 1; /* Lnorm = I - sqrtDegMatrix */
        for (j = i + 1; j < arraySize; j++) {

            value = -(*diagDegMatrix)[0][i] * (*weightedAdjMatrix)[i][j] *
                    (*diagDegMatrix)[0][j];
            (normLaplacianMatrix)[i][j] = value;
            (normLaplacianMatrix)[j][i] = value;
        }
    }

    return normLaplacianMatrix;
}

/**
 * Takes a symmetrical matrix, and returns the result of running the
 * Jacobi algorithm on it.
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

    A = (*inputMatrix);
    offAtag = calcOff(arraySize, &A); /* calculate off(A)^2 */

    /** run until convergence */
    do {
        offA = offAtag;

        /** calculate P by finding the indices of the max absolute value
         * and c,s */
        findmatrixP(arraySize, &A, &c, &s, &row, &col);

        /** calculate Atag: A' = P^T * A * P and update off(A') accordingly */
        updateAtag(arraySize, &A, c, s, row, col, &offAtag);

        /** calculate V : V *= P */
        updateV(arraySize, &V, c, s, row, col);

        /** check convergence due to delta which is smaller/equal to epsilon */
        converge = (offA - offAtag) <= epsilon ? 1 : 0;
        iterations++;

        /* as long as (delta > epsilon) or number of iterations is under 100 */
     } while (!converge && (iterations < 100));

    /** create the jacobiMatrix as mix of two matrices :
     * A has the eigenvalues & V has the eigenvectors */
    jacobiMatrix = copyJacoby(arraySize, &A, &V);

    freeMatrix(&V);
    /* A is either normLaplacian or datapoints, which are freed later on
     * according to needs */

    return jacobiMatrix;
}


/**
 * Copies the eigenValues and eigenVectors into the final jacobiMatrix
 * and returns it
 * @param arraySize - amount of datapoints
 * @param Atag - pointer to 2-D array containing the eigenValues on the diag
 * @param V - pointer to 2-D array containing the eigenVectors as columns
 * @return The jacobi matrix
 */
double **copyJacoby(int arraySize, double ***Atag, double ***V) {

    int i, j;
    double **jacobiMatrix;

    jacobiMatrix = createMatrix(arraySize+1, arraySize);

    /** copy the EigenValues */
    for(j = 0; j < arraySize; j++) {
        jacobiMatrix[0][j] = (*Atag)[j][j];
    }
    /** copy the EigenVectors */
    for(i = 1; i < arraySize + 1; i++) {
        for(j = 0; j < arraySize; j++){
            jacobiMatrix[i][j] = (*V)[i-1][j];
        }
    }

    return jacobiMatrix;
}

/**
 * Find the sum of squares of all off-diagonal elements of the given matrix
 * @param arraySize - amount of datapoints
 * @param matrix - 2-D symmetric array
 * @return the calculated off(*matrix)^2
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
 * Calculate the paramaters c and s in the matrix P using the matrix A
 * @param arraySize - amount of datapoints
 * @param A - 2-D symmetric array
 * @param c - pointer to the value in the index (row,row) in matrix P
 * @param s - pointer to the value in the index (row,col) in matrix P
 * @param row - pointer to the row index of the max value in A (off diag),
 *              to be calculated
 * @param col - pointer to the col index of the max value in A (off diag),
 *              to be calculated
 */
void findmatrixP(int arraySize, double ***A, double *c, double *s,
                 int *row, int *col) {

    double teta, t, sign = -1;

    /* find indices of the max absolute value in the matrix A */
    findMaxOffDiag(arraySize, A, row, col);

    /** update *c,*s using given formulas */
    teta = ((*A)[*col][*col] - (*A)[*row][*row]) / (2 * (*A)[*row][*col]);

    if (teta >= 0) {
        sign = 1;
    }
    t = sign / (fabs(teta) + sqrt(pow(teta, 2) + 1));
    (*c) = 1 / sqrt(pow(t, 2) + 1);
    (*s) = t * (*c);
}

/**
 * Find the index of the off-diag maximum absolute value in A.
 * @param arraySize - amount of datapoints
 * @param A - pointer to 2-D symmetric array
 * @param row - pointer to the row index of the max value in A (off diag),
 *              to be calculated
 * @param col - pointer to the col index of the max value in A (off diag),
 *              to be calculated
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
 * Calculate the matrix Atag using the current A and P (based on c&s) :
 * (A' = P^T * A * P)  & update the offAtag accordingly to the changes made
 * @param arraySize - amount of datapoints
 * @param A - pointer 2-D array to be calculated based on the given one
 * @param c - the value of P[row][row]
 * @param s - the value of P[row][col]
 * @param row - row index of the max value in Atag (off diag)
 * @param col - col index of the max value in Atag (off diag)
 * @param offAtag - pointer to the current off(A) value
 */
void updateAtag(int arraySize, double ***A, double c, double s,
                int row, int col, double *offAtag) {

    int r;
    double ri, rj, ij, ii, jj;

    /** for (r != i,j) update A'[r][i] & A'[r][j]
     * can cut computation in half because the matrix is symmetric */
    for (r = 0; r < arraySize; r++) {

        if ((r != row) && (r != col)) {
            ri = (*A)[r][row];
            rj = (*A)[r][col];
            (*A)[r][row] = (c * ri) - (s * rj);
            (*A)[row][r] = (*A)[r][row];
            (*A)[r][col] = (c * rj) + (s * ri);
            (*A)[col][r] = (*A)[r][col];
            *offAtag += 2 * ((*A)[r][row] * (*A)[r][row] +
                            (*A)[r][col] * (*A)[r][col]) -
                        2 * (ri * ri + rj * rj);
        }
    }
    /*  A[i][i] / A[j][j]/ A[i][j] == A[j][i] */
    ij = (*A)[row][col];
    ii = (*A)[row][row];
    jj = (*A)[col][col];

    /* update A' using given formulas */
    (*A)[row][row] = (c * c * ii) + (s * s * jj) - (2 * s * c * ij);
    (*A)[col][col] = (s * s * ii) + (c * c * jj) + (2 * s * c * ij);
    (*A)[row][col] = 0;
    (*A)[col][row] = 0;

    *offAtag -= 2 * (ij * ij);
}

/**
 *Calculate the current matrix V, which is the Eigenvectors matrix.
 * V *= P (based on c & s)
 * @param arraySize - amount of datapoints
 * @param V - pointer 2-D array to be calculated
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

        /* V'[r][row] = (V[r][row] * P[row][row]) + (V[r][col] * P[col][row]) ==
         * V'[r][i] = (V[r][i] * c) + (V[r][j] * (-s)) */
        (*V)[r][row] = (ri * c) + (rj * -s);

        /* V'[r][col] = (V[r][row] * P[row][col]) + (V[r][col] * P[col][col]) ==
         * V'[r][j] = (V[r][i] * s) + (V[r][j] * c) */
        (*V)[r][col] = (ri * s) + (rj * c);
    }
}

/** Calculate the spk matrix using KMeans & returns it if !isCAPI,
 * else - returns the matrix T for KMeans++
* @param k - pointer to the number of clusters
* @param arraySize - amount of datapoints
* @param isCAPI - 1 if running in Python, 0 if running in C
* @param jacobiMatrix - pointer to the JacobiMatrix calculated before
* @return if (isCAPI == 1) - The spk matrix after running KMeans in C,
*         else - returns the matrix T for running KMeans++ in python */
double **calcSpectralClusters(int *k, int arraySize, int isCAPI,
                              double ***jacobiMatrix) {

    double **U, **T, **spkMatrix, **combined;
    int i, *init_centroids;

    /** create (2 * n) matrix, which includes the eigenValues in the 1'st line,
      * and the original index of each eigenValue in the line below */
    combined = createMatrix(2, arraySize);

    /** sort jacobi's eigenValues*/
    sortEigenValues(arraySize, jacobiMatrix, &combined);

    /** determine k (if k==0) AND change k variable accordingly,
     * important for CAPI */
    if(*k == 0){
        *k = getK(arraySize, &combined);
    }

    /** create U with the first k sorted jacobi's eigenVectors & normalize it */
    U = getMatrixU(*k, arraySize, jacobiMatrix, combined);
    freeMatrix(&combined);

    normalizeU(*k, arraySize, &U);

    T = U;

    if (isCAPI){ /* return T for KMeans++ algorithm */
        return T;
    }

    /** Initialization of the k clusters as the first k points Respectively */
    init_centroids = calloc((*k), sizeof(int));
    ASSERT_ERROR(init_centroids != NULL)
    for (i=0 ; i < (*k); i++){
        init_centroids[i] = i;
    }

    /** run KMeans to find spkMatrix */
    spkMatrix = KMeansAlgorithm(*k, arraySize, &T, &init_centroids);

    /* free memory and return results */
    free(init_centroids);
    freeMatrix(&T);

    return spkMatrix;
}

/**
 * Takes the jacobi matrix, and update the combined matrix to possess
 * the eigenValues in sorted order.
 * @param arraySize - size of the nXn matrix
 * @param jacobiMatrix - pointer to the jacobi matrix
 * @param combined - pointer to 2D array, represents pairs of eigenValue &
 *                  original index
 */
void sortEigenValues(int arraySize, double ***jacobiMatrix, double ***combined){
    int i;
    double **tmp;
    /* Auxiliary array to use while sorting,
     * because mergeSort does not sort in place */
    tmp = createMatrix(2, arraySize);

    /** initialization of combined with the initial values */
    for (i = 0; i < arraySize; i++){
        (*combined)[0][i] = (*jacobiMatrix)[0][i]; /* the eigenValue */
        (*combined)[1][i] = i; /* the index */
    }

    /** sort the eigenValues */
    mergeSortEigenValues(combined, &tmp, 0, arraySize - 1);

    freeMatrix(&tmp);
}


/**
 * Takes the jacobi matrix, sort it's eigenValues, and update the
 * combined matrix.
 * @param combined - pointer to 2D array, includes pairs of eigenValue &
 *                  originalValue
 * @param tmp - pointer to 2D Auxiliary array to combined
 * @param low - the lower index of the subArray to be sorted
 * @param high - the higher index of the subArray to be sorted
 */
void mergeSortEigenValues(double ***combined, double ***tmp, int low, int high){
    int mid;

    if(low < high){
        mid = (low + high) / 2;
        mergeSortEigenValues(combined, tmp, low, mid);
        mergeSortEigenValues(combined, tmp, mid + 1, high);
        mergeCombined(combined, tmp, low, mid, high);
    }
    else{
        return;
    }
}


/**
 * Takes the combined 2D array, and merge 2 sorted subArrays of it,
 * based on the given indices
 * @param combined - pointer to 2D array, includes pairs of eigenValue &
 *                  originalValue
 * @param tmp - pointer to 2D Auxiliary array to combined
 * @param low - the lower index of the 1'st subArray to be merged
 * @param mid - the higher index of the 1'st subArray to be merged
 * @param high - the higher index of the 2'nd subArray to be merged
 */
void mergeCombined(double ***combined, double ***tmp, int low,
                   int mid, int high){
    int l1, l2, i;

    for(l1 = low, l2 = mid+1, i = low; l1 <= mid && l2 <= high; i++) {
        /* left value is smaller\equal to the right one */
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

    /** copy the rest of the "1'st array" */
    while(l1 <= mid){
        (*tmp)[0][i] = (*combined)[0][l1];
        (*tmp)[1][i] = (*combined)[1][l1];

        i++;
        l1++;
    }

    /** copy the rest of the "2'nd array" */
    while(l2 <= high){
        (*tmp)[0][i] = (*combined)[0][l2];
        (*tmp)[1][i] = (*combined)[1][l2];
        i++;
        l2++;
    }

    /** copy the sorted eigenValues from the tmp array, back to the
     * original combined array */
    for(i = low; i <= high; i++){
        (*combined)[0][i] = (*tmp)[0][i];
        (*combined)[1][i] = (*tmp)[1][i];
    }
}


/**
 * Applies the Eigengap Heuristic and returns k based on the largest
 * difference (in absolute value) between two adjacent eigenvalues
 * @param arraySize - amount of datapoints
 * @param combined - pointer to 2-D array containing pairs of
 *                   sorted eigenValues & original indices
 * @return k - the number of clusters
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
 * Returns matrix U (n*k), which is the first K columns (eigenVectors)
 * in the sorted jacobiMatrix, based on combined
 * @param k - number of initial clusters
 * @param arraySize - amount of datapoints
 * @param jacobiMatrix - pointer to the jacobi matrix
 * @param combined - 2D array contains pairs of sorted eigenValues &
 *                  original index
 * @return the matrix U
 */
double **getMatrixU(int k, int arraySize, double ***jacobiMatrix,
                    double **combined){

    int i, j, col;
    double **U = createMatrix(arraySize, k);

    /** take the first K eigenVectors of jacobiMatrix assuming
     * it is sorted by ascending order of the eigenValues */
    for(j = 0; j < k; j++){
        /* The original index of the j eigenValue is sorted order */
        col = (int) combined[1][j];
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
 * @param U - pointer to the matrix to be normalized by rows
 */
void normalizeU(int k, int arraySize, double ***U){
    int i, j;
    double valueOfRow;

    /** loop over each lime, calculate the root of the sums, and normalize each
     * element in the row */
    for(i = 0; i < arraySize; i++){

        valueOfRow = 0;
        for(j = 0; j < k; j++){
            valueOfRow += pow((*U)[i][j], 2);
        }

        valueOfRow = sqrt(valueOfRow);

        if(valueOfRow != 0){
            for (j = 0; j < k; j++) {
                (*U)[i][j] /= valueOfRow; /* T[i][j] = U[i][j] / valueOfRow */
            }
        }
    }
}

/**
 * Calculates the KMeans algorithm on T, mostly copied from HW1 / HW2
 * @param k - pointer to the number of clusters
 * @param arraySize - amount of datapoints
 * @param T - datapoints on which we will run KMeans algorithm on
 * @param init_centroids - indices of initial centroids from K-Means++
 *                         (or 1,2,...,K if C)
 * @return k X k matrix containing the spectral clusters for the datapoints
 */
double **KMeansAlgorithm(int k, int arraySize,
                         double ***T, int **init_centroids) {

    /** Initialization*/
    int max_iter, update, *datap_cluster_assignment, *countArray;
    double **t_centroids, **sumArrayHead;

    max_iter = MAX_ITER;
    update = 1;

    /* gets centroids into matrix, according to list of indices */
    t_centroids = centroidsFromList(T, init_centroids, k);

    datap_cluster_assignment = calloc( arraySize , sizeof (int));
    ASSERT_ERROR(datap_cluster_assignment != NULL)

    sumArrayHead = createMatrix(k, k);

    countArray = calloc(k, sizeof(int));
    ASSERT_ERROR(countArray != NULL)

    /* until no update is seen in cluster assignment or up to 300 iterations */
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

/**
 * Calculates distances of datapoints from each cluster, and assigning them
 * to their appropriate clusters according to shortest distance
 * Copied from HW1
 * @param datapoint - matrix containing datapoints points that are to be
 *                    assigned to clusters
 * @param centroid - matrix containing initial centroids of the clusters
 * @param datap_cluster_assignment - array that remembers the assignment
 *                                   of each datapoint
 * @param d - dimension of datapoints
 * @param k - amount of clusters to be
 * @param arraySize - amount of datapoints
 * @param sumArrayHead - matrix used internally for summation purposes
 * @param countArray - array used internally for counting purposes
 * @return 1 if the cluster assignment was updated, 0 otherwise
 */
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
                dist += ((*datapoint)[i][v] - (*centroid)[j][v]) *
                        ((*datapoint)[i][v] - (*centroid)[j][v]);
            }

            if (min_dist > dist){
                min_dist = dist;
                min_cluster = j;
            }
        }
        if((*datap_cluster_assignment)[i] != min_cluster){
            /* there is a change in one or more of the data points
             * cluster assignment */
            update = 1;
        }
        (*datap_cluster_assignment)[i] = min_cluster;
    }

    /* loop to initialize sum/counter */

    for(i = 0; i < arraySize; i++){ /* count and sum up all the sizes */
        currCluster = (*datap_cluster_assignment)[i];
        (*countArray)[currCluster]++;
        for(v = 0; v < d; v++){
            (*sumArrayHead)[currCluster][v] += (*datapoint)[i][v];
        }
    }

    /** update the new clusters and initialize to 0 */
    for(j = 0; j < k; j++) {
        /* each loop for different cluster */
        for (v = 0; v < d; v++) {
            /* each loop for feature of the current cluster*/

            new_value = (*sumArrayHead)[j][v] / (*countArray)[j];
            (*centroid)[j][v] = new_value;
            (*sumArrayHead)[j][v] = 0;
        }
        (*countArray)[j] = 0;
    }

    return update;
}