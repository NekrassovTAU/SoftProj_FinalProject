#include "spkmeans.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <math.h>

#define FLT_MAX 3.402823466e+38F
#define FLT_MIN 1.1754943508e-38F
#define BUFFER_SIZE 1000
#define INITIAL_ARRAY_SIZE 1000

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

enum goalEnum {
    spk, wam, ddg, lnorm, jacobi, INVALID
};

void checkArgs(int argc, char **origArgv);

double **
initProcess(int k, int d, int arraySize, enum goalEnum goal,
            double ***datapoint, int **init_centroids, int isCAPI);

double **
goalBasedProcess(int d, int arraySize, double ***datapoint, enum goalEnum goal,
                 int **init_centroids);

void printResult(double ***retMatrix, enum goalEnum goal);

static enum goalEnum checkGoal(char *string);

static double **makeTwoDimArray(double **array, int n, int m);

int *processInitCentroidFile(char *initCentrFilename);

void processDatapoints(char *filename, double **datap_array, double ***datapoint,
                       int **dArraySizeInfo);

/****************************/
double **calcWeightedAdjMatrix(int d, int arraySize, double ***datapoint);

double calcWeight(int d, double ***datapoint, int i, int j);

double **calcDiagDegMatrix(int d, int arraySize, double ***weightedAdjMatrix);

double **calcNormLaplacian(int d, int arraySize, double ***weightedAdjMatrix, double ***diagDegMatrix);

double **twoDinitialization(int arraySize, int identity);

double **calcJacobi(int d, int arraySize, double ***normLaplacian);

void copymatrix(int arraySize, double ***matrix1, double ***matrix2);

void findmatrixP(int arraySize, double ***A, double ***P, int *row, int *col);

void findMaxOffDiag(int arraySize, double ***A, int *row, int *col);

void updateAtag(int arraySize, double ***Atag, double ***P, int row, int col);

void updateV(int arraySize, double ***V, double ***P, int row, int col);

/**
 * main, a shell function for the spectral clustering algorithm implementation
 */
int main(int argc, char *argv[]) {
    checkArgs(argc, argv);
    return 0;
}

/**
 * Checks the validity of arguments passed, reads and parses info in files
 * and starts the implementation of spectral k-means clustering
 * @param argc - # of arguments passed to main
 * @param origArgv - the arguments passed to main
 */
void checkArgs(int argc, char **origArgv) {
    int k, *init_centroids, d, arraySize, *dArraySizeInfo;
    double *datap_array, **datapoint;
    enum goalEnum goal;
    char *ptr;
    ASSERT_ARGS(argc > 3)

    k = strtol(origArgv[1], &ptr, 10);

    /** checking whether k>=0,
     * strcmp with "0" used for the case which strtol fails and returns 0*/
    ASSERT_ARGS(k > 0 || !strcmp(origArgv[1], "0"))

    goal = checkGoal(origArgv[2]);
    ASSERT_ARGS(goal != INVALID)

    dArraySizeInfo = (int *) calloc(2, sizeof(int));

    /** deals with optional init_centroids argument, process files accordingly
     * and initialize the process to achieve the provided goal*/
    if (argc > 4) {
        init_centroids = processInitCentroidFile(origArgv[3]);
        processDatapoints(origArgv[4], &datap_array, &datapoint, &dArraySizeInfo);

        d = dArraySizeInfo[0];
        arraySize = dArraySizeInfo[1];
        free(dArraySizeInfo);

        ASSERT_ARGS(k < arraySize)

        initProcess(k, d, arraySize, goal, &datapoint, &init_centroids, 0);

        free(init_centroids);
    } else {
        /** if init_centroids is not provided, spk cannot be executed */
        ASSERT_ARGS(goal != spk)
        processDatapoints(origArgv[3], &datap_array, &datapoint, &dArraySizeInfo);

        d = dArraySizeInfo[0];
        arraySize = dArraySizeInfo[1];
        free(dArraySizeInfo);

        ASSERT_ARGS(k < arraySize)

        initProcess(k, d, arraySize, goal, &datapoint, NULL, 0);
    }
    free(datapoint);
    free(datap_array);
}

/**
* Initializes datapoint-related arrays on a row-by-row basis
* Mostly copied from HW1
* TODO: DEAL WITH FREEING MEMORY ON ASSERT_ERROR
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

    for (i = 0, d = 0; currRow[i]; i++) {
        d += (currRow[i] == ',');
    }

    arraySize = INITIAL_ARRAY_SIZE;

    *datap_array = (double *) calloc(d * arraySize, sizeof(double));
    ASSERT_ERROR(*datap_array != NULL)

    /** Populates datap_array with info from the file */
    /**TODO: ?Move datap_array population to different method?*/
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

    *datapoint = makeTwoDimArray(datap_array, arraySize, d);

    *dArraySizeInfo[0] = d;
    *dArraySizeInfo[1] = arraySize;

}

/**
 * Read the initial centroid file and creates an int-array with the info in it
 * @param initCentrFilename - path incl. filename containing initial centroid
 * @return int array with initial centroid indexes
 */
int *processInitCentroidFile(char *initCentrFilename) {
    return NULL;
}

/**
 * Initializes datapoints and related arrays. Mostly copied from HW1.
 * @param k - amount of clusters the data needs to divide into.
 *            if 0, need to apply Eigengap Heuristic
 * @param goal - goal of call to main (spk / wam / ddg / lnorm / jacobi)
 * @param datapoint - path to file with datapoints
 * @param init_centroids - path to file with initial centroids positions.
 *                      used only in C since K-means++ impl. was in Python
 * @param isCAPI - used in order to distinguish between
 *                 C-based and C-API-based call to this function
 */
double **
initProcess(int k, int d, int arraySize, enum goalEnum goal,
            double ***datapoint, int **init_centroids, int isCAPI) {

    double **ret_matrix;

    ret_matrix = goalBasedProcess(d, arraySize, datapoint, goal,
                                  init_centroids);

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
 * Takes the datapoint info, and returns relevant matrix according to the goal.
 * @param d - dimension of datapoints
 * @param arraySize - amount of datapoints
 * @param datapoint - 2D-array containing all datapoints
 * @param goal - the goal of the computation
 * @return relevant matrix according to goal
 */
double **
goalBasedProcess(int d, int arraySize, double ***datapoint, enum goalEnum goal,
                 int **init_centroids) {
    double **weightedAdjMatrix, **diagDegMatrix, **normLaplacian,
            **jacobiMatrix, **spkMatrix;
    /** TODO: implement calc functions */
    weightedAdjMatrix = calcWeightedAdjMatrix(d, arraySize, datapoint);
    if (goal == wam) {
        return weightedAdjMatrix;
    }
    diagDegMatrix = calcDiagDegMatrix(d, arraySize, &weightedAdjMatrix);
    if (goal == ddg) {
        free(weightedAdjMatrix);
        return diagDegMatrix;
    }
    normLaplacian = calcNormLaplacian(d, arraySize, &weightedAdjMatrix,
                                      &diagDegMatrix);
    if (goal == lnorm) {
        free(weightedAdjMatrix);
        free(diagDegMatrix);
        return normLaplacian;
    }
    jacobiMatrix = calcJacobi(d, arraySize, &normLaplacian);
    if (goal == jacobi) {
        free(weightedAdjMatrix);
        free(diagDegMatrix);
        free(normLaplacian);
        return jacobiMatrix;
    }
    /** TODO: if Python, we should return to implement K-Means++*/
    spkMatrix = calcSpectralClusters();
    free(weightedAdjMatrix);
    free(diagDegMatrix);
    free(normLaplacian);
    free(jacobiMatrix);
    return spkMatrix;
}

/**
 * Prints the result of the algorithm according to goal
 * wam / ddg / lnorm - print a nXn matrix
 * jacobi - first line is eigenvalues,
 *          next are nXn eigenvecor matrix
 * spk - first line is k-means++ init_centroid indices
 *       next are the k centroids after k-means algorithm
 * @param retMatrix - matrix to be printed to user
 * @param goal - distinguishes between spk-printing and matrix-printing
 */
void printResult(double ***retMatrix, enum goalEnum goal) {
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
 * @param array - a n*m-sized array
 * @param n - # of rows
 * @param m - # of columns
 * @return a 2D array referencing this array
 */
static double **makeTwoDimArray(double **array, int n, int m) {
    double **twoDimArray;
    int i;

    twoDimArray = calloc(n, sizeof(double *));
    if (twoDimArray == NULL) {
        return NULL;
    }

    for (i = 0; i < n; i++) {
        twoDimArray[i] = (*array) + i * m;
    }

    return twoDimArray;
}

/**
 * Takes the datapoint info, and returns the weighted adj matrix.
 * @param d - dimension of datapoints
 * @param arraySize - amount of datapoints
 * @param datapoint - 2D-array containing all datapoints
 * @return weighted adj matrix
 */
double **calcWeightedAdjMatrix(int d, int arraySize, double ***datapoint) {
    double value;

    /** creating 2-D array */
    double **weightedAdjMatrix = (double **) calloc(arraySize, sizeof(double *));
    for (int i = 0; i < arraySize; i++) {
        weightedAdjMatrix[i] = (double *) calloc(arraySize, sizeof(double));
    }
    /** try this line instead
    *  double **diagDegMatrix = twoDinitialization(int arraySize, int identity);
    */

    /** calculating the weights */
    for (i = 0; i < arraySize; i++) {
        for (int j = i + 1; j < arraySize; j++) {
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
    double value = 0;

    /** calculating the sum of the distances */
    for (int t = 0; t < d; t++) {
        value += pow(((*datapoint)[i][t] - (*datapoint)[j][t]), 2);
    }

    value = exp(sqrt(value) / (-2));
    return value;
}

/**
 * Takes the weightedAdjMatrix, and returns the sum of weights for each point.
 * @param d - dimension of datapoints
 *   do you need it ?
 * @param arraySize - amount of datapoints
 * @param weightedAdjMatrix - 2-D array containing all the weights
 * @return 2-D array with sum of weights for each point
 */
double **calcDiagDegMatrix(int d, int arraySize, double ***weightedAdjMatrix) {

    double sumofweights;

    double **diagDegMatrix = (double **) calloc(arraySize, sizeof(double *));
    for (int i = 0; i < arraySize; i++) {
        diagDegMatrix[i] = (double *) calloc(arraySize, sizeof(double));
    }

/** try this line instead
 *  double **diagDegMatrix = twoDinitialization(int arraySize, 0);
 */

    for (i = 0; i < arraySize; i++) {
        sumofweights = 0;
        for (int j = 0; j < arraySize; j++) {
            sumofweights += (*weightedAdjMatrix)[i][j];
        }
        diagDegMatrix[i][i] = sumofweights;
    }

    return diagDegMatrix;
}

/**
 * Takes the weighted and the diagonal matrices  info, and returns the normLaplacian matrix.
 * @param d - dimension of datapoints
 *  do you need it ?
 * @param arraySize - amount of datapoints
 * @param weightedAdjMatrix - 2-D array containing all the weights
 * @param diagDegMatrix - 2-D array containing the sum of weights for each datapoint
 * @return NormLaplacian matrix
 */
double **calcNormLaplacian(int d, int arraySize, double ***weightedAdjMatrix, double ***diagDegMatrix) {

    double value;
    /** 1-D array represents sqrt of diagDegMatrix */
    double *sqrtDegMatrix = (double *) calloc(arraySize, sizeof(double));
    /** TODO: ASSERT */

    for (int i = 0; i < arraySize; i++) {
        sqrtDegMatrix[i] = 1 / sqrt((*diagDegMatrix)[i][i]);
    }

    /**  sqrtDegMatrix[i][j] = sqrtDeg[i] * weight[i][j] * sqrtDeg[j] = sqrtDegMatrix[j][i]
     * convince your self its true ;)   */
    double **normLaplacian = twoDinitialization(arraySize, 0);
    for (i = 0; i < arraySize; i++) {
        normLaplacian[i][i] = 1;
        for (int j = i + 1; j < arraySize; j++) {
            value = -(sqrtDegMatrix[i] * (*weightedAdjMatrix)[i][j] * sqrtDegMatrix[j]);
            normLaplacian[i][j] = value;
            normLaplacian[j][i] = value;
        }
    }
    free(sqrtDegMatrix);
    return normLaplacian;
}

/**
 * Gets size (N), and returns a matrix of size (N * N).
 * @param arraySize - amount of datapoints
 * @param identity - if identity = 1, initialize identity matrix. else, initialize matrix of '0'
 * @return 2-D initialized array
 */
double **twoDinitialization(int arraySize, int identity) {

    double **array = (double **) calloc(arraySize, sizeof(double *));
    /** TODO: ASSERT */
    if (identity == 1) {
        for (int i = 0; i < arraySize; i++) {
            array[i] = (double *) calloc(arraySize, sizeof(double));
            /** TODO: ASSERT */
            array[i][i] = 1;
        }
    } else {
        for (int i = 0; i < arraySize; i++) {
            array[i] = (double *) calloc(arraySize, sizeof(double));
            /** TODO: ASSERT */
        }
    }
    return array;
}

/**
 * Takes the Lnorm matrix, and  ?? calculate Eigenvalues and Eigenvectors ??
 * @param d - dimension of datapoints
 *  do you need it ?
 * @param arraySize - amount of datapoints
 * @param normLaplacian - 2-D array
 * @return  ??
 */
double **calcJacobi(int d, int arraySize, double ***normLaplacian) {
    int row, col;
    double **jacobiMatrix, **A, **Atag, **V, **P;

    /** initialization */
    A = twoDinitialization(arraySize, 0);
    Atag = twoDinitialization(arraySize, 0);
    V = twoDinitialization(arraySize, 1);
    P = twoDinitialization(arraySize, 1);
    jacobiMatrix = twoDinitialization(arraySize, 0);

    /** TODO: ASSERT after each line ?  */

    // probably left here by mistake. make sure before you delete it.  copymatrix(arraySize, normLaplacian, &curA); // curA = normLaplacian

    /** run until convergence -  */
    copymatrix(arraySize, normLaplacian, &Atag);
    do {
        copymatrix(arraySize, &Atag, &A);                   // A = Atag
        findmatrixP(arraySize, &A, &P, &row, &col);  // P
        updateAtag(arraySize, &Atag, &P, row, col);         // A' = P^T * A * P
        updateV(arraySize, &V, &P, row, col);                                // V *= P
    } while (!converge(&A, &Atag));                            // as long as delta > epsilon

    return NULL; // has to be changed
}

/**
 * Gets 2 matrices, and cupy the 1'st to the 2'nd a matrix of size (N * N).
 * @param arraySize - amount of datapoints
 * @param matrix1 - source matrix
 * @param matrix2 - matrix to be updated
 * TODO if ill copy " symmetric matrix, will that be better ?
 */
void copymatrix(int arraySize, double ***matrix1, double ***matrix2) {
    for (int i = 0; i < arraySize; i++) {
        for (int j = i; j < arraySize; j++) {
            (*matrix2)[i][j] = (*matrix1)[i][j];
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
void findmatrixP(int arraySize, double ***A, double ***P, int *row, int *col) {

    double teta, t, c, s, sign = -1;
    findMaxOffDiag(arraySize, A, row, col);

    /** update P values */
    teta = ((*A)[*col][*col] - (*A)[*row][*row]) / (2 * (*A)[*row][*col]);
    if (teta >= 0) {
        sign = 1;
    }
    t = sign / (fabs(teta) + sqrt(pow(teta, 2) + 1));
    c = 1 / sqrt(pow(t, 2) + 1);
    s = t * c;
    (*P)[*row][*row] = c;
    (*P)[*col][*col] = c;
    (*P)[*row][*col] = s;
    (*P)[*col][*row] = -s;

}

/**
 * Find the index of the maximum odd-diag value in A.
 * @param arraySize - amount of datapoints
 * @param A - 2-D symmetric array
 * @param row - pointer to the row index of the max value in A (off diag), to be calculated
 * @param col - pointer to the col index of the max value in A (off diag), to be calculated
 */
void findMaxOffDiag(int arraySize, double ***A, int *row, int *col) {
    double curmax = FLT_MIN;

    for (int i = 0; i < arraySize; i++) {
        for (int j = i + 1; j < arraySize; j++) {
            if ((*A)[i][j] > curmax) {
                curmax = (*A)[i][j];
                *row = i;
                *col = j;
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
void updateAtag(int arraySize, double ***Atag, double ***P, int row, int col) {

    double c, s, ri, rj, ij;
    c = (*Atag)[row][row];
    s = (*Atag)[row][col];

    /** update each Atag[r,i]/Atag[r,j], for : r != i,j
     * can cut computation in half because the matrix is symmetric */
    for (int r = 0; r < arraySize; r++) {
        ri = (*Atag)[r][row];
        rj = (*Atag)[r][col];
        if ((r != row) && (r != col)) {
            (*Atag)[r][row] = (c * ri) - (s * rj);
            (*Atag)[r][col] = (c * rj) + (s * ri);
        }
    }
    /**  Atag[i,i] / Atag[j,j]/ Atag[i,j] / Atag[j,i] */
    ij = (*Atag)[row][col];
    (*Atag)[row][row] = (c * c * (*Atag)[row][row]) + (s * s * (*Atag)[col][col]) - (2 * s * c * ij);
    (*Atag)[col][col] = (s * s * (*Atag)[row][row]) + (c * c * (*Atag)[col][col]) - (2 *s * c* ij);
    (*Atag)[row][col] = 0;
    (*Atag)[col][row] = 0;

}

void updateV(int arraySize, double ***V, double ***P, int row, int col) {
    double ri, rj, c, s;
    c = (*P)[row][row];
    s = (*P)[row][col];

    for (int r = 0; r < arraySize; r++) {
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