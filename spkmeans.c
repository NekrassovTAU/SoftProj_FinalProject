#include "spkmeans.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define FLT_MAX 3.402823466e+38F
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

    *datap_array = (double*)calloc(d * arraySize, sizeof(double));
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
    jacobiMatrix = calcJacobi(d, arraySize, &weightedAdjMatrix, &diagDegMatrix);
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
    if(twoDimArray == NULL){
        return NULL;
    }

    for (i = 0; i < n; i++) {
        twoDimArray[i] = (*array) + i * m;
    }

    return twoDimArray;
}