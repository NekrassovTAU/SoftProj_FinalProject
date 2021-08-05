#ifndef SOFTPROJ_FINAL_SPKMEANS_H
#define SOFTPROJ_FINAL_SPKMEANS_H

enum goalEnum {
    spk, wam, ddg, lnorm, jacobi, INVALID
};

double **checkArgs(int argc, char **origArgv, int isCAPI, int *returnRowCount,
                   int *returnColCount);

double **
initProcess(int *k, int d, int arraySize, enum goalEnum goal,
            double ***datapoint, int isCAPI);

double **
goalBasedProcess(int *k, int d, int arraySize, double ***datapoint,
                 enum goalEnum goal, int isCAPI);

void printResult(double ***retArray, enum goalEnum goal);

static enum goalEnum checkGoal(char *string);

static double **arrayToTwoDimArray(double **array, int n, int m);

static double ** createMatrix(int n, int m);

void processDatapoints(char *filename, double **datap_array, double ***datapoint,
                       int **dArraySizeInfo);

/****************************/
double ** calcWeightedAdjMatrix(int d, int arraySize, double ***datapoint);

double calcWeight(int d, double ***datapoint, int i, int j);

double ** calcDiagDegMatrix(int arraySize, double ***weightedAdjMatrix);

double ** calcNormLaplacian(int arraySize, double ***weightedAdjMatrix,
                            double ***diagDegMatrix);

double ** calcJacobi(int arraySize, double ***inputMatrix);

void copymatrix(int arraySize, double ***matrix1, double ***matrix2);

void findmatrixP(int arraySize, double ***A, double *c, double *s, int *row, int *col);

void findMaxOffDiag(int arraySize, double ***A, int *row, int *col);

void updateAtag(int arraySize, double ***Atag, double ***A, double c, double s, int row, int col);

void updateV(int arraySize, double ***V, double c, double s, int row, int col);

int convergenceCheck(int arraySize, double ***A, double ***Atag);

double calcOff(int arraySize, double ***matrix);

double** copyJacoby(int arraySize, double ***Atag, double ***V);


double **
calcSpectralClusters(int *k, int arraySize, int isCAPI, double ***jacobiMatrix);

void printTest(double **matrix, int n, int m);

void makeIntoIdentityMatrix(double ***emptyMatrix, int matrixSize);

void freeMatrix(double ***matrix);

double ** KMeansAlgorithm(int k, int *d, int arraySize, double ***datapoints,
                          int isCAPI, int **init_centroids);

#endif //SOFTPROJ_FINAL_SPKMEANS_H
