#ifndef SOFTPROJ_FINAL_SPKMEANS_H
#define SOFTPROJ_FINAL_SPKMEANS_H

enum goalEnum {
    spk, wam, ddg, lnorm, jacobi, INVALID
};

double **startCSPKMeansExec(int argc, char **origArgv, int isCAPI,
                            int *returnRowCount, int *returnColCount);

double **
initProcess(int *k, int d, int arraySize, enum goalEnum goal,
            double ***datapoint, int isCAPI);

double **
goalBasedProcess(int *k, int d, int arraySize, double ***datapoint,
                 enum goalEnum goal, int isCAPI);

void printResult(double ***ret_matrix, enum goalEnum goal,
        int arraySize, int k);

void printRegular(double ***ret_matrix, int rows, int cols);

void printJacobi(double ***retArray, int arraySize);

void printDiagonal(double ***ret_matrix , int arraySize);


enum goalEnum checkGoal(char *string);

double ** createMatrix(int n, int m);

void
processDatapoints(char *filename, double ***datapoint, int *d, int *arraySize);

/****************************/
double ** calcWeightedAdjMatrix(int d, int arraySize, double ***datapoint);

double calcWeight(int d, double ***datapoint, int i, int j);

double ** calcDiagDegMatrix(int arraySize, double ***weightedAdjMatrix);

double ** calcNormLaplacian(int arraySize, double ***weightedAdjMatrix,
                            double ***diagDegMatrix);

double ** calcJacobi(int arraySize, double ***inputMatrix);

void findmatrixP(int arraySize, double ***A, double *c, double *s,
                 int *row, int *col);

void findMaxOffDiag(int arraySize, double ***A, int *row, int *col);

void updateAtag(int arraySize, double ***A, double c, double s, int row,
                int col, double *offAtag);

void updateV(int arraySize, double ***V, double c, double s, int row, int col);

double calcOff(int arraySize, double ***matrix);

double **copyJacoby(int arraySize, double ***Atag, double ***V);


double **calcSpectralClusters(int *k, int arraySize, int isCAPI,
                              double ***jacobiMatrix);

/***/
void sortEigenValues(int arraySize, double ***jacobiMatrix, double ***combined);
void mergeSortEigenValues(double ***combined, double ***tmp, int low, int high);
void mergeCombined(double ***combined, double ***tmp, int low,
                   int mid, int high);

/***/

int getK(int arraySize, double ***combined);

double **getMatrixU(int k, int arraySize, double ***jacobiMatrix,
                    double **combined);

void normalizeU(int k, int arraySize, double ***U);




void makeIntoIdentityMatrix(double ***emptyMatrix, int matrixSize);

void freeMatrix(double ***matrix);

double **
KMeansAlgorithm(int k, int arraySize, double ***T, int **init_centroids);


/** Testers */
void printTest(double **matrix, int n, int m);

void TesterToSortJacobi();

double **TesterToWeight();

double **randomMatrix(int n, int m);

void print1Darray(int *array, int arraySize);

void determineRowAndCol(enum goalEnum goal, int k, int arraySize, int *rowCount,
                        int *colCount, int isCAPI);

double **centroidsFromList(double ***pointList, int **init_centroids, int k);

int updateCentroidsPerDatap(double ***datapoint, double ***centroid,
                            int **datap_cluster_assignment, int d, int k,
                            int arraySize, double ***sumArrayHead,
                            int **countArray);

void
copyDatapoint(double ***datapoint, double **datap_array, int arraySize, int d);

void fixZeros(double ***matrix, int rows, int cols);

#endif
