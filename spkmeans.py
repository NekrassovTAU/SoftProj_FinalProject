import sys
from math import inf
import numpy as np
from spkmeansmodule import initializeCProcess, KMeansPlusPlusIntegration


def main():
    arguments = sys.argv.copy()

    # run C program through module
    ret_matrix = initializeCProcess(arguments)

    # in case of goal = spk, we need to run KMeans++ on T and
    # return initial centroids to KMeans in C
    if sys.argv[2] == "spk":
        init_centroids = k_means_plus_plus(ret_matrix)
        print_init_centroids(init_centroids)
        ret_matrix = KMeansPlusPlusIntegration(ret_matrix,
                                               init_centroids.copy())

    print_results(ret_matrix, sys.argv[2])


# code copied from HW2
# calculates KMeans++ on T
def k_means_plus_plus(datapoints):
    array_size = len(datapoints)
    k = len(datapoints[0])
    d = len(datapoints[0])
    np.random.seed(0)

    # construct array of clusters - the first cluster is randomly selected
    clusters = np.zeros((k, d))
    initial_clusters = [0 for _ in range(k)]

    # initialize first cluster
    first_cluster_idx = np.random.choice(array_size)
    clusters[0] = datapoints[first_cluster_idx]
    initial_clusters[0] = first_cluster_idx

    # loop over and update k initial clusters
    z = 1
    while z < k:

        probability = [float(inf) for _ in
                       range(array_size)]  # initialize the probability array

        # loop over all combinations of datapoints and centroids
        # and find for each datapoint a minimum distance to a centroid
        for i in range(array_size):
            for j in range(z):
                cur_norm = np.inner(datapoints[i] - clusters[j],
                                    datapoints[i] - clusters[j])
                if cur_norm < probability[i]:
                    probability[i] = cur_norm

        # create probability array based on distances
        # (saved in probability array apriori)
        sum_of_probs = sum(probability)
        for i in range(len(probability)):
            probability[i] = probability[i] / sum_of_probs

        # use 'probability' to randomly choose point to be cluster
        next_cluster_idx = np.random.choice(array_size, p=probability)
        clusters[z] = datapoints[next_cluster_idx]  # initialize another cluster
        initial_clusters[z] = next_cluster_idx

        z += 1

    return initial_clusters


def print_init_centroids(init_centroids):
    k = len(init_centroids)
    for i in range(k):
        print(init_centroids[i], end="")
        if i != k - 1:
            print(",", end="")
        else:
            print("\n", end="")


def print_results(ret_matrix, goal):

    ret_matrix = fix_zeros(ret_matrix)

    if goal == "ddg":
        ret_matrix = ddg_printable(ret_matrix)

    if goal == "jacobi":
        ret_matrix = invert_jacobi(ret_matrix)

    row = len(ret_matrix)
    col = len(ret_matrix[0])

    for i in range(row):
        for j in range(col):
            print("{:.4f}".format(ret_matrix[i][j]), end="")
            if j != col - 1:
                print(',', end="")
            elif i != row - 1:
                print("\n", end="")


def fix_zeros(ret_matrix):
    row = len(ret_matrix)
    col = len(ret_matrix[0])

    for i in range(row):
        for j in range(col):
            if -0.00005 < ret_matrix[i][j] < 0:
                ret_matrix[i][j] = 0

    return ret_matrix


def ddg_printable(ret_matrix):
    n = len(ret_matrix[0])

    matrix = [[0.0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        matrix[i][i] = ret_matrix[0][i]

    return matrix


def invert_jacobi(ret_matrix):
    transposed = [ret_matrix[0]]
    transposed.extend(np.transpose(ret_matrix[1:]))

    return transposed


if __name__ == "__main__":
    main()
