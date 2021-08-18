import sys
from math import inf
import numpy as np
from spkmeansmodule import initializeCProcess, KMeansPlusPlusIntegration


def main():
    arguments = sys.argv.copy()

    ret_matrix = initializeCProcess(arguments)

    if sys.argv[2] == "spk":
        init_centroids = k_means_plus_plus(ret_matrix)
        print(init_centroids)
        ret_matrix = KMeansPlusPlusIntegration(ret_matrix,
                                               init_centroids.copy())

    #print_results(ret_matrix)


# code copied from HW2
def k_means_plus_plus(datapoints):
    array_size, d, k = len(datapoints), len(datapoints[0]), len(datapoints[0])
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


def print_results(ret_matrix):
    print(ret_matrix)


if __name__ == "__main__":
    main()
