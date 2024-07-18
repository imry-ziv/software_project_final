# Imports.
import argparse
import sys
from typing import List, Tuple
import numpy as np
from sklearn.metrics import silhouette_score
from kmeans import kmeans, Point
from symnmf import compute_symnmf, file_to_lists

# Environment variables.
MAX_ITER = 300
EPSILON = 1e-4
ERROR_MSG = 'An Error Has Occurred'

# Helper functions.

def flatten_3d_to_2d(
        points_list: List[List[Point]],
) -> List[Point]:
    """
    Given a 3d representation of the clusters,
    return a flat representation of all points in the clusters.
    """
    return [point for sublist in points_list for point in sublist]

def file_to_points(
        file_name: str,
) -> List[Point]:
    """
    Read data from .txt file, such that each point is represented as a list of Point objects.
    """
    data_matrix_list = []
    try:
        with open(file_name, 'r') as file:
            for line in file:
                row = line.strip().split(',')
                data_matrix_list.append(Point([float(value) for value in row]))

    except FileNotFoundError:
        print(ERROR_MSG)
        sys.exit(1)
    return data_matrix_list

def compute_sscore(
        nmf_clusters: List[List[Point]],
        kmeans_clusters: List[List[Point]],
) -> Tuple[float, float]:
    """
    Given the representation of nmf_clusters and kmeans_clusters,
    represented by lists corresponding to clusters,
    compute the sscore for each clustering scheme,
    and return a tuple with both scores: (nmf_score, kmeans_score)
    """

    return (
        compute_sscore_for_method(nmf_clusters),
        compute_sscore_for_method(kmeans_clusters),
    )


def compute_sscore_for_method(
        clusters: List[List[Point]]
) -> float:
    """
    Given a list of clusters (Point representation),
    derive a List representation of the points using point.coordinates,
    and calculate silhouette score for clustering solution with sklearn.metrics.silhouette.score.
    """

    # Flatten the list of clusters to a list of points.
    flat_clusters = flatten_3d_to_2d(clusters)
    n = len(flat_clusters)

    # Create a 2D numpy array for data and a 1D numpy array for labels.
    data_array = np.array([point.coordinates for point in flat_clusters])
    labels_array = np.empty(n, dtype=int)

    index = 0
    for cluster_index, cluster in enumerate(clusters):
        for datapoint_index, datapoint in enumerate(cluster):
            labels_array[index] = cluster_index
            index += 1

    # Calculate the silhouette score.
    score = silhouette_score(data_array, labels_array)
    return score



def derive_clustering_solution_symnmf(
        data_matrix: List[Point],
        best_H: List[List[float]],
        k: int,
) -> List[List[Point]]:
    """
    Given the association matrix best_H and the data_matrix,
    derive list of lists, each inner list corresponding to a cluster,
    and containing Point objects.
    """
    clusters = [[] for _ in range(k)] # Initialize clusters as empty lists.
    for i in range(len(data_matrix)):
        datapoint = data_matrix[i]
        association_row = best_H[i]
        # Cluster index is the maximum value in the corresponding row
        cluster_index = association_row.index(max(association_row))
        clusters[cluster_index].append(datapoint)
    return clusters





if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument(
        'k',
        type=int,
        help='An integer representing the number of clusters'
    )
    parser.add_argument(
        'input_file',
        type=str,
        help='The path to the input .txt file containing the points'
    )
    try:
        args = parser.parse_args()
    except SystemExit as e:
        print(ERROR_MSG)
        sys.exit(1)

    k = args.k
    input_file = args.input_file

    # Parse data in both ways: as Point objects and as lists of floats.
    data_matrix_as_points = file_to_points(input_file)
    data_matrix_as_lists = file_to_lists(input_file)

    n = len(data_matrix_as_points)
    try:
        d = len(data_matrix_as_points[0]) # We assume that all points have the same d.
    except:
        print(ERROR_MSG)
        sys.exit(1)

    # Compute nmf clustering solution.
    best_H = compute_symnmf(n, k, d, data_matrix_as_lists)
    nmf_clusters = derive_clustering_solution_symnmf(data_matrix_as_points, best_H, k)

    # Compute kmeans clustering solution.
    _, kmeans_clusters = kmeans(n, k, MAX_ITER, data_matrix_as_points)

    # Final sscores.
    sscore_nmf, sscore_kmeans = compute_sscore(nmf_clusters, kmeans_clusters)

    print(f'nmf: {sscore_nmf:.4f}')
    print(f'kmeans: {sscore_kmeans:.4f}')
