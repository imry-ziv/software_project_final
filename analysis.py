import argparse
from typing import List, Tuple

from sklearn.metrics import silhouette_score

from kmeans import kmeans, euclidean_distance, Point
from symnmf import get_points_from_file, compute_symnmf, derive_clustering_solution_symnmf
import numpy as np

MAX_ITER = 300
EPSILON = 1e-4

def compute_sscore(
        nmf_clusters: List[List[List[float]]],
        kmeans_clusters: List[List[List[float]]],
) -> Tuple[float, float]:

    return (
        compute_sscore_for_method(nmf_clusters),
        compute_sscore_for_method(kmeans_clusters),
    )


def compute_sscore_for_method(clusters: List[List[Point]]) -> float:
    # Flatten the list of clusters to a list of points
    flat_clusters = flatten_3d_to_2d(clusters)
    n = len(flat_clusters)

    # Create a 2D numpy array for data and a 1D numpy array for labels
    data_array = np.array([point.coordinates for point in flat_clusters])
    labels_array = np.empty(n, dtype=int)

    index = 0
    for cluster_index, cluster in enumerate(clusters):
        for datapoint_index, datapoint in enumerate(cluster):
            labels_array[index] = cluster_index
            index += 1

    # Calculate the silhouette score
    score = silhouette_score(data_array, labels_array)
    return score









#
# def compute_within_cluster_distance_for_point(
#         point: Tuple[float, ...],
#         cluster_of_point: List[Tuple[float, ...]],
#
# ) -> float: # Based on assumption that all data points are different!
#     avg = 0
#     for other_point in cluster_of_point:
#         if other_point == point:
#             continue
#         d = euclidean_distance(point, other_point)
#         avg += d
#     avg /= len(cluster_of_point) - 1 # Minus one because we are not counting point itself.
#     return avg
#
# def compute_between_cluster_distance_for_point(
#         point: Tuple[float, ...],
#         all_clusters: List[List[Tuple[float, ...]]],
#         index_of_point_cluster: int, # To be skipped
# ) -> float:
#     min_dist = float('inf')
#     for i in range(len(all_clusters)):
#         other_cluster = all_clusters[i]
#         avg = 0
#         if i == index_of_point_cluster: # Skip point's own cluster
#             continue
#         for other_point in other_cluster:
#             d = euclidean_distance(point, other_point)
#             avg += d
#         avg /= len(other_cluster)
#
#         if avg < min_dist:
#             min_dist = avg
#
#     return min_dist
#
#

def flatten_3d_to_2d(points_list: List[List[Point]]) -> List[Point]:
    return [point for sublist in points_list for point in sublist]



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
    args = parser.parse_args()
    k = args.k
    input_file = args.input_file

    data_matrix = get_points_from_file(input_file)
    n = len(data_matrix)
    try:
        d = len(data_matrix[0]) # We assume that all points have the same d.
    except:
        print('tbd')
    # For both methods,
    # resulting "clusters" is List[List[float]]
    # where the cluster is denoted by its index in list.

    # Compute symnmnf clusters
    best_H = compute_symnmf(n, k, d, data_matrix)
    nmf_clusters = derive_clustering_solution_symnmf(data_matrix, best_H, k)

    # Compute kmeans clusters
    _, kmeans_clusters = kmeans(n, k, MAX_ITER, data_matrix)

    sscore_nmf, sscore_kmeans = compute_sscore(nmf_clusters, kmeans_clusters)
    print(f'nmf: {sscore_nmf}')
    print(f'kmeans: {sscore_kmeans}')
