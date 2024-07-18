import argparse
import sys
from typing import List, Tuple

from sklearn.metrics import silhouette_score

from kmeans import kmeans, euclidean_distance, Point
from symnmf import compute_symnmf, derive_clustering_solution_symnmf, file_to_lists
import numpy as np

MAX_ITER = 300
EPSILON = 1e-4


YEET_MATRIX = [[0.0824,0.0924,0.0965,0.0701,0.0434,0.1544],
[0.1134,0.0967,0.1654,0.0706,0.1077,0.0673],
[0.1508,0.1582,0.0609,0.0707,0.0250,0.1771],
[0.1121,0.1165,0.1393,0.1045,0.0811,0.1333],
[0.0552,0.1284,0.0744,0.1239,0.1214,0.1398],
[0.0959,0.1273,0.1376,0.0986,0.0101,0.1500],
[0.1077,0.1144,0.1255,0.0975,0.0709,0.0842],
[0.1667,0.0350,0.1768,0.0970,0.0802,0.0476],
[0.1064,0.1125,0.1370,0.1092,0.1224,0.0547],
[0.0821,0.0968,0.1463,0.1009,0.1010,0.1457],
[0.0699,0.0673,0.1615,0.0683,0.0798,0.1526],
[0.1679,0.0374,0.1606,0.0361,0.0984,0.0908],
[0.1694,0.0918,0.1770,0.0196,0.0815,0.0466],
[0.1208,0.0574,0.1281,0.0983,0.0431,0.1532],
[0.0912,0.0951,0.1178,0.0479,0.0928,0.1756]]

def file_to_points(
        file_name: str,
) -> List[List[float]]:


    matrix = []
    with open(file_name, 'r') as file:
        for line in file:
            row = line.strip().split(',')
            matrix.append(Point([float(value) for value in row]))
    return matrix

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

    # input_file = 'test_inputs/test_input9.txt'
    # k=6

    data_matrix_as_points = file_to_points(input_file)
    data_matrix_as_lists = file_to_lists(input_file)

    n = len(data_matrix_as_points)
    try:
        d = len(data_matrix_as_points[0]) # We assume that all points have the same d.
    except:
        print('tbd')
    # For both methods,
    # resulting "clusters" is List[List[float]]
    # where the cluster is denoted by its index in list.

    # Compute symnmnf clusters
    best_H = compute_symnmf(n, k, d, data_matrix_as_lists)
    #best_H = YEET_MATRIX
    nmf_clusters = derive_clustering_solution_symnmf(data_matrix_as_points, best_H, k)

    # Compute kmeans clusters
    _, kmeans_clusters = kmeans(n, k, MAX_ITER, data_matrix_as_points)

    sscore_nmf, sscore_kmeans = compute_sscore(nmf_clusters, kmeans_clusters)

    print(f'nmf: {sscore_nmf:.4f}')
    print(f'kmeans: {sscore_kmeans:.4f}')
