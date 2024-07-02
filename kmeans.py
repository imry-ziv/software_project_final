import sys
from typing import List, Tuple

EPSILON = 0.001
def read_file(
        file_path:str
) -> List[Tuple[float, ...]]:
    """
    Reads .txt, removes newline characters and returns list of tuples of floats.
    Each tuple in the list corresponds to a datapoint.
    """
    with open(file_path, 'r') as file:
        points = []
        for line in file:
            parts = tuple([float(val) for val in line.strip().split(',')])
            points.append(parts)
    return points


def euclidean_distance(
        v1: Tuple[float, ...],
        v2: Tuple[float, ...]
) -> float:

    return sum((x - y) ** 2 for x, y in zip(v1, v2)) ** 0.5

def find_closest_cluster(
        point: Tuple[float, ...],
        centroids: List[Tuple[float, ...]],
):

    closest_cluster_index = None
    min_distance = float('inf')
    for i in range(len(centroids)):
        cent = centroids[i]
        dist = euclidean_distance(point, cent)
        if dist < min_distance:
            min_distance = dist
            closest_cluster_index = i

    return closest_cluster_index

def average(
        cluster: List[Tuple[float, ...]]
) -> Tuple[float]:

    result = [0.0 for _ in range(len(cluster[0]))]

    for t in cluster:
        for i in range(len(t)):
            result[i] += t[i]
    result_divided = [i / len(cluster) for i in result]
    return tuple(result_divided)
def update_centroids(
        centroids: List[Tuple[float, ...]],
        clusters: List[List[Tuple[float, ...]]],
) -> Tuple[List[Tuple[float, ...]], bool]:

    stop_flag = False
    for i in range(len(centroids)):
        corr_cluster = clusters[i] # List of points corresponding to cluster
        average_of_cluster = average(corr_cluster)
        delta = euclidean_distance(average_of_cluster, centroids[i])
        if delta < EPSILON:
            stop_flag = True
        # Update value of centroid
        centroids[i] = average_of_cluster
    return centroids, stop_flag



def kmeans(
        N: int,
        k: int,
        iter: int,
        points: List[Tuple[float, ...]],
) -> Tuple[List[Tuple[float, ...]], List[List[Tuple[float, ...]]]]:

    assert len(points) == N
    centroids = points[:k]
    clusters = [[] for _ in range(k)] # Initialize clusters as empty lists.

    for i in range(iter):
        for point in points:
            closest_cluster_index = find_closest_cluster(point, centroids)
            clusters[closest_cluster_index].append(point)

    # Updating centroids
        centroids, stop_flag = update_centroids(centroids, clusters)
        if stop_flag:
            break
    return centroids, clusters

def show_centroids(
        centroids: List[Tuple[float, ...]],
) -> None:
    for cent in centroids:
        formatted_cent = [f"{value:.4f}" for value in cent]
        print(formatted_cent)
    return

if __name__ == "__main__":

    # Validating input requirements.
    file_path = sys.argv[-1] # We assume that this is valid and provided.
    points = read_file(file_path)
    N = len(points)

    if len(sys.argv) != 4: # General error: too many / too little arguments
        print("An Error Has Occured")
        sys.exit(1)

    try:
        k = int(sys.argv[1])
        assert k > 1 and k < N
    except:
        print("Invalid number of clusters!")
        sys.exit(1)

    try:
        iter = int(sys.argv[2])
        assert 1 < iter and iter < 1000
    except:
        print("Invalid maximum iteration!")
        sys.exit(1)

    try:
        centroids, _ = kmeans(N, k, iter, points)
        show_centroids(centroids)
    except:
        print("An Error Has Occured")
        sys.exit(1)
    sys.exit(0)