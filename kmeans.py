import sys
from typing import List, Tuple
import math
EPSILON = 0.001


class Point:
    def __init__(self, coordinates: List[float]):
        self.coordinates = coordinates
        self.d = len(self.coordinates)
    def __repr__(self):
        return ','.join(['{:.4f}'.format(coord) for coord in self.coordinates])

    def __eq__(self, other: 'Point') -> bool:
        return self.coordinates == other.coordinates

    def __len__(self):
        return len(self.coordinates)

    def __getitem__(self, index: int) -> float:
        return self.coordinates[index]


def read_file(
        file_path:str
) -> List[Point]:
    """
    Reads .txt, removes newline characters and returns list of tuples of floats.
    Each tuple in the list corresponds to a datapoint.
    """
    with open(file_path, 'r') as file:
        points = []
        for line in file:
            parts = Point([float(val) for val in line.strip().split(',')])
            points.append(parts)
    return points


def euclidean_distance(
        v1: Point,
        v2: Point,
) -> float:

    squared_diff = [(a - b) ** 2 for a, b in zip(v1.coordinates, v2.coordinates)]
    distance = math.sqrt(sum(squared_diff))
    return distance

def find_closest_cluster(
        point: Point,
        centroids: List[Point],
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
        cluster: List[Point]
) -> Point:

    result = [0.0 for _ in range(len(cluster[0]))]

    for t in cluster:
        for i in range(len(t)):
            result[i] += t[i]
    result_divided = [i / float(len(cluster)) for i in result]
    result_point = Point(result_divided)
    return result_point
def update_centroids(
        centroids: List[Point],
        clusters: List[List[Point]],
) -> Tuple[List[Point], bool]:

    stop_flag = True
    for i in range(len(centroids)):
        corr_cluster = clusters[i] # List of points corresponding to cluster
        average_of_cluster = average(corr_cluster) # Returns Point
        delta = euclidean_distance(average_of_cluster, centroids[i])
        if stop_flag == True and delta >= EPSILON:
            stop_flag = False
        # Update value of centroid
        centroids[i] = average_of_cluster

    return centroids, stop_flag



def kmeans(
        N: int,
        k: int,
        iter: int,
        points: List[Point]
) -> Tuple[List[Point], List[List[Point]]]:

    assert len(points) == N
    centroids = points[:k]

    for i in range(iter):
        clusters = [[] for _ in range(k)]  # Initialize clusters as empty lists.
        for point in points:
            closest_cluster_index = find_closest_cluster(point, centroids)
            clusters[closest_cluster_index].append(point)

    # Updating centroids
        centroids, stop_flag = update_centroids(centroids, clusters)
        if stop_flag or i == 50:
            break
    return centroids, clusters

def show_centroids(
        centroids: List[Tuple[float, ...]],
) -> None:
    for cent in centroids:
        print(cent)

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