import argparse
import numpy as np
from typing import List
import symnmf_c_api as sym
from kmeans import Point

# Set random seed for reproducibility
np.random.seed(0)


# Helper functions.
def get_points_from_file(
        file_name: str
) -> List[Point]:
    """
    Opens .txt file, returns matrix (List[Point) containing data.
    """

    matrix = []
    with open(file_name, 'r') as file:
        for line in file:
            row = line.strip().split(',')
            matrix.append(Point([float(value) for value in row]))

    return matrix

def show_matrix(matrix: List[List[float]]) -> None:
    for row in matrix:
        formatted_row = ['{:.4f}'.format(num) for num in row]
        print(','.join(formatted_row))

def initialize_H_matrix(
        n:int,
        k:int,
        m:float,
) -> List[List[float]]:
    upper_bound = 2 * np.sqrt(m / k)
    return np.random.uniform(0, upper_bound, (n, k)).tolist()


def average_value_over_matrix(
        matrix: List[List[float]]
) -> np.floating:
    """
    Return average value across all matrix entries.
    """
    np_matrix = np.array(matrix)
    return np.mean(np_matrix)

def compute_symnmf(
        n:int,
        k:int,
        d:int,
        data_matrix: List[Point],
) -> List[List[float]]:

    data_matrix_lists = [point.coordinates for point in data_matrix]

    W = sym.norm(
        n,
        d,
        data_matrix_lists,
    )

    m = average_value_over_matrix(W)
    initial_H = initialize_H_matrix(n, k, m)  # Returns List[List[float]]
    x = sym.symnmf(
        n,
        k,
        W,
        initial_H,
    ) # Returns best H
    return x

def derive_clustering_solution_symnmf(
        data_matrix: List[Point],
        best_H: List[List[float]],
        k: int,
) -> List[List[float]]:
    """
    Derives list of lists, each inner list corresponding to a cluster, based on
    the association matrix H.
    To be used in analysis.
    """
    clusters = [[] for _ in range(k)] # Initialize clusters as empty lists.

    for i in range(len(data_matrix)):
        datapoint = data_matrix[i]
        association_row = best_H[i]
        cluster_index = association_row.index(max(association_row))
        clusters[cluster_index].append(datapoint)
    return clusters




if __name__ == '__main__': # We currently assume inputs are valid.

    parser = argparse.ArgumentParser()

    parser.add_argument(
        'k',
        type=int,
        help='An integer representing the number of clusters'
    )
    parser.add_argument(
        'goal',
        choices=['symnmf','sym','ddg','norm'],
        help='Indicates which function to run.'
    )
    parser.add_argument(
        'file_name',
        type=str,
        help='The path to the .txt file assumed to contain N valid points.'
    )
    args = parser.parse_args()
    k = args.k
    goal = args.goal
    file_name = args.file_name


    data_matrix = get_points_from_file(file_name) # List of Point values
    n = len(data_matrix)
    try:
        d = data_matrix[0].d # All values are supposed to have the same d.
    except:
        print('An Error Has Occured')

    if goal == 'symnmf':
        required_matrix = compute_symnmf(n, k, d, data_matrix)
    else:
        data_matrix_lists = [point.coordinates for point in data_matrix]
        if goal == 'sym':
            required_matrix = sym.sym(n, d, data_matrix_lists)
        elif goal == 'ddg':
            required_matrix = sym.ddg(n, d, data_matrix_lists)
        else: # goal == 'norm' - inputs considered valid
            required_matrix = sym.norm(n,d,data_matrix_lists)

    show_matrix(required_matrix)
