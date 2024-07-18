# Imports.
import argparse
import sys
import numpy as np
from typing import List, Any
import symnmf_c_api as sym
from numpy import floating

# Set seeds and env vars.
np.random.seed(0)
ERROR_MSG = 'An Error Has Occurred'


# Helper functions.
def file_to_lists(
        file_name: str,
) -> List[List[float]]:
    """
    Read data from .txt file, such that each point is represented as a list of floats.
    """

    data_matrix_list = []
    try:
        with open(file_name, 'r') as file:
            for line in file:
                row = line.strip().split(',')
                data_matrix_list.append([float(value) for value in row])

    except FileNotFoundError:
        print(ERROR_MSG)
        sys.exit(1)
    return data_matrix_list

def show_matrix(
        matrix: List[List[float]]
) -> None:
    """
    Prints matrix with point entries delimited by commas,
    points delimited by newlines,
    and point values truncated to .4f.
    """
    for row in matrix:
        formatted_list = ["{:.4f}".format(num) for num in row]
        print(','.join(map(str, formatted_list)))

def initialize_H_matrix(
        n: int,
        k:int,
        m:float,
) -> List[List[float]]:
    """
    Creates a n by k initial H matrix, populated by random values from [0,2*sqrt(m/k)],
    where m is the average value of the W matrix.
    """
    initial_h = []
    for i in range(n):
        initial_h.append([])
        for j in range(k):
            initial_h[i].append(2 * np.sqrt(m / k) * np.random.uniform())
    return initial_h


def average_value_over_matrix(
        matrix: List[List[float]]
):
    """
    Return average value across all matrix entries.
    """
    np_matrix = np.array(matrix)
    return np.mean(np_matrix)

def compute_symnmf(
        n:int,
        k:int,
        d:int,
        data_matrix: List[List[float]],
) -> List[List[float]]:
    """
    Given a data matrix of size n by k, with each point of dimension d:
        1. Computes W, the normalized similarity matrix
        2. Initializes H based on m, the average entry value of W
        3. Returns the best H by running symnmf optimization starting from the initial H
    """

    W = sym.norm(n, d, data_matrix)
    m = average_value_over_matrix(W)
    initial_H = initialize_H_matrix(n, k, m)

    return sym.symnmf(n, k, W, initial_H) # Best H following optimization.


if __name__ == '__main__':
    """
    As per instructions, we assume argument inputs to be valid.
    """

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

    # Parsing user arguments.
    try:
        args = parser.parse_args()
    except SystemExit as e:
        print(ERROR_MSG)
        sys.exit(1)

    k = args.k
    goal = args.goal
    file_name = args.file_name
    data_matrix = file_to_lists(file_name)
    n = len(data_matrix)

    try:
        d = len(data_matrix[0]) # We assume that all points have the same d.
    except:
        print(ERROR_MSG)
        sys.exit(1)

    if goal == 'symnmf':
        required_matrix = compute_symnmf(n, k, d, data_matrix)
    elif goal == 'sym':
        required_matrix = sym.sym(n, d, data_matrix)
    elif goal == 'ddg':
        required_matrix = sym.ddg(n, d, data_matrix)
    else: # Implies goal == 'norm', thanks to previous try/except block.
        required_matrix = sym.norm(n, d, data_matrix)

    show_matrix(required_matrix)
