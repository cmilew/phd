import numpy as np


def generate_shifted_matrix(input_array):
    """Function generating a matrix where the given array is shifted at each row and
    missing values are filled with zeros"""
    n = len(input_array)
    result_matrix = np.zeros((n, n), dtype=int)
    if n % 2 == 0:
        center_index = n // 2 - 1
        shift = 1
    else:
        center_index = n // 2
        shift = 0
    index_result_matrix = 0

    for i in range(center_index + 1):
        result_matrix[index_result_matrix][: center_index + i + 1 + shift] = (
            input_array[center_index - i :]
        )
        index_result_matrix += 1
    for i in range(center_index + shift):
        result_matrix[index_result_matrix][i + 1 :] = input_array[0 : n - i - 1]
        index_result_matrix += 1
    return result_matrix


# Test the function with the provided input array
input_array = [2, 3, 2, 1]
result = generate_shifted_matrix(input_array)
print(result)
