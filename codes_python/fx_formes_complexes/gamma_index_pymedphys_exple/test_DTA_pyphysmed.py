import numpy as np
import pymedphys
import sys


def calc_DTA(ref_coords, comp_coords, distance_tolerance):
    # reformat the coordinates, z = 1 because 2D comparison
    axes_reference = (
        ref_coords[:, 1],
        ref_coords[:, 0],
    )
    axes_evaluation = (comp_coords[:, 1], comp_coords[:, 0])

    # Create dose grids filled with ones to simulate uniform doses (for distance-only gamma)
    dose_reference = np.ones((len(ref_coords[:, 1]), len(ref_coords[:, 0])))
    dose_evaluation = np.ones((len(comp_coords[:, 1]), len(comp_coords[:, 0])))

    gamma_options = {
        "dose_percent_threshold": 0,
        "distance_mm_threshold": 1,
        "lower_percent_dose_cutoff": 0,
        "interp_fraction": 10,  # Should be 10 or more for more accurate results
        "max_gamma": 2,
        "random_subset": None,
        "local_gamma": True,
        "ram_available": 2**29,  # 1/2 GB
    }

    gamma = pymedphys.gamma(
        axes_reference,
        dose_reference,
        axes_evaluation,
        dose_evaluation,
        **gamma_options,
    )
    return gamma


# Example usage
reference_array = np.array([[0, 0], [0, 1], [1, 0], [1, 1]])
comparison_array = np.array([[0, 0], [0, 1], [1, 0.5], [1, 1.5]])

# sort arrays so that y coord is ascending
sorted_reference_array = reference_array[np.argsort(reference_array[:, 1])]
sorted_comparison_array = comparison_array[np.argsort(comparison_array[:, 1])]


distance_tolerance = 3  # Distance-to-agreement tolerance in mm

gamma_result = calc_DTA(
    sorted_reference_array, sorted_comparison_array, distance_tolerance
)
print("Gamma Index result (distance-only):", gamma_result)
