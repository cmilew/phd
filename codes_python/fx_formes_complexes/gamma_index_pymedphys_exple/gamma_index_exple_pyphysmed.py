import numpy as np
import matplotlib.pyplot as plt
import pydicom
import pymedphys


reference_filepath = r"C:\Users\milewski\Desktop\these\phd\codes_python\fx_formes_complexes\original_dose_beam_4.dcm"
evaluation_filepath = r"C:\Users\milewski\Desktop\these\phd\codes_python\fx_formes_complexes\logfile_dose_beam_4.dcm"

reference = pydicom.read_file(str(reference_filepath), force=True)
evaluation = pydicom.read_file(str(evaluation_filepath), force=True)

axes_reference, dose_reference = pymedphys.dicom.zyx_and_dose_from_dataset(reference)
axes_evaluation, dose_evaluation = pymedphys.dicom.zyx_and_dose_from_dataset(evaluation)

(z_ref, y_ref, x_ref) = axes_reference
(z_eval, y_eval, x_eval) = axes_evaluation

gamma_options = {
    "dose_percent_threshold": 1,
    "distance_mm_threshold": 1,
    "lower_percent_dose_cutoff": 20,
    "interp_fraction": 10,  # Should be 10 or more for more accurate results
    "max_gamma": 2,
    "random_subset": None,
    "local_gamma": True,
    "ram_available": 2**29,  # 1/2 GB
}

gamma = pymedphys.gamma(
    axes_reference, dose_reference, axes_evaluation, dose_evaluation, **gamma_options
)

print(z_ref)
