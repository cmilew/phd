import os
import SimpleITK as sitk
import numpy as np
from openpyxl import load_workbook
import sys
import matplotlib.pyplot as plt


def get_strip_edep(edep_array, edep_uncertain_array):
    """Dose actor pixellised : 23 pixels/strip and 8 pixels/interstrip zones
    This function sum the dose in 23 adjacent pixels every 8 pixels starting by 8
    8 pixels (side )"""

    # mhd x y z => numpy z y x
    edep_x = edep_array[0][0]
    uncertain_x = edep_uncertain_array[0][0]
    strip_edep = []
    strip_uncertain = []
    n_strip_pix = 23
    n_interstrip_pix = 8
    total_pix = len(edep_x)
    pixel = 0

    while pixel < total_pix:
        if pixel + n_strip_pix <= total_pix:
            strip_edep.append(np.sum(edep_x[pixel : pixel + n_strip_pix]))
            strip_uncertain.append(
                np.sqrt(np.sum(uncertain_x[pixel : pixel + n_strip_pix] ** 2))
            )
        else:
            break
        pixel += n_strip_pix + n_interstrip_pix

    return np.array(strip_edep), np.array(strip_uncertain)


def fill_excel(excel_path, ws_name, data_to_fill, start_l, start_col):
    wb_res = load_workbook(excel_path)
    ws = wb_res[ws_name]
    for i in range(len(data_to_fill)):
        ws.cell(row=start_l + i, column=start_col, value=data_to_fill[i])
    wb_res.save(excel_path)


### DATA TO FILL ######
slab_thickness = 0.0
slab_material = "RW3"
ws_name = "ESRF_edep_136_strips"
fill_excel_bool = False
plot_bool = True
col_to_fill = 4

# simulation files
# path_edep = os.path.join(
#     os.path.dirname(__file__),
#     rf"{slab_thickness}mm_{slab_material}_detector_edep_tot.mhd",
# )
# path_edep_uncertain = os.path.join(
#     os.path.dirname(__file__),
#     rf"{slab_thickness}mm_{slab_material}_detector_edep_uncertainty_tot.mhd",
# )


# gets dose actor = img
# edep = sitk.ReadImage(path_edep)
# edep_uncertain = sitk.ReadImage(path_edep_uncertain)

# converts dose actor img in numpy array
# edep_array = sitk.GetArrayFromImage(edep)
# edep_uncertain_array = sitk.GetArrayFromImage(edep_uncertain)
# strip_edep, strip_uncertain = get_strip_edep(edep_array, edep_uncertain_array)
array_test = np.zeros((1, 1, 93))  # Initialize with zeros
array_test[0, 0, 0:23] = 2  # Set the first block of 23 ones
array_test[0, 0, 31:54] = 3  # Set the second block of 23 ones
array_test[0, 0, 62:85] = 4  # Set the third block of 23 ones
test_array, test_uncertain_array = get_strip_edep(array_test, array_test)
print(array_test)
print(test_uncertain_array)

sys.exit()


if plot_bool:
    slab_thickness = int(slab_thickness)
    plt.bar(range(1, 137), strip_edep)
    # plt.yscale("log")
    plt.xlabel("Strip number")
    plt.ylabel("Edep (MeV)")
    plt.title("Edep measured by strips ESRF")
    plt.show()
