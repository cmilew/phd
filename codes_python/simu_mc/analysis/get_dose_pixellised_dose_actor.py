import os
import SimpleITK as sitk
import numpy as np
from openpyxl import load_workbook
import sys
import matplotlib.pyplot as plt


def get_strip_edep(edep_array):
    """Dose actor pixellised : 23 pixels/strip and 8 pixels/interstrip zones
    This function sum the dose in 23 adjacent pixels every 8 pixels starting by 8
    8 pixels (side )"""

    # mhd x y z => numpy z y x
    edep_array_x = edep_array[0][0]
    strip_edep = []
    n_strip_pix = 23
    n_interstrip_pix = 8
    total_pix = len(edep_array_x)
    pixel = 0

    while pixel < len(edep_array_x):
        print(f"pixel = {pixel}")
        if pixel % (n_strip_pix + n_interstrip_pix) != 0 or pixel == 0:
            n_pix_to_sum = n_interstrip_pix
            print(f"n_pix_to_sum = {n_pix_to_sum}")
        else:
            n_pix_to_sum = n_strip_pix
            print(f"n_pix_to_sum = {n_pix_to_sum}")
            print(f"sum = {np.sum(edep_array_x[pixel : pixel + n_pix_to_sum])}")
            strip_edep.append(np.sum(edep_array_x[pixel : pixel + n_pix_to_sum]))
        pixel += n_pix_to_sum

    return np.array(strip_edep)


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
path_edep = rf"C:\Users\milewski\Desktop\these\phd\codes_python\simu_mc\output\{slab_thickness}mm_{slab_material}_detector_edep.mhd"

# path_edep_uncertain = rf"C:\Users\milewski\Desktop\these\phd\codes_python\simu_mc\output\{slab_thickness}mm_{slab_material}_detector-edep-uncertainty.mhd"


# excel file for results
excel_file = os.path.join(
    os.path.dirname(__file__),
    r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\simulations_MC\results_step_phantom.xlsx",
)


# gets dose actor = img
edep = sitk.ReadImage(path_edep)
# edep_uncertain = sitk.ReadImage(path_edep_uncertain)

# converts dose actor img in numpy array
edep_array = sitk.GetArrayFromImage(edep)


array_test = np.zeros((1, 1, 93))  # Initialize with zeros
array_test[0, 0, 8:31] = 1  # Set the first block of 23 ones
array_test[0, 0, 39:62] = 1  # Set the second block of 23 ones
array_test[0, 0, 70:93] = 1  # Set the third block of 23 ones


strip_edep = get_strip_edep(edep_array)


strip_num = range(1, 137)

if plot_bool:
    slab_thickness = int(slab_thickness)
    plt.bar(strip_num, strip_edep)
    plt.xlabel("Strip number")
    plt.ylabel("Edep (MeV)")
    plt.title(f"Edep measured by strips for mono energy 123keV beam ESRF")
    plt.show()

# fill excel
if fill_excel_bool:
    fill_excel(excel_file, ws_name, strip_edep, 15, col_to_fill)
