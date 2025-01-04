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

    while pixel < total_pix:
        if pixel + n_strip_pix <= total_pix:
            strip_edep.append(np.sum(edep_array_x[pixel : pixel + n_strip_pix]))
        else:
            break
        pixel += n_strip_pix + n_interstrip_pix

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
path_edep = os.path.join(
    os.path.dirname(__file__),
    rf"{slab_thickness}mm_{slab_material}_detector_edep.mhd",
)

# gets dose actor = img
edep = sitk.ReadImage(path_edep)
# edep_uncertain = sitk.ReadImage(path_edep_uncertain)

# converts dose actor img in numpy array
edep_array = sitk.GetArrayFromImage(edep)
strip_edep = get_strip_edep(edep_array)


if plot_bool:
    slab_thickness = int(slab_thickness)
    plt.bar(range(1, 137), strip_edep)
    # plt.yscale("log")
    plt.xlabel("Strip number")
    plt.ylabel("Edep (MeV)")
    plt.title("Edep measured by strips ESRF")
    plt.show()
