import os
import SimpleITK as sitk
from openpyxl import load_workbook

# DATA TO FILL IN ######################
rw3_thickness = 0
output_path = r"C:\Users\milewski\Desktop\these\papiers\caracterisation_detecteur_153_voies\simulations_MC\test_det_135_strips"
excel_path = r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\simulations_MC\results_step_phantom.xlsx"
ws_name = "edep_135_strips"
start_l = 6
########################################

wb_res = load_workbook(excel_path)
ws = wb_res[ws_name]

for filename in os.listdir(output_path):
    if filename.endswith(".mhd"):
        # gets strip number
        strip_number = int("".join(filter(str.isdigit, filename)))
        print(f"strip number: {strip_number}")

        # reads dose actor file
        file_path = os.path.join(output_path, filename)
        print(f"file path: {file_path}")
        img = sitk.ReadImage(file_path)

        # converts img of dose actor in array (only one pixel)
        strip_value = sitk.GetArrayFromImage(img)
        if "uncertainty" in filename.lower():
            ws.cell(row=strip_number + start_l, column=6, value=strip_value[0][0][0])
            print(f"uncertainty strip value: {strip_value}")
        else:
            ws.cell(row=strip_number + start_l, column=5, value=strip_value[0][0][0])
            print(f"edep strip value: {strip_value}")

wb_res.save(excel_path)
# path_edep = os.path.join(os.path.dirname(__file__), f"{rw3_thickness}cm_w-edep.mhd")
# path_edep_uncertain = os.path.join(
#     os.path.dirname(__file__), f"{rw3_thickness}cm_w-edep-uncertainty.mhd"
# )

# # gets dose actor = img
# img_edep = sitk.ReadImage(path_edep)
# img_edep_uncertain = sitk.ReadImage(path_edep_uncertain)

# # converts img dose actor img in numpy array (only one pixel)
# edep = sitk.GetArrayFromImage(img_edep)
# edep_unertain = sitk.GetArrayFromImage(img_edep_uncertain)
