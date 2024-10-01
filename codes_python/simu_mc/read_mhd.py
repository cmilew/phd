import os
import SimpleITK as sitk

rw3_thickness = 0
strip_number = 8

path_edep = os.path.join(
    os.path.dirname(__file__),
    f"output/{rw3_thickness}cm_rw3_strip_{strip_number}-edep.mhd",
)

path_edep_uncertain = os.path.join(
    os.path.dirname(__file__),
    f"output/{rw3_thickness}cm_rw3_strip_{strip_number}-edep-uncertainty.mhd",
)

# gets dose actor = img
edep = sitk.ReadImage(path_edep)
edep_uncertain = sitk.ReadImage(path_edep_uncertain)

# converts dose actor img in numpy array
edep_array = sitk.GetArrayFromImage(edep)
edep_unertain_array = sitk.GetArrayFromImage(edep_uncertain)

# only one pixel in array
print(f"edep strip {strip_number} = {edep_array}")
print(f"edep_uncertain {strip_number} = {edep_unertain_array}")
