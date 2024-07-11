import os

# import opengate as gate
# import numpy as np
# import pandas as pd
import uproot
import SimpleITK as sitk
import matplotlib.pyplot as plt

# f = uproot.open(os.path.join(os.path.dirname(__file__), "central_strip_hits.root"))
# energy_deposited = f['hits_central_strip;1']['TotalEnergyDeposit'].array().to_numpy()
# plt.hist(f['hits_central_strip;1']['TotalEnergyDeposit'].array(), bins=100, range=(0, 1.4))
# plt.show()
water_thickness = 3

path_edep = (
    f"/home/candice/Documents/phd/diamond_detector_{water_thickness}cm_water-edep.mhd"
)
path_edep_uncertain = f"/home/candice/Documents/phd/diamond_detector_{water_thickness}cm_water-edep-uncertainty.mhd"
# path_edep = "/home/candice/Documents/phd/central_strip_dose-dose.mhd"
# path_edep_uncertain = "/home/candice/Documents/phd/central_strip_dose-dose-uncertainty.mhd"

# gets dose actor = img
edep = sitk.ReadImage(path_edep)
edep_uncertain = sitk.ReadImage(path_edep_uncertain)

# converts dose actor img in numpy array
edep_array = sitk.GetArrayFromImage(edep)
edep_unertain_array = sitk.GetArrayFromImage(edep_uncertain)
print(edep_array)
print(edep_unertain_array)
