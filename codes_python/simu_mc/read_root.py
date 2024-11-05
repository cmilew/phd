import uproot
import glob
import numpy as np
import os


# path = "./rootfile."
# folder_list = glob.glob(path + "*")
# name_file = "/*.root"
# for i in range(len(folder_list)):
#     file_path = glob.glob(folder_list[i] + name_file)
#     assert len(file_path) == 1, "More than one root file in rootfile folder"
#     file_path = file_path[0]
#     print(f"Num entries for {file_path} :")
#     with uproot.open(file_path) as file:
#         print(
#             f"Entries {file['hits_central_strip;1']['TotalEnergyDeposit'].num_entries}"
#         )
#         print(
#             f"Edep {np.sum(file['hits_central_strip;1']['TotalEnergyDeposit'].array(library='np'))}"
#         )
#         print(
#             f"Std {np.std(file['hits_central_strip;1']['TotalEnergyDeposit'].array(library='np'))}"
#         )

file_name = "phsp_esrf_line_1.000E+2_events"

path_root_file = os.path.join(os.path.dirname(__file__), f"output/{file_name}.root")
with uproot.open(path_root_file) as file:
    print(f"Entries {file['PhaseSpace;1']['KineticEnergy'].num_entries}")
    print(f"Edep {np.mean(file['PhaseSpace;1']['KineticEnergy'].array(library='np'))}")
    print(f"Std {np.std(file['PhaseSpace;1']['KineticEnergy'].array(library='np'))}")
