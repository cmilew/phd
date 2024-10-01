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

rw3_thickness = 1
strip_number = 13

path_root_file = os.path.join(
    os.path.dirname(__file__), f"output/{rw3_thickness}cm_rw3_strip_{strip_number}.root"
)
with uproot.open(path_root_file) as file:
    hits_name = f"hits_{strip_number};1"
    print(f"Entries {file[hits_name]['TotalEnergyDeposit'].num_entries}")
    print(f"Edep {np.sum(file[hits_name]['TotalEnergyDeposit'].array(library='np'))}")
    print(f"Std {np.std(file[hits_name]['TotalEnergyDeposit'].array(library='np'))}")
