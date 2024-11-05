import uproot
import glob
import numpy as np
import os

# get path of root files
path = "./rootfile."
folder_list = glob.glob(path + "*")
name_file = "/*.root"

for i in range(len(folder_list)):
    file_path = glob.glob(folder_list[i] + name_file)
    assert len(file_path) == 1, "More than one root file in rootfile folder"
    file_path = file_path[0]
    with uproot.open(file_path) as file:
        print(f"Entries {file['PhaseSpace;1']['KineticEnergy'].num_entries}")
        print(
            f"Edep {np.mean(file['PhaseSpace;1']['KineticEnergy'].array(library='np'))}"
        )
        print(
            f"Std {np.std(file['PhaseSpace;1']['KineticEnergy'].array(library='np'))}"
        )
