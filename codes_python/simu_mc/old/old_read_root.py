import uproot
import glob
import numpy as np
import os
import matplotlib.pyplot as plt

# # get path of root files
# path = "./rootfile."
# folder_list = glob.glob(path + "*")
# name_file = "/*.root"

# for i in range(len(folder_list)):
#     file_path = glob.glob(folder_list[i] + name_file)
#     assert len(file_path) == 1, "More than one root file in rootfile folder"
#     file_path = file_path[0]

phsp_name_1 = "phsp_end_ESRF_line1.000E+5_events"
# phsp_name_2 = "phsp_test_post_col"
phsp_name_3 = "phsp_test_post_msc"
file_path = os.path.join(os.path.dirname(__file__), f"output/{phsp_name_1}.root")
# with uproot.open(file_path) as file:
#     # print(f"Entries {file['PhaseSpace;1']['KineticEnergy'].num_entries}")
#     # print(
#     #     f"Edep {np.mean(file['PhaseSpace;1']['KineticEnergy'].array(library='np'))}"
#     # )
#     # print(
#     #     f"Std {np.std(file['PhaseSpace;1']['KineticEnergy'].array(library='np'))}"
#     # )
#     phsp_name = f"{phsp_name_1};1"
#     edep = file[phsp_name]["KineticEnergy"].array(library="np")
#     pos_x = file[phsp_name]["PrePositionLocal_X"].array(library="np")
#     pos_y = file[phsp_name]["PrePositionLocal_Y"].array(library="np")

# # plt edep in function of x and y
# plt.scatter(pos_x, pos_y, s=1, c=edep, cmap="viridis")
# cbar = plt.colorbar()
# cbar.set_label("Energy (MeV)", rotation=270, labelpad=15)
# plt.xlabel("x position (mm)")
# plt.ylabel("y position (mm)")
# # plt.xlim(-20, 20)
# plt.title(phsp_name_1)
# plt.show()

# print(f"edep mean phsp end line = {np.mean(edep)}")
# print(f"nb particules phps end line {len(edep)}")

# plt.hist(edep, bins=100, color="skyblue", edgecolor="black")
# plt.xlabel("Energy (MeV)")
# plt.ylabel("Counts")
# plt.title("Spectrum phsp pre col")
# plt.show()

# file_path = os.path.join(os.path.dirname(__file__), f"output/{phsp_name_2}.root")
# with uproot.open(file_path) as file:
#     phsp_name = f"{phsp_name_2};1"
#     edep = file[phsp_name]["KineticEnergy"].array(library="np")
#     pos_x = file[phsp_name]["PrePositionLocal_X"].array(library="np")
#     pos_y = file[phsp_name]["PrePositionLocal_Y"].array(library="np")

# print(f"edep mean phsp post col = {np.mean(edep)}")
# print(f"nb particules phps {len(edep)}")

# # plt edep in function of x and y
# plt.scatter(pos_x, pos_y, c=edep, cmap="viridis")
# cbar = plt.colorbar()
# cbar.set_label("Energy (MeV)", rotation=270, labelpad=15)
# plt.xlabel("x position (mm)")
# plt.ylabel("y position (mm)")
# plt.title(phsp_name_2)
# plt.show()

# plt.hist(edep, bins=100, color="skyblue", edgecolor="black")
# plt.xlabel("Energy (MeV)")
# plt.ylabel("Counts")
# plt.title("Spectrum phsp post col")
# plt.show()


file_path = os.path.join(os.path.dirname(__file__), f"output/{phsp_name_3}.root")
with uproot.open(file_path) as file:
    phsp_name = f"{phsp_name_3};1"
    edep = file[phsp_name]["KineticEnergy"].array(library="np")
    pos_x = file[phsp_name]["PrePositionLocal_X"].array(library="np")
    pos_y = file[phsp_name]["PrePositionLocal_Y"].array(library="np")

print(f"edep mean phsp post source = {np.mean(edep)}")
print(f"nb particules phps {len(edep)}")

# plt edep in function of x and y
plt.scatter(pos_x, pos_y, s=0.5, c=edep, cmap="viridis")
cbar = plt.colorbar()
cbar.set_label("Energy (MeV)", rotation=270, labelpad=15)
plt.ylim(-0.5, 0.5)
plt.xlim(-18, 18)
plt.xlabel("x position (mm)")
plt.ylabel("y position (mm)")
plt.title(phsp_name_3)
plt.show()

# plt.hist(edep, bins=100, color="skyblue", edgecolor="black")
# plt.xlabel("Energy (MeV)")
# plt.ylabel("Counts")
# plt.title("Spectrum phsp post msc")
# plt.show()
