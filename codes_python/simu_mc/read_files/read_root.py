import uproot
import sys
import numpy as np
import os
import matplotlib.pyplot as plt


phsp_name = "phsp_ANSTO_poly_CuCu_1.000E+4_events"
# phsp_name = "phsp_beh_det"
# phsp_name = "phsp_end_ANSTO_line1.000E+4_events"
beam_width = 35  # mm
beam_height = 1.059  # mm
file_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    f"{phsp_name}.root",
)


with uproot.open(file_path) as file:
    phsp_name = f"{phsp_name};1"
    edep = file[phsp_name]["KineticEnergy"].array(library="np")
    pos_x = file[phsp_name]["PrePositionLocal_X"].array(library="np")
    pos_y = file[phsp_name]["PrePositionLocal_Y"].array(library="np")

print(f"nb of particles: {len(edep)}")
print(f"mean energy: {np.mean(edep):.6f} MeV")
plt.hist(edep, bins=100, color="skyblue", edgecolor="black")
plt.xlabel("Energy (MeV)")
plt.ylabel("Counts")
plt.title("Spectrum phsp post msc")
plt.show()

# # plt edep in function of x and y
# # plt.scatter(pos_x, pos_y, s=0.5, c=edep, cmap="RdPu")
# plt.scatter(pos_x, pos_y, s=0.5, c=edep, cmap="viridis")
# cbar = plt.colorbar()
# cbar.set_label("Energy (MeV)", rotation=270, labelpad=15)
# plt.ylim(-beam_height / 2 - 1, beam_height / 2 + 1)
# plt.xlim(-beam_width - 1, beam_width + 1)
# plt.xlabel("x position (mm)")
# plt.ylabel("y position (mm)")
# plt.title(phsp_name)
# plt.show()
