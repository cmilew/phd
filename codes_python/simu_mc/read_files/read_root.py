import uproot
import glob
import numpy as np
import os
import matplotlib.pyplot as plt


phsp_name = "phsp_end_ESRF_line1.000E+7_events"
file_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    f"{phsp_name}.root",
)


with uproot.open(file_path) as file:
    phsp_name = f"{phsp_name};1"
    edep = file[phsp_name]["KineticEnergy"].array(library="np")
    # pos_x = file[phsp_name]["PrePositionLocal_X"].array(library="np")
    # pos_y = file[phsp_name]["PrePositionLocal_Y"].array(library="np")

print(f"nb of particles: {len(edep)}")
plt.hist(edep, bins=100, color="skyblue", edgecolor="black")
plt.xlabel("Energy (MeV)")
plt.ylabel("Counts")
plt.title("Spectrum phsp post msc")
plt.show()

# plt edep in function of x and y
# plt.scatter(pos_x, pos_y, s=0.5, c=edep, cmap="viridis")
# cbar = plt.colorbar()
# cbar.set_label("Energy (MeV)", rotation=270, labelpad=15)
# plt.ylim(-0.5, 0.5)
# plt.xlim(-18, 18)
# plt.xlabel("x position (mm)")
# plt.ylabel("y position (mm)")
# plt.title(phsp_name)
# plt.show()
