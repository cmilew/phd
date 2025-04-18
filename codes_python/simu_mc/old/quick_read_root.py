import uproot
import numpy as np
import os
import matplotlib.pyplot as plt

phsp_name = "phsp_end_ESRF_line1.000E+5_events"
file_path = os.path.join(os.path.dirname(__file__), f"data/{phsp_name}.root")

# get root data
with uproot.open(file_path) as file:
    phsp_name = f"{phsp_name};1"
    edep = file[phsp_name]["KineticEnergy"].array(library="np")
    pos_x = file[phsp_name]["PrePositionLocal_X"].array(library="np")
    pos_y = file[phsp_name]["PrePositionLocal_Y"].array(library="np")

# plt edep in function of x and y
plt.scatter(pos_x, pos_y, s=1, c=edep, cmap="viridis")
cbar = plt.colorbar()
cbar.set_label("Energy (MeV)", rotation=270, labelpad=15)
plt.xlabel("x position (mm)")
plt.ylabel("y position (mm)")
plt.title(phsp_name)
plt.show()

print(f"edep mean phsp end line = {np.mean(edep)}")
print(f"nb particules phsp end line {len(edep)}")
