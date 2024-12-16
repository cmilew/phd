import numpy as np
import uproot
import os

m = np.float64(1)
cm = np.float64(1e-2)
mm = np.float64(1e-3)
um = np.float64(1e-6)
nm = np.float64(1e-9)

input_phsp_path = os.path.join(
    os.path.dirname(__file__), "parallel_beam_phsp_1.000E+2.root"
)
phsp_name = "parallel_beam_phsp_1.000E+2"
output_phsp_name = os.path.join(
    os.path.dirname(__file__), "divergent_beam_phsp_1.000E+2.root"
)

source_X = np.float64(0.0 * m)
source_Y = np.float64(0.0 * m)
source_Z = np.float64(-40.2 * m)


with uproot.open(input_phsp_path) as input_phsp_file:
    input_phsp_tree = input_phsp_file[phsp_name]
    phsp_keys = input_phsp_tree.keys()

input_phsp_branches = input_phsp_tree.arrays()

input_phsp_df = input_phsp_tree.arrays([*phsp_keys], library="pd")


norm = np.sqrt(
    (np.float64(input_phsp_df["PrePositionLocal_X"]) * mm - source_X) ** 2
    + (np.float64(input_phsp_df["PrePositionLocal_Y"]) * mm - source_Y) ** 2
    + (np.float64(input_phsp_df["PrePositionLocal_Z"]) * mm - source_Z) ** 2
)
dX = (np.float64(input_phsp_df["PrePositionLocal_X"]) * mm - source_X) / norm
dY = (np.float64(input_phsp_df["PrePositionLocal_Y"]) * mm - source_Y) / norm
dZ = (np.float64(input_phsp_df["PrePositionLocal_Z"]) * mm - source_Z) / norm

input_phsp_df["PreDirectionLocal_X"] = dX
input_phsp_df["PreDirectionLocal_Y"] = dY
input_phsp_df["PreDirectionLocal_Z"] = dZ

with uproot.recreate(output_phsp_name) as output_phsp_file:
    output_phsp_file.mktree(
        phsp_name,
        {
            phsp_keys[0]: np.float64,
            phsp_keys[1]: np.float64,
            phsp_keys[2]: np.float64,
            phsp_keys[3]: np.float64,
            phsp_keys[4]: np.float64,
            phsp_keys[5]: np.float64,
            phsp_keys[6]: np.float64,
            phsp_keys[7]: np.float64,
        },
    )

    output_phsp_file[phsp_name].extend(
        {
            phsp_keys[0]: input_phsp_df[phsp_keys[0]],
            phsp_keys[1]: input_phsp_df[phsp_keys[1]],
            phsp_keys[2]: input_phsp_df[phsp_keys[2]],
            phsp_keys[3]: input_phsp_df[phsp_keys[3]],
            phsp_keys[4]: input_phsp_df[phsp_keys[4]],
            phsp_keys[5]: input_phsp_df[phsp_keys[5]],
            phsp_keys[6]: input_phsp_df[phsp_keys[6]],
            phsp_keys[7]: input_phsp_df[phsp_keys[7]],
        }
    )
