import numpy as np
import itk
import glob
import os
import sys

# get path of root files
path = "./run."
folder_list = glob.glob(path + "*")

# folder_path = glob.glob(path, recursive=True)
assert len(folder_list) == 1
# only one run folder
folder_path = folder_list[0]
# dict containing name of mhd files in keys and pixel sums in values
output = dict()

for mhd_path in glob.glob(os.path.join(folder_path, "output.*/*.mhd"), recursive=True):
    # converts img of dose actor in array
    mhd_array = itk.GetArrayFromImage(itk.imread(mhd_path))
    mhd_name = mhd_path.split("/")[-1]

    if mhd_name.endswith("uncertainty.mhd"):
        # if uncertainty sum the square of the values
        if mhd_name in output.keys():
            output[mhd_name] += mhd_array**2
        else:
            output[mhd_name] = mhd_array**2
    else:
        if mhd_name in output.keys():
            output[mhd_name] += mhd_array
        else:
            output[mhd_name] = mhd_array
print(f"output.keys() = {output.keys()}")
for k, v in output.items():
    if k.endswith("uncertainty"):
        # propagation of errors
        output[k] = np.sqrt(v)

    # creates new output mhd files with pixel sums
    img_tot = itk.GetImageFromArray(output[k])

    itk.imwrite(img_tot, os.path.join(folder_path, k[:-4] + "_tot.mhd"))
