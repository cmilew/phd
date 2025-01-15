import numpy as np
import matplotlib.pyplot as plt
import uproot
import pandas

m = np.float64(1)
cm = np.float64(1e-2)
mm = np.float64(1e-3)
um = np.float64(1e-6)
nm = np.float64(1e-9)

input_phsp_name = "output/parallel_source.root"
output_phsp_name = "output/divergent_source.root"

source_X = np.float64(0.*m)
source_Y = np.float64(0.*m)
source_Z = np.float64(-40.*m)

with uproot.open(input_phsp_name) as input_phsp_file:
        input_phsp_tree = input_phsp_file['PhaseSpace']
        phsp_keys = input_phsp_tree.keys()

print(phsp_keys)
input_phsp_branches = input_phsp_tree.arrays()

for key in phsp_keys:
        print(type(input_phsp_branches[key][0]))

input_phsp_df = input_phsp_tree.arrays([*phsp_keys], library='pd')


norm = np.sqrt((np.float64(input_phsp_df['X'])*mm - source_X)**2 + (np.float64(input_phsp_df['Y'])*mm - source_Y)**2 + (np.float64(input_phsp_df['Z'])*mm - source_Z)**2)
dX = (np.float64(input_phsp_df['X'])*mm - source_X) / norm
dY = (np.float64(input_phsp_df['Y'])*mm - source_Y) / norm
dZ = (np.float64(input_phsp_df['Z'])*mm - source_Z) / norm

input_phsp_df['dX'] = dX
input_phsp_df['dY'] = dY
input_phsp_df['dZ'] = dZ

with uproot.recreate(output_phsp_name) as output_phsp_file:
        output_phsp_file.mktree("PhaseSpace", {phsp_keys[0]: np.int32,
                                               phsp_keys[1]: np.float64,
                                               phsp_keys[2]: np.float32,
                                               phsp_keys[3]: np.float32,
                                               phsp_keys[4]: np.int32,
                                               phsp_keys[5]: np.float32,
                                               phsp_keys[6]: np.float32,
                                               phsp_keys[7]: np.float32,
                                               phsp_keys[8]: np.float32,
                                               phsp_keys[9]: np.float32,
                                               phsp_keys[10]: np.float32,
                                               phsp_keys[11]: str,
                                               phsp_keys[12]: np.int32,
                                               phsp_keys[13]: np.int32,
                                               phsp_keys[14]: np.int32,
                                               phsp_keys[15]: np.int32 })

        output_phsp_file['PhaseSpace'].extend({phsp_keys[0]: input_phsp_df[phsp_keys[0]],
                                               phsp_keys[1]: input_phsp_df[phsp_keys[1]],
                                               phsp_keys[2]: input_phsp_df[phsp_keys[2]],
                                               phsp_keys[3]: input_phsp_df[phsp_keys[3]],
                                               phsp_keys[4]: input_phsp_df[phsp_keys[4]],
                                               phsp_keys[5]: input_phsp_df[phsp_keys[5]],
                                               phsp_keys[6]: input_phsp_df[phsp_keys[6]],
                                               phsp_keys[7]: input_phsp_df[phsp_keys[7]],
                                               phsp_keys[8]: input_phsp_df[phsp_keys[8]],
                                               phsp_keys[9]: input_phsp_df[phsp_keys[9]],
                                               phsp_keys[10]: input_phsp_df[phsp_keys[10]],
                                               phsp_keys[11]: input_phsp_df[phsp_keys[11]],
                                               phsp_keys[12]: input_phsp_df[phsp_keys[12]],
                                               phsp_keys[13]: input_phsp_df[phsp_keys[13]],
                                               phsp_keys[14]: input_phsp_df[phsp_keys[14]],
                                               phsp_keys[15]: input_phsp_df[phsp_keys[15]] })
                                               
                                               




