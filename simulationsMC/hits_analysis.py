import os
# import opengate as gate
# import numpy as np
# import pandas as pd
import uproot
import matplotlib.pyplot as plt 

f = uproot.open(os.path.join(os.path.dirname(__file__), "central_strip_hits.root"))
# energy_deposited = f['central_strip_hits;1']['TotalEnergyDeposit'].array().to_numpy()
# plt.hist(f['central_strip_hits;1']['TotalEnergyDeposit'].array(), bins=100, range=(0, 1.4))
# plt.show()
energy_deposited = f['hits_central_strip;1']['TotalEnergyDeposit'].array().to_numpy()
plt.hist(f['hits_central_strip;1']['TotalEnergyDeposit'].array(), bins=100, range=(0, 1.4))
plt.s