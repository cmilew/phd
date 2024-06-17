import os
# import opengate as gate
# import numpy as np
# import pandas as pd
import uproot
import matplotlib.pyplot as plt 

f = uproot.open(os.path.join(os.path.dirname(__file__), "detector_hits.root"))
energy_deposited = f['hits_detector;1']['TotalEnergyDeposit'].array().to_numpy()
plt.hist(f['hits_detector;1']['TotalEnergyDeposit'].array(), bins=100, range=(0, 1.4))
plt.show()