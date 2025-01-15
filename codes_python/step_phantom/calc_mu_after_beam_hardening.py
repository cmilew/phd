import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score
import pandas as pd
import sys
from decimal import Decimal
from scipy.interpolate import interp1d, CubicSpline
import pandas as pd

font_size = 15

# Synchrotron spectrum
spectrum_path = r"C:\Users\milewski\Desktop\these\papiers\caracterisation_detecteur_153_voies\simulations_MC\ESRF_spectrum.xlsx"
df = pd.read_excel(spectrum_path, skiprows=6)
spectrum_energy = df["Energy [eV]"].to_numpy()
spectrum_energy = spectrum_energy * 10**-6  # eV to MeV
spectrum_weight = df["Flux_3_intensity"].to_numpy()

# # plot ESRF spectrum
# plt.figure()
# plt.plot(spectrum_energy, spectrum_weight, linewidth=5)
# plt.xlabel("Energy (MeV)", fontsize=font_size)
# plt.tick_params(axis="both", which="major", labelsize=font_size)
# plt.ylabel("Intensity (%)", fontsize=font_size)
# plt.show()

# x-ray mass attenuation coeff in water (NIST data)
df = pd.read_excel(spectrum_path, skiprows=17)
nist_energy = df["E (MeV)"][:25].to_numpy()
nist_mu_rw3 = df["mu(RW3) cm-1"][:25].to_numpy()
spectrum_weight_harden_6cm = []

# calculation of spectrum weight for each energy after exponential attenuation in 6 cm of RW3
lin_interpol_mu = interp1d(
    nist_energy, nist_mu_rw3, kind="linear", fill_value="extrapolate"
)
# lin_intepol_mu = CubicSpline(
#     nist_energy,
#     nist_mu_rw3,
#     axis=0,
#     bc_type="not-a-knot",
#     extrapolate=True,
# )
# for iweight, energy in enumerate(spectrum_energy):
#     interp_mu = lin_interpol_mu(energy)
#     spectrum_weight_harden_6cm.append(spectrum_weight[iweight] * np.exp(-interp_mu * 6))

spectrum_weight_harden_6cm = (
    np.exp(-lin_interpol_mu(spectrum_energy) * 0) * spectrum_weight
)

# new mean energy of spectrum hardened by 6 cm of water
mean_energy_6cm_hardening = (
    spectrum_energy * spectrum_weight_harden_6cm
).sum() / spectrum_weight_harden_6cm.sum()


print(f"Energie moyenne faisceau durci 0cm RW3 :{mean_energy_6cm_hardening} MeV")
