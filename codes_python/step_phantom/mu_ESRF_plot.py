import pandas as pd
import matplotlib.pyplot as plt
from decimal import Decimal

# Data to fill in
file_path = r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\step_phantom\param_et_resultats.xlsx"
ws = "mu_ESRF"
df = pd.read_excel(file_path, sheet_name=ws, skiprows=2)
mu_theo = 0.1617
mu_theo_harden_spec = 0.1606


strip_num = df["Strip num"].values
mu = df["Mu (cm-1)"].values

# Sort even and odd strips
even_strip_num = []
mu_even_strip = []
odd_strip_num = []
mu_odd_strip = []

for strip_num, mu in zip(strip_num, mu):
    if strip_num % 2 == 0:
        even_strip_num.append(strip_num)
        mu_even_strip.append(mu)
    else:
        odd_strip_num.append(strip_num)
        mu_odd_strip.append(mu)


plt.scatter(even_strip_num, mu_even_strip, label="even strips", color="b", marker="o")
plt.scatter(odd_strip_num, mu_odd_strip, label="odd strips", color="y", marker="o")
plt.axhline(
    y=mu_theo,
    color="r",
    linestyle="--",
    label=f"$\mu$(RW3) = {Decimal(mu_theo):.3E} cm\u207b\u00b9 for 121.1 keV",
)
plt.axhline(
    y=mu_theo_harden_spec,
    color="k",
    linestyle="--",
    label=f"$\mu$(RW3) = {Decimal(mu_theo_harden_spec):.3E} cm\u207b\u00b9 for 123.5 keV",
)

plt.xlabel("Strip number")
plt.ylabel("Linear attenuation coefficient (cm-1)")
plt.legend()
plt.show()
