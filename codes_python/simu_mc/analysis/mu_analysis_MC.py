from openpyxl import load_workbook
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from decimal import Decimal
import sys


def get_excel_data(mc_data_file, mc_ws_name, start_l, start_col, len_data):
    """gets normalization values from excel file calculated on step 0 of step phantom
    measurements"""
    wb_res = load_workbook(mc_data_file)
    ws = wb_res[mc_ws_name]
    excel_data = [
        ws.cell(row=start_l + i, column=start_col).value for i in range(len_data)
    ]
    return np.asarray(excel_data)


def exponential_model(x, a, b):
    return a * np.exp(-b * x)


def fill_excel(excel_path, mc_ws_name, data_to_fill, start_l, start_col):
    wb_res = load_workbook(excel_path)
    ws = wb_res[mc_ws_name]
    for i in range(len(data_to_fill)):
        ws.cell(row=start_l + i, column=start_col, value=data_to_fill[i])
    wb_res.save(excel_path)


def generate_normal_matrix(matrix):
    """function normalizing every line of a matrix by the first element of the line"""
    normal_matrix = np.zeros(matrix.shape)
    for i in range(matrix.shape[0]):
        normal_matrix[i] = matrix[i] / matrix[i, 0]
    return normal_matrix


def get_fit_expo_data(x_data, y_data):
    popt, pcov = curve_fit(
        exponential_model,
        x_data,
        y_data,
        p0=(1, 1),
    )
    a_opt, b_opt = popt
    mu = b_opt
    perr = np.sqrt(np.diag(pcov))
    a_err, b_err = perr
    mu_uncertain = b_err

    x_fit = np.linspace(min(step_thickness), max(step_thickness), len(step_thickness))
    y_fit = exponential_model(x_fit, *popt)

    return mu, mu_uncertain, x_fit, y_fit


# DATA TO FILL IN ######################
n_strips = 136
n_steps = 7
step_thickness = [0, 1, 2, 2.5, 3, 4, 5]
mc_data_file = r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\simulations_MC\results_step_phantom.xlsx"
mes_data_file = r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\step_phantom\param_et_resultats.xlsx"
mc_ws_name = "sum_136_strips"
mes_ws_name = "mean_results"
start_l = 8
start_c = 10
start_l_mes = 47
start_c_mes = 6
fill_excel_bool = False
mu_val_by_strip = False
mu_NIST = 0.170456
# ########################################

# get sum edep of strips obtained from MC simu from excel file
df = pd.read_excel(mc_data_file, sheet_name=mc_ws_name)
subset_df = df.iloc[7 : 7 + n_strips, 2 : 2 + n_steps]
subset_df_uncertain = df.iloc[7 : 7 + n_strips, 13 : 13 + n_steps]
strip_sum_edep = subset_df.to_numpy()
strip_sum_edep_uncertain = subset_df_uncertain.to_numpy()

# normalization of sum edep by the first step and propagates uncertainties
normal_sum_edep = generate_normal_matrix(strip_sum_edep)
prop_uncertain = normal_sum_edep * np.sqrt(
    (
        (strip_sum_edep_uncertain / strip_sum_edep) ** 2
        + (strip_sum_edep_uncertain[0] / strip_sum_edep[0]) ** 2
    ).astype(float)
)


# get measurements
mes = get_excel_data(mes_data_file, mes_ws_name, start_l_mes, start_c_mes, 7)
mes_uncertain = get_excel_data(
    mes_data_file, mes_ws_name, start_l_mes, start_c_mes + 1, 7
)

# exponential fit to retrieve linear att coeff by strips
if mu_val_by_strip:
    mu_values = np.zeros(n_strips)
    mu_uncertain = np.zeros(n_strips)
    for strip in range(n_strips):
        popt, pcov = curve_fit(
            exponential_model,
            step_thickness,
            normal_sum_edep[strip],
            p0=(1, 1),
        )
        a_opt, b_opt = popt
        mu_values[strip] = b_opt

        # uncertainties on exponential fit
        perr = np.sqrt(np.diag(pcov))
        a_err, b_err = perr
        mu_uncertain[strip] = b_err

if fill_excel_bool:
    fill_excel(mc_data_file, mc_ws_name, mu_values, start_l, start_c)
    fill_excel(mc_data_file, mc_ws_name, mu_uncertain, start_l, start_c + 1)


in_beam_normal_sum_edep = normal_sum_edep[::2]
in_beam_uncertain = prop_uncertain[::2]
in_beam_mean_mc = np.mean(in_beam_normal_sum_edep, axis=0)
in_beam_mean_uncertain_mc = (1 / len(in_beam_uncertain)) * np.sqrt(
    np.sum(in_beam_uncertain**2)
)

# theoratical curve
x_theo = np.arange(0, 5, 0.05, dtype=float)
y_theo = np.exp(-mu_NIST * x_theo)

# fit exponential model to mean of in beam sum edep
mean_mu_mc, mean_mu_uncertain_mc, x_fit_mc, y_fit_mc = get_fit_expo_data(
    step_thickness, in_beam_mean_mc
)
mu_mes, mu_mes_uncertain, x_fit_mes, y_fit_mes = get_fit_expo_data(step_thickness, mes)

plt.plot(x_theo, y_theo, label="theory (NIST)", linestyle="--", color="black")
plt.scatter(step_thickness, in_beam_mean_mc, label="MC simulations")
plt.errorbar(
    step_thickness,
    in_beam_mean_mc,
    yerr=in_beam_mean_uncertain_mc,
    fmt="o",
    capsize=5,
    color="blue",
)
plt.scatter(step_thickness, mes, marker="^", label="experimental")
plt.errorbar(
    step_thickness, mes, yerr=mes_uncertain, fmt="^", capsize=5, color="orange"
)
plt.plot(
    [],
    [],
    color="white",
    label=f"$\mu$ MC = {mean_mu_mc:.3f} $\pm$ {Decimal(mean_mu_uncertain_mc):.2E} cm\u207b\u00b9",
)
plt.plot(
    [],
    [],
    color="white",
    label=f"$\mu$ exp = {mu_mes:.3f} $\pm$ {Decimal(mu_mes_uncertain):.2E} cm\u207b\u00b9",
)
plt.plot(
    [],
    [],
    color="white",
    label=f"$\mu$ NIST = {mu_NIST:.3f} cm\u207b\u00b9",
)
# plt.title(
#     "Comparison between attenuation coefficient $\mu$ measured at ANSTO and simulated by MC",
#     fontsize=15,
# )
plt.xlabel("SolidHE water thickness (cm)", fontsize=15)
plt.ylabel("Response / Response step 0", fontsize=15)
plt.legend()
plt.show()
