from openpyxl import load_workbook
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
from decimal import Decimal
import sys


def get_excel_data(mc_data_file, mc_ws_name, start_l, start_col, n_l, n_c):
    """gets normalization values from excel file calculated on step 0 of step phantom
    measurements"""
    wb_res = load_workbook(mc_data_file)
    ws = wb_res[mc_ws_name]
    excel_data = []
    for line in range(n_l):
        for col in range(n_c):
            excel_data.append(ws.cell(row=start_l + line, column=start_col + col).value)
    return np.asarray(excel_data)


def exponential_model(x, a, b):
    return a * np.exp(-b * x)


def calc_mu(step_thickness, edep):
    popt, pcov = curve_fit(
        exponential_model,
        step_thickness,
        edep,
        p0=(1, 1),
    )
    a_opt, b_opt = popt
    mu = b_opt

    # uncertainties on exponential fit
    perr = np.sqrt(np.diag(pcov))
    a_err, b_err = perr
    mu_uncertain = b_err

    return mu, mu_uncertain


def fill_excel(excel_path, mc_ws_name, data_to_fill, start_l, start_col):
    wb_res = load_workbook(excel_path)
    ws = wb_res[mc_ws_name]
    for i in range(len(data_to_fill)):
        ws.cell(row=start_l + i, column=start_col, value=data_to_fill[i])
    wb_res.save(excel_path)


def generate_normal_matrix(matrix, matrix_uncertain):
    """function normalizing every line of a matrix by the first element of the line and
    propagates uncertainties"""

    if matrix.shape != matrix_uncertain.shape:
        raise ValueError("matrix and matrix_uncertain must have the same shape")
    if matrix.ndim != 2:
        raise ValueError("matrix and matrix_uncertain must be 2D")

    normal_matrix = matrix / matrix[:, :1]
    normal_matrix_uncertain = normal_matrix * np.sqrt(
        (matrix_uncertain[:, :1] / matrix[:, :1]) ** 2
        + (matrix_uncertain / matrix) ** 2
    )
    return normal_matrix, normal_matrix_uncertain


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
step_thickness = np.array([0, 1, 2, 3, 4, 6])
step_phant_mat = "RW3"
mc_data_file = r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\simulations_MC\results_step_phantom.xlsx"
mes_data_file = r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\step_phantom\param_et_resultats.xlsx"
mc_ws_name = "ESRF_force_coll_test"
mes_ws_name = "mean_results"
start_l_mc = 8
start_c_mc = 2
start_l_mes = 8
start_c_mes = 6
mu_NIST = 0.1616615
# ########################################

# get edep of strips obtained from MC simu from excel file
df = pd.read_excel(mc_data_file, sheet_name=mc_ws_name)
subset_df = df.iloc[
    start_l_mc : start_l_mc + n_strips, start_c_mc : start_c_mc + len(step_thickness)
]
subset_df_uncertain = df.iloc[
    start_l_mc : start_l_mc + n_strips,
    start_c_mc + 7 : start_c_mc + 7 + len(step_thickness) * 2,
]
all_edep_mc = subset_df.to_numpy()
all_uncertain_mc = subset_df_uncertain.to_numpy()

# retrieves only statistical uncertainties
stat_uncertain_mc = all_uncertain_mc[:, 0::2]


# retrieves edep only for the strips facing microbeams and make sure values are float
edep_mc = np.array(all_edep_mc[1::2], dtype=float)
uncertain_mc = np.array(stat_uncertain_mc[1::2], dtype=float)

# normalizes edep by the value at step 0 and propagates uncertainties
normal_edep_mc, normal_uncertain_mc = generate_normal_matrix(edep_mc, uncertain_mc)

# calc mean value and propagates uncertainties
mean_edep_mc = np.mean(normal_edep_mc, axis=0)
mean_uncertain_mc = (1 / len(normal_uncertain_mc)) * np.sqrt(
    np.array(np.sum(normal_uncertain_mc**2, axis=0), dtype=np.float64)
)

# get measurements from excel file
edep_mes = get_excel_data(
    mes_data_file, mes_ws_name, start_l_mes, start_c_mes, len(step_thickness), 1
)
uncertain_mes = get_excel_data(
    mes_data_file, mes_ws_name, start_l_mes, start_c_mes + 1, len(step_thickness), 1
)

# exponential fit to retrieve linear att coeff
mu_mc, uncertain_mu_mc = calc_mu(step_thickness, mean_edep_mc)
mu_mes, uncertain_mu_mes = calc_mu(step_thickness, edep_mes)

# print(f"mu_mc {mu_mc:.3f} +/- {Decimal(uncertain_mu_mc):.2E}")
# sys.exit()

# theoratical curve
x_theo = np.arange(0, 5, 0.05, dtype=float)
y_theo = np.exp(-mu_NIST * x_theo)

# # fit exponential model to mean of in beam sum edep
# mean_mu_mc, mean_mu_uncertain_mc, x_fit_mc, y_fit_mc = get_fit_expo_data(
#     step_thickness, in_beam_mean_mc
# )
# mu_mes, mu_mes_uncertain, x_fit_mes, y_fit_mes = get_fit_expo_data(step_thickness, mes)

# plt.plot(x_theo, y_theo, label="theory (NIST)", linestyle="--", color="black")
plt.scatter(step_thickness, mean_edep_mc, label="MC simulations")
# plt.errorbar(
#     step_thickness,
#     mean_edep_mc,
#     yerr=mean_uncertain_mc,
#     fmt="o",
#     capsize=5,
#     color="blue",
# )
plt.scatter(step_thickness, edep_mes, marker="^", label="experimental")
plt.errorbar(
    step_thickness, edep_mes, yerr=uncertain_mes, fmt="^", capsize=5, color="orange"
)
plt.plot(
    [],
    [],
    color="white",
    label=f"$\mu$ MC = {mu_mc:.3f} $\pm$ {Decimal(uncertain_mu_mc):.2E} cm\u207b\u00b9",
)
plt.plot(
    [],
    [],
    color="white",
    label=f"$\mu$ exp = {mu_mes:.3f} $\pm$ {Decimal(uncertain_mu_mes):.2E} cm\u207b\u00b9",
)
# plt.plot(
#     [],
#     [],
#     color="white",
#     label=f"$\mu$ NIST = {mu_NIST:.3f} cm\u207b\u00b9",
# )
# plt.title(
#     "Comparison between attenuation coefficient $\mu$ measured at ANSTO and simulated by MC",
#     fontsize=15,
# )
x_label = step_phant_mat + " thickness (cm)"
plt.xlabel(x_label, fontsize=15)
plt.ylabel("Response / Response step 0", fontsize=15)
plt.legend()
plt.show()
