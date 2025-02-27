import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math

# from shapely.geometry import Polygon
import pydicom
import sys
from openpyxl import load_workbook
from cmasher import get_sub_cmap
from scipy.interpolate import interp1d
from scipy.signal import find_peaks, peak_widths


def read_data(data_file, correspondence_table_file):
    """Function to read data from measurement .bin file"""
    # .bin measurements file
    dt = np.dtype("uint16")
    f = open(data_file, "rb")
    data = f.read()
    zdata = np.frombuffer(data, dt)

    # correspondence of QDC number and strip number file
    pf = open(correspondence_table_file, "r")
    correspondence_table = pf.readlines()

    # number of measurements
    nb_events = np.size(zdata) // 309

    # time conversion in seconds (integration time = 0.01s + 0.0005 s of dead time)
    time_values = [event * 0.0105 for event in range(nb_events)]

    # strips responses matrix (line = strips, columns = strip responses)
    raw_strip_resp = np.zeros((153, len(time_values)))

    # 17 first strips on the missing diamond => 0 response
    for strip_num in range(18, 153):
        corresponding_QDC_num = int(correspondence_table[strip_num])
        for event in range(nb_events):
            raw_strip_resp[strip_num, event] = np.uint32(
                ((zdata[3 + corresponding_QDC_num * 2 + event * 309]) << 16)
                + (zdata[4 + corresponding_QDC_num * 2 + event * 309])
                >> 6
            )

    return time_values, raw_strip_resp


def mask_hull(tuple_points):
    """Function to compute the convex hull of a set of points"""
    points = []
    for [x, y] in tuple_points:
        middle = math.ceil(len(points) / 2)
        if len(points) >= 2 and points[middle - 1][0] == x and points[middle][0] == x:
            points.pop(middle)
        points.insert(middle, (x, y))

    return points


def get_excel_data(excel_file, ws_name, start_l, start_col, len_data):
    """gets  noize and normalization values from excel file calculated from lateral
    microbeam scan"""
    wb_res = load_workbook(excel_file)
    ws = wb_res[ws_name]
    excel_data = [
        ws.cell(row=start_l + i, column=start_col).value for i in range(len_data)
    ]
    np.concatenate((np.zeros(18), excel_data))
    return np.array(excel_data)


def find_closest_index(list_of_values, target_value):
    return min(
        range(len(list_of_values)), key=lambda x: abs(list_of_values[x] - target_value)
    )


def find_x_closest_peak(x_data, y_data, x_target):
    x_peaks_index, _ = find_peaks(y_data)
    x_peaks = x_data[x_peaks_index]
    closest_x_peak = x_peaks[np.abs(x_peaks - x_target).argmin()]
    return closest_x_peak


def get_x_coord_theo_circle(x_val, v_val):
    """Function to calculate the x coordinate of the center of the theoretical beam
    circle shape. Integrates an X profile of measurements taken at middle of Y axis.
    Calculates abscisse of the middle of the integral = x coordinate of the center of
    theoratical circle"""

    # find all peaks
    peaks, _ = find_peaks(v_val, height=0.5)
    if len(peaks) > 0:
        # get index of central peak
        central_peak_index = peaks[len(peaks) // 2 - 1]
        # retrieves its x position = x pos of center of theo circle
        x_center_circ = x_val[central_peak_index]

    return x_center_circ


def get_and_format_noize_normal_val(
    excel_file, ws_name, start_l, start_col, time_values
):
    noize = get_excel_data(excel_file, ws_name, start_l, start_col, 135)
    normal_val = get_excel_data(excel_file, ws_name, start_l, start_col + 1, 135)
    noize = np.repeat(noize, len(time_values)).reshape(-1, len(time_values))
    normal_val = np.repeat(normal_val, len(time_values)).reshape(-1, len(time_values))

    return noize, normal_val


def get_y_coord_of_center_y_profile(mid_strip, couch_shift, strip_resp):
    """Fonction retrieving the Y profile at a given strip (should be middle strip),
    calc the FWHM of this profile and return the y coord of the center of the FWHM"""

    # get mid strip resp
    mid_strip_resp = strip_resp[mid_strip, :]

    # calc FWHM middle strip
    max_index = mid_strip_resp.argmax()
    fwhm = peak_widths(mid_strip_resp, [max_index], rel_height=0.5)
    fwhm_left = np.interp(fwhm[2][0], np.arange(len(couch_shift)), couch_shift)
    fwhm_right = np.interp(fwhm[3][0], np.arange(len(couch_shift)), couch_shift)

    return (fwhm_right - fwhm_left) / 2 + fwhm_left


def plot_resp(time, resp, title, y_label):
    for strip_num in range(18, 153):
        plt.plot(time, resp[strip_num, :])
    plt.xlabel("Time (s)")
    plt.ylabel(y_label)
    plt.legend(loc="best")
    plt.title(title, fontsize=18, fontstyle="oblique")
    plt.show()


## TO FILL #######
plot_raw_resp = False
plot_strip_resp = False
plot_x_profile = True
fontsize_value = 20
CORRESPONDANCE_FILE = r"C:\Users\milewski\Desktop\these\mesures\analyse_data\codes_python\150_voies\add_piste.txt"
DATA_FILE = r"C:\Users\milewski\Desktop\these\papiers\caracterisation_detecteur_153_voies\fx_formes_complexes\fx_circulaire\mesure\zData_150V_150ubeam_795mu_24p8v0_40collim17mmvitesse10.bin"
DCM_FILE = "RP.1.3.6.1.4.1.33868.20201021131037.169743"
COUCH_SPEED = 10
THRESHOLD = 0.15
NORMAL_VAL_FILE = r"C:\Users\milewski\Desktop\these\papiers\caracterisation_detecteur_153_voies\step_phantom\param_et_resultats.xlsx"

# strip exact positioning
strip_pos_file = r"C:\Users\milewski\OneDrive - Université Grenoble Alpes\these\papiers\caracterisation_detecteur_153_voies\microbeam_scan_analysis\strip_exact_positioning.xlsx"
strip_pos = get_excel_data(strip_pos_file, "strip_pos", 5, 7, 153)


# strip positionning at mask
pitch_strip = 0.2325
pitch_at_mask = 0.2075
dim_fact = pitch_at_mask / pitch_strip
strip_pos_at_mask = strip_pos * dim_fact
dic_strip_pos = {i: strip_pos_at_mask[i] for i in range(153)}
central_strip = 73 - 1  # centered on 9 diamonds because index starts at 0

time_values, raw_strip_resp = read_data(DATA_FILE, CORRESPONDANCE_FILE)
couch_shift = np.array(time_values) * COUCH_SPEED
noize_val = get_excel_data(NORMAL_VAL_FILE, "mb_scan_ESRF", 4, 4, 135)
normal_val = get_excel_data(NORMAL_VAL_FILE, "mb_scan_ESRF", 4, 5, 135)

# raw strips response plot
if plot_raw_resp:
    plot_resp(
        time_values,
        raw_strip_resp,
        "Raw responses circular beam ESRF",
        "QDC response (AU)",
    )

# cut noize and normalize strip resp
noize_val, normal_val = get_and_format_noize_normal_val(
    NORMAL_VAL_FILE, "mb_scan_ESRF", 4, 4, time_values
)
strip_resp = (raw_strip_resp[18:153] - noize_val) / normal_val * 100
strip_resp = np.concatenate((np.zeros((18, len(time_values))), strip_resp))

# strips response plot
if plot_strip_resp:
    plot_resp(
        time_values,
        strip_resp,
        "Normalized strip responses circular beam ESRF",
        "Normalized resp (%)",
    )


# retrieves measurements whose resp > threshold and stores them in "points" array
points = np.empty((0, 3))
for istrip, strip in enumerate(strip_resp):
    for iresp, resp in enumerate(strip):
        if resp >= THRESHOLD:
            couch_mvt = time_values[iresp] * COUCH_SPEED
            # dist_between_strips = istrip * STRIP_PITCH
            points = np.append(
                points, [(dic_strip_pos[istrip], couch_mvt, resp)], axis=0
            )
x, y, v = points[:, 0], points[:, 1], points[:, 2]


# starts Y axis at 0
min_y = min(y)
y = [y[i] - min_y for i in range(len(y))]
points[:, 1] = points[:, 1] - min_y
couch_shift = couch_shift - min_y

# Calc X profile at middle of Y axis
theo_y_mid = (max(points[:, 1]) - min(points[:, 1])) / 2
mid_y = y[find_closest_index(points[:, 1], theo_y_mid)]
x_profile_x = points[np.where(points[:, 1] == mid_y)][:, 0]
x_profile_v = points[np.where(points[:, 1] == mid_y)][:, 2]

# find central peak of x profile = center of theo circle
x_center_circ = get_x_coord_theo_circle(x_profile_x, x_profile_v)

if plot_x_profile:
    plt.plot(x_profile_x, x_profile_v)
    plt.scatter(
        x_center_circ,
        x_profile_v[np.where(x_profile_x == x_center_circ)],
        color="r",
        label="Center of circle",
    )
    plt.xlabel("Distance between strips at mask position (mm)", fontsize=fontsize_value)
    plt.ylabel("Normalized response (%)", fontsize=fontsize_value)
    plt.title("X profile at middle of Y axis", fontsize=fontsize_value)
    plt.tick_params(axis="both", which="major", labelsize=fontsize_value)
    plt.show()


# print(mid_y)
# sys.exit()
x_center_circ = get_x_coord_theo_circle(points)
# x_center_circ = dic_strip_pos[central_strip]

# get y coord of theo circle beam shape = middle of FWHM of central strip but central
# strip not in front of a microbeam so take average of middle of FWHM of neighboring
# strips

y_center_circ_right = get_y_coord_of_center_y_profile(
    central_strip - 1, couch_shift, strip_resp
)
y_center_circ_left = get_y_coord_of_center_y_profile(
    central_strip + 1, couch_shift, strip_resp
)
y_center_circ = (y_center_circ_right + y_center_circ_left) / 2

# create theoratical circle shape of 17 mm centered at calculated points
rad = 17 / 2
theta = np.linspace(0, 2 * np.pi, 100)
x_circ = x_center_circ + rad * np.cos(theta)
y_circ = y_center_circ + rad * np.sin(theta)


# Plot theoratical and measured beam shapes
fig, ax = plt.subplots()
# ax.plot(x_circ, y_circ, "k-", label="17 mm circle beam shape")
y_value = (max(points[:, 1]) - min(points[:, 1])) / 2
ax.axhline(y=mid_y, color="r", linestyle="--", linewidth=2, label=f"y = {y_value}")

cmap = get_sub_cmap("YlGn", 0.2, 0.8)
norm = colors.Normalize(vmin=np.min(v), vmax=np.max(v))
im = ax.scatter(x, y, s=30, marker="s", cmap=cmap, vmin=v.min(), vmax=v.max(), c=v)
cbar = fig.colorbar(im, label="Normalized response (%)")
cbar.ax.tick_params(labelsize=fontsize_value)
cbar.set_label("Normalized response (%)", fontsize=fontsize_value)
ax.set_xlabel("Distance between strips at mask position (mm)", fontsize=fontsize_value)


# To get same scale on both axis
ticks = np.arange(-12, 30, 2)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_ylabel("Couch shift (mm)", fontsize=fontsize_value)
ax.set_ylim(-2, 20)
ax.set_xlim(-12, 12)
plt.tick_params(axis="both", which="major", labelsize=fontsize_value)
plt.show()
