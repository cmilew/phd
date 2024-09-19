import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math
from shapely.geometry import Polygon
import pydicom
import sys
from openpyxl import load_workbook
from scipy.interpolate import interp1d
from cmasher import get_sub_cmap
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


def get_beam_contour_from_dcm(plan_name):
    """Function to read dcm file containing beam shape (from TPS) and returning beam
    contours"""
    plan = pydicom.read_file(plan_name, force=True)
    dummy_X, dummy_Y, centroid = [], [], []
    beam_seq = plan[(0x300A, 0x00B0)].value
    beam_contours = []
    for i in range(len(beam_seq)):
        block_seq = beam_seq[i][(0x300A, 0x00F4)].value
        block_data = block_seq[0][(0x300A, 0x0106)].value
        dummy = []
        for entry in block_data:
            dummy.append(entry)
        dummy_array = np.array(dummy).astype(float)
        dummy_X.append(dummy_array[::2])
        dummy_X[i] = dummy_X[i].reshape(dummy_X[i].shape[0], 1)
        dummy_Y.append(dummy_array[1::2])
        dummy_Y[i] = dummy_Y[i].reshape(dummy_Y[i].shape[0], 1)
        dummy_X = dummy_X
        dummy_Y = dummy_Y
        centroid.append(
            (
                sum(dummy_X[i][:, 0]) / len(dummy_X[i][:, 0]),
                sum(dummy_Y[i][:, 0]) / len(dummy_Y[i][:, 0]),
            )
        )
        beam_contours.append(np.concatenate((dummy_X[i], dummy_Y[i]), axis=1))
    return beam_contours


def get_excel_data(excel_file, ws_name, start_l, start_col):
    """gets normalization values from excel file calculated on step 0 of step phantom
    measurements"""
    wb_res = load_workbook(excel_file)
    ws = wb_res[ws_name]
    excel_data = [ws.cell(row=start_l + i, column=start_col).value for i in range(135)]
    np.concatenate((np.zeros(18), excel_data))
    return np.array(excel_data)


def find_closest_index(list_of_values, target_value):
    return min(
        range(len(list_of_values)), key=lambda x: abs(list_of_values[x] - target_value)
    )


def get_and_format_noize_normal_val(
    excel_file, ws_name, start_l, start_col, time_values
):
    noize = get_excel_data(excel_file, ws_name, start_l, start_col)
    normal_val = get_excel_data(excel_file, ws_name, start_l, start_col + 1)
    noize = np.repeat(noize, len(time_values)).reshape(-1, len(time_values))
    noize_strips = np.concatenate((np.zeros((18, len(time_values))), noize))
    normal_val = np.repeat(normal_val, len(time_values)).reshape(-1, len(time_values))
    normal_val_strips = np.concatenate((np.zeros((18, len(time_values))), normal_val))

    return noize_strips, normal_val_strips


def plot_resp(time, resp, title, y_label):
    for strip_num in range(18, 153):
        plt.plot(time, resp[strip_num, :], label=f"Strip {strip_num}")
    plt.xlabel("Time (s)")
    plt.ylabel(y_label)
    plt.legend(loc="best")
    plt.title(title, fontsize=18, fontstyle="oblique")
    plt.show()


def get_central_x_val_fwhm(x, y):
    """Function to return the central value of the FWHM of a profile"""

    max_index = y.argmax()
    fwhm = peak_widths(y, [max_index], rel_height=0.5)
    fwhm_left = np.interp(fwhm[2][0], np.arange(len(x)), x)
    fwhm_right = np.interp(fwhm[3][0], np.arange(len(x)), x)

    return (fwhm_right - fwhm_left) / 2 + fwhm_left


# TO FILL #############
plot_raw_resp = False
plot_strip_resp = False
fontsize_value = 10
# strip measuring the microbeam on the furthest right of the mask
extreme_right_strip = 146

CORRESPONDANCE_FILE = r"C:\Users\milewski\Desktop\these\mesures\analyse_data\codes_python\150_voies\add_piste.txt"
DATA_FILE = r"C:\Users\milewski\Desktop\these\papiers\caracterisation_detecteur_153_voies\fx_formes_complexes\fx_cochons\mesures\zData_150V_150ubeam_795mu_24p8v0_40_collimCochon_vitesse10.bin"
DCM_FILE = r"C:\Users\milewski\Desktop\these\phd\codes_python\fx_formes_complexes\RP.1.3.6.1.4.1.33868.20210210162030.730598"
NORMAL_VAL_FILE = r"C:\Users\milewski\Desktop\these\papiers\caracterisation_detecteur_153_voies\step_phantom\param_et_resultats.xlsx"

# strip pitch 0.2075 mm calculated at mask position
STRIP_PITCH = 0.2075
COUCH_SPEED = 10
THRESHOLD = 0.15


time_values, raw_strip_resp = read_data(DATA_FILE, CORRESPONDANCE_FILE)
couch_shift = np.array(time_values) * COUCH_SPEED
noize_val = get_excel_data(NORMAL_VAL_FILE, "mb_scan_ESRF", 4, 4)
normal_val = get_excel_data(NORMAL_VAL_FILE, "mb_scan_ESRF", 4, 5)
strip_width = 0.1725  # mm
inter_strip_width = 0.06  # mm


# raw strips response plot
if plot_raw_resp:
    plot_resp(
        time_values,
        raw_strip_resp,
        "Raw responses pig beam ESRF",
        "QDC response (AU)",
    )


# cut noize and normalize strip resp
noize_val, normal_val = get_and_format_noize_normal_val(
    NORMAL_VAL_FILE, "mb_scan_ESRF", 4, 4, time_values
)
strip_resp = (raw_strip_resp - noize_val) / normal_val * 100

# strips response plot
if plot_strip_resp:
    plot_resp(
        time_values,
        strip_resp,
        "Normalized strip responses pig beam ESRF",
        "Normalized resp (%)",
    )


# Compute points
points = np.empty((0, 3))
for istrip, strip in enumerate(strip_resp):
    for iresp, resp in enumerate(strip):
        if resp >= THRESHOLD:
            dist_between_strips = istrip * STRIP_PITCH
            points = np.append(
                points, [(dist_between_strips, couch_shift[iresp], resp)], axis=0
            )

x, y, v = points[:, 0], points[:, 1], points[:, 2]

# to start Y axis at 0
min_y = min(y)
y = [y[i] - min_y for i in range(len(y))]
points[:, 1] = points[:, 1] - min_y
couch_shift = couch_shift - min_y

# beam contours from dcm
beam_cont_dcm = get_beam_contour_from_dcm(DCM_FILE)
beam_contours = beam_cont_dcm[1][:-1]

# add first point of beam contours at the end to complete the shape
beam_contours = np.vstack([beam_contours, beam_contours[0]])

# rotate shape to retrieve same orientation as strip measurements
beam_contours[:, 0] = -beam_contours[:, 0]

# x position of strip 146 (most extreme right strip reading microbeam)
x_extreme_r_strip = extreme_right_strip * STRIP_PITCH

# right side x position of mask shape
x_r_side_mask = extreme_right_strip * STRIP_PITCH + strip_width / 2 + inter_strip_width

# translate mask shape to match measurements in x direction
x_shift = abs(np.max(beam_contours[:, 0]) - x_r_side_mask)
x_dcm = [x + x_shift for x in x_dcm]
beam_contours[:, 0] = x_dcm

# Strip furthest on the right of all the strips facing a microbeam
max_strip_resp = np.max(strip_resp, axis=1)
r_mb_strip = (np.where(max_strip_resp > 90)[0])[-1]

# get central value of FWHM of the profile of the furthest right strip = y center
# position of the .dcm beam shape
y_center_r_mb_strip = get_central_x_val_fwhm(couch_shift, strip_resp[r_mb_strip, :])

# get max x coord of beam contours
max_x_dcm = np.max(beam_contours[:, 0])

# get y dcm coord of furthest right side of TPS beam shape
y_r_side_coord = beam_contours[beam_contours[:, 0] == max_x_dcm][:, 1]

# compare center of right side TPS beam shame and center right mb strip = y shift
y_r_side_center = (y_r_side_coord[1] - y_r_side_coord[0]) / 2 + y_r_side_coord[0]
print(f"y_r_side_coord = {y_r_side_coord}")
print(f"y_r_side_center = {y_r_side_center}")
print(f"y_center_r_mb_strip = {y_center_r_mb_strip}")
y_shift = abs(y_r_side_center - y_center_r_mb_strip)

y_dcm = [y + y_shift for y in y_dcm]

# x coord swhere to cut TPS beam shape
x_cut = 18 * STRIP_PITCH - strip_width / 2
interpolation = interp1d(x_dcm, y_dcm, kind="linear")
y_cut = interpolation(x_cut)
x_dcm = [x for x in x_dcm if x > x_cut]
y_dcm = [y for y in y_dcm if y > y_cut]

print(f"beam_contours = {beam_contours}")
print(f"x_dcm = {x_dcm}")
print(f"y_dcm = {y_dcm}")

sys.exit()
# Plot beam shape
fontsize_value = 15
fig, ax = plt.subplots()
ax.plot(x_dcm, y_dcm, "k-", label=theo_shape_label)

# ax.axhline(y=mid_y, color="k", linestyle="--", label="Middle of Y axis")
cmap = get_sub_cmap("YlGn", 0.2, 0.8)
norm = colors.Normalize(vmin=np.min(v), vmax=np.max(v))

im = ax.scatter(x, y, s=5, cmap=cmap, vmin=v.min(), vmax=v.max(), c=v)
cbar = fig.colorbar(im, label="Normalized response (%)")
cbar.ax.tick_params(labelsize=fontsize_value)
cbar.set_label("Normalized response (%)", fontsize=fontsize_value)
ax.set_xlabel("Distance between strips at mask position (mm)", fontsize=fontsize_value)
ax.legend(fontsize=fontsize_value)

# To get same scale on both axis
ticks = np.arange(0, 38, 2)
ax.set_xticks(ticks)
ax.set_yticks(ticks)
ax.set_aspect("equal", adjustable="box")
ax.set_ylabel("Couch shift (mm)", fontsize=fontsize_value)
ax.set_ylim(-2, 24)
ax.set_xlim(0, 34)
plt.show()
