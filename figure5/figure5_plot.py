################################################################################
# FILE: figure5_plot.py
#
# DESCRIPTION:
#   This script generates Figure 5, illustrating the simultaneous effects of
#   three stressors (pmax, temp, ma) on the temporal lag between seagrass
#   biomass (S) and stable carbon (CB) reaching new equilibria. The first
#   column shows pairwise stressor interactions (no third stressor), and the
#   remaining columns present scenarios in which all three stressors vary.
#
#   Three discrete rates of change are highlighted (slow, medium, fast). A
#   pastel rainbow colormap denotes the lag magnitude, with red representing
#   longer lags and blue shorter lags. This script loads data from:
#     - df1_pattern_lag.csv, df2_pattern_lag.csv, df3_pattern_lag.csv
#       (pairwise combinations)
#     - lag_rate_3D.csv (full 3-stressor simulations)
#
# DEPENDENCIES:
#   Requires the following Python libraries:
#     - os
#     - pandas
#     - numpy
#     - matplotlib
################################################################################

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation

# ------------------------------------------------------------------------------
# 1) SET WORKING DIRECTORY
# ------------------------------------------------------------------------------
os.chdir("")

# ------------------------------------------------------------------------------
# 2) LOAD PAIRWISE STRESSOR DATA FOR FIRST COLUMN
# ------------------------------------------------------------------------------
data_pmax_temp = pd.read_csv("df1_pattern_lag.csv")
data_pmax_ma   = pd.read_csv("df2_pattern_lag.csv")
data_temp_ma   = pd.read_csv("df3_pattern_lag.csv")

# ------------------------------------------------------------------------------
# 3) LOAD FULL 3D INTERACTION DATA FOR THE NEXT THREE COLUMNS
# ------------------------------------------------------------------------------
data = pd.read_csv("lag_rate_3D.csv")

# Convert 'temp', 'ma', and 'pmax' rates from days to years
data['temp'] /= 365
data['ma']   /= 365
data['pmax'] /= 365

# ------------------------------------------------------------------------------
# 4) COMBINE LAG VALUES FOR A SINGLE COLOR SCALE
# ------------------------------------------------------------------------------
all_lags = pd.concat([
    data_pmax_ma['lag'],
    data_pmax_temp['lag'],
    data_temp_ma['lag'],
    data['lag']
], ignore_index=True)
lag_min, lag_max = all_lags.min(), all_lags.max()

# ------------------------------------------------------------------------------
# 5) DEFINE A PASTEL RAINBOW COLORMAP FOR VISUALIZING LAG
# ------------------------------------------------------------------------------
pastel_lag_colors = plt.cm.rainbow(np.linspace(0, 1, 256)) * 0.7 + 0.3
pastel_cmap       = plt.cm.colors.ListedColormap(pastel_lag_colors)

# ------------------------------------------------------------------------------
# 6) INITIALIZE THE 4×3 FIGURE (4 columns, 3 rows)
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(nrows=3, ncols=4, figsize=(20, 14))

# ------------------------------------------------------------------------------
# HELPER FUNCTION: 2D LAG PLOT WITH OPTIONAL CONTOUR
# ------------------------------------------------------------------------------
def plot_2d_lag(ax, xdata, ydata, zdata, xlab, ylab, title_str):
    """
    Plots the temporal lag (zdata) vs. two stressor axes (xdata, ydata).
    If insufficient data points exist for a triangulated contour, it falls
    back to a scatter plot.
    """
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)

    if len(xdata) < 3:
        # Fall back to scatter for very few points
        ax.scatter(
            xdata, ydata, c=zdata,
            cmap=pastel_cmap, vmin=lag_min, vmax=lag_max
        )
        ax.set_title(f"{title_str} (scatter)")
    else:
        # Use triangulation-based contour
        tri = Triangulation(xdata, ydata)
        ax.tricontourf(
            tri, zdata, levels=20,
            cmap=pastel_cmap, vmin=lag_min, vmax=lag_max
        )
        ax.set_title(title_str)

# ------------------------------------------------------------------------------
# 7) FIRST COLUMN: PAIRWISE STRESSORS (pmax–ma, pmax–temp, temp–ma)
# ------------------------------------------------------------------------------
axes[0, 0].scatter(
    data_pmax_ma["pmax"], data_pmax_ma["ma"],
    c=data_pmax_ma["lag"], cmap=pastel_cmap,
    s=200, marker='s', vmin=lag_min, vmax=lag_max
)
axes[0, 0].set_xlabel("pmax")
axes[0, 0].set_ylabel("ma")
axes[0, 0].set_title("pmax vs ma (Lag)")

axes[1, 0].scatter(
    data_pmax_temp["pmax"], data_pmax_temp["temp"],
    c=data_pmax_temp["lag"], cmap=pastel_cmap,
    s=200, marker='s', vmin=lag_min, vmax=lag_max
)
axes[1, 0].set_xlabel("pmax")
axes[1, 0].set_ylabel("temp")
axes[1, 0].set_title("pmax vs temp (Lag)")

axes[2, 0].scatter(
    data_temp_ma["temp"], data_temp_ma["ma"],
    c=data_temp_ma["lag"], cmap=pastel_cmap,
    s=200, marker='s', vmin=lag_min, vmax=lag_max
)
axes[2, 0].set_xlabel("temp")
axes[2, 0].set_ylabel("ma")
axes[2, 0].set_title("temp vs ma (Lag)")

# ------------------------------------------------------------------------------
# 8) NEXT THREE COLUMNS: FULL 3D INTERACTION (lag_rate_3D.csv)
# ------------------------------------------------------------------------------
# Discrete rate values (in years) for the third stressor
val_high   = 1000
val_medium = 556.0
val_low    = 1
param_vals = [val_high, val_medium, val_low]

# Row 1: Fix temp → (x=pmax, y=ma)
for col, val in enumerate(param_vals):
    subset = data[np.isclose(data['temp'], val, atol=1.0)]
    x = subset['pmax'].values
    y = subset['ma'].values
    z = subset['lag'].values
    ax = axes[0, col + 1]
    plot_2d_lag(ax, x, y, z, "pmax", "ma", f"Fix temp ≈ {val:.1f}")

# Row 2: Fix ma → (x=pmax, y=temp)
for col, val in enumerate(param_vals):
    subset = data[np.isclose(data['ma'], val, atol=1.0)]
    x = subset['pmax'].values
    y = subset['temp'].values
    z = subset['lag'].values
    ax = axes[1, col + 1]
    plot_2d_lag(ax, x, y, z, "pmax", "temp", f"Fix ma ≈ {val:.1f}")

# Row 3: Fix pmax → (x=temp, y=ma)
for col, val in enumerate(param_vals):
    subset = data[np.isclose(data['pmax'], val, atol=1.0)]
    x = subset['temp'].values
    y = subset['ma'].values
    z = subset['lag'].values
    ax = axes[2, col + 1]
    plot_2d_lag(ax, x, y, z, "temp", "ma", f"Fix pmax ≈ {val:.1f}")

# ------------------------------------------------------------------------------
# 9) SAVE MAIN FIGURE (NO LEGEND)
# ------------------------------------------------------------------------------
plt.tight_layout()
plt.savefig("figure5.png", dpi=300, bbox_inches='tight')
plt.show()

# ------------------------------------------------------------------------------
# 10) CREATE AND SAVE A SEPARATE FIGURE FOR THE VERTICAL COLORBAR
# ------------------------------------------------------------------------------
fig_legend = plt.figure(figsize=(2, 6))  # Narrow figure, taller for vertical bar
norm = plt.Normalize(lag_min, lag_max)
scalar_mappable = plt.cm.ScalarMappable(norm=norm, cmap=pastel_cmap)
scalar_mappable.set_array([])

cbar = plt.colorbar(
    mappable=scalar_mappable,
    orientation='vertical',
    fraction=0.3,
    pad=0.1
)
# Position label/ticks on the right side
cbar.ax.yaxis.set_label_position("right")
cbar.ax.yaxis.tick_right()
cbar.set_label("Lag", labelpad=15)

plt.savefig("figure5_legend_bar.png", dpi=300, bbox_inches='tight')
plt.show()
