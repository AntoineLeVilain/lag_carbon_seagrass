################################################################################
# FILE: figure4_plot.py
#
# DESCRIPTION:
#   This script creates a two-row, three-column set of panels illustrating the
#   combined pairwise effects of stressor type and rate of change on transient
#   patterns and temporal lag duration in a seagrass–soil model. It loads three
#   CSV files (df1_pattern_lag.csv, df2_pattern_lag.csv, df3_pattern_lag.csv)
#   that contain pre-computed data for each pairwise stressor scenario.
#
#   The top row of panels visualizes the "lag" in response to each pairwise
#   interaction, using a color gradient. The bottom row depicts pattern types
#   (e.g., monotonic increase/decrease of seagrass biomass vs. carbon). A shared
#   legend appears in the bottom-right panel.
#
# DEPENDENCIES:
#   - pandas
#   - matplotlib
#   - numpy
################################################################################

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# ------------------------------------------------------------------------------
# Set working directory
# ------------------------------------------------------------------------------
os.chdir("")

# ------------------------------------------------------------------------------
# Load precomputed datasets for each pair of stressors
# ------------------------------------------------------------------------------
data_pmax_temp = pd.read_csv("df1_pattern_lag.csv")
data_pmax_ma   = pd.read_csv("df2_pattern_lag.csv")
data_temp_ma   = pd.read_csv("df3_pattern_lag.csv")

# ------------------------------------------------------------------------------
# Set up the figure layout and custom color maps
# ------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(18, 12))

# Create a pastel-like colormap for visualizing lag
pastel_lag_colors = plt.cm.rainbow(np.linspace(0, 1, 256)) * 0.7 + 0.3
pastel_cmap = plt.cm.colors.ListedColormap(pastel_lag_colors)

# Identify all unique patterns from the three datasets
unique_patterns = sorted(
    set(data_pmax_ma['pattern']) 
    | set(data_pmax_temp['pattern']) 
    | set(data_temp_ma['pattern'])
)

# Assign colors for each unique pattern
pattern_colors = plt.cm.tab10(np.linspace(0, 1, len(unique_patterns)))
pattern_color_dict = {
    pat: pattern_colors[i] for i, pat in enumerate(unique_patterns)
}

# ------------------------------------------------------------------------------
# TOP ROW: Plot lag for each stressor combination
# ------------------------------------------------------------------------------
# 1) pmax vs. ma (Lag)
scatter_1 = axes[0, 0].scatter(
    data_pmax_ma['pmax'],
    data_pmax_ma['ma'],
    c=data_pmax_ma['lag'],
    cmap=pastel_cmap,
    s=300,
    marker='s'
)
axes[0, 0].set_xlabel("pmax")
axes[0, 0].set_ylabel("ma")
axes[0, 0].set_title("pmax vs ma (Lag)")
cbar1 = fig.colorbar(scatter_1, ax=axes[0, 0], orientation='horizontal',
                     pad=0.1)
cbar1.set_label("Lag")

# 2) pmax vs. temp (Lag)
scatter_2 = axes[0, 1].scatter(
    data_pmax_temp['pmax'],
    data_pmax_temp['temp'],
    c=data_pmax_temp['lag'],
    cmap=pastel_cmap,
    s=300,
    marker='s'
)
axes[0, 1].set_xlabel("pmax")
axes[0, 1].set_ylabel("temp")
axes[0, 1].set_title("pmax vs temp (Lag)")
cbar2 = fig.colorbar(scatter_2, ax=axes[0, 1], orientation='horizontal',
                     pad=0.1)
cbar2.set_label("Lag")

# 3) temp vs. ma (Lag)
scatter_3 = axes[0, 2].scatter(
    data_temp_ma['temp'],
    data_temp_ma['ma'],
    c=data_temp_ma['lag'],
    cmap=pastel_cmap,
    s=300,
    marker='s'
)
axes[0, 2].set_xlabel("temp")
axes[0, 2].set_ylabel("ma")
axes[0, 2].set_title("temp vs ma (Lag)")
cbar3 = fig.colorbar(scatter_3, ax=axes[0, 2], orientation='horizontal',
                     pad=0.1)
cbar3.set_label("Lag")

# ------------------------------------------------------------------------------
# BOTTOM ROW: Plot pattern type for each stressor combination
# ------------------------------------------------------------------------------
# Prepare iteration over (axis, dataset, x-col, y-col, title)
plot_specs = zip(
    axes[1],
    [data_pmax_ma, data_pmax_temp, data_temp_ma],
    ["pmax",        "pmax",        "temp"],
    ["ma",          "temp",        "ma"],
    [
        "pmax vs ma (Pattern)",
        "pmax vs temp (Pattern)",
        "temp vs ma (Pattern)"
    ]
)

# For each panel in the bottom row, plot by pattern
for ax, df, xcol, ycol, panel_title in plot_specs:
    for pat in unique_patterns:
        subset = df[df['pattern'] == pat]
        ax.scatter(
            subset[xcol],
            subset[ycol],
            color=pattern_color_dict[pat],
            label=pat,
            s=300,
            marker='s'
        )
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    ax.set_title(panel_title)

# ------------------------------------------------------------------------------
# Add a shared legend for the pattern types
# ------------------------------------------------------------------------------
handles = [
    plt.Line2D(
        [0], [0],
        marker='s',
        color=color,
        linestyle='None',
        label=pat,
        markersize=10
    )
    for pat, color in pattern_color_dict.items()
]

# Place the legend in the bottom-right panel, anchored outside the axes
axes[1, 2].legend(
    handles=handles,
    title="Pattern",
    bbox_to_anchor=(1.05, 1),
    loc='upper left'
)

# ------------------------------------------------------------------------------
# Final layout adjustments and figure export
# ------------------------------------------------------------------------------
plt.tight_layout()
plt.savefig("figure4.png", dpi=300)
plt.show()
