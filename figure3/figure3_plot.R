################################################################################
# FILE: figure3_plot.R
#
# DESCRIPTION:
#   This script plots the rate-dependent lag between seagrass biomass (S) and stable
#   carbon (CB) dynamics, as computed by "figure3_lag.R" (stored in "lag_rate.csv").
#   Three subplots correspond to different stressors (ma, pmax, temp). Each plot
#   shows the lag (years) vs. the time before high stress (years). The figure is
#   saved as "figure3.png".
#
# DEPENDENCIES:
#   Requires the following R packages:
#     - ggplot2
#     - dplyr
#     - gridExtra
################################################################################

# ------------------------------------------------------------------------------
# Load Required Packages
# ------------------------------------------------------------------------------
package_list <- c("ggplot2", "dplyr", "gridExtra")

for (pkg in package_list) {
  library(pkg, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# Set Working Directory (optional; adjust as needed)
# ------------------------------------------------------------------------------
setwd("")  # set your own working directory

# ------------------------------------------------------------------------------
# Load and Prepare Data
# ------------------------------------------------------------------------------
data <- read.csv("lag_rate.csv", header = TRUE)

# Convert rate from days to years
data$rate <- data$rate / 365

# Subset data for each stressor
data_ma <- data %>%
  filter(scenario == "ma") %>%
  mutate(lag = round(lag, 0))

data_pmax <- data %>%
  filter(scenario == "pmax")

data_temp <- data %>%
  filter(scenario == "temp")

# ------------------------------------------------------------------------------
# Build Individual Plots
# ------------------------------------------------------------------------------
plot_ma <- ggplot(data_ma, aes(x = rate, y = lag)) +
  geom_line() +
  ylab("Lag (years)") +
  xlab("Time before high stress (years)") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_reverse() +
  scale_y_continuous(limits = c(300, 490))

plot_pmax <- ggplot(data_pmax, aes(x = rate, y = lag)) +
  geom_line() +
  ylab("Lag (years)") +
  xlab("Time before high stress (years)") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_reverse()

plot_temp <- ggplot(data_temp, aes(x = rate, y = lag)) +
  geom_line() +
  ylab("Lag (years)") +
  xlab("Time before high stress (years)") +
  theme_classic() +
  theme(legend.position = "none") +
  scale_x_reverse()

# ------------------------------------------------------------------------------
# Arrange Plots into a Single Figure
# ------------------------------------------------------------------------------
combined_plot <- grid.arrange(
  plot_ma,
  plot_pmax,
  plot_temp,
  ncol = 3
)

# ------------------------------------------------------------------------------
# Save Figure
# ------------------------------------------------------------------------------
ggsave(
  "figure3.png",
  combined_plot,
  width = 13,
  height = 3,
  dpi = 600
)
