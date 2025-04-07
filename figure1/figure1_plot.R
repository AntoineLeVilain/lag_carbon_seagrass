################################################################################
# FILE: figure1_plot.R
#
# DESCRIPTION:
#   This script generates a set of plots (Figure 1) using our seagrass–soil 
#   model, including bifurcation diagrams and transient dynamics for 
#   different rates of environmental change. It saves the final figure as 
#   "figure1.png" and also creates a separate legend file 
#   "figure1_legend_bar.png".
#
# DEPENDENCIES:
#   Requires the following R packages:
#     - ggplot2
#     - readxl
#     - dplyr
#     - scales
#     - deSolve
#     - nleqslv
#     - rootSolve
#     - gridExtra
#     - ggnewscale
#     - grid
#     - cowplot
################################################################################

# ------------------------------------------------------------------------------
# Set working directory (adjust as needed)
# ------------------------------------------------------------------------------
setwd("")  # set your own working directory

# ------------------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------------------
package_list <- c(
  "ggplot2", "readxl", "dplyr", "scales", 
  "deSolve", "nleqslv", "rootSolve", "gridExtra", 
  "ggnewscale", "grid", "cowplot"
)

for (package in package_list) {
  library(package, character.only = TRUE)
}

# ------------------------------------------------------------------------------
# multisolve2()
#
# DESCRIPTION:
#   Finds equilibrium solutions for a given model. It uses `multiroot` to
#   search for roots from multiple initial guesses (between lower_limit and
#   upper_limit), and compiles any unique solutions found.
#
# ARGUMENTS:
#   model        - function; accepts a numeric vector and returns a list of
#                  dynamics (one element per dimension).
#   lower_limit  - numeric vector; lower bounds for the initial guesses in each
#                  dimension.
#   upper_limit  - numeric vector; upper bounds for the initial guesses in each
#                  dimension.
#   iter         - numeric; number of iteration steps per dimension (total
#                  combinations = iter^dimensions).
#
# RETURNS:
#   A data frame of unique equilibrium solutions, rounded to two decimal places.
# ------------------------------------------------------------------------------
multisolve2 <- function(model, lower_limit, upper_limit, iter) {
  
  # Initialize an empty data frame to store equilibria
  equilibrium <- data.frame(
    matrix(0, ncol = length(model(0)), nrow = 1)
  )
  colnames(equilibrium) <- names(model(0))
  
  # Generate sequences for each dimension
  lower_upper_list <- lapply(
    1:length(lower_limit),
    function(i) seq(lower_limit[i], upper_limit[i], length.out = iter)
  )
  
  # Generate all possible combinations (including zeros row)
  combinations <- expand.grid(lower_upper_list)
  colnames(combinations) <- names(model(0))
  combinations <- rbind(rep(0, length(lower_limit)), combinations)
  
  # Iterate over each combination
  for (i in 1:dim(combinations)[1]) {
    
    solutions_found <- FALSE
    
    # Attempt to find a root with multiroot
    tryCatch({
      solutions <- multiroot(
        f = model,
        start = as.numeric(combinations[i, ]),
        positive = TRUE,
        useFortran = TRUE
      )
      solutions_found <- TRUE
    },
    warning = function(w) {
      # Silently handle warnings
    },
    error = function(e) {
      # Silently handle errors
    })
    
    # If no solution found, move on
    if (!solutions_found) {
      next
    }
    
    # If first solution, save directly
    if (i == 1) {
      equilibrium[1, ] <- solutions$root
    } else {
      # Otherwise, check if already found
      for (j in 1:dim(equilibrium)[2]) {
        if (j == 1) {
          current <- equilibrium %>%
            filter(
              (solutions$root[j] * 0.99 - max(upper_limit) * 0.0001)
              <= (!!sym(colnames(equilibrium)[j]))
            ) %>%
            filter(
              (!!sym(colnames(equilibrium)[j]))
              <= (solutions$root[j] * 1.01 + max(upper_limit) * 0.0001)
            )
        } else {
          current <- current %>%
            filter(
              (solutions$root[j] * 0.99 - max(upper_limit) * 0.0001)
              <= (!!sym(colnames(equilibrium)[j]))
            ) %>%
            filter(
              (!!sym(colnames(equilibrium)[j]))
              <= (solutions$root[j] * 1.01 + max(upper_limit) * 0.0001)
            )
        }
      }
      # If not in our list, add it
      if (dim(current)[1] == 0) {
        equilibrium <- rbind(equilibrium, solutions$root)
      }
    }
  }
  
  # Return solutions, rounded to five decimals
  return(round(equilibrium, 5))
}

# ------------------------------------------------------------------------------
# multisolve()
#
# DESCRIPTION:
#   Main function to find equilibrium solutions. It attempts to use searchZeros
#   first (from 'nleqslv'), and if that fails to return a solution, it calls
#   'multisolve2'.
#
# ARGUMENTS:
#   model        - function; accepts a numeric vector and returns a list of
#                  dynamics (one element per dimension).
#   lower_limit  - numeric vector; lower bounds for initial guesses in each
#                  dimension.
#   upper_limit  - numeric vector; upper bounds for initial guesses in each
#                  dimension.
#   iter         - numeric; number of iteration steps per dimension (total
#                  combinations = iter^dimensions).
#
# RETURNS:
#   A data frame of unique equilibrium solutions, all >= 0, and rounded to
#   two decimal places.
# ------------------------------------------------------------------------------
multisolve <- function(model, lower_limit, upper_limit, iter) {
  
  # Generate sequences for each dimension
  lower_upper_list <- lapply(
    1:length(lower_limit),
    function(i) seq(lower_limit[i], upper_limit[i], length.out = iter)
  )
  
  # Generate combinations (include zeros row)
  combinations <- expand.grid(lower_upper_list)
  colnames(combinations) <- names(model(0))
  combinations <- rbind(rep(0, length(lower_limit)), combinations)
  
  # Solve using searchZeros
  equilibrium <- searchZeros(
    as.matrix(combinations), 
    fn = model,
    control = list(xtol = 10e-12, ftol = 10e-12, allowSingular = TRUE)
  )$x
  
  # If no solution from searchZeros, fall back to brute force
  if (is.null(equilibrium)) {
    equilibrium <- multisolve2(model, lower_limit, upper_limit, iter)
  }
  
  # Clean up solutions
  equilibrium <- round(as.data.frame(equilibrium), 5)
  equilibrium <- equilibrium %>%
    select_if(is.numeric) %>%
    filter(apply(., 1, function(x) all(x >= 0)))
  
  equilibrium <- distinct(equilibrium)
  
  return(equilibrium)
}

# ------------------------------------------------------------------------------
# model()
#
# DESCRIPTION:
#   Defines the seagrass–soil model. Takes a numeric vector x = [S, CA, CS,
#   CB, N] and returns a named vector of the system's dynamics 
#   (dS/dt, dCA/dt, dCS/dt, dCB/dt, dN/dt).
#
# ARGUMENTS:
#   x - numeric vector of state variables:
#       x[1] = S   (Seagrass biomass)
#       x[2] = CA  (Active carbon)
#       x[3] = CS  (Slow carbon)
#       x[4] = CB  (Buried carbon)
#       x[5] = N   (Nutrients)
#
# RETURNS:
#   A named numeric vector with rates of change for each state variable
#   (S, CA, CS, CB, N).
# ------------------------------------------------------------------------------
model <- function(x) {
  
  S <- rmax * ((tmax - temp) / (tmax - topt)) *
    ((temp / topt)^(topt / (tmax - topt))) *
    x[5] / (x[5] + nr) *
    (i0 * exp(-a * pmax * z) / (i0 * exp(-a * pmax * z) + ir)) *
    x[1] -
    x[1] * (ms + ma + mh * hmax * sh / (sh + x[1]))
  
  CA <- alpha * h * pmax *
    (1 - (hmax * sh / (sh + x[1])) /
       (hmax * sh / (sh + x[1]) + hp)) -
    x[2] *
    (phida * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
       (hmax * sh / (sh + x[1])) /
       (hmax * sh / (sh + x[1]) + hea) +
       1 / (1 + exp(-fb1 * (x[2] + x[3] - fb2))))
  
  CS <- beta * x[1] * (ms + mh * hmax * sh / (sh + x[1])) -
    x[3] *
    (phids * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
       (hmax * sh / (sh + x[1])) /
       (hmax * sh / (sh + x[1]) + hes) +
       1 / (1 + exp(-fb1 * (x[2] + x[3] - fb2))))
  
  CB <- x[2] * 1 / (1 + exp(-fb1 * (x[2] + x[3] - fb2))) +
    x[3] * 1 / (1 + exp(-fb1 * (x[2] + x[3] - fb2))) -
    (1 - ma * da) * x[4] * phidb *
    1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
    ma * da * x[4] * (phida + phids) / 2 *
    1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
  
  N <- i + gamma * (
    x[2] * phida * 1 / (exp(-ea / (8.314 * td))) *
      exp(-ea / (8.314 * temp)) +
      x[3] * phids * 1 / (exp(-ea / (8.314 * td))) *
      exp(-ea / (8.314 * temp)) +
      (1 - ma * da) * x[4] * phidb *
      1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
      ma * da * x[4] * (phida + phids) / 2 *
      1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
  ) -
    delta *
    (rmax * ((tmax - temp) / (tmax - topt)) *
       ((temp / topt)^(topt / (tmax - topt))) *
       x[5] / (x[5] + nr) *
       (i0 * exp(-a * pmax * z) /
          (i0 * exp(-a * pmax * z) + ir))
    ) * x[1] -
    l * x[5]
  
  c(S = S, CA = CA, CS = CS, CB = CB, N = N)
}

# ------------------------------------------------------------------------------
# Parameter assignments and initial equilibrium calculations
# ------------------------------------------------------------------------------
rmax  <- 0.011
ms    <- 0.0014
topt  <- 298.85
tmax  <- 307.05
nr    <- 0.015
i0    <- 709
ir    <- 296
sh    <- 22.7
hp    <- 0.023
a     <- 0.07
mh    <- 0.028
alpha <- 0.194
h     <- 0.68
phida <- 0.024
phids <- 0.00024
phidb <- 0.00005
ea    <- 58000
r     <- 8.314
td    <- 303.65
hea   <- 12.53
hes   <- 12.53
fb1   <- 0.0022
fb2   <- 1335
beta  <- 0.292
gamma <- 0.137
delta <- 0.052
l     <- 0.01
da    <- 0.125 / 6

z      <- 10
temp   <- 298.85
hmax   <- 0.15
pmax   <- 0.5
ma     <- 0
i      <- 0.005865335

iv1_max <- 5000
iv2_max <- 5
iv3_max <- 5
iv4_max <- 5
iv5_max <- 10000

upper_limit <- c(iv1_max, iv2_max, iv3_max, iv4_max, iv5_max)
lower_limit <- c(10, 0, 0, 0, 0)
iter        <- 4

# Initial equilibrium (max seagrass)
initial_eq <- multisolve(model, lower_limit, upper_limit, iter) %>%
  filter(S == max(S))

# ------------------------------------------------------------------------------
# Parameter values for transient solving
# ------------------------------------------------------------------------------
params <- c(
  rmax = 0.011, ms = 0.0014, topt = 298.85, tmax = 307.05, nr = 0.015,
  i0 = 709, ir = 296, sh = 22.7, hp = 0.023, a = 0.07, mh = 0.028,
  alpha = 0.194, h = 0.68, phida = 0.024, phids = 0.00024, phidb = 0.00005,
  ea = 58000, r = 8.314, td = 303.65, hea = 12.53, hes = 12.53,
  fb1 = 0.0022, fb2 = 1335, beta = 0.292, gamma = 0.137, delta = 0.052,
  l = 0.01, da = 0.125 / 6, z = 10, hmax = 0.15, i = 0.005865335
)

pmax_0 <- 0.5    # pmax initial
pmax_f <- 4.14   # pmax final

temp_0 <- 298.85 # temp initial
temp_f <- 308    # temp final

ma_0 <- 0        # ma initial
ma_f <- 0.00330  # ma final

# Calculate final equilibria for each stress scenario
pmax <- pmax_f
final_eq_pmax <- multisolve(model, lower_limit, upper_limit, iter) %>% filter(S == max(S))
pmax <- pmax_0

temp <- temp_f
final_eq_temp <- multisolve(model, lower_limit, upper_limit, iter) %>% filter(S == max(S))
temp <- temp_0

ma <- ma_f
final_eq_ma <- multisolve(model, lower_limit, upper_limit, iter) %>% filter(S == max(S))
ma <- ma_0

precision <- 0.05

rate_list <- c(10 * 365, 1000 * 365)
max_rate  <- max(rate_list)

# ------------------------------------------------------------------------------
# Loop over parameters for bifurcation data and transient simulations
# ------------------------------------------------------------------------------
for (par in c("ma", "pmax", "temp")) {
  
  par_f <- get(paste0(par, "_f"))
  
  # --------------------------------------------------------------------------
  # Read & filter equilibrium data for seagrass, bistable, bare states
  # --------------------------------------------------------------------------
  assign(
    paste0("data_", par, "_seagrass"),
    read.csv(paste0(par, "_equilibrium.csv")) %>%
      filter(stability == "seagrass" | stability == "bistable") %>%
      filter(!!sym(par) <= par_f * 1.01)
  )
  
  assign(
    paste0("data_", par, "_bistable"),
    read.csv(paste0(par, "_equilibrium.csv")) %>%
      filter(stability == "bistable" & !!sym(par) <= par_f * 1.01)
  )
  
  assign(
    paste0("data_", par, "_bare"),
    read.csv(paste0(par, "_equilibrium.csv")) %>%
      filter(stability == "bistable" | stability == "bare") %>%
      filter(!!sym(par) <= par_f * 1.01)
  )
  
  # --------------------------------------------------------------------------
  # Build S (Seagrass) bifurcation diagram
  # --------------------------------------------------------------------------
  assign(
    paste0("plot_", par, "_S"),
    ggplot() +
      geom_line(
        data = get(paste0("data_", par, "_seagrass")),
        aes(x = !!sym(par), y = S_seagrass),
        color = "black", size = 1, linetype = "solid"
      ) +
      xlab(par) + ylab("S") +
      theme_classic() +
      theme(legend.position = "none")
  )
  
  assign(
    paste0("plot_", par, "_S"),
    get(paste0("plot_", par, "_S")) +
      geom_line(
        data = get(paste0("data_", par, "_bistable")),
        aes(x = !!sym(par), y = S_bistable),
        color = "black", size = 1, linetype = "dashed"
      )
  )
  
  assign(
    paste0("plot_", par, "_S"),
    get(paste0("plot_", par, "_S")) +
      geom_line(
        data = get(paste0("data_", par, "_bare")),
        aes(x = !!sym(par), y = S_bare),
        color = "black", size = 1, linetype = "solid"
      )
  )
  
  # --------------------------------------------------------------------------
  # Build CB (Stable Carbon) bifurcation diagram
  # --------------------------------------------------------------------------
  assign(
    paste0("data_", par, "_seagrass"),
    read.csv(paste0(par, "_equilibrium.csv")) %>%
      filter(stability == "seagrass" | stability == "bistable") %>%
      filter(!!sym(par) <= par_f * 1.01)
  )
  
  assign(
    paste0("data_", par, "_bistable"),
    read.csv(paste0(par, "_equilibrium.csv")) %>%
      filter(stability == "bistable" & !!sym(par) <= par_f * 1.01)
  )
  
  assign(
    paste0("data_", par, "_bare"),
    read.csv(paste0(par, "_equilibrium.csv")) %>%
      filter(stability == "bistable" | stability == "bare") %>%
      filter(!!sym(par) <= par_f * 1.01)
  )
  
  assign(
    paste0("plot_", par, "_CB"),
    ggplot() +
      geom_line(
        data = get(paste0("data_", par, "_seagrass")),
        aes(x = !!sym(par), y = CB_seagrass),
        color = "black", size = 1, linetype = "solid"
      ) +
      xlab(par) + ylab("CB") +
      theme_classic() +
      theme(legend.position = "none")
  )
  
  assign(
    paste0("plot_", par, "_CB"),
    get(paste0("plot_", par, "_CB")) +
      geom_line(
        data = get(paste0("data_", par, "_bistable")),
        aes(x = !!sym(par), y = CB_bistable),
        color = "black", size = 1, linetype = "dashed"
      )
  )
  
  assign(
    paste0("plot_", par, "_CB"),
    get(paste0("plot_", par, "_CB")) +
      geom_line(
        data = get(paste0("data_", par, "_bare")),
        aes(x = !!sym(par), y = CB_bare),
        color = "black", size = 1, linetype = "solid"
      )
  )
  
  # --------------------------------------------------------------------------
  # Transient simulations for rates in rate_list
  # --------------------------------------------------------------------------
  for (rate in rate_list) {
    
    times_1 <- seq(0, rate, 1)
    start_1 <- initial_eq
    
    assign("final_eq", get(paste0("final_eq_", par)))
    
    # ------------------------------------------------------------------------
    # Define time-varying seagrass-soil function based on parameter
    # ------------------------------------------------------------------------
    if (par == "pmax") {
      seagrass <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          
          temp <- 298.85
          pmax <- 0.5 + (pmax_f - pmax_0) * t * 1 / rate
          ma   <- 0
          
          dS <- rmax * ((tmax - temp) / (tmax - topt)) *
            ((temp / topt)^(topt / (tmax - topt))) *
            N / (N + nr) *
            (i0 * exp(-a * pmax * z) /
               (i0 * exp(-a * pmax * z) + ir)) *
            S -
            S * (ms + ma + mh * hmax * sh / (sh + S))
          
          dCA <- alpha * h * pmax *
            (1 - (hmax * sh / (sh + S)) /
               (hmax * sh / (sh + S) + hp)) -
            CA * (
              phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hea) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
            CS * (
              phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hes) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCB <- CA * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) +
            CS * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) -
            (1 - ma * da) * CB * phidb *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
            ma * da * CB * (phida + phids) / 2 *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
          
          dN <- i +
            gamma * (
              CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (1 - ma * da) * CB * phidb *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                ma * da * CB * (phida + phids) / 2 *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
            ) -
            delta * (
              rmax * ((tmax - temp) / (tmax - topt)) *
                ((temp / topt)^(topt / (tmax - topt))) *
                N / (N + nr) *
                (i0 * exp(-a * pmax * z) /
                   (i0 * exp(-a * pmax * z) + ir))
            ) * S -
            l * N
          
          list(c(dS, dCA, dCS, dCB, dN))
        })
      }
      
      seagrass2 <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          
          temp <- 298.85
          pmax <- pmax_f
          ma   <- 0
          
          dS <- rmax * ((tmax - temp) / (tmax - topt)) *
            ((temp / topt)^(topt / (tmax - topt))) *
            N / (N + nr) *
            (i0 * exp(-a * pmax * z) /
               (i0 * exp(-a * pmax * z) + ir)) *
            S -
            S * (ms + ma + mh * hmax * sh / (sh + S))
          
          dCA <- alpha * h * pmax *
            (1 - (hmax * sh / (sh + S)) /
               (hmax * sh / (sh + S) + hp)) -
            CA * (
              phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hea) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
            CS * (
              phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hes) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCB <- CA * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) +
            CS * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) -
            (1 - ma * da) * CB * phidb *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
            ma * da * CB * (phida + phids) / 2 *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
          
          dN <- i +
            gamma * (
              CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (1 - ma * da) * CB * phidb *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                ma * da * CB * (phida + phids) / 2 *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
            ) -
            delta * (
              rmax * ((tmax - temp) / (tmax - topt)) *
                ((temp / topt)^(topt / (tmax - topt))) *
                N / (N + nr) *
                (i0 * exp(-a * pmax * z) /
                   (i0 * exp(-a * pmax * z) + ir))
            ) * S -
            l * N
          
          list(c(dS, dCA, dCS, dCB, dN))
        })
      }
      
    } else if (par == "temp") {
      
      seagrass <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          
          temp <- 298.85 + (temp_f - temp_0) * t * 1 / rate
          pmax <- 0.5
          ma   <- 0
          
          dS <- rmax * ((tmax - temp) / (tmax - topt)) *
            ((temp / topt)^(topt / (tmax - topt))) *
            N / (N + nr) *
            (i0 * exp(-a * pmax * z) /
               (i0 * exp(-a * pmax * z) + ir)) *
            S -
            S * (ms + ma + mh * hmax * sh / (sh + S))
          
          dCA <- alpha * h * pmax *
            (1 - (hmax * sh / (sh + S)) /
               (hmax * sh / (sh + S) + hp)) -
            CA * (
              phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hea) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
            CS * (
              phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hes) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCB <- CA * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) +
            CS * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) -
            (1 - ma * da) * CB * phidb *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
            ma * da * CB * (phida + phids) / 2 *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
          
          dN <- i +
            gamma * (
              CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (1 - ma * da) * CB * phidb *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                ma * da * CB * (phida + phids) / 2 *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
            ) -
            delta * (
              rmax * ((tmax - temp) / (tmax - topt)) *
                ((temp / topt)^(topt / (tmax - topt))) *
                N / (N + nr) *
                (i0 * exp(-a * pmax * z) /
                   (i0 * exp(-a * pmax * z) + ir))
            ) * S -
            l * N
          
          list(c(dS, dCA, dCS, dCB, dN))
        })
      }
      
      seagrass2 <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          
          temp <- temp_f
          pmax <- 0.5
          ma   <- 0
          
          dS <- rmax * ((tmax - temp) / (tmax - topt)) *
            ((temp / topt)^(topt / (tmax - topt))) *
            N / (N + nr) *
            (i0 * exp(-a * pmax * z) /
               (i0 * exp(-a * pmax * z) + ir)) *
            S -
            S * (ms + ma + mh * hmax * sh / (sh + S))
          
          dCA <- alpha * h * pmax *
            (1 - (hmax * sh / (sh + S)) /
               (hmax * sh / (sh + S) + hp)) -
            CA * (
              phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hea) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
            CS * (
              phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hes) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCB <- CA * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) +
            CS * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) -
            (1 - ma * da) * CB * phidb *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
            ma * da * CB * (phida + phids) / 2 *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
          
          dN <- i +
            gamma * (
              CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (1 - ma * da) * CB * phidb *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                ma * da * CB * (phida + phids) / 2 *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
            ) -
            delta * (
              rmax * ((tmax - temp) / (tmax - topt)) *
                ((temp / topt)^(topt / (tmax - topt))) *
                N / (N + nr) *
                (i0 * exp(-a * pmax * z) /
                   (i0 * exp(-a * pmax * z) + ir))
            ) * S -
            l * N
          
          list(c(dS, dCA, dCS, dCB, dN))
        })
      }
      
    } else {
      
      seagrass <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          
          temp <- 298.85
          pmax <- 0.5
          ma   <- 0 + (ma_f - ma_0) * t * 1 / rate
          
          dS <- rmax * ((tmax - temp) / (tmax - topt)) *
            ((temp / topt)^(topt / (tmax - topt))) *
            N / (N + nr) *
            (i0 * exp(-a * pmax * z) /
               (i0 * exp(-a * pmax * z) + ir)) *
            S -
            S * (ms + ma + mh * hmax * sh / (sh + S))
          
          dCA <- alpha * h * pmax *
            (1 - (hmax * sh / (sh + S)) /
               (hmax * sh / (sh + S) + hp)) -
            CA * (
              phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hea) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
            CS * (
              phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hes) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCB <- CA * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) +
            CS * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) -
            (1 - ma * da) * CB * phidb *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
            ma * da * CB * (phida + phids) / 2 *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
          
          dN <- i +
            gamma * (
              CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (1 - ma * da) * CB * phidb *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                ma * da * CB * (phida + phids) / 2 *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
            ) -
            delta * (
              rmax * ((tmax - temp) / (tmax - topt)) *
                ((temp / topt)^(topt / (tmax - topt))) *
                N / (N + nr) *
                (i0 * exp(-a * pmax * z) /
                   (i0 * exp(-a * pmax * z) + ir))
            ) * S -
            l * N
          
          list(c(dS, dCA, dCS, dCB, dN))
        })
      }
      
      seagrass2 <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          
          temp <- 298.85
          pmax <- 0.5
          ma   <- ma_f
          
          dS <- rmax * ((tmax - temp) / (tmax - topt)) *
            ((temp / topt)^(topt / (tmax - topt))) *
            N / (N + nr) *
            (i0 * exp(-a * pmax * z) /
               (i0 * exp(-a * pmax * z) + ir)) *
            S -
            S * (ms + ma + mh * hmax * sh / (sh + S))
          
          dCA <- alpha * h * pmax *
            (1 - (hmax * sh / (sh + S)) /
               (hmax * sh / (sh + S) + hp)) -
            CA * (
              phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hea) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
            CS * (
              phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hes) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCB <- CA * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) +
            CS * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) -
            (1 - ma * da) * CB * phidb *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
            ma * da * CB * (phida + phids) / 2 *
            1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
          
          dN <- i +
            gamma * (
              CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (1 - ma * da) * CB * phidb *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                ma * da * CB * (phida + phids) / 2 *
                1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
            ) -
            delta * (
              rmax * ((tmax - temp) / (tmax - topt)) *
                ((temp / topt)^(topt / (tmax - topt))) *
                N / (N + nr) *
                (i0 * exp(-a * pmax * z) /
                   (i0 * exp(-a * pmax * z) + ir))
            ) * S -
            l * N
          
          list(c(dS, dCA, dCS, dCB, dN))
        })
      }
    }
    
    start_1_vector <- as.numeric(start_1[1, ])
    names(start_1_vector) <- names(start_1)
    
    # Stress increase
    out_1 <- lsoda(
      y = start_1_vector,
      times = times_1,
      func = seagrass,
      parms = params
    )
    
    times_2 <- seq(rate, max_rate * 2, 1)
    start_2 <- out_1[dim(out_1)[1], 2:6]
    
    out_1 <- as.data.frame(out_1) %>%
      mutate(
        !!par := get(paste0(par, "_0")) +
          (get(paste0(par, "_f")) - get(paste0(par, "_0"))) *
          time * 1 / rate
      )
    
    # Then stress remains at its final value
    out_2 <- lsoda(
      y = start_2,
      times = times_2,
      func = seagrass2,
      parms = params
    )
    
    out_2 <- as.data.frame(out_2) %>%
      mutate(
        !!par := get(paste0(par, "_f"))
      )
    
    assign(
      paste0("data_", par),
      as.data.frame(rbind(out_1, out_2[-1, ]))
    )
    
    rm(out_1)
    rm(out_2)
    
    s_asymptotic <- get(paste0("data_", par)) %>%
      mutate(across(everything(), ~round(., digits = 2))) %>%
      filter(
        round(final_eq[, 1] * (1 - precision), 0) <= S &
          S <= round(final_eq[, 1] * (1 + precision), 0)
      )
    
    if (nrow(s_asymptotic) > 0) {
      s_asymptotic <- s_asymptotic %>%
        filter(time == min(time)) %>%
        select(time)
    }
    
    c_asymptotic <- get(paste0("data_", par)) %>%
      mutate(across(everything(), ~round(., digits = 2))) %>%
      filter(
        round(final_eq[, 4] * (1 - precision), 0) <= CB &
          CB <= round(final_eq[, 4] * (1 + precision), 0)
      )
    
    if (nrow(c_asymptotic) > 0) {
      c_asymptotic <- c_asymptotic %>%
        filter(time == min(time)) %>%
        select(time)
    }
    
    index <- 1
    
    # Keep extending time until equilibrium is reached
    while (nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0) {
      
      start_n <- as.numeric(
        get(paste0("data_", par)) %>% 
          filter(time == max(time)) %>%
          select(S, CA, CS, CB, N)
      )
      names(start_n) <- names(
        get(paste0("data_", par)) %>% 
          filter(time == max(time)) %>% 
          select(S, CA, CS, CB, N)
      )
      
      times_n <- seq(max_rate * (1 + index), max_rate * (2 + index), 1)
      out_n <- lsoda(
        y = start_n, 
        times = times_n, 
        func = seagrass2, 
        parms = params
      )
      
      out_n <- as.data.frame(out_n) %>%
        mutate(
          !!par := get(paste0(par, "_f"))
        )
      
      assign(
        paste0("data_", par),
        as.data.frame(rbind(get(paste0("data_", par)), out_n[-1, ]))
      )
      
      rm(out_n)
      
      s_asymptotic <- get(paste0("data_", par)) %>%
        mutate(across(everything(), ~round(., digits = 2))) %>%
        filter(
          round(final_eq[, 1] * (1 - precision), 0) <= S &
            S <= round(final_eq[, 1] * (1 + precision), 0)
        )
      
      if (nrow(s_asymptotic) > 0) {
        s_asymptotic <- s_asymptotic %>%
          filter(time == min(time)) %>%
          select(time)
      }
      
      c_asymptotic <- get(paste0("data_", par)) %>%
        mutate(across(everything(), ~round(., digits = 2))) %>%
        filter(
          round(final_eq[, 4] * (1 - precision), 0) <= CB &
            CB <= round(final_eq[, 4] * (1 + precision), 0)
        )
      
      if (nrow(c_asymptotic) > 0) {
        c_asymptotic <- c_asymptotic %>%
          filter(time == min(time)) %>%
          select(time)
      }
      
      index <- index + 1
    }
    
    # ------------------------------------------------------------------------
    # Create color gradient
    # ------------------------------------------------------------------------
    n_segments <- 50
    data_par <- get(paste0("data_", par)) %>%
      filter(time <= max_rate * 2)
    
    data_par$segment <- as.integer(
      cut(data_par[["time"]],
          breaks = n_segments,
          labels = FALSE,
          include.lowest = TRUE
      )
    )
    
    # Jet color palette (white→red)
    jet_colors <- colorRampPalette(c("#FCBBA1", "red"))
    
    # HSV-based blue palette (white→blue)
    hsv_blue <- colorRampPalette(c("#C7E9B4", "blue"))
    
    if (rate == min(rate_list)) {
      palette_colors <- jet_colors(n_segments)
    } else {
      palette_colors <- hsv_blue(n_segments)
    }
    
    assign(
      paste0("plot_", par, "_S"),
      get(paste0("plot_", par, "_S")) +
        geom_line(
          data = data_par,
          aes(x = !!sym(par), y = S, color = factor(segment)),
          size = 1, linetype = "solid"
        ) +
        scale_color_manual(values = palette_colors) +
        new_scale_color()
    )
    
    assign(
      paste0("plot_", par, "_CB"),
      get(paste0("plot_", par, "_CB")) +
        geom_line(
          data = data_par,
          aes(x = !!sym(par), y = CB, color = factor(segment)),
          size = 1, linetype = "solid"
        ) +
        scale_color_manual(values = palette_colors) +
        new_scale_color()
    )
    
    if (rate == min(rate_list)) {
      
      assign(
        paste0("plot_", par),
        ggplot(get(paste0("data_", par)),
               aes(x = time / 365, y = !!sym(par))) +
          geom_line(color = "red", size = 1) +
          scale_x_continuous(
            name = "time (years)",
            limits = c(0, (max_rate * 1.1) / 365)
          ) +
          theme_classic() +
          theme(legend.position = "none")
      )
      
    } else {
      
      assign(
        paste0("plot_", par),
        get(paste0("plot_", par)) +
          geom_line(
            data = get(paste0("data_", par)),
            aes(x = time / 365, y = !!sym(par)),
            color = "blue",
            size = 1,
            linetype = "solid"
          )
      )
      
    }
  }
}

# ------------------------------------------------------------------------------
# Arrange final plots into a grid
# ------------------------------------------------------------------------------
plot <- grid.arrange(
  arrangeGrob(plot_ma, plot_pmax, plot_temp, ncol = 3),
  arrangeGrob(plot_ma_S, plot_pmax_S, plot_temp_S, ncol = 3),
  arrangeGrob(plot_ma_CB, plot_pmax_CB, plot_temp_CB, ncol = 3),
  heights = c(0.7, 1, 1)  # first row is shorter
)

ggsave(
  "figure1.png",
  plot,
  width = 13,
  height = 7,
  dpi = 2300
)

# ------------------------------------------------------------------------------
# Build color bars for legend (white→red, white→blue)
# ------------------------------------------------------------------------------
red_50  <- jet_colors(50)
blue_50 <- hsv_blue(50)

# Dummy data spanning from 0 to 2 * max_rate
df <- data.frame(z = seq(0, max_rate * 2, length.out = 101))

# White→red color bar
legend_bar_red <- ggplot(df, aes(x = 1, y = z, fill = z)) +
  geom_raster() +
  scale_fill_stepsn(
    colours = red_50,
    breaks  = seq(0, max_rate * 2, length.out = 51),
    labels  = c("0", rep("", 49), as.character(max_rate * 2)),
    limits  = c(0, max_rate * 2),
    name    = "Value"
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.key.height = unit(3, "cm"),
    legend.key.width  = unit(0.5, "cm")
  )

# White→blue color bar
legend_bar_blue <- ggplot(df, aes(x = 1, y = z, fill = z)) +
  geom_raster() +
  scale_fill_stepsn(
    colours = blue_50,
    breaks  = seq(0, max_rate * 2, length.out = 51),
    labels  = c("0", rep("", 49), as.character(max_rate * 2)),
    limits  = c(0, max_rate * 2),
    name    = "Value"
  ) +
  theme_void() +
  theme(
    legend.position   = "right",
    legend.key.height = unit(3, "cm"),
    legend.key.width  = unit(0.5, "cm")
  )

# Combine both legends
legend_combined <- plot_grid(legend_bar_red, legend_bar_blue, ncol = 2)

# Save the combined legend figure
ggsave(
  "figure1_legend_bar.png",
  legend_combined,
  width = 6,
  height = 7,
  dpi = 600
)
