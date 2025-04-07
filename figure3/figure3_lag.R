################################################################################
# FILE: figure3_lag.R
#
# DESCRIPTION:
#   This script runs a time-lag analysis for the transient seagrass-soil dynamics.
#   It estimates how long it takes for seagrass biomass (S) and stable carbon
#   (CB) to reach new equilibria in response to gradually increasing stressors
#   (pmax, temp, ma). Different rates of environmental change (rate_values) are
#   tested in parallel, and for each scenario, a "lag" (the delay between
#   seagrass collapse and carbon loss) is computed and saved in CSV format.
#
# DEPENDENCIES:
#   Requires the following R packages:
#     - deSolve   (for ODE solving)
#     - dplyr
#     - nleqslv
#     - rootSolve
#     - doMPI     (for parallelization)
################################################################################

# ------------------------------------------------------------------------------
# Load Required Packages
# ------------------------------------------------------------------------------
package_list <- c("deSolve", "dplyr", "nleqslv", "rootSolve")
for (package in package_list) {
  library(package, character.only = TRUE)
}

library(doMPI)

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

rate_values <- round(seq(1000 * 365, 1 * 365, length.out = 64), 0)
max_rate    <- max(rate_values)

pmax_0 <- 0.5
pmax_f <- 4.14

temp_0 <- 298.85
temp_f <- 308

ma_0   <- 0
ma_f   <- 0.00330

# ------------------------------------------------------------------------------
# Final Equilibria for Each Stressor
# ------------------------------------------------------------------------------
pmax <- pmax_f
final_eq_pmax <- multisolve(model, lower_limit, upper_limit, iter) %>%
  filter(S == max(S))
pmax <- pmax_0

temp <- temp_f
final_eq_temp <- multisolve(model, lower_limit, upper_limit, iter) %>%
  filter(S == max(S))
temp <- temp_0

ma <- ma_f
final_eq_ma <- multisolve(model, lower_limit, upper_limit, iter) %>%
  filter(S == max(S))
ma <- ma_0

precision <- 0.05

# ------------------------------------------------------------------------------
# Start Parallelization
# ------------------------------------------------------------------------------
cl <- startMPIcluster()
registerDoMPI(cl)

# ------------------------------------------------------------------------------
# Time-Lag Estimation for Each Stressor and Rate
# ------------------------------------------------------------------------------
for (par1 in c("pmax", "temp", "ma")) {
  
  # Parallel loop over rate_values
  foreach(rate = rate_values,
          .packages = c("rootSolve", "dplyr", "rlang", "deSolve")) %dopar% {
            
            # Simulation from 0..rate (stress increase)
            times_1 <- seq(0, rate, 1)
            start_1 <- initial_eq
            
            # Retrieve final equilibrium object (already computed)
            assign("final_eq", get(paste0("final_eq_", par1)))
            
            # --------------------------------------------------------------------------
            # Define forcing functions for each stressor
            #   seagrass()  - ramp up the stress
            #   seagrass2() - hold final stress
            # --------------------------------------------------------------------------
            if (par1 == "pmax") {
              
              # Ramps pmax from pmax_0 to pmax_f
              seagrass <- function(t, state, parameters) {
                with(as.list(c(state, parameters)), {
                  
                  temp  <- 298.85
                  pmax  <- 0.5 + (pmax_f - pmax_0) * t * 1 / rate
                  ma    <- 0
                  
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
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp)) -
                    ma * da * CB * (phida + phids) / 2 *
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp))
                  
                  dN <- i +
                    gamma * (
                      CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        (1 - ma * da) * CB * phidb *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        ma * da * CB * (phida + phids) / 2 *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp))
                    ) -
                    delta * (
                      rmax * ((tmax - temp) / (tmax - topt)) *
                        ((temp / topt)^(topt / (tmax - topt))) *
                        N / (N + nr) *
                        (i0 * exp(-a * pmax * z) /
                           (i0 * exp(-a * pmax * z) + ir)
                        )
                    ) * S -
                    l * N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                })
              }
              
              # Holds pmax at pmax_f
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
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp)) -
                    ma * da * CB * (phida + phids) / 2 *
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp))
                  
                  dN <- i +
                    gamma * (
                      CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        (1 - ma * da) * CB * phidb *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        ma * da * CB * (phida + phids) / 2 *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp))
                    ) -
                    delta * (
                      rmax * ((tmax - temp) / (tmax - topt)) *
                        ((temp / topt)^(topt / (tmax - topt))) *
                        N / (N + nr) *
                        (i0 * exp(-a * pmax * z) /
                           (i0 * exp(-a * pmax * z) + ir)
                        )
                    ) * S -
                    l * N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                })
              }
              
            } else if (par1 == "temp") {
              
              # Ramps temp from temp_0 to temp_f
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
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp)) -
                    ma * da * CB * (phida + phids) / 2 *
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp))
                  
                  dN <- i +
                    gamma * (
                      CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        (1 - ma * da) * CB * phidb *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        ma * da * CB * (phida + phids) / 2 *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp))
                    ) -
                    delta * (
                      rmax * ((tmax - temp) / (tmax - topt)) *
                        ((temp / topt)^(topt / (tmax - topt))) *
                        N / (N + nr) *
                        (i0 * exp(-a * pmax * z) /
                           (i0 * exp(-a * pmax * z) + ir)
                        )
                    ) * S -
                    l * N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                })
              }
              
              # Holds temp at temp_f
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
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp)) -
                    ma * da * CB * (phida + phids) / 2 *
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp))
                  
                  dN <- i +
                    gamma * (
                      CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        (1 - ma * da) * CB * phidb *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        ma * da * CB * (phida + phids) / 2 *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp))
                    ) -
                    delta * (
                      rmax * ((tmax - temp) / (tmax - topt)) *
                        ((temp / topt)^(topt / (tmax - topt))) *
                        N / (N + nr) *
                        (i0 * exp(-a * pmax * z) /
                           (i0 * exp(-a * pmax * z) + ir)
                        )
                    ) * S -
                    l * N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                })
              }
              
            } else {
              
              # Ramps ma from ma_0 to ma_f
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
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp)) -
                    ma * da * CB * (phida + phids) / 2 *
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp))
                  
                  dN <- i +
                    gamma * (
                      CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        (1 - ma * da) * CB * phidb *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        ma * da * CB * (phida + phids) / 2 *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp))
                    ) -
                    delta * (
                      rmax * ((tmax - temp) / (tmax - topt)) *
                        ((temp / topt)^(topt / (tmax - topt))) *
                        N / (N + nr) *
                        (i0 * exp(-a * pmax * z) /
                           (i0 * exp(-a * pmax * z) + ir)
                        )
                    ) * S -
                    l * N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                })
              }
              
              # Holds ma at ma_f
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
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp)) -
                    ma * da * CB * (phida + phids) / 2 *
                    1 / (exp(-ea / (8.314 * td))) *
                    exp(-ea / (8.314 * temp))
                  
                  dN <- i +
                    gamma * (
                      CA * phida * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        CS * phids * 1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        (1 - ma * da) * CB * phidb *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp)) +
                        ma * da * CB * (phida + phids) / 2 *
                        1 / (exp(-ea / (8.314 * td))) *
                        exp(-ea / (8.314 * temp))
                    ) -
                    delta * (
                      rmax * ((tmax - temp) / (tmax - topt)) *
                        ((temp / topt)^(topt / (tmax - topt))) *
                        N / (N + nr) *
                        (i0 * exp(-a * pmax * z) /
                           (i0 * exp(-a * pmax * z) + ir)
                        )
                    ) * S -
                    l * N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                })
              }
            }
            
            # --------------------------------------------------------------------------
            # Run Phase 1: Stress from _0 to _f
            # --------------------------------------------------------------------------
            start_1_vector <- as.numeric(start_1[1, ])
            names(start_1_vector) <- names(start_1)
            
            out_1 <- lsoda(
              y = start_1_vector,
              times = times_1,
              func = seagrass,
              parms = params
            )
            
            # --------------------------------------------------------------------------
            # Run Phase 2: Hold final stress
            # --------------------------------------------------------------------------
            times_2 <- seq(rate, max_rate * 2, 1)
            start_2 <- out_1[nrow(out_1), 2:6]
            
            out_1 <- as.data.frame(out_1) %>%
              mutate(
                !!par1 :=
                  get(paste0(par1, "_0")) +
                  (get(paste0(par1, "_f")) - get(paste0(par1, "_0"))) * time * 1 / rate
              )
            
            out_2 <- lsoda(
              y = start_2,
              times = times_2,
              func = seagrass2,
              parms = params
            )
            
            out_2 <- as.data.frame(out_2) %>%
              mutate(
                !!par1 := get(paste0(par1, "_f"))
              )
            
            data <- rbind(out_1, out_2[-1, ])
            rm(out_1, out_2)
            
            # --------------------------------------------------------------------------
            # Check times at which S and CB approach final equilibrium
            # --------------------------------------------------------------------------
            s_asymptotic <- data %>%
              mutate(across(everything(), ~ round(., digits = 2))) %>%
              filter(
                round(final_eq[, 1] * (1 - precision), 0) <= S &
                  S <= round(final_eq[, 1] * (1 + precision), 0)
              )
            
            if (nrow(s_asymptotic) > 0) {
              s_asymptotic <- s_asymptotic %>%
                filter(time == min(time)) %>%
                select(time)
            }
            
            c_asymptotic <- data %>%
              mutate(across(everything(), ~ round(., digits = 2))) %>%
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
            
            # Extend simulation if system hasn't converged
            while (nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0) {
              
              start_n <- as.numeric(data[nrow(data), 2:6])
              names(start_n) <- names(data[, 2:6])
              
              times_n <- seq(max_rate * (1 + index), max_rate * (2 + index), 1)
              out_n <- lsoda(y = start_n, times = times_n, func = seagrass2, parms = params)
              
              out_n <- as.data.frame(out_n) %>%
                mutate(
                  !!par1 := get(paste0(par1, "_f"))
                )
              
              data <- rbind(data, out_n[-1, ])
              rm(out_n)
              
              s_asymptotic <- data %>%
                mutate(across(everything(), ~ round(., digits = 2))) %>%
                filter(
                  round(final_eq[, 1] * (1 - precision), 0) <= S &
                    S <= round(final_eq[, 1] * (1 + precision), 0)
                )
              
              if (nrow(s_asymptotic) > 0) {
                s_asymptotic <- s_asymptotic %>%
                  filter(time == min(time)) %>%
                  select(time)
              }
              
              c_asymptotic <- data %>%
                mutate(across(everything(), ~ round(., digits = 2))) %>%
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
            
            s_asymptotic <- as.numeric(s_asymptotic)
            c_asymptotic <- as.numeric(c_asymptotic)
            
            rm(data)
            
            # Compute lag (delay) between S and CB equilibrium times
            lag <- (c_asymptotic - s_asymptotic) / 365
            
            # Write scenario results to CSV
            write.csv(
              data.frame(scenario = par1, rate = rate, lag = lag),
              file = paste0(par1, rate, ".csv"),
              row.names = FALSE
            )
          }
}

# ------------------------------------------------------------------------------
# Gather and Combine Results
# ------------------------------------------------------------------------------
assign("results", read.csv(paste0("pmax", rate_values[1], ".csv"), header = TRUE))
unlink(paste0("pmax", rate_values[1], ".csv"))

for (par1 in c("pmax", "temp", "ma")) {
  for (rate in rate_values) {
    
    if (par1 == "temp" || par1 == "ma" || rate != rate_values[1]) {
      assign(paste0(par1, rate),
             read.csv(paste0(par1, rate, ".csv"), header = TRUE))
      
      results <- rbind(results, get(paste0(par1, rate)))
      
      unlink(paste0(par1, rate, ".csv"))
      rm(list = paste0(par1, rate))
    }
  }
}

write.csv(results, file = "lag_rate.csv", row.names = FALSE)

# ------------------------------------------------------------------------------
# End Parallelization
# ------------------------------------------------------------------------------
closeCluster(cl)
mpi.quit()
