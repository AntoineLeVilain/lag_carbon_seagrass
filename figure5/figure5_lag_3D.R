################################################################################
# FILE: figure5_lag_3D.R
#
# DESCRIPTION:
#     This script runs a time-lag analysis in three dimensions (pmax, temp, ma)
#     for a seagrass–soil model, accounting for different rates of
#     parameter change. It uses a multi-phase ODE approach, checking whether
#     seagrass biomass (S) and stable carbon (CB) reach equilibrium within each
#     incremental phase of stressor increase. If not, it continues iterating
#     until the system is asymptotic or until a certain limit. The final lag
#     time is computed as the difference between CB and S reaching their new
#     stable states. Results are exported as CSV files (one for each iteration),
#     then combined into "lag_rate_3D.csv".
#
# DEPENDENCIES:
#   Requires the following R packages:
#     - deSolve      (for ODE solving)
#     - dplyr
#     - nleqslv
#     - rootSolve
#     - doMPI        (for parallelization)
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

pmax_0 <- 0.5    # pmax initial value
pmax_f <- 4.14   # pmax final value
temp_0 <- 298.85 # temp initial value
temp_f <- 308    # temp final value
ma_0   <- 0      # ma initial value
ma_f   <- 0.00330 # ma final value

pmax <- pmax_f
temp <- temp_f
ma   <- ma_f

final_eq <- multisolve(model, lower_limit, upper_limit, iter) %>%
  filter(S == max(S))

pmax <- pmax_0
temp <- temp_0
ma   <- ma_0

# For "slow" vs. "fast" parameter changes (in days)
slow <- 365000  # e.g., 1000 years
fast <- 365     # e.g., 1 year

# --------------------------------------------------------------------------------
# Define the number of values per parameter (e.g., 10)
# --------------------------------------------------------------------------------

number_values <- 10
pmax_values <- as.integer(seq(slow, fast, length.out = number_values))
temp_values <- as.integer(seq(slow, fast, length.out = number_values))
ma_values   <- as.integer(seq(slow, fast, length.out = number_values))

# Create the dataframe of all possible combinations
values <- expand.grid(pmax = pmax_values, temp = temp_values, ma = ma_values)

precision   <- 0.05
precision_c <- 10
precision_s <- 10

# --------------------------------------------------------------------------------
# Start parallelization
# --------------------------------------------------------------------------------

cl <- startMPIcluster()
registerDoMPI(cl)

# ------------------------------------------------------------------------------
# Analyze Rate Combinations in Parallel
# ------------------------------------------------------------------------------

foreach(
  rate_index = 1:dim(values)[1],
  .packages = c("rootSolve", "dplyr", "rlang", "deSolve")
) %dopar% {
  
  ma   <- ma_0
  pmax <- pmax_0
  temp <- temp_0
  
  rates <- values[rate_index, ]
  
  # Transpose the rates df to sort columns by the first row
  rates <- as.data.frame(t(rates))
  # Sort by ascending numeric values
  rates <- rates[order(rates[, 1], decreasing = FALSE), , drop = FALSE]
  # Transpose back
  rates <- as.data.frame(t(rates))
  
  rate_pmax <- as.numeric(rates %>% select(pmax))
  rate_temp <- as.numeric(rates %>% select(temp))
  rate_ma   <- as.numeric(rates %>% select(ma))
  
  min_par <- colnames(rates[1])
  mid_rate <- as.numeric(rates[2])
  mid_par <- colnames(rates[2])
  max_rate <- as.numeric(rates[3])
  max_par <- colnames(rates[3])
  
  # ----------------------------------------------------------------------------
  # Inner model function for ODE solver (Phase 1)
  # ----------------------------------------------------------------------------
  seagrass <- function(t, state, parameters) {
    with(as.list(c(state, parameters)), {
      
      # Parameters changing linearly with time in the first phase
      temp <- temp_0 + (temp_f - temp_0) * t * 1 / rate_temp
      pmax <- pmax_0 + (pmax_f - pmax_0) * t * 1 / rate_pmax
      ma   <- ma_0 + (ma_f - ma_0) * t * 1 / rate_ma
      
      dS <- rmax * ((tmax - temp) / (tmax - topt)) *
        ((temp / topt)^(topt / (tmax - topt))) *
        N / (N + nr) *
        (i0 * exp(-a * pmax * z) / (i0 * exp(-a * pmax * z) + ir)) *
        S - S * (ms + ma + mh * hmax * sh / (sh + S))
      
      dCA <- alpha * h * pmax *
        (1 - (hmax * sh / (sh + S)) /
           (hmax * sh / (sh + S) + hp)) -
        CA * (
          phida *
            1 / (exp(-ea / (8.314 * td))) *
            exp(-ea / (8.314 * temp)) +
            (hmax * sh / (sh + S)) /
            (hmax * sh / (sh + S) + hea) +
            1 / (1 + exp(-fb1 * (CA + CS - fb2)))
        )
      
      dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
        CS * (
          phids *
            1 / (exp(-ea / (8.314 * td))) *
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
      
      dN <- i + gamma * (
        CA * phida *
          1 / (exp(-ea / (8.314 * td))) *
          exp(-ea / (8.314 * temp)) +
          CS * phids *
          1 / (exp(-ea / (8.314 * td))) *
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
        ) * S - l * N
      
      list(c(dS, dCA, dCS, dCB, dN))
    })
  }
  
  # ----------------------------------------------------------------------------
  # 1st integration phase
  # ----------------------------------------------------------------------------
  min_rate_val <- min(rates)
  times_1 <- seq(0, min_rate_val, 1)
  
  start_1_vector <- as.numeric(as.numeric(initial_eq))
  names(start_1_vector) <- names(initial_eq)
  
  data <- as.data.frame(lsoda(
    y = start_1_vector,
    times = times_1,
    func = seagrass,
    parms = params
  ))
  
  data <- data %>%
    mutate(pmax = pmax_0 + (pmax_f - pmax_0) * time * 1 / rate_pmax) %>%
    mutate(temp = temp_0 + (temp_f - temp_0) * time * 1 / rate_temp) %>%
    mutate(ma   = ma_0 + (ma_f - ma_0) * time * 1 / rate_ma)
  
  # Check if S, CB are near final equilibrium
  s_asymptotic <- data %>%
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
  
  c_asymptotic <- data %>%
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
  
  # ----------------------------------------------------------------------------
  # 2nd integration phase if not asymptotic after phase 1
  # ----------------------------------------------------------------------------
  if (nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0) {
    
    seagrass <- function(t, state, parameters) {
      with(as.list(c(state, parameters)), {
        
        # In phase 2, min_par is set to its final value
        assign(min_par, get(paste0(min_par, "_f")))
        assign(mid_par,
               get(paste0(mid_par, "_0")) +
                 (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0")))
               * t * 1 / mid_rate
        )
        assign(max_par,
               get(paste0(max_par, "_0")) +
                 (get(paste0(max_par, "_f")) - get(paste0(max_par, "_0")))
               * t * 1 / max_rate
        )
        
        dS <- rmax * ((tmax - temp) / (tmax - topt)) *
          ((temp / topt)^(topt / (tmax - topt))) *
          N / (N + nr) *
          (i0 * exp(-a * pmax * z) /
             (i0 * exp(-a * pmax * z) + ir)) * S -
          S * (ms + ma + mh * hmax * sh / (sh + S))
        
        dCA <- alpha * h * pmax *
          (1 - (hmax * sh / (sh + S)) /
             (hmax * sh / (sh + S) + hp)) -
          CA * (
            phida *
              1 / (exp(-ea / (8.314 * td))) *
              exp(-ea / (8.314 * temp)) +
              (hmax * sh / (sh + S)) /
              (hmax * sh / (sh + S) + hea) +
              1 / (1 + exp(-fb1 * (CA + CS - fb2)))
          )
        
        dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
          CS * (
            phids *
              1 / (exp(-ea / (8.314 * td))) *
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
        
        dN <- i + gamma * (
          CA * phida *
            1 / (exp(-ea / (8.314 * td))) *
            exp(-ea / (8.314 * temp)) +
            CS * phids *
            1 / (exp(-ea / (8.314 * td))) *
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
                 (i0 * exp(-a * pmax * z) + ir))
          ) * S - l * N
        
        list(c(dS, dCA, dCS, dCB, dN))
      })
    }
    
    if (length(seq(min_rate_val, mid_rate, 1)) != 1) {
      
      times_2 <- seq(min_rate_val, mid_rate, 1)
      start_2_vector <- as.numeric(
        data %>% filter(time == max(time)) %>% select(S, CA, CS, CB, N)
      )
      names(start_2_vector) <- names(initial_eq)
      
      out_2 <- as.data.frame(
        lsoda(y = start_2_vector, times = times_2,
              func = seagrass, parms = params)
      )
      
      out_2 <- out_2 %>%
        mutate(!!min_par := rep(get(paste0(min_par, "_f")),
                                length(times_2))) %>%
        mutate(!!mid_par :=
                 get(paste0(mid_par, "_0")) +
                 (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0")))
               * time * 1 / mid_rate) %>%
        mutate(!!max_par :=
                 get(paste0(max_par, "_0")) +
                 (get(paste0(max_par, "_f")) - get(paste0(max_par, "_0")))
               * time * 1 / max_rate)
      
      data <- rbind(data, out_2[-1, ])
      rm(out_2)
      
      s_asymptotic <- data %>%
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
      
      c_asymptotic <- data %>%
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
    }
    
    # ------------------------------------------------------------------------
    # 3rd integration phase if still not asymptotic
    # ------------------------------------------------------------------------
    if (nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0) {
      
      seagrass <- function(t, state, parameters) {
        with(as.list(c(state, parameters)), {
          
          # Now min_par and mid_par remain at final,
          # only max_par changes
          assign(min_par, get(paste0(min_par, "_f")))
          assign(mid_par, get(paste0(mid_par, "_f")))
          assign(max_par,
                 get(paste0(max_par, "_0")) +
                   (get(paste0(max_par, "_f")) -
                      get(paste0(max_par, "_0"))) *
                   t * 1 / max_rate
          )
          
          dS <- rmax * ((tmax - temp) / (tmax - topt)) *
            ((temp / topt)^(topt / (tmax - topt))) *
            N / (N + nr) *
            (i0 * exp(-a * pmax * z) /
               (i0 * exp(-a * pmax * z) + ir)) * S -
            S * (ms + ma + mh * hmax * sh / (sh + S))
          
          dCA <- alpha * h * pmax *
            (1 - (hmax * sh / (sh + S)) /
               (hmax * sh / (sh + S) + hp)) -
            CA * (
              phida *
                1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                (hmax * sh / (sh + S)) /
                (hmax * sh / (sh + S) + hea) +
                1 / (1 + exp(-fb1 * (CA + CS - fb2)))
            )
          
          dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
            CS * (
              phids *
                1 / (exp(-ea / (8.314 * td))) *
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
          
          dN <- i + gamma * (
            CA * phida *
              1 / (exp(-ea / (8.314 * td))) *
              exp(-ea / (8.314 * temp)) +
              CS * phids *
              1 / (exp(-ea / (8.314 * td))) *
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
            ) * S - l * N
          
          list(c(dS, dCA, dCS, dCB, dN))
        })
      }
      
      if (length(seq(mid_rate, max_rate, 1)) != 1) {
        
        times_3 <- seq(mid_rate, max_rate, 1)
        start_3_vector <- as.numeric(
          data %>% filter(time == max(time)) %>%
            select(S, CA, CS, CB, N)
        )
        names(start_3_vector) <- names(initial_eq)
        
        out_3 <- as.data.frame(lsoda(
          y = start_3_vector,
          times = times_3,
          func = seagrass,
          parms = params
        ))
        
        out_3 <- out_3 %>%
          mutate(!!min_par := rep(get(paste0(min_par, "_f")),
                                  length(times_3))) %>%
          mutate(!!mid_par :=
                   rep(get(paste0(mid_par, "_f")), length(times_3))) %>%
          mutate(!!max_par :=
                   get(paste0(max_par, "_0")) +
                   (get(paste0(max_par, "_f")) - get(paste0(max_par, "_0")))
                 * time * 1 / max_rate)
        
        data <- rbind(data, out_3[-1, ])
        rm(out_3)
        
        s_asymptotic <- data %>%
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
        
        c_asymptotic <- data %>%
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
      }
      
      # --------------------------------------------------------------------
      # Further extension if STILL not asymptotic
      # --------------------------------------------------------------------
      index <- 1
      while (nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0) {
        
        seagrass <- function(t, state, parameters) {
          with(as.list(c(state, parameters)), {
            
            # In further extensions, all parameters remain final
            assign(min_par, get(paste0(min_par, "_f")))
            assign(mid_par, get(paste0(mid_par, "_f")))
            assign(max_par, get(paste0(max_par, "_f")))
            
            dS <- rmax * ((tmax - temp) / (tmax - topt)) *
              ((temp / topt)^(topt / (tmax - topt))) *
              N / (N + nr) *
              (i0 * exp(-a * pmax * z) /
                 (i0 * exp(-a * pmax * z) + ir)) * S -
              S * (ms + ma + mh * hmax * sh / (sh + S))
            
            dCA <- alpha * h * pmax *
              (1 - (hmax * sh / (sh + S)) /
                 (hmax * sh / (sh + S) + hp)) -
              CA * (
                phida *
                  1 / (exp(-ea / (8.314 * td))) *
                  exp(-ea / (8.314 * temp)) +
                  (hmax * sh / (sh + S)) /
                  (hmax * sh / (sh + S) + hea) +
                  1 / (1 + exp(-fb1 * (CA + CS - fb2)))
              )
            
            dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
              CS * (
                phids *
                  1 / (exp(-ea / (8.314 * td))) *
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
            
            dN <- i + gamma * (
              CA * phida *
                1 / (exp(-ea / (8.314 * td))) *
                exp(-ea / (8.314 * temp)) +
                CS * phids *
                1 / (exp(-ea / (8.314 * td))) *
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
              ) * S - l * N
            
            list(c(dS, dCA, dCS, dCB, dN))
          })
        }
        
        times_n <- seq(max_rate * (index), max_rate * (index + 1), 1)
        start_n_vector <- as.numeric(
          data %>% filter(time == max(time)) %>%
            select(S, CA, CS, CB, N)
        )
        names(start_n_vector) <- names(initial_eq)
        
        out_n <- as.data.frame(lsoda(
          y = start_n_vector,
          times = times_n,
          func = seagrass,
          parms = params
        ))
        
        out_n <- out_n %>%
          mutate(!!min_par := rep(get(paste0(min_par, "_f")),
                                  length(times_n))) %>%
          mutate(!!mid_par := rep(get(paste0(mid_par, "_f")),
                                  length(times_n))) %>%
          mutate(!!max_par := rep(get(paste0(max_par, "_f")),
                                  length(times_n)))
        
        data <- rbind(data, out_n[-1, ])
        rm(out_n)
        
        s_asymptotic <- data %>%
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
        
        c_asymptotic <- data %>%
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
    }
  }
  
  # Compute lag time (in years) as difference in times
  lag <- (c_asymptotic - s_asymptotic) / 365
  lag <- lag %>% rename(lag = time)
  
  rm(data)
  
  write.csv(
    cbind(rates, lag),
    file = paste0(rate_index, ".csv"),
    row.names = FALSE
  )
}

# Collect results into "lag_rate_3D.csv"
assign("data", read.csv("1.csv", header = TRUE))
unlink("1.csv")

for (i in 2:dim(values)[1]) {
  tryCatch({
    assign(paste0("data_", i), read.csv(paste0(i, ".csv"), header = TRUE))
    data <- rbind(data, get(paste0("data_", i)))
    unlink(paste0(i, ".csv"))
    rm(list = paste0("data_", i))
  }, error = function(e) {
    message("Error encountered for file ", i, ". Skipping to next iteration.")
  })
}

write.csv(data, file = "lag_rate_3D.csv", row.names = FALSE)

# --------------------------------------------------------------------------------
# End parallelisation
# --------------------------------------------------------------------------------
closeCluster(cl)
mpi.quit()
