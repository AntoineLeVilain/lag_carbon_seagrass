################################################################################
# FILE: figure4_pattern_lag_2D.R
#
# DESCRIPTION:
#   This script explores pairwise parameter change scenarios for three stressors
#   (pmax, temp, ma), potentially changing them at different rates. It records
#   lag and patterns of S and CB changes (e.g., “++” for both S and CB increasing,
#   “-+” for S decreasing and CB increasing, etc.). Results are saved as separate
#   CSV files (e.g., "df1_pattern_lag.csv") for each stressor combination
#   scenario.
#
# DEPENDENCIES:
#   Requires the following R packages:
#     - deSolve
#     - dplyr
#     - nleqslv
#     - rootSolve
#     - doMPI
################################################################################

# ------------------------------------------------------------------------------
# Load Required Packages
# ------------------------------------------------------------------------------
package_list <- c("deSolve", "dplyr", "nleqslv", "rootSolve")
for (pkg in package_list) {
  library(pkg, character.only = TRUE)
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

pmax_0 <- 0.5
pmax_f <- 4.14

temp_0 <- 298.85
temp_f <- 308

ma_0 <- 0
ma_f <- 0.00330

slow       <- 1000 * 365  # "slow" rate
fast       <- 0           # "fast" rate
very_high  <- 1000000     # extremely high rate
number_values <- 20

# Generate the values for the combinations
param_1 <- seq(slow, fast, length.out = number_values)
param_2 <- seq(slow, fast, length.out = number_values)

df_1 <- expand.grid(param_1, param_2)
df_1 <- cbind(df_1, rep(very_high, ncol(df_1)))
df_1 <- round(df_1[!(df_1[[1]] == 0 & df_1[[2]] == 0), ])

colnames(df_1) <- c("pmax", "temp", "ma")
rownames(df_1) <- NULL

df_2 <- df_1
colnames(df_2) <- c("pmax", "ma", "temp")

df_3 <- df_1
colnames(df_3) <- c("temp", "ma", "pmax")

precision    <- 0.05
precision_c  <- 10
precision_s  <- 10

dimo <- seq_len(nrow(df_1))

# ------------------------------------------------------------------------------
# Start Parallelization
# ------------------------------------------------------------------------------
cl <- startMPIcluster()
registerDoMPI(cl)

# ------------------------------------------------------------------------------
# Analyze Rate Combinations in Parallel
# ------------------------------------------------------------------------------
for (df_index in 1:3) {
  
  values <- get(paste0("df_", df_index))
  
  foreach(rate_index = dimo,
          .packages = c("rootSolve", "dplyr", "rlang", "deSolve", "nleqslv")) %dopar% {
            
            # Reset parameter variables
            ma   <- ma_0
            pmax <- pmax_0
            temp <- temp_0
            
            rates <- values[rate_index, ]
            
            # Transpose the row (like pivot) to sort columns by the first row
            rates <- as.data.frame(t(rates))
            # Sort the transposed dataframe by ascending numeric values
            rates <- rates[order(rates[, 1], decreasing = FALSE), , drop = FALSE]
            # Transpose back
            rates <- as.data.frame(t(rates))
            
            rate_pmax <- as.numeric(rates %>% dplyr::select(pmax))
            rate_temp <- as.numeric(rates %>% dplyr::select(temp))
            rate_ma   <- as.numeric(rates %>% dplyr::select(ma))
            
            min_rate <- min(rates)
            min_par  <- colnames(rates[1])
            
            mid_rate <- as.numeric(rates[2])
            mid_par  <- colnames(rates[2])
            
            max_rate <- as.numeric(rates[3])
            max_par  <- colnames(rates[3])
            
            # If the "lowest" rate parameter is nonzero, set that param to final
            if (rates[1] != 0) {
              assign(colnames(rates)[1], get(paste0(colnames(rates)[1], "_f")))
            }
            
            # If the "middle" rate param is nonzero, set that param to final
            if (rates[2] != 0) {
              assign(colnames(rates)[2], get(paste0(colnames(rates)[2], "_f")))
            }
            
            # The "highest" param we keep at initial
            assign(colnames(rates)[3], get(paste0(colnames(rates)[3], "_0")))
            
            # Compute new final equilibrium with these param setups
            final_eq <- multisolve(model, lower_limit, upper_limit, iter) %>%
              dplyr::filter(S == max(S))
            
            # -------------------------------
            # If min_rate == 0, only mid param changes from 0..mid_rate
            # -------------------------------
            if (min_rate == 0) {
              
              # seagrass() => only mid param changes from init..final
              seagrass <- function(t, state, parameters) {
                with(as.list(c(state, parameters)), {
                  
                  assign(min_par, get(paste0(min_par, "_0")))
                  assign(mid_par, get(paste0(mid_par, "_0")) +
                           (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0"))) *
                           t * 1 / mid_rate
                  )
                  assign(max_par, get(paste0(max_par, "_0")))
                  
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
                      phida * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                        (hmax * sh / (sh + S)) / (hmax * sh / (sh + S) + hea) +
                        1 / (1 + exp(-fb1 * (CA + CS - fb2)))
                    )
                  
                  dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
                    CS * (
                      phids * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                        (hmax * sh / (sh + S)) / (hmax * sh / (sh + S) + hes) +
                        1 / (1 + exp(-fb1 * (CA + CS - fb2)))
                    )
                  
                  dCB <- CA * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) +
                    CS * 1 / (1 + exp(-fb1 * (CA + CS - fb2))) -
                    (1 - ma * da) * CB * phidb *
                    1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) -
                    ma * da * CB * (phida + phids) / 2 *
                    1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp))
                  
                  dN <- i + gamma * (
                    CA * phida * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                      CS * phids * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
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
              
              times_1 <- seq(0, mid_rate, 1)
              start_1_vector <- as.numeric(as.numeric(initial_eq))
              names(start_1_vector) <- names(initial_eq)
              
              data <- as.data.frame(lsoda(
                y = start_1_vector,
                times = times_1,
                func = seagrass,
                parms = params
              ))
              
              data <- data %>%
                dplyr::mutate(
                  !!min_par := get(paste0(min_par, "_0")),
                  !!mid_par := get(paste0(mid_par, "_0")) +
                    (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0"))) *
                    time / mid_rate,
                  !!max_par := get(paste0(max_par, "_0"))
                )
              
              s_asymptotic <- data %>%
                dplyr::mutate(across(everything(), ~ round(., 2))) %>%
                dplyr::filter(
                  round(final_eq[, 1] * (1 - precision), 0) <= S &
                    S <= round(final_eq[, 1] * (1 + precision), 0)
                )
              
              if (nrow(s_asymptotic) > 0) {
                s_asymptotic <- s_asymptotic %>%
                  dplyr::filter(time == min(time)) %>%
                  dplyr::select(time)
              }
              
              c_asymptotic <- data %>%
                dplyr::mutate(across(everything(), ~ round(., 2))) %>%
                dplyr::filter(
                  round(final_eq[, 4] * (1 - precision), 0) <= CB &
                    CB <= round(final_eq[, 4] * (1 + precision), 0)
                )
              
              if (nrow(c_asymptotic) > 0) {
                c_asymptotic <- c_asymptotic %>%
                  dplyr::filter(time == min(time)) %>%
                  dplyr::select(time)
              }
              
              index    <- 1
              max_time <- max(data$time)
              
              # If not converged, hold mid_par at final
              while (nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0) {
                
                seagrass <- function(t, state, parameters) {
                  with(as.list(c(state, parameters)), {
                    
                    assign(min_par, get(paste0(min_par, "_0")))
                    assign(mid_par, get(paste0(mid_par, "_f")))
                    assign(max_par, get(paste0(max_par, "_0")))
                    
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
                        phida * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                          (hmax * sh / (sh + S)) /
                          (hmax * sh / (sh + S) + hea) +
                          1 / (1 + exp(-fb1 * (CA + CS - fb2)))
                      )
                    
                    dCS <- beta * S * (ms + mh * hmax * sh / (sh + S)) -
                      CS * (
                        phids * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
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
                    
                    dN <- i + gamma * (
                      CA * phida * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
                        CS * phids * 1 / (exp(-ea / (8.314 * td))) * exp(-ea / (8.314 * temp)) +
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
                
                times_n <- seq(max_time * index, max_time * (index + 1), 1)
                
                start_n_vector <- as.numeric(data %>%
                                               dplyr::filter(time == max(time)) %>%
                                               dplyr::select(S, CA, CS, CB, N))
                
                names(start_n_vector) <- names(initial_eq)
                
                out_n <- as.data.frame(lsoda(
                  y = start_n_vector,
                  times = times_n,
                  func = seagrass,
                  parms = params
                ))
                
                out_n <- out_n %>%
                  dplyr::mutate(
                    !!min_par := rep(get(paste0(min_par, "_0")), length(times_n)),
                    !!mid_par := rep(get(paste0(mid_par, "_f")), length(times_n)),
                    !!max_par := get(paste0(max_par, "_0"))
                  )
                
                data <- rbind(data, out_n[-1, ])
                
                rm(out_n)
                
                s_asymptotic <- data %>%
                  dplyr::mutate(across(everything(), ~ round(., 2))) %>%
                  dplyr::filter(
                    round(final_eq[, 1] * (1 - precision), 0) <= S &
                      S <= round(final_eq[, 1] * (1 + precision), 0)
                  )
                
                if (nrow(s_asymptotic) > 0) {
                  s_asymptotic <- s_asymptotic %>%
                    dplyr::filter(time == min(time)) %>%
                    dplyr::select(time)
                }
                
                c_asymptotic <- data %>%
                  dplyr::mutate(across(everything(), ~ round(., 2))) %>%
                  dplyr::filter(
                    round(final_eq[, 4] * (1 - precision), 0) <= CB &
                      CB <= round(final_eq[, 4] * (1 + precision), 0)
                  )
                
                if (nrow(c_asymptotic) > 0) {
                  c_asymptotic <- c_asymptotic %>%
                    dplyr::filter(time == min(time)) %>%
                    dplyr::select(time)
                }
                
                index <- index + 1
              }
              
            } else if (min_rate == mid_rate) {
              # ------------------------------------------------------------------
              # If min_rate == mid_rate, both parameters ramp 0..final together
              # (Code logic identical to above pattern)
              # ------------------------------------------------------------------
              
              seagrass <- function(t, state, parameters){
                
                with(as.list(c(state, parameters)),{
                  
                  assign(min_par, get(paste0(min_par, "_0")) + (get(paste0(min_par, "_f")) - get(paste0(min_par, "_0")))*t*1/min_rate)
                  assign(mid_par, get(paste0(mid_par, "_0")) + (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0")))*t*1/mid_rate)
                  assign(max_par, get(paste0(max_par, "_0")))
                  
                  dS <- rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z)/(i0*exp(-a*pmax*z) + ir))*S - S*(ms + ma + mh*hmax*sh/(sh + S))
                  
                  dCA <- alpha*h*pmax*(1 - (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hp)) - CA*(phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hea) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                  
                  dCS <- beta*S*(ms + mh*hmax*sh/(sh + S)) - CS*(phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hes) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                  
                  dCB <- CA*1/(1+exp(-fb1*(CA + CS - fb2))) + CS*1/(1+exp(-fb1*(CA + CS - fb2))) - (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) - ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))
                  
                  dN <- i + gamma*(CA*phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + CS*phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))) - delta*(rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z))/(i0*exp(-a*pmax*z) + ir))*S-l*N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                  
                })
                
              }
              
              times_1 <- seq(0, mid_rate, 1)
              
              start_1_vector <- as.numeric(as.numeric(initial_eq))
              names(start_1_vector) <- names(initial_eq)
              
              data <- as.data.frame(lsoda(y = start_1_vector, times = times_1, func = seagrass, parms = params))
              
              data <- data %>%
                mutate(
                  !!min_par := get(paste0(min_par, "_0")) + (get(paste0(min_par, "_f")) - get(paste0(min_par, "_0"))) * time / min_rate,
                  !!mid_par := get(paste0(mid_par, "_0")) + (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0"))) * time / mid_rate,
                  !!max_par := get(paste0(max_par, "_0"))  # Repeat the initial value for max_par
                )
              
              s_asymptotic <- data %>%
                mutate(across(everything(), ~round(., digits = 2))) %>%
                filter(round(final_eq[, 1]*(1 - precision), 0) <= S & S <= round(final_eq[, 1]*(1 + precision), 0))
              
              if(nrow(s_asymptotic) > 0){
                
                s_asymptotic <- s_asymptotic %>%
                  filter(time == min(time)) %>%
                  select(time)
                
              }
              
              c_asymptotic <- data %>%
                mutate(across(everything(), ~round(., digits = 2))) %>%
                filter(round(final_eq[, 4]*(1 - precision), 0) <= CB & CB <= round(final_eq[, 4]*(1 + precision), 0))
              
              if (nrow(c_asymptotic) > 0){
                
                c_asymptotic <- c_asymptotic %>%
                  filter(time == min(time)) %>%
                  select(time)
                
              }
              
              index <- 1
              max_time <- max(data$time)
              
              while(nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0){
                
                seagrass <- function(t, state, parameters){
                  
                  with(as.list(c(state, parameters)),{
                    
                    assign(min_par, get(paste0(min_par, "_f")))
                    assign(mid_par, get(paste0(mid_par, "_f")))
                    assign(max_par, get(paste0(max_par, "_0")))
                    
                    dS <- rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z)/(i0*exp(-a*pmax*z) + ir))*S - S*(ms + ma + mh*hmax*sh/(sh + S))
                    
                    dCA <- alpha*h*pmax*(1 - (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hp)) - CA*(phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hea) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                    
                    dCS <- beta*S*(ms + mh*hmax*sh/(sh + S)) - CS*(phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hes) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                    
                    dCB <- CA*1/(1+exp(-fb1*(CA + CS - fb2))) + CS*1/(1+exp(-fb1*(CA + CS - fb2))) - (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) - ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))
                    
                    dN <- i + gamma*(CA*phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + CS*phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))) - delta*(rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z))/(i0*exp(-a*pmax*z) + ir))*S-l*N
                    
                    list(c(dS, dCA, dCS, dCB, dN))
                    
                  })
                  
                }
                
                times_n <- seq(max_time*index, max_time*(index + 1), 1)
                
                start_n_vector <- as.numeric(data %>% 
                                               filter(time == max(time)) %>%
                                               select(S, CA, CS, CB, N))
                
                names(start_n_vector) <- names(initial_eq)
                
                out_n <- as.data.frame(lsoda(y = start_n_vector, times = times_n, func = seagrass, parms = params))
                
                out_n <- out_n %>%
                  mutate(!!min_par := rep(get(paste0(min_par, "_f")), length(times_n))) %>%
                  mutate(!!mid_par := rep(get(paste0(mid_par, "_f")), length(times_n))) %>%
                  mutate(!!max_par := get(paste0(max_par, "_0")))
                
                data <- rbind(data, out_n[-1, ])
                
                rm(out_n)
                
                s_asymptotic <- data %>%
                  mutate(across(everything(), ~round(., digits = 2))) %>%
                  filter(round(final_eq[, 1]*(1 - precision), 0) <= S & S <= round(final_eq[, 1]*(1 + precision), 0))
                
                if(nrow(s_asymptotic) > 0){
                  
                  s_asymptotic <- s_asymptotic %>%
                    filter(time == min(time)) %>%
                    select(time)
                  
                }
                
                c_asymptotic <- data %>%
                  mutate(across(everything(), ~round(., digits = 2))) %>%
                  filter(round(final_eq[, 4]*(1 - precision), 0) <= CB & CB <= round(final_eq[, 4]*(1 + precision), 0))
                
                if (nrow(c_asymptotic) > 0){
                  
                  c_asymptotic <- c_asymptotic %>%
                    filter(time == min(time)) %>%
                    select(time)
                  
                }
                
                index <- index + 1
                
              }
              
            } else {
              # ------------------------------------------------------------------
              # Otherwise, the min param changes at a certain rate, then the mid
              # param changes, while the max param remains at initial. 
              # *** Again, identical logic repeated ***
              # ------------------------------------------------------------------
              
              seagrass <- function(t, state, parameters){
                
                with(as.list(c(state, parameters)),{
                  
                  assign(min_par, get(paste0(min_par, "_0")) + (get(paste0(min_par, "_f")) - get(paste0(min_par, "_0")))*t*1/min_rate)
                  assign(mid_par, get(paste0(mid_par, "_0")) + (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0")))*t*1/mid_rate)
                  assign(max_par, get(paste0(max_par, "_0")))
                  
                  dS <- rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z)/(i0*exp(-a*pmax*z) + ir))*S - S*(ms + ma + mh*hmax*sh/(sh + S))
                  
                  dCA <- alpha*h*pmax*(1 - (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hp)) - CA*(phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hea) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                  
                  dCS <- beta*S*(ms + mh*hmax*sh/(sh + S)) - CS*(phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hes) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                  
                  dCB <- CA*1/(1+exp(-fb1*(CA + CS - fb2))) + CS*1/(1+exp(-fb1*(CA + CS - fb2))) - (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) - ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))
                  
                  dN <- i + gamma*(CA*phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + CS*phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))) - delta*(rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z))/(i0*exp(-a*pmax*z) + ir))*S-l*N
                  
                  list(c(dS, dCA, dCS, dCB, dN))
                  
                })
                
              }
              
              times_1 <- seq(0, min_rate, 1)
              
              start_1_vector <- as.numeric(as.numeric(initial_eq))
              names(start_1_vector) <- names(initial_eq)
              
              data <- as.data.frame(lsoda(y = start_1_vector, times = times_1, func = seagrass, parms = params))
              
              data <- data %>%
                mutate(
                  !!min_par := get(paste0(min_par, "_0")) + (get(paste0(min_par, "_f")) - get(paste0(min_par, "_0"))) * time / min_rate,
                  !!mid_par := get(paste0(mid_par, "_0")) + (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0"))) * time / mid_rate,
                  !!max_par := get(paste0(max_par, "_0"))  # Repeat the initial value for max_par
                )
              
              s_asymptotic <- data %>%
                mutate(across(everything(), ~round(., digits = 2))) %>%
                filter(round(final_eq[, 1]*(1 - precision), 0) <= S & S <= round(final_eq[, 1]*(1 + precision), 0))
              
              if(nrow(s_asymptotic) > 0){
                
                s_asymptotic <- s_asymptotic %>%
                  filter(time == min(time)) %>%
                  select(time)
                
              }
              
              c_asymptotic <- data %>%
                mutate(across(everything(), ~round(., digits = 2))) %>%
                filter(round(final_eq[, 4]*(1 - precision), 0) <= CB & CB <= round(final_eq[, 4]*(1 + precision), 0))
              
              if (nrow(c_asymptotic) > 0){
                
                c_asymptotic <- c_asymptotic %>%
                  filter(time == min(time)) %>%
                  select(time)
                
              }
              
              if(nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0){
                
                seagrass <- function(t, state, parameters){
                  
                  with(as.list(c(state, parameters)),{
                    
                    assign(min_par, get(paste0(min_par, "_f")))
                    assign(mid_par, get(paste0(mid_par, "_0")) + (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0")))*t*1/mid_rate)
                    assign(max_par, get(paste0(max_par, "_0")))
                    
                    dS <- rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z)/(i0*exp(-a*pmax*z) + ir))*S - S*(ms + ma + mh*hmax*sh/(sh + S))
                    
                    dCA <- alpha*h*pmax*(1 - (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hp)) - CA*(phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hea) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                    
                    dCS <- beta*S*(ms + mh*hmax*sh/(sh + S)) - CS*(phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hes) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                    
                    dCB <- CA*1/(1+exp(-fb1*(CA + CS - fb2))) + CS*1/(1+exp(-fb1*(CA + CS - fb2))) - (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) - ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))
                    
                    dN <- i + gamma*(CA*phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + CS*phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))) - delta*(rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z))/(i0*exp(-a*pmax*z) + ir))*S-l*N
                    
                    list(c(dS, dCA, dCS, dCB, dN))
                    
                  })
                  
                }
                
                times_2 <- seq(min_rate, mid_rate, 1)
                
                start_2_vector <- as.numeric(data %>% filter(time == max(time)) %>% select(S, CA, CS, CB, N))
                
                names(start_2_vector) <- names(initial_eq)
                
                out_2 <- as.data.frame(lsoda(y = start_2_vector, times = times_2, func = seagrass, parms = params))
                
                out_2 <- out_2 %>%
                  mutate(
                    !!min_par := get(paste0(min_par, "_f")),
                    !!mid_par := get(paste0(mid_par, "_0")) + (get(paste0(mid_par, "_f")) - get(paste0(mid_par, "_0"))) * time / mid_rate,
                    !!max_par := get(paste0(max_par, "_0"))
                  )
                
                data <- rbind(data, out_2[-1, ])
                
                rm(out_2)
                
                s_asymptotic <- data %>%
                  mutate(across(everything(), ~round(., digits = 2))) %>%
                  filter(round(final_eq[, 1]*(1 - precision), 0) <= S & S <= round(final_eq[, 1]*(1 + precision), 0))
                
                if(nrow(s_asymptotic) > 0){
                  
                  s_asymptotic <- s_asymptotic %>%
                    filter(time == min(time)) %>%
                    select(time)
                  
                }
                
                c_asymptotic <- data %>%
                  mutate(across(everything(), ~round(., digits = 2))) %>%
                  filter(round(final_eq[, 4]*(1 - precision), 0) <= CB & CB <= round(final_eq[, 4]*(1 + precision), 0))
                
                if (nrow(c_asymptotic) > 0){
                  
                  c_asymptotic <- c_asymptotic %>%
                    filter(time == min(time)) %>%
                    select(time)
                  
                }
                
                index <- 1
                max_time <- max(data$time)
                
                while(nrow(s_asymptotic) == 0 | nrow(c_asymptotic) == 0){
                  
                  seagrass <- function(t, state, parameters){
                    
                    with(as.list(c(state, parameters)),{
                      
                      assign(min_par, get(paste0(min_par, "_f")))
                      assign(mid_par, get(paste0(mid_par, "_f")))
                      assign(max_par, get(paste0(max_par, "_0")))
                      
                      dS <- rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z)/(i0*exp(-a*pmax*z) + ir))*S - S*(ms + ma + mh*hmax*sh/(sh + S))
                      
                      dCA <- alpha*h*pmax*(1 - (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hp)) - CA*(phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hea) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                      
                      dCS <- beta*S*(ms + mh*hmax*sh/(sh + S)) - CS*(phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (hmax*sh/(sh + S))/(hmax*sh/(sh + S) + hes) + 1/(1+exp(-fb1*(CA + CS - fb2))))
                      
                      dCB <- CA*1/(1+exp(-fb1*(CA + CS - fb2))) + CS*1/(1+exp(-fb1*(CA + CS - fb2))) - (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) - ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))
                      
                      dN <- i + gamma*(CA*phida*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + CS*phids*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + (1-ma*da)*CB*phidb*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp)) + ma*da*CB*(phida+phids)/2*1/(exp(-ea/(8.314*td)))*exp(-ea/(8.314*temp))) - delta*(rmax*((tmax - temp)/(tmax - topt))*((temp/topt)^(topt/(tmax - topt)))*N/(N + nr)*(i0*exp(-a*pmax*z))/(i0*exp(-a*pmax*z) + ir))*S-l*N
                      
                      list(c(dS, dCA, dCS, dCB, dN))
                      
                    })
                    
                  }
                  
                  times_n <- seq(max_time*index, max_time*(index + 1), 1)
                  
                  start_n_vector <- as.numeric(data %>% 
                                                 filter(time == max(time)) %>%
                                                 select(S, CA, CS, CB, N))
                  
                  names(start_n_vector) <- names(initial_eq)
                  
                  out_n <- as.data.frame(lsoda(y = start_n_vector, times = times_n, func = seagrass, parms = params))
                  
                  out_n <- out_n %>%
                    mutate(!!min_par := rep(get(paste0(min_par, "_f")), length(times_n))) %>%
                    mutate(!!mid_par := rep(get(paste0(mid_par, "_f")), length(times_n))) %>%
                    mutate(!!max_par := get(paste0(max_par, "_0")))
                  
                  data <- rbind(data, out_n[-1, ])
                  
                  rm(out_n)
                  
                  s_asymptotic <- data %>%
                    mutate(across(everything(), ~round(., digits = 2))) %>%
                    filter(round(final_eq[, 1]*(1 - precision), 0) <= S & S <= round(final_eq[, 1]*(1 + precision), 0))
                  
                  if(nrow(s_asymptotic) > 0){
                    
                    s_asymptotic <- s_asymptotic %>%
                      filter(time == min(time)) %>%
                      select(time)
                    
                  }
                  
                  c_asymptotic <- data %>%
                    mutate(across(everything(), ~round(., digits = 2))) %>%
                    filter(round(final_eq[, 4]*(1 - precision), 0) <= CB & CB <= round(final_eq[, 4]*(1 + precision), 0))
                  
                  if (nrow(c_asymptotic) > 0){
                    
                    c_asymptotic <- c_asymptotic %>%
                      filter(time == min(time)) %>%
                      select(time)
                    
                  }
                  
                  index <- index + 1
                  
                }
                
              }
              
            }
            
            # Once we have data, compute the final lag, pattern, etc.
            lag <- (c_asymptotic - s_asymptotic) / 365
            lag <- lag %>% dplyr::rename(lag = time)
            
            low_time <- 0
            high_time <- dim(data)[1]
            bif_time <- round((low_time + high_time)/2, 0)
            
            assign("pmax", as.numeric(data %>% filter(time == bif_time) %>% select(pmax)))
            assign("ma", as.numeric(data %>% filter(time == bif_time) %>% select(ma)))
            assign("temp", as.numeric(data %>% filter(time == bif_time) %>% select(temp)))
            
            bifurcation_1 <- dim(multisolve(model, lower_limit, upper_limit, iter))[1]
            
            assign("pmax", as.numeric(data %>% filter(time == bif_time + 1) %>% select(pmax)))
            assign("ma", as.numeric(data %>% filter(time == bif_time + 1) %>% select(ma)))
            assign("temp", as.numeric(data %>% filter(time == bif_time + 1) %>% select(temp)))
            
            bifurcation_2 <- dim(multisolve(model, lower_limit, upper_limit, iter))[1]
            
            iter_count <- 0
            
            while((bifurcation_2 != 1 | bifurcation_1 < 2) && (iter_count < 50)){
              
              iter_count <- iter_count + 1
              
              if(bifurcation_2 == 1 & bifurcation_1 == 1){
                
                high_time <- bif_time
                bif_time <- round((high_time + low_time)/2, 0)
                
              } else if(bifurcation_2 > 1) {
                
                low_time <- bif_time
                bif_time <- round((high_time + low_time)/2, 0)
                
              }
              
              assign("pmax", as.numeric(data %>% filter(time == bif_time) %>% select(pmax)))
              assign("ma", as.numeric(data %>% filter(time == bif_time) %>% select(ma)))
              assign("temp", as.numeric(data %>% filter(time == bif_time) %>% select(temp)))
              
              bifurcation_1 <- dim(multisolve(model, lower_limit, upper_limit, iter))[1]
              
              assign("pmax", as.numeric(data %>% filter(time == bif_time + 1) %>% select(pmax)))
              assign("ma", as.numeric(data %>% filter(time == bif_time + 1) %>% select(ma)))
              assign("temp", as.numeric(data %>% filter(time == bif_time + 1) %>% select(temp)))
              
              bifurcation_2 <- dim(multisolve(model, lower_limit, upper_limit, iter))[1]
              
              print(bif_time)
              
            }
            
            data <- data %>% filter(time <= bif_time)
            
            max_CB <- data %>%
              filter(CB == max(CB)) %>%
              select(time)
            
            max_S <- data %>%
              filter(S == max(S)) %>%
              select(time)
            
            mismatch <- as.numeric(max_CB - max_S)
            
            s_pattern <- (as.numeric(data %>% filter(time == max(time)) %>% select(S)) - as.numeric(data %>% filter(time == 0) %>% select(S)))/as.numeric(data %>% filter(time == 0) %>% select(S))*100
            
            c_pattern <- (as.numeric(data %>% filter(time == max(time)) %>% select(CB)) - as.numeric(data %>% filter(time == 0) %>% select(CB)))/as.numeric(data %>% filter(time == 0) %>% select(CB))*100
            
            if (identical(colnames(rates)[1:2], c("pmax", "temp")) | identical(colnames(rates)[1:2], c("temp", "pmax"))){
              
              if((s_pattern > precision_s) & (c_pattern < -precision_c)){
                
                pattern <- "+-"
                
              } else if ((s_pattern > precision_s) & (c_pattern < precision_c) & (c_pattern > -precision_c)){
                
                pattern <- "+="
                
              } else if ((s_pattern > precision_s) & (c_pattern > precision_c)){
                
                pattern <- "++"
                
              } else if ((s_pattern < precision_s) & (s_pattern > -precision_s) & (c_pattern < -precision_c)){
                
                pattern <- "=-"
                
              } else if ((s_pattern < precision_s) & (s_pattern > -precision_s) & (c_pattern < precision_c) & (c_pattern > -precision_c)){
                
                pattern <- "=="
                
              } else if ((s_pattern < precision_s) & (s_pattern > -precision_s) & (c_pattern > precision_c)){
                
                pattern <- "=+"
                
              } else if ((s_pattern < -precision_s) & (c_pattern < -precision_c)){
                
                pattern <- "--"
                
              } else if ((s_pattern < -precision_s) & (c_pattern < precision_c) & (c_pattern > -precision_c)){
                
                pattern <- "-="
                
              } else if ((s_pattern < -precision_s) & (c_pattern > precision_c)){
                
                pattern <- "-+"
                
              }
              
            } else {
              
              if((s_pattern >= 0) & (c_pattern < 0)){
                
                pattern <- "+-"
                
              } else if ((s_pattern >= 0) & (c_pattern >= 0)){
                
                pattern <- "++"
                
              } else if ((s_pattern < 0) & (c_pattern < 0)){
                
                pattern <- "--"
                
              } else if ((s_pattern < 0) & (c_pattern >= 0)){
                
                pattern <- "-+"
                
              }
              
            }
            
            rm(data)
            
            rates[3] <- 0
            
            write.csv(
              cbind(rates / 365, data.frame(pattern), lag),
              file = paste0("df", df_index, "_", rate_index, ".csv"),
              row.names = FALSE
            )
            
          }
}

# ------------------------------------------------------------------------------
# Combine & Clean Up
# ------------------------------------------------------------------------------
for (df_index in 1:3) {
  
  assign("data", read.csv(paste0("df", df_index, "_", "1.csv"), header = TRUE))
  unlink(paste0("df", df_index, "_", "1.csv"))
  
  for (i in 2:length(dimo)) {
    assign(paste0("data_", i), read.csv(paste0("df", df_index, "_", i, ".csv"), header = TRUE))
    
    data <- rbind(data, get(paste0("data_", i)))
    
    unlink(paste0("df", df_index, "_", i, ".csv"))
    rm(list = paste0("data_", i))
  }
  
  write.csv(
    data,
    file = paste0("df", df_index, "_", "pattern_lag.csv"),
    row.names = FALSE
  )
}

# ------------------------------------------------------------------------------
# End Parallelization
# ------------------------------------------------------------------------------
closeCluster(cl)
mpi.quit()
