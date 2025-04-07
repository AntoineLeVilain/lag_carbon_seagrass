################################################################################
# FILE: figure1_equilibrium.R
#
# DESCRIPTION:
#   This script generates equilibrium solutions of our seagrass–soil model.
#   It includes parallelized routines to scan across stressor ranges, 
#   find equilibria, and write results to CSV files for Figure 1 plotting.
#
# DEPENDENCIES:
#   Requires the following R packages:
#     - readxl
#     - rootSolve
#     - dplyr
#     - rlang
#     - pracma
#     - matrixcalc
#     - data.table
#     - nleqslv
#     - doMPI
################################################################################

# ------------------------- #
#      Load Packages        #
# ------------------------- #

package_list <- c(
  "readxl", "rootSolve", "dplyr", "rlang", 
  "pracma", "matrixcalc", "data.table", "nleqslv"
)

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
  
  # Initialize an empty dataframe to store equilibria
  equilibrium <- data.frame(
    matrix(0, ncol = length(model(0)), nrow = 1)
  )
  colnames(equilibrium) <- names(model(0))
  
  # Generate sequences for each dimension
  lower_upper_list <- lapply(
    1:length(lower_limit),
    function(i) seq(lower_limit[i], upper_limit[i], length.out = iter)
  )
  
  # Generate all possible combinations (including a row of zeros)
  combinations <- expand.grid(lower_upper_list)
  colnames(combinations) <- names(model(0))
  combinations <- rbind(rep(0, length(lower_limit)), combinations)
  
  # Iterate over each combination of initial values
  for (i in 1:dim(combinations)[1]) {
    
    solutions_found <- FALSE  # Flag to track if solution was found
    
    # Attempt to find a root using multiroot
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
    
    # If no solution was found, move to the next combination
    if (!solutions_found) {
      next
    }
    
    # If this is the first successful solution, save it directly
    if (i == 1) {
      equilibrium[1, ] <- solutions$root
    } else {
      # For subsequent solutions, check if we've already found it
      
      # Start by filtering from all known solutions
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
      
      # If solution is not in the list, add it
      if (dim(current)[1] == 0) {
        equilibrium <- rbind(equilibrium, solutions$root)
      }
    }
  }
  
  # Return equilibrium solutions, rounded to two decimals
  return(round(equilibrium, 2))
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
  
  # Generate all possible combinations (including a row of zeros)
  combinations <- expand.grid(lower_upper_list)
  colnames(combinations) <- names(model(0))
  combinations <- rbind(rep(0, length(lower_limit)), combinations)
  
  # Solve using searchZeros
  equilibrium <- searchZeros(
    as.matrix(combinations), fn = model,
    control = list(
      xtol = 10e-12, 
      ftol = 10e-12, 
      allowSingular = TRUE
    )
  )$x
  
  # If searchZeros returns no solution, fall back to multisolve2
  if (is.null(equilibrium)) {
    equilibrium <- multisolve2(model, lower_limit, upper_limit, iter)
  }
  
  # Round, filter non-numeric columns, keep only non-negative solutions
  equilibrium <- round(as.data.frame(equilibrium), 2)
  equilibrium <- equilibrium %>%
    select_if(is.numeric) %>%
    filter(apply(., 1, function(x) all(x >= 0)))
  
  # Remove duplicate solutions
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
  
  S <- rmax *
    ((tmax - temp) / (tmax - topt)) *
    ((temp / topt)^(topt / (tmax - topt))) *
    x[5] / (x[5] + nr) *
    (i0 * exp(-a * pmax * z) /
       (i0 * exp(-a * pmax * z) + ir)) *
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
    (rmax *
       ((tmax - temp) / (tmax - topt)) *
       ((temp / topt)^(topt / (tmax - topt))) *
       x[5] / (x[5] + nr) *
       (i0 * exp(-a * pmax * z) /
          (i0 * exp(-a * pmax * z) + ir))) * x[1] -
    l * x[5]
  
  c(S = S, CA = CA, CS = CS, CB = CB, N = N)
}

# ------------------------------- #
# Bifurcation Analysis (MPI)   #
# ------------------------------- #

# Start parallelization
cl <- startMPIcluster()
registerDoMPI(cl)

# For each stressor, run a bifurcation analysis along its gradient
for (par1 in c("ma", "pmax", "temp")) {
  
  # ----------------------------- #
  #        Parameter Reset        #
  # ----------------------------- #
  
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
  sh    <- 22.7
  
  z     <- 10
  temp  <- 298.85
  hmax  <- 0.15
  pmax  <- 0.5
  ma    <- 0
  i     <- 0.005865335
  
  # ----------------------------- #
  #   multisolve() Parameters     #
  # ----------------------------- #
  
  iv1_max <- 5000
  iv2_max <- 5
  iv3_max <- 5
  iv4_max <- 5
  iv5_max <- 10000
  
  upper_limit <- c(iv1_max, iv2_max, iv3_max, iv4_max, iv5_max)
  lower_limit <- c(10, 0, 0, 0, 0)
  iter        <- 4
  
  # Define a range for the current stressor
  par1_ranges <- list(
    ma   = seq(0, 0.006, length.out = 150),
    pmax = seq(0.5, 5, length.out = 150),
    temp = seq(298.85, 308.15, length.out = 150)
  )
  par1_range <- par1_ranges[[par1]]
  
  # --------------------------------- #
  #   Bifurcation Analysis (Parallel)  #
  # --------------------------------- #
  
  foreach(par1_index = par1_range, .packages = c("rootSolve", "dplyr", "rlang")) %dopar% {
    
    # Assign the new parameter value to the relevant variable
    assign(par1, par1_index)
    
    # Initialize storage for parameter values and states
    assign(paste("liste", par1, sep = "_"), c())
    liste_stability <- c()
    
    # Initialize storage for equilibrium solutions for each state
    for (index in 1:length(model(0))) {
      assign(paste("liste", names(model(0))[index], "seagrass", sep = "_"), c())
      assign(paste("liste", names(model(0))[index], "bistable", sep = "_"), c())
      assign(paste("liste", names(model(0))[index], "bare", sep = "_"), c())
    }
    
    # Get the equilibrium solutions
    equilibrium <- multisolve(model, lower_limit, upper_limit, iter)
    
    # Determine system state from number of equilibrium solutions found
    if (dim(equilibrium)[1] == 3) {
      
      # Bistable scenario: three distinct equilibria
      assign(paste("liste", par1, sep = "_"), 
             c(get(paste("liste", par1, sep = "_")), get(par1)))
      liste_stability <- c(liste_stability, "bistable")
      
      max_save    <- equilibrium[which(equilibrium[, 1] == max(equilibrium[, 1])), ]
      medium_save <- equilibrium[which(
        equilibrium[, 1] != max(equilibrium[, 1]) & 
          equilibrium[, 1] != min(equilibrium[, 1])
      ), ]
      min_save    <- equilibrium[which(equilibrium[, 1] == min(equilibrium[, 1])), ]
      
      for (index2 in 1:length(model(0))) {
        
        if (dim(equilibrium)[2] > 1) {
          assign(
            paste("liste", names(model(0))[index2], "seagrass", sep = "_"),
            c(get(paste("liste", names(model(0))[index2], "seagrass", sep = "_")),
              max_save[1, index2])
          )
          assign(
            paste("liste", names(model(0))[index2], "bistable", sep = "_"),
            c(get(paste("liste", names(model(0))[index2], "bistable", sep = "_")),
              medium_save[1, index2])
          )
          assign(
            paste("liste", names(model(0))[index2], "bare", sep = "_"),
            c(get(paste("liste", names(model(0))[index2], "bare", sep = "_")),
              min_save[1, index2])
          )
        } else {
          assign(
            paste("liste", names(model(0))[index2], "seagrass", sep = "_"),
            c(get(paste("liste", names(model(0))[index2], "seagrass", sep = "_")),
              max_save[1])
          )
          assign(
            paste("liste", names(model(0))[index2], "bistable", sep = "_"),
            c(get(paste("liste", names(model(0))[index2], "bistable", sep = "_")),
              medium_save[1])
          )
          assign(
            paste("liste", names(model(0))[index2], "bare", sep = "_"),
            c(get(paste("liste", names(model(0))[index2], "bare", sep = "_")),
              min_save[1])
          )
        }
      }
      
    } else if (dim(equilibrium)[1] == 2) {
      
      # Seagrass scenario: two distinct equilibria
      assign(paste("liste", par1, sep = "_"), 
             c(get(paste("liste", par1, sep = "_")), get(par1)))
      liste_stability <- c(liste_stability, "seagrass")
      
      max_save <- equilibrium[which(equilibrium[, 1] == max(equilibrium[, 1])), ]
      min_save <- equilibrium[which(equilibrium[, 1] == min(equilibrium[, 1])), ]
      
      for (index3 in 1:length(model(0))) {
        
        if (dim(equilibrium)[2] > 1) {
          assign(
            paste("liste", names(model(0))[index3], "seagrass", sep = "_"),
            c(get(paste("liste", names(model(0))[index3], "seagrass", sep = "_")),
              max_save[1, index3])
          )
          assign(
            paste("liste", names(model(0))[index3], "bistable", sep = "_"),
            c(get(paste("liste", names(model(0))[index3], "bistable", sep = "_")),
              NA)
          )
          assign(
            paste("liste", names(model(0))[index3], "bare", sep = "_"),
            c(get(paste("liste", names(model(0))[index3], "bare", sep = "_")),
              min_save[1, index3])
          )
        } else {
          assign(
            paste("liste", names(model(0))[index3], "seagrass", sep = "_"),
            c(get(paste("liste", names(model(0))[index3], "seagrass", sep = "_")),
              max_save[1])
          )
          assign(
            paste("liste", names(model(0))[index3], "bistable", sep = "_"),
            c(get(paste("liste", names(model(0))[index3], "bistable", sep = "_")),
              NA)
          )
          assign(
            paste("liste", names(model(0))[index3], "bare", sep = "_"),
            c(get(paste("liste", names(model(0))[index3], "bare", sep = "_")),
              min_save[1])
          )
        }
      }
      
    } else {
      
      # Bare scenario: single equilibrium
      assign(paste("liste", par1, sep = "_"), 
             c(get(paste("liste", par1, sep = "_")), get(par1)))
      liste_stability <- c(liste_stability, "bare")
      
      for (index4 in 1:length(model(0))) {
        
        if (dim(equilibrium)[2] > 1) {
          assign(
            paste("liste", names(model(0))[index4], "seagrass", sep = "_"),
            c(get(paste("liste", names(model(0))[index4], "seagrass", sep = "_")),
              NA)
          )
          assign(
            paste("liste", names(model(0))[index4], "bistable", sep = "_"),
            c(get(paste("liste", names(model(0))[index4], "bistable", sep = "_")),
              NA)
          )
          assign(
            paste("liste", names(model(0))[index4], "bare", sep = "_"),
            c(get(paste("liste", names(model(0))[index4], "bare", sep = "_")),
              equilibrium[1, index4])
          )
        } else {
          assign(
            paste("liste", names(model(0))[index4], "seagrass", sep = "_"),
            c(get(paste("liste", names(model(0))[index4], "seagrass", sep = "_")),
              NA)
          )
          assign(
            paste("liste", names(model(0))[index4], "bistable", sep = "_"),
            c(get(paste("liste", names(model(0))[index4], "bistable", sep = "_")),
              NA)
          )
          assign(
            paste("liste", names(model(0))[index4], "bare", sep = "_"),
            c(get(paste("liste", names(model(0))[index4], "bare", sep = "_")),
              equilibrium[1, 1])
          )
        }
      }
    }
    
    # --------------------------- #
    #  Store Data in CSV Format  #
    # --------------------------- #
    
    data <- data.frame(
      par1 = get(paste("liste", par1, sep = "_")), 
      stability = liste_stability
    )
    list_names <- c(par1, "stability")
    
    for (index5 in 1:length(model(0))) {
      data <- cbind(
        data, 
        get(paste("liste", names(model(0))[index5], "seagrass", sep = "_"))
      )
      data <- cbind(
        data, 
        get(paste("liste", names(model(0))[index5], "bistable", sep = "_"))
      )
      data <- cbind(
        data, 
        get(paste("liste", names(model(0))[index5], "bare", sep = "_"))
      )
      
      list_names <- c(
        list_names,
        paste(names(model(0))[index5], "seagrass", sep = "_"),
        paste(names(model(0))[index5], "bistable", sep = "_"),
        paste(names(model(0))[index5], "bare", sep = "_")
      )
    }
    
    colnames(data) <- list_names
    
    # Write a CSV for each parameter value
    write.csv(
      data,
      file = paste0(par1, par1_index, ".csv"),
      row.names = FALSE
    )
  }
  
  # ----------------------- #
  #  Combine & Clean Up     #
  # ----------------------- #
  
  assign('data', read.csv(paste0(par1, par1_range[1], ".csv"), header = TRUE))
  unlink(paste0(par1, par1_range[1], ".csv"))
  
  for (par1_index in par1_range[-1]) {
    assign(
      paste0(par1, par1_index),
      read.csv(paste0(par1, par1_index, ".csv"), header = TRUE)
    )
    data <- rbindlist(list(data, get(paste0(par1, par1_index))))
    unlink(paste0(par1, par1_index, ".csv"))
    rm(list = (paste0(par1, par1_index)))
  }
  
  # Write final combined data for current stressor
  write.csv(data, paste0(par1, "_equilibrium.csv"))
  rm(data)
}

# End parallelization
closeCluster(cl)
mpi.quit()
