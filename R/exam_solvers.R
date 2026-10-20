#' THE ULTIMATE EXAM SOLVERS (TOPICS 1, 2, 3)
#' @name exam_solvers
NULL

# ==============================================================================
# TOPIC 1: DESCRIPTIVE STATISTICS (Handles both raw data and frequency tables)
# ==============================================================================
#' @export
solve_topic1 <- function(x, n = NULL) {
  if(!is.null(n)) x <- rep(x, n) 
  
  cat("\n========== TOPIC 1: DESCRIPTIVE STATISTICS ==========\n")
  cat("1. CENTRAL TENDENCY\n")
  cat("   Mean:             ", mean(x), "\n")
  cat("   Median:           ", median(x), "\n")
  freq_table <- table(x)
  cat("   Mode:             ", as.numeric(names(freq_table)[which.max(freq_table)]), "\n\n")
  
  cat("2. DISPERSION\n")
  pop_var <- sum((x - mean(x))^2) / length(x)
  pop_sd <- sqrt(pop_var)
  cat("   Pop. Variance:    ", pop_var, "\n")
  cat("   Pop. Std Dev:     ", pop_sd, "\n")
  cat("   Coef. Variation:  ", (pop_sd / abs(mean(x))) * 100, "% \n\n")
  
  cat("3. SHAPE (Fisher)\n")
  m3 <- sum((x - mean(x))^3) / length(x)
  m4 <- sum((x - mean(x))^4) / length(x)
  cat("   Skewness (g1):    ", m3 / (pop_sd^3), "\n")
  cat("   Kurtosis (g2):    ", (m4 / (pop_sd^4)) - 3, "\n\n")
  
  cat("4. QUARTILES & OUTLIERS\n")
  q <- quantile(x, probs = c(0.25, 0.5, 0.75), type = 2)
  cat("   Q1 =", q[1], "| Q2 =", q[2], "| Q3 =", q[3], "\n")
  iqr <- q[3] - q[1]
  cat("   Moderate Outliers Limits: [", q[1] - 1.5*iqr, ",", q[3] + 1.5*iqr, "]\n")
  cat("   Extreme Outliers Limits:  [", q[1] - 3*iqr, ",", q[3] + 3*iqr, "]\n")
  cat("=====================================================\n")
}

# ==============================================================================
# TOPIC 2: STATISTICAL MODELLING (Linear, Exponential and Hyperbolic)
# ==============================================================================
#' @export
solve_topic2 <- function(x, y) {
  n <- length(x)
  cov_xy <- sum((x - mean(x)) * (y - mean(y))) / n
  var_x <- sum((x - mean(x))^2) / n
  var_y <- sum((y - mean(y))^2) / n
  r <- cov_xy / sqrt(var_x * var_y)
  
  cat("\n========== TOPIC 2: STATISTICAL MODELLING ==========\n")
  cat("1. BASIC STATS\n")
  cat("   Covariance:       ", cov_xy, "\n")
  cat("   Correlation (r):  ", r, "\n")
  cat("   Var(X):", var_x, " | Var(Y):", var_y, "\n\n")
  
  cat("2. LINEAR FIT ( y = a + bx )\n")
  b_lin <- cov_xy / var_x
  a_lin <- mean(y) - b_lin * mean(x)
  r2_lin <- r^2
  cat("   Equation: y =", a_lin, "+", b_lin, "* x\n")
  cat("   R-squared: ", r2_lin, "\n\n")
  
  cat("3. EXPONENTIAL FIT ( y = a * e^(bx) )\n")
  Y_ln <- log(y)
  cov_x_Yln <- sum((x - mean(x)) * (Y_ln - mean(Y_ln))) / n
  var_Yln <- sum((Y_ln - mean(Y_ln))^2) / n
  b_exp <- cov_x_Yln / var_x
  A_exp <- mean(Y_ln) - b_exp * mean(x)
  a_exp <- exp(A_exp)
  r2_exp <- (cov_x_Yln / sqrt(var_x * var_Yln))^2
  cat("   Equation: y =", a_exp, "* e^(", b_exp, "* x )\n")
  cat("   R-squared: ", r2_exp, "\n\n")
  
  cat("4. HYPERBOLIC FIT ( y = 1 / (a + bx) )\n")
  Y_inv <- 1/y
  cov_x_Yinv <- sum((x - mean(x)) * (Y_inv - mean(Y_inv))) / n
  var_Yinv <- sum((Y_inv - mean(Y_inv))^2) / n
  b_hyp <- cov_x_Yinv / var_x
  a_hyp <- mean(Y_inv) - b_hyp * mean(x)
  r2_hyp <- (cov_x_Yinv / sqrt(var_x * var_Yinv))^2
  cat("   Equation: y = 1 / (", a_hyp, "+", b_hyp, "* x )\n")
  cat("   R-squared: ", r2_hyp, "\n")
  cat("====================================================\n")
}

# ==============================================================================
# TOPIC 3: TIME SERIES (Trend, Seasonality, Deseasonalization, Deflation)
# ==============================================================================
#' @export
solve_topic3_decomposition <- function(data_vector, periods = 4) {
  cat("\n========== TOPIC 3: TIME SERIES DECOMPOSITION ==========\n")
  
  if(periods %% 2 != 0) {
    trend <- stats::filter(data_vector, rep(1/periods, periods), sides = 2)
  } else {
    weights <- c(0.5, rep(1, periods - 1), 0.5) / periods
    trend <- stats::filter(data_vector, weights, sides = 2)
  }
  
  est_ale <- data_vector / trend
  est_matrix <- matrix(est_ale, ncol = periods, byrow = TRUE)
  uncorrected_indices <- colMeans(est_matrix, na.rm = TRUE)
  correction_factor <- periods / sum(uncorrected_indices)
  corrected_indices <- uncorrected_indices * correction_factor
  
  full_indices <- rep(corrected_indices, length.out = length(data_vector))
  deseasonalized <- data_vector / full_indices
  
  cat("1. Corrected Seasonal Indices (E):\n")
  print(corrected_indices * 100) 
  
  cat("\n2. Deseasonalized Series (X / E):\n")
  print(deseasonalized)
  
  cat("\n3. True Trend (Regression of Deseasonalized Series):\n")
  t_index <- seq_along(deseasonalized)
  model <- lm(deseasonalized ~ t_index)
  cat("   Equation: X =", coef(model)[1], "+", coef(model)[2], "* t\n")
  cat("========================================================\n")
}

#' @export
solve_topic3_deflation <- function(nominal_values, cpi_percentages) {
  cat("\n========== TOPIC 3: DEFLATION (REAL VALUES) ==========\n")
  real_values <- nominal_values / (cpi_percentages / 100)
  cat("Nominal Values: ", nominal_values, "\n")
  cat("CPI (%):        ", cpi_percentages, "\n")
  cat("Real Values:    ", real_values, "\n")
  cat("======================================================\n")
}