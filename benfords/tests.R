##### Mantissa #####
mantissa_distance <- function(vector) {
  # create a distribution table
  distribution_table <- check_mantissae(vector)
  tabulated <- table(distribution_table$mantissa)
  probabilities <- tabulated / sum(tabulated) 
  DT_probabilities <- as.data.table(probabilities)[
    j = actCumFreq:=cumsum(N)
  ][
    j = expCumFreq:=.I / .N
  ][
    j = distance:=.(actCumFreq-expCumFreq)
  ]
  return(max(DT_probabilities$distance))
}

check_mantissae <- function(vector) {
  value <- log10(vector)
  mantissa <- value - floor(value)
  DT <- data.table(mantissa = sort(mantissa))[, id:=.I]
  return(DT)
}


##### SSD #####
sum.sqr.dev <- function(data_frequencies, expected_frequencies) {
  sum((data_frequencies*100 - expected_frequencies*100)^2)
}

SSD.conformity <- 	function(SSD = NULL,
                            digits.used = c("First Digit",
                                            "First-Two Digits")){ 
  # from Kossovsky, A. E. (2021)
  Conformity.Levels <- c("Perfectly Benford", 
                         "Acceptably close", 
                         "Marginally Benford", 
                         "Non Benford")
  
  ssd.intervals <- switch(digits.used,
                          "First Digit" = c(0, 2, 25, 100),
                          "First-Two Digits" = c(0, 2, 10, 50)
  )
  
  conformity <- Conformity.Levels[findInterval(SSD, ssd.intervals)]
  
  return(conformity)
}

##### MAD #####
mean.abs.dev <- function(data_frequencies, expected_frequencies) {
  sum(abs(data_frequencies - expected_frequencies)/(length(expected_frequencies)))
}

MAD.conformity <- 	function(MAD = NULL,
                            digits.used = c("First Digit",
                                            "Second Digit",
                                            "First-Two Digits")){
  Conformity.Levels <- c("Close conformity", 
                         "Acceptable conformity", 
                         "Marginally acceptable conformity", 
                         "Nonconformity")
  mad.intervals <- switch(digits.used,
                          "First Digit" = c(0.000, 0.006, 0.012, 0.015),
                          "Second Digit" = c(0.000, 0.008, 0.010, 0.011),
                          "First-Two Digits" = c(0.000, 0.0012, 0.0018, 0.0022)
  )
  
  conformity <- Conformity.Levels[findInterval(MAD, mad.intervals)]
  
  out <- list(MAD = MAD, digits.used = digits.used, conformity = conformity)
  
  return(out)
}

##### overBenford #####
overBenford <- function(vector, digits) {
  vector <- vector[vector!=0]
  n <- length(vector)
  observed_distribution <- distribution_table(vector, digits)$N
  expected_distribution <- get_expected_frequency(digits) * n
  actual_chi_sq <- chisq.distance(observed_distribution, expected_distribution)
  expected_benford_distribution <- monte_carlo_overBenfords_chisq(digits, R = repetitions, n = n)    
  
  over <- data.table(expected_benford_distribution)
  over <- over[j = exceeds:= fifelse(expected_benford_distribution > actual_chi_sq, T, F)]
  if (mean(as.double(over$exceeds)) < 0.01) {
    p.value <- paste0(format(round(mean(as.double(over$exceeds)),3)), " **")
  } else if (mean(as.double(over$exceeds)) < 0.05) {
    p.value <- paste0(format(round(mean(as.double(over$exceeds)),3)), " *")
  } else {
    p.value <- format(round(mean(as.double(over$exceeds)),3))
  }
  if (mean(as.double(over$exceeds))>=0.05) {
    result <- "Accepted"   
  } else {
    result <- "Rejected"
  }
  
  return(list(result = result, p.value = p.value, observed = actual_chi_sq, boot_table = over, N = n, R = repetitions))
}

chisq.distance <- function(observed, expected){
  # chi-squared distance calculated the same way as chi-squared
  # seee https://www.researchgate.net/post/What-is-chi-squared-distance-I-need-help-with-the-source-code
  # https://www.geeksforgeeks.org/chi-square-distance-in-python/
  
  # chi-squared test statistic
  # chisq.test <- sum(((observed - expected)^2)/expected)
  
  # chi-squared distance (for comparing two histograms)
  chisq.distance <- sum(((observed - expected)^2)/ (observed + expected)) / 2
  return(chisq.distance)
}

monte_carlo_overBenfords_chisq <- function(digits, R, n) {
  results <- c(chi_sq = numeric())
  expected <- get_expected_frequency(digits) * n
  for (r in 1:R) {
    simulated <- create_benfords_distribution(n, digits)
    observed <- distribution_table(simulated, digits)$N
    chi_sq <- chisq.distance(observed, expected)
    results <- c(results, chi_sq)
  }
  return(results)
}

##### BREG #####
library(boot)
run_BREG <- function(vector) {
  vector <- vector[vector!=0]
  if (max(floor(log10(vector)), na.rm = T) == 0) {
    # insufficient range (there need to be some two digit numbers)
    BREG_output <- NA_character_
    BREG_ideal <- NA_character_
    BREG <- evaluate_BREG(BREG_output, BREG_ideal)
  } else {
    BREG_output <- boot(vector, bootstrap_regression, R=repetitions)
    BREG_ideal <- boot(create_benfords_distribution(length(vector),2), bootstrap_regression, R = repetitions)
    BREG <- evaluate_BREG(BREG_output, BREG_ideal)   
  }
  
  return(BREG)
}

benford_regression_model <- function (actual_first_two_digits) {
  # the first-two digits distribution is derived from the distributions of the first and the second digits
  # an actual first-two digits can therefore be compared with the one expected from the ideal first and second digits 
  # the regression intercept (alpha_zero) and slope (alpha_one) are calculated in both directions for 
  # actual first two vs ideal first and second
  
  mean_d1 <- sum(seq(from = 1, to = 9) * first_digit_frequency())
  mean_d2 <- sum(seq(from = 0, to = 9) * second_digit_frequency())
  sd_d1 <- (sum(seq(from = 1, to = 9)^2 * first_digit_frequency()) - mean_d1^2)^0.5
  sd_d2 <- (sum(seq(from = 0, to = 9)^2 * second_digit_frequency()) - mean_d2^2)^0.5
  
  # turn the data vector into a frequency table
  digit_frequency <- frequency_table(actual_first_two_digits, 2)
  
  first_two_distribution <- data.frame(
    first_two_digits = seq(10,99), 
    first_two_freqs = as.numeric(digit_frequency$frequency)
  )
  
  first_two_distribution$first_digit <- get_leading_digits(first_two_distribution$first_two_digits)
  first_two_distribution$second_digit <-first_two_distribution$first_two_digits%%10
  
  actual_EV_first_two <- sum(
    first_two_distribution$first_digit * first_two_distribution$second_digit * first_two_distribution$first_two_freqs
  )
  
  pearsons_correlation <- (actual_EV_first_two - (mean_d1 * mean_d2)) / (sd_d1 * sd_d2)
  
  slope_alpha_one <- pearsons_correlation * sd_d1 / sd_d2
  intercept_alpha_zero <- mean_d1 - slope_alpha_one * mean_d2
  
  slope_beta_one <-pearsons_correlation * sd_d2 / sd_d1
  intercept_beta_zero <- mean_d2 - slope_beta_one * mean_d1
  
  regression_model <- data.frame(intercept_alpha_zero, slope_alpha_one, intercept_beta_zero, slope_beta_one)
  return(regression_model)
}

bootstrap_regression <- function(data, indices){ 
  data <- as.data.frame(data)
  dt<-data[indices,]
  
  first_two_digits <- get_leading_digits(dt, 2)
  
  # calculate the regression model between first-two and first and second digits
  regression_model <- benford_regression_model(first_two_digits)
  
  c(regression_model$intercept_alpha_zero,
    regression_model$slope_alpha_one,
    regression_model$intercept_beta_zero,
    regression_model$slope_beta_one
  )    
}
evaluate_BREG <- function(BREG_output, BREG_ideal) { 
  if (length(BREG_output)==11) {
    output_parameters <- BREG_output$t0
    ideal_parameters <- c(3.239170, 0.048018, 3.962264, 0.065439)
    standard_errors_output <- apply(BREG_output$t,2,sd)
    standard_errors_ideal <- apply(BREG_ideal$t,2,sd)
    
    # 95% two-tailed confidence limit
    # subject to Bonferroni correction as there are two estimates involved
    CI_output <- standard_errors_output * qnorm(p=1 - (0.025 / 2 ))
    CI_ideal <- standard_errors_ideal * qnorm(p=1 - (0.025 / 2 ))
    
    # does the output CI include the ideal parameter? this is the two star result
    two_stars <- ideal_parameters - CI_ideal <= output_parameters & 
      ideal_parameters + CI_ideal >= output_parameters
    
    # does the output CI overlap with the CI of the ideal parameters? this is the one star result  
    one_star <- ideal_parameters - CI_ideal >= output_parameters + CI_output | 
      ideal_parameters + CI_ideal >= output_parameters - CI_output
    stars <- str_replace(two_stars, "TRUE", "**")
    stars <- str_replace(stars, "FALSE", "")
    
    for (i in 1:4) {
      if (stars[i] != "**") {
        if(one_star[i] == TRUE) {stars[i] <- "*"} else {stars[i] <- ""}
      }
    }
    
    bootstrap_CIs <-c()
    benford_CIs <-c()
    
    bootstrap_CIs[1] <- str_trim(paste(round(output_parameters[1] - CI_output[1],3), "<=", "a-intercept", "<=", round(output_parameters[1] + CI_output[1],3), stars[1], " "))  
    bootstrap_CIs[2] <- str_trim(paste(round(output_parameters[2] - CI_output[2],3), "<=", "a-slope", "<=", round(output_parameters[2] + CI_output[2],3), stars[2], " ")) 
    bootstrap_CIs[3] <- str_trim(paste(round(output_parameters[3] - CI_output[3],3), "<=", "b-intercept", "<=", round(output_parameters[3] + CI_output[3],3), stars[3], " "))
    bootstrap_CIs[4] <- str_trim(paste(round(output_parameters[4] - CI_output[4],3), "<=", "b-slope", "<=", round(output_parameters[4] + CI_output[4],3), stars[4], " ")) 
    
    benford_CIs[1] <- paste(round(ideal_parameters[1] - CI_ideal[1],3), "<=", "a-intercept", "<=", round(ideal_parameters[1] + CI_ideal[1],3), " ")  
    benford_CIs[2] <- paste(round(ideal_parameters[2] - CI_ideal[2],3), "<=", "a-slope", "<=", round(ideal_parameters[2] + CI_ideal[2],3), " ") 
    benford_CIs[3] <- paste(round(ideal_parameters[3] - CI_ideal[3],3), "<=", "b-intercept", "<=", round(ideal_parameters[3] + CI_ideal[3],3), " ") 
    benford_CIs[4] <- paste(round(ideal_parameters[4] - CI_ideal[4],3), "<=", "b-slope", "<=", round(ideal_parameters[4] + CI_ideal[4],3), " ") 
    
    score <- 0
    for (i in 1:4) {
      if (sum(two_stars[i] + one_star[i]) >= 1) {
        score <- score + 1
      }
    }
    
    if (score == 4) {
      verdict <- "Conformant"
    } else  {
      verdict <- "Not conformant"
    }
  } else {
    verdict = NA
    output_parameters = NA
    ideal_parameters = NA 
    bootstrap_CIs=NA 
    benford_CIs= NA
  }
  
  return(list(
    verdict = verdict, 
    output_parameters = output_parameters, 
    ideal_parameters = ideal_parameters, 
    bootstrap_CIs=bootstrap_CIs, 
    benford_CIs=benford_CIs)
  )
}

##### Kuiper #####
kuiper_test <- function(vector, digits) {
  # strip out any zeros
  vector <- vector[vector!=0]
  
  # determine the type of distribution
  mantissa_difference = mantissa_distance(vector)
  
  data <- get_leading_digits(vector, digits)
  # bootstrap some random samples that are the same size as the data
  # to determine how much variation may be expected 
  noise_factor <- 4
  
  results <- c(kuiper_simulated = numeric())
  # if (mantissa_difference <= 0.3) {
  n <- ceiling(length(vector) / noise_factor)
  expected <-  data.table(freq=get_expected_frequency(digits))[j = cum_freq:=cumsum(freq)]$cum_freq
  for (rep in 1:repetitions) {
    simulated <- create_benfords_distribution(n, digits)
    simulated <- frequency_table(simulated, digits)$cum_freq
    kuiper_simulated <- kuiper.statistic(simulated, expected)
    results <- c(results, kuiper_simulated)
  }    
  
  #  } else {
  # n <- length(vector) 
  # expected <- data.table(freq=get_extreme_frequency(digits))[j = cum_freq:=cumsum(freq)]$cum_freq
  # for (rep in 1:repetitions) {
  #   simulated <- create_extreme_distribution(n, digits)
  #   simulated <- frequency_table(simulated, digits)$cum_freq
  #   kuiper_simulated <- kuiper.statistic(simulated, expected)
  #   results <- c(results, kuiper_simulated)
  # }      
  # }
  
  # compare actual kuiper with the bootstrapped one
  kuiper_actual <- kuiper.statistic(frequency_table(data, digits)$cum_freq, expected)
  
  critical_table <- data.table(results)
  critical_table <- critical_table[j = exceeds:= fifelse(results > kuiper_actual, T, F)]
  if (mean(as.double(critical_table$exceeds)) < 0.01) {
    p.value <- paste0(format(round(mean(as.double(critical_table$exceeds)),3)), " **")
  } else if (mean(as.double(critical_table$exceeds)) < 0.05) {
    p.value <- paste0(format(round(mean(as.double(critical_table$exceeds)),3)), " *")
  } else {
    p.value <- format(round(mean(as.double(critical_table$exceeds)),3))
  }
  if (mean(as.double(critical_table$exceeds))>=0.05) {
    result <- "Accepted"   
  } else {
    result <- "Rejected"
  } 
  return(list(
    result = result, p.value = p.value, observed = kuiper_actual, boot_table = critical_table, N = n, R = repetitions))
}

kuiper.distance <- function(observed, expected){
  # equation 4 in Lanzante 2021
  # kuiper difference V = |D+| + |D−|
  
  max_positive <- max(observed - expected)
  max_negative <- max(expected - observed)
  
  return(max_positive + max_negative)
}

kuiper.statistic <- function(observed, expected){
  # from Koch and Okamura 2020
  # Koch and Okamura use the cumulative frequencies of the first digit d in the observed and the Benfords distribution
  
  # calculate Kuiper distance
  V <- kuiper.distance(observed, expected)
  
  # calculate sample size adjustment
  N_adj <- length(observed)^0.5 + 0.155 + 0.24 / length(observed)^0.5
  
  # return test statistic  
  return(V * N_adj)
}


kuiper_extreme <- function(vector_of_data, digits) {
  # strip out any zeros
  vector_of_data <- vector_of_data[vector_of_data!=0]
  
  # determine the type of distribution
  mantissa_difference = mantissa_distance(vector_of_data)
  
  data <- get_leading_digits(vector_of_data, digits)
  # bootstrap some random samples that are the same size as the data
  # to determine how much variation may be expected 
  noise_factor <- 4
  
  results <- c(kuiper_simulated = numeric())
  if (mantissa_difference <= 0.30) {
    n <- ceiling(length(vector_of_data) / noise_factor)
    expected <-  data.table(freq=get_expected_frequency(digits))[j = cum_freq:=cumsum(freq)]$cum_freq
    for (rep in 1:repetitions) {
      simulated <- create_benfords_distribution(n, digits)
      simulated <- frequency_table(simulated, digits)$cum_freq
      kuiper_simulated <- kuiper.statistic(simulated, expected)
      results <- c(results, kuiper_simulated)
    }    
    
  } else {
    n <- length(vector_of_data) 
    expected <- data.table(freq=get_extreme_frequency(digits))[j = cum_freq:=cumsum(freq)]$cum_freq
    for (rep in 1:repetitions) {
      simulated <- create_extreme_distribution(n, digits)
      simulated <- frequency_table(simulated, digits)$cum_freq
      kuiper_simulated <- kuiper.statistic(simulated, expected)
      results <- c(results, kuiper_simulated)
    }      
  }
  
  # compare actual kuiper with the bootstrapped one
  kuiper_actual <- kuiper.statistic(frequency_table(data, digits)$cum_freq, expected)
  critical_table <- data.table(results)
  critical_table <- critical_table[j = exceeds:= fifelse(results > kuiper_actual, T, F)]
  if (mean(as.double(critical_table$exceeds)) < 0.01) {
    p.value <- paste0(format(round(mean(as.double(critical_table$exceeds)),3)), " **")
  } else if (mean(as.double(critical_table$exceeds)) < 0.05) {
    p.value <- paste0(format(round(mean(as.double(critical_table$exceeds)),3)), " *")
  } else {
    p.value <- format(round(mean(as.double(critical_table$exceeds)),3))
  }
  if (mean(as.double(critical_table$exceeds))>=0.05) {
    result <- "Accepted"   
  } else {
    result <- "Rejected"
  } 
  return(list(
    result = result, p.value = p.value, observed = kuiper_actual, boot_table = critical_table, N = n, R = repetitions))
}
