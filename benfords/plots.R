
plot_mantissae <- function(vector, title = "") {
  # remove zeros
  vector <- vector[vector!=0]
  
  # find the p value for the mantissa arc line
  benfords_1_D <- benford(vector, number.of.digits  = 1)
  mac <- paste(
    "p-value = ", as.character(signif(benfords_1_D$stats$mantissa.arc.test$p.value,3)))
  
  # plot the mantissae 
  DT_mantissa <- check_mantissae(vector)
  mantissa_plot <- ggplot(DT_mantissa, aes(x = id, y = mantissa)) + 
    geom_line() + 
    annotate("text", x = nrow(DT_mantissa) / 1.4, y = 0.3, label = "Mantissa Arc Test") +
    annotate("text", x = nrow(DT_mantissa) / 1.4, y = 0.2, label = mac) +
    annotate("text", x = nrow(DT_mantissa) / 1.4, y = 0.1, label = paste0(
      "N = ", format(nrow(DT_mantissa), big.mark = ","))) +
    geom_abline(colour = "red", slope = 1 / nrow(DT_mantissa)) + 
    ggtitle(title) + xlab("Record")+ coord_cartesian(ylim=c(0,1))
  return(mantissa_plot)
}

plot_benford_leading <- function(vector, 
                                 test.type = "First Digit", 
                                 title = "", subtitle = "", 
                                 variable = "", 
                                 test_overBenford = TRUE, 
                                 test_BREG = TRUE, 
                                 test_Kuiper = TRUE) {
  # remove 0's
  vector <- vector[vector!=0]
  
  # check magnitudes
  oom_10 <- robust_magnitude_10(vector)
  oom_5 <- robust_magnitude_5(vector)
  oom_naive <- naive_magnitude(vector)
  
  # check overSkewness
  ES12 <- ES12(vector)
  
  if (test.type == "First Digit") {
    offset <- 0
    benford_result <- benford(vector,1)
    # OverBenford with N replications
    if (test_overBenford == T) {
      oB <- overBenford(vector, 1)
      overBenford.result <- oB$result
      overBenford.p <- oB$p.value}

    # Bootstrap with N replications
    if (test_BREG == T) {
      BREG <- run_BREG(vector)}
    
    # Kuiper with N replications
    if (test_Kuiper ==T) {
      kuiper_result <- kuiper_test(vector, 1)
      kuiper.result <- kuiper_result$result
      kuiper.p <- kuiper_result$p.value}
  } 
  
  if (test.type == "First-Two Digits") {
    offset <- 10
    # check magnitude
    magnitude <- ifelse(vector > 0, floor(log10(vector)), NA_integer_)
    if (max(magnitude, na.rm = T) == 0) {
      # there are not enough for a full first-two digit
      benford_result <- benford(vector,2)
      # OverBenford with N replications
      if (test_overBenford == T) {
        overBenford.result <- NA_character_
        overBenford.p <- NA}
      
      # Bootstrap with N replications
      if (test_BREG == T) {
        BREG <- run_BREG(vector)}
      
      # Kuiper with N replications
      if (test_Kuiper ==T) {
        kuiper.result <- NA_character_
        kuiper.p <- NA}
      
    } else {
      # remove single digits (they are already covered by the first digit test)
      benford_result <- benford(vector[vector>=10],2)
      # OverBenford with N replications
      if (test_overBenford == T) {
        oB <- overBenford(vector[vector>=10], 2)
        overBenford.result <- oB$result
        overBenford.p <- oB$p.value}
      
      # Bootstrap with N replications
      if (test_BREG == T) {
        BREG <- run_BREG(vector)}
      
      # Kuiper with N replications
      if (test_Kuiper ==T) {
        kuiper_result <- kuiper_test(vector[vector>=10], 2)
        kuiper.result <- kuiper_result$result
        kuiper.p <- kuiper_result$p.value}
    }
  }

  plot_data <- benford_result$bfd    
  max_y <- max(max(plot_data$data.dist), plot_data$benford.dist.freq[1] / benford_result$info$n)
  max_x <- max(plot_data$digits)
  x_digits <- plot_data$digits
  min_x <- min(plot_data$digits)
  expected_dist <- plot_data$benford.dist
  nigrini <- benford_result$MAD.conformity
  mad <- benford_result$MAD
  kossovsky <- SSD.conformity(sum((plot_data$data.dist*100 - plot_data$benford.dist*100)^2), test.type)
  ssd <- sum((plot_data$data.dist*100 - plot_data$benford.dist*100)^2)
  chisq <- benford_result$stats$chisq$statistic
  N <- benford_result$info$n
  df <- benford_result$stats$chisq$parameter
  critical_value <- qchisq(p=0.05, df = benford_result$stats$chisq$parameter, lower.tail = F)
  p_value <- benford_result$stats$chisq$p.value
  lines <-function(x){return(0.1 +1 - x/10)}
  
  p <- ggplot(plot_data, aes(x = as.factor(digits), y = data.dist)) + geom_col() + 
    annotate("text", x = max_x-offset, y =max_y *lines(1), label = paste0("Records: ", benford_result$info$n), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(2), label = paste0("MAD: ", round(benford_result$MAD, 4), " (", benford_result$MAD.conformity, ")"), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(3), label = paste0("SSD: ",sprintf("%0.2f", sum((benford_result$bfd$data.dist*100 - benford_result$bfd$benford.dist*100)^2)), " (", SSD.conformity(sum((benford_result$bfd$data.dist*100 - benford_result$bfd$benford.dist*100)^2), test.type), ")"), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(4), label = paste0("Chi-Squared: ", round(benford_result$stats$chisq$statistic, 1), ", p: ", as.character(signif(benford_result$stats$chisq$p.value,2)), ", cv: ", round(qchisq(p=0.05, df = benford_result$stats$chisq$parameter, lower.tail = F),1)),size = 3, hjust = 1) +
    geom_line(aes(x=as.factor(x_digits), y = expected_dist, group = 1, col = "red"), show.legend = FALSE) +
    xlab("First digit")+ ylab("Proportion") + ggtitle(title, subtitle = subtitle)
  
  if (test.type == "First-Two Digits") {
    p <- p +
      scale_x_discrete(breaks = seq(from = min_x, to = max_x, by = min_x), labels = seq(from = min_x, to = max_x, by = min_x)) +
      xlab("First-two digits")
  }
  
  line_count <- 4
  if (test_overBenford == T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("OverBenford: ", overBenford.result, " (", overBenford.p, ")"), size = 3, hjust = 1)
  }
  if (test_BREG == T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("BREG: ", BREG$verdict), size = 3, hjust = 1)
  }
  if (test_Kuiper ==T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("Kuiper: ", kuiper.result, " (", kuiper.p, ")"), size = 3, hjust = 1)
  }  
  
  if (test_overBenford == F) {
    overBenford.result = NA_character_
    overBenford.p = NA}
  
  if (test_BREG == F) {
    BREG <- data.frame(verdict = NA_character_, output_parameters = NA, ideal_parameters = NA, bootstrap_CIs=NA_character_, benford_CIs=NA_character_)
  }
  
  if (test_Kuiper ==F) {
    kuiper.result <- NA_character_
    kuiper.p <- NA}
    
  stats <- data.table(
    region = subtitle, 
    dataset = variable,
    test = test.type,
    total = sum(benford_result$data$data.used),
    range = max(fifelse(benford_result$data$data.used > 0, ceiling(log10(benford_result$data$data.used)), 0)),
    mantissa_difference = mantissa_distance(vector),
    oom_10 = oom_10,
    oom_5 = oom_5,
    oom_naive= oom_naive,
    ES12 = ES12,
    Nigrini = nigrini, 
    MAD = mad,
    Kossovsky = kossovsky,
    SSD = round(ssd, 3),
    Chi.squared = round(chisq, 3),
    n = N,
    df = df,
    p = signif(p_value,3),
    kuiper.result = kuiper.result,
    kuiper.p = kuiper.p,
    overBenford.result = overBenford.result,
    overBenford.p = overBenford.p,
    BREG_verdict = BREG$verdict,
    BREG_output = BREG$output_parameters[1], 
    BREG_ideal = BREG$ideal_parameters[1], 
    BREG_bootstrap_CI = BREG$bootstrap_CIs[1], 
    BREG_benford_CI = BREG$benford_CIs[1]
  )
  
  return(list(plot = p, stats = stats))
}


plot_benford_second <- function(tbl_distribution, title = "", subtitle = "") {
  
  plot_data <- tbl_distribution   
  max_y <- max(plot_data$data.prop)*1.3
  max_x <- max(plot_data$second_digit)+1
  mad <- mean.abs.dev(plot_data$data.prop, second_digit_frequency())
  ssd <- sum.sqr.dev(plot_data$data.prop, second_digit_frequency())
  chisq <- chisq.test(plot_data[j = .(data.dist)][j = exp:=sum(data.dist) * second_digit_frequency()])$statistic
  N <- sum(plot_data$data.dist)
  p_value <- chisq.test(plot_data[j = .(data.dist)][j = exp:=sum(data.dist) * second_digit_frequency()])$p.value  
  
  ggplot(plot_data, aes(x = as.factor(second_digit), y = data.prop)) + geom_bar(stat = 'identity') + 
    geom_line(aes(x=as.factor(second_digit), y = second_digit_frequency(), group = 1, col = "red")) +
    annotate("text", x = max_x, y = max_y * 1, label = paste0("Records: ", as.character(N)), size = 3, hjust = 1) +
    annotate("text", x = max_x, y = max_y *0.9, label = paste0("MAD: ", round(mad,4)), size = 3, hjust = 1) + 
    annotate("text", x = max_x, y = max_y *0.8, label = paste0("SSD: ", round(ssd,4)), size = 3, hjust = 1) +
    annotate("text", x = max_x, y = max_y *0.7, label = paste0("chi-squared: ", round(chisq, 3), ", p: ", as.character(signif(p_value,3))), size = 3, hjust = 1) +
    scale_x_discrete(drop = FALSE, breaks = seq(from = 0, to = 9, by =1))+
    xlab("Digit") + ylab("Proportion") + ggtitle(title, subtitle = subtitle)
  
}

plot_last_one <- function(vector, variable, magnitude_range, title = "", subtitle = "") {
  if (length(vector) > 0) {
    # make a template of zeros
    template <- data.table(last_digit = seq(from = 0, to = 9), data.dist = 0)
    # find the last digit
    last_one <- data.table(last_digit = vector %% 10)[
      j = .(data.dist=.N), by = .(last_digit)][j = data.prop:=data.dist / sum(data.dist)]
    # remove the last_one items from the template
    missing_numbers <- template[!last_one, on=.(last_digit)]
    # combine the two tables
    last_one <- rbindlist(list(last_one, missing_numbers), fill=T)[j = data.prop:=fifelse(is.na(data.prop), 0, data.prop)]
    # make the digit a factor
    last_one <- last_one[order(last_digit)][j = last_digit:=factor(last_digit, labels = seq(from = 0, to = 9))]
    
    if (magnitude_range > 1) {
      # stats
      last_MAD <- mean.abs.dev(last_one$data.prop, rep(0.1,10))
      last_SSD <- sum.sqr.dev(last_one$data.prop, rep(0.1,10))
      last_chi <- chisq.test(last_one[j = .(data.dist)])
      last_chi_stat <- last_chi$statistic
      last_chi_df <- last_chi$parameter
      last_chi_p <- last_chi$p.value     
    } else {
      # the last digit is actually the second digit
      last_MAD <- mean.abs.dev(last_one$data.prop, second_digit_frequency())
      last_SSD <- sum.sqr.dev(last_one$data.prop, second_digit_frequency())
      last_chi <- chisq.test(last_one[j = .(data.dist)], p = second_digit_frequency())   
      last_chi_stat <- last_chi$statistic
      last_chi_df <- last_chi$parameter
      last_chi_p <- last_chi$p.value
    }
  } else {
    last_one <- data.table(last_digit=seq(from = 0, to = 9), data.dist = 0, data.prop = 0)
    last_one <- last_one[order(last_digit)][j = last_digit:=factor(last_digit, labels = seq(from = 0, to = 9))]
    last_MAD <- NA
    last_SSD <- NA
    last_chi_stat <- NA
    last_chi_df <- NA
    last_chi_p <- NA
  }
  # plot
  last_one_plot <- ggplot(last_one, aes(x = last_digit, y = data.prop)) + 
    geom_bar(stat='identity') + 
    annotate("text", x = 10, y =max(last_one$data.prop) *1.15, 
             label = paste0("MAD: ",round(last_MAD,4)), 
             size = 2, hjust = 1) +
    annotate("text", x = 10, y =max(last_one$data.prop) *1.05, 
             label = paste0("SSD: ",round(last_SSD,4)), 
             size = 2, hjust = 1) +
    annotate("text", x = 4, y =max(last_one$data.prop) *1.15, 
             label = paste0(
               "X-sq: ", round(last_chi_stat, 3)), 
             size = 2, hjust = 1) +
    annotate("text", x = 4, y =max(last_one$data.prop) *1.05, 
             label = paste0(
               ", p = ", as.character(signif(last_chi_p,3))), 
             size = 2, hjust = 1) 
  # if (magnitude_range == 1) {
  #   last_one_plot <- last_one_plot + 
  #     geom_line(aes(x=last_one$last_digit, y = first_digit_frequency(), group = 1, col = "red"), show.legend = FALSE) +
  #     xlab("last_digit")
  # } else {
    last_one_plot <- last_one_plot + geom_hline(yintercept=1/10, color = "red") # + ggtitle(subtitle)  
  # }  

  return(last_one_plot)
}

plot_last_two <- function(vector, variable, title = "", subtitle = "") {
  if (length(vector) > 0) {
    # make a template of zeros
    template <- data.table(last_2digits = seq(from = 0, to = 99), data.dist = 0)
    # remove values less than 10
    vector <- vector[vector >= 10]
    # find the last two digits
    last_two <- data.table(last_2digits = vector %% 100)[
      j = .(data.dist=.N), by = .(last_2digits)][j = data.prop:=data.dist / sum(data.dist)]
    # remove the last_two items from the template
    missing_numbers <- template[!last_two, on=.(last_2digits)]
    # combine the two tables
    last_two <- rbindlist(list(last_two, missing_numbers), fill=T)[j = data.prop:=fifelse(is.na(data.prop), 0, data.prop)]
    # make the digit a factor
    last_two <- last_two[order(last_2digits)][j = last_2digits:=factor(last_2digits, labels = seq(from = 0, to = 99))]
    
    # stats
    last_MAD <- mean.abs.dev(last_two$data.prop, rep(0.1,10))
    last_SSD <- sum.sqr.dev(last_two$data.prop, rep(0.1,10))
    last_chi_stat <- chisq.test(last_two[j = .(data.dist)])$statistic
    last_chi_df <- chisq.test(last_two[j = .(data.dist)])$parameter
    last_chi_p <- chisq.test(last_two[j = .(data.dist)])$p.value
    
  } else {
    last_two <- data.table(last_2digits=seq(from = 0, to = 99), data.dist = 0, data.prop = 0)
    last_MAD <- NA
    last_SSD <- NA
    last_chi_stat <- NA
    last_chi_df <- NA
    last_chi_p <- NA
  }
  # plot
  last_two_plot <- ggplot(last_two, aes(x = last_2digits, y = data.prop)) + 
    geom_col(stat = 'identity') + scale_x_discrete(breaks=seq(from = 0, to = 90, by =10)) +
    annotate("text", x = 90, y =max(0.02, last_two$data.prop) *1.15, 
             label = paste0("MAD: ",round(last_MAD,4)), 
             size = 2, hjust = 1) +
    annotate("text", x = 90, y =max(0.02, last_two$data.prop) *1.05, 
             label = paste0("SSD: ",round(last_SSD,4)), 
             size = 2, hjust = 1) +
    annotate("text", x = 40, y =max(0.02, last_two$data.prop) *1.15,
             label = paste0(
               "X-sq: ", round(last_chi_stat, 3)),
             size = 2, hjust = 1) +
    annotate("text", x = 40, y =max(0.02, last_two$data.prop) *1.05,
             label = paste0(
               ", p = ", as.character(signif(last_chi_p,3))),
             size = 2, hjust = 1) +
    geom_hline(yintercept=1/100, color = "red") # + ggtitle(subtitle)
  
 
  return(last_two_plot)
}


plot_benford_leading_timed <- function(vector, 
                                 test.type = "First Digit", 
                                 title = "", subtitle = "", 
                                 variable = "", 
                                 test_overBenford = TRUE, 
                                 test_BREG = TRUE, 
                                 test_Kuiper = TRUE) {
  # remove 0's
  vector <- vector[vector!=0]
  
  # check magnitudes
  oom_10 <- robust_magnitude_10(vector)
  oom_5 <- robust_magnitude_5(vector)
  oom_naive <- naive_magnitude(vector)
  
  # check overSkewness
  ES12 <- ES12(vector)
  
  if (test.type == "First Digit") {
    offset <- 0
    start <- Sys.time()
    benford_result <- benford(vector,1)
    end <- Sys.time()
    benford_time <- end - start
    
    # OverBenford with N replications
    start <- Sys.time()  
    if (test_overBenford == T) {
      oB <- overBenford(vector, 1)
      overBenford.result <- oB$result
      overBenford.p <- oB$p.value
      
      end <- Sys.time()
      overBenford_time <- end - start}
    
    # Bootstrap with N replications
    start <- Sys.time()
    if (test_BREG == T) {
      BREG <- run_BREG(vector)
      end <- Sys.time()
      BREG_time <- end - start}
    
    # Kuiper with N replications
    start <- Sys.time()
    if (test_Kuiper ==T) {
      kuiper_result <- kuiper_test(vector, 1)
      kuiper.result <- kuiper_result$result
      kuiper.p <- kuiper_result$p.value
      end <- Sys.time()
      kuiper_time <- end - start}
  } 
  
  if (test.type == "First-Two Digits") {
    offset <- 10
    # check magnitude
    magnitude <- ifelse(vector > 0, floor(log10(vector)), NA_integer_)
    if (max(magnitude, na.rm = T) == 0) {
      # there are not enough for a full first-two digit
      start <- Sys.time()
      benford_result <- benford(vector,2)
      end <- Sys.time()
      benford_time <- end - start
      
      # OverBenford with N replications
      if (test_overBenford == T) {
        start <- Sys.time()
        overBenford.result <- NA_character_
        overBenford.p <- NA
        end <- Sys.time()
        overBenford_time <- end - start
        }
      
      # Bootstrap with N replications
      if (test_BREG == T) {
        start <- Sys.time()
        BREG <- run_BREG(vector)
        end <- Sys.time()
        BREG_time <- end - start}
      
      # Kuiper with N replications
      if (test_Kuiper ==T) {
        start <- Sys.time()
        kuiper.result <- NA_character_
        kuiper.p <- NA
        end <- Sys.time()
        kuiper_time <- end - start
        }
      
    } else {
      # remove single digits (they are already covered by the first digit test)
      start <- Sys.time()
      benford_result <- benford(vector[vector>=10],2)
      end <- Sys.time()
      benford_time <- end - start
      
      # OverBenford with N replications
      if (test_overBenford == T) {
        start <- Sys.time()
        oB <- overBenford(vector[vector>=10], 2)
        overBenford.result <- oB$result
        overBenford.p <- oB$p.value
        end <- Sys.time()
        overBenford_time <- end - start
        }
      
      # Bootstrap with N replications
      if (test_BREG == T) {
        start <- Sys.time()
        BREG <- run_BREG(vector)
        end <- Sys.time()
        BREG_time <- end - start
        }
      
      # Kuiper with N replications
      if (test_Kuiper ==T) {
        start <- Sys.time()
        kuiper_result <- kuiper_test(vector[vector>=10], 2)
        kuiper.result <- kuiper_result$result
        kuiper.p <- kuiper_result$p.value
        end <- Sys.time()
        kuiper_time <- end - start
        }
    }
  }
  
  plot_data <- benford_result$bfd    
  max_y <- max(max(plot_data$data.dist), plot_data$benford.dist.freq[1] / benford_result$info$n)
  max_x <- max(plot_data$digits)
  x_digits <- plot_data$digits
  min_x <- min(plot_data$digits)
  expected_dist <- plot_data$benford.dist
  nigrini <- benford_result$MAD.conformity
  mad <- benford_result$MAD
  kossovsky <- SSD.conformity(sum((plot_data$data.dist*100 - plot_data$benford.dist*100)^2), test.type)
  ssd <- sum((plot_data$data.dist*100 - plot_data$benford.dist*100)^2)
  chisq <- benford_result$stats$chisq$statistic
  N <- benford_result$info$n
  df <- benford_result$stats$chisq$parameter
  critical_value <- qchisq(p=0.05, df = benford_result$stats$chisq$parameter, lower.tail = F)
  p_value <- benford_result$stats$chisq$p.value
  lines <-function(x){return(0.1 +1 - x/10)}
  
  p <- ggplot(plot_data, aes(x = as.factor(digits), y = data.dist)) + geom_col() + 
    annotate("text", x = max_x-offset, y =max_y *lines(1), label = paste0("Records: ", benford_result$info$n), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(2), label = paste0("MAD: ", round(benford_result$MAD, 4), " (", benford_result$MAD.conformity, ")"), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(3), label = paste0("SSD: ",sprintf("%0.2f", sum((benford_result$bfd$data.dist*100 - benford_result$bfd$benford.dist*100)^2)), " (", SSD.conformity(sum((benford_result$bfd$data.dist*100 - benford_result$bfd$benford.dist*100)^2), test.type), ")"), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(4), label = paste0("Chi-Squared: ", round(benford_result$stats$chisq$statistic, 1), ", p: ", as.character(signif(benford_result$stats$chisq$p.value,2)), ", cv: ", round(qchisq(p=0.05, df = benford_result$stats$chisq$parameter, lower.tail = F),1)),size = 3, hjust = 1) +
    geom_line(aes(x=as.factor(x_digits), y = expected_dist, group = 1, col = "red"), show.legend = FALSE) +
    xlab("First digit")+ ylab("Proportion") + ggtitle(title, subtitle = subtitle)
  
  if (test.type == "First-Two Digits") {
    p <- p +
      scale_x_discrete(breaks = seq(from = min_x, to = max_x, by = min_x), labels = seq(from = min_x, to = max_x, by = min_x)) +
      xlab("First-two digits")
  }
  
  line_count <- 4
  if (test_overBenford == T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("OverBenford: ", overBenford.result, " (", overBenford.p, ")"), size = 3, hjust = 1)
  }
  if (test_BREG == T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("BREG: ", BREG$verdict), size = 3, hjust = 1)
  }
  if (test_Kuiper ==T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("Kuiper: ", kuiper.result, " (", kuiper.p, ")"), size = 3, hjust = 1)
  }  
  
  if (test_overBenford == F) {
    overBenford.result = NA_character_
    overBenford.p = NA}
  
  if (test_BREG == F) {
    BREG <- data.frame(verdict = NA_character_, output_parameters = NA, ideal_parameters = NA, bootstrap_CIs=NA_character_, benford_CIs=NA_character_)
  }
  
  if (test_Kuiper ==F) {
    kuiper.result <- NA_character_
    kuiper.p <- NA}
  
  stats <- data.table(
    region = subtitle, 
    dataset = variable,
    test = test.type,
    total = sum(benford_result$data$data.used),
    range = max(fifelse(benford_result$data$data.used > 0, ceiling(log10(benford_result$data$data.used)), 0)),
    mantissa_difference = mantissa_distance(vector),
    oom_10 = oom_10,
    oom_5 = oom_5,
    oom_naive= oom_naive,
    ES12 = ES12,
    Nigrini = nigrini, 
    MAD = mad,
    Kossovsky = kossovsky,
    SSD = round(ssd, 3),
    Chi.squared = round(chisq, 3),
    n = N,
    df = df,
    p = signif(p_value,3),
    kuiper.result = kuiper.result,
    kuiper.p = kuiper.p,
    overBenford.result = overBenford.result,
    overBenford.p = overBenford.p,
    BREG_verdict = BREG$verdict,
    BREG_output = BREG$output_parameters[1], 
    BREG_ideal = BREG$ideal_parameters[1], 
    BREG_bootstrap_CI = BREG$bootstrap_CIs[1], 
    BREG_benford_CI = BREG$benford_CIs[1]
  )
  
  timings <- data.table(state = subtitle, 
                        variable = variable,
                        MAD = as.numeric(benford_time),
                        BREG = as.numeric(BREG_time),
                        overBenford = as.numeric(overBenford_time),
                        kuiper = as.numeric(kuiper_time))
  
  return(list(plot = p, stats = stats, timings = timings))
}

plot_benford_leading_extreme <- function(vector, 
                                 test.type = "First Digit", 
                                 title = "", subtitle = "", 
                                 variable = "", 
                                 test_overBenford = TRUE, 
                                 test_BREG = TRUE, 
                                 test_Kuiper = TRUE) {
  # remove 0's
  vector <- vector[vector!=0]
  
  # check magnitudes
  oom_10 <- robust_magnitude_10(vector)
  oom_5 <- robust_magnitude_5(vector)
  oom_naive <- naive_magnitude(vector)
  
  # check overSkewness
  ES12 <- ES12(vector)
  
  if (test.type == "First Digit") {
    offset <- 0
    benford_result <- benford(vector,1)
    # OverBenford with N replications
    if (test_overBenford == T) {
      oB <- overBenford(vector, 1)
      overBenford.result <- oB$result
      overBenford.p <- oB$p.value}
    
    # Bootstrap with N replications
    if (test_BREG == T) {
      BREG <- run_BREG(vector)}
    
    # Kuiper with N replications
    if (test_Kuiper ==T) {
      kuiper_result <- kuiper_extreme(vector, 1)
      kuiper.result <- kuiper_result$result
      kuiper.p <- kuiper_result$p.value}
  } 
  
  if (test.type == "First-Two Digits") {
    offset <- 10
    # check magnitude
    magnitude <- ifelse(vector > 0, floor(log10(vector)), NA_integer_)
    if (max(magnitude, na.rm = T) == 0) {
      # there are not enough for a full first-two digit
      benford_result <- benford(vector,2)
      # OverBenford with N replications
      if (test_overBenford == T) {
        overBenford.result <- NA_character_
        overBenford.p <- NA}
      
      # Bootstrap with N replications
      if (test_BREG == T) {
        BREG <- run_BREG(vector)}
      
      # Kuiper with N replications
      if (test_Kuiper ==T) {
        kuiper.result <- NA_character_
        kuiper.p <- NA}
      
    } else {
      # remove single digits (they are already covered by the first digit test)
      benford_result <- benford(vector[vector>=10],2)
      # OverBenford with N replications
      if (test_overBenford == T) {
        oB <- overBenford(vector[vector>=10], 2)
        overBenford.result <- oB$result
        overBenford.p <- oB$p.value}
      
      # Bootstrap with N replications
      if (test_BREG == T) {
        BREG <- run_BREG(vector)}
      
      # Kuiper with N replications
      if (test_Kuiper ==T) {
        kuiper_result <- kuiper_extreme(vector[vector>=10], 2)
        kuiper.result <- kuiper_result$result
        kuiper.p <- kuiper_result$p.value}
    }
  }
  
  plot_data <- benford_result$bfd    
  max_y <- max(max(plot_data$data.dist), plot_data$benford.dist.freq[1] / benford_result$info$n)
  max_x <- max(plot_data$digits)
  x_digits <- plot_data$digits
  min_x <- min(plot_data$digits)
  expected_dist <- plot_data$benford.dist
  nigrini <- benford_result$MAD.conformity
  mad <- benford_result$MAD
  kossovsky <- SSD.conformity(sum((plot_data$data.dist*100 - plot_data$benford.dist*100)^2), test.type)
  ssd <- sum((plot_data$data.dist*100 - plot_data$benford.dist*100)^2)
  chisq <- benford_result$stats$chisq$statistic
  N <- benford_result$info$n
  df <- benford_result$stats$chisq$parameter
  critical_value <- qchisq(p=0.05, df = benford_result$stats$chisq$parameter, lower.tail = F)
  p_value <- benford_result$stats$chisq$p.value
  lines <-function(x){return(0.1 +1 - x/10)}
  
  p <- ggplot(plot_data, aes(x = as.factor(digits), y = data.dist)) + geom_col() + 
    annotate("text", x = max_x-offset, y =max_y *lines(1), label = paste0("Records: ", benford_result$info$n), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(2), label = paste0("MAD: ", round(benford_result$MAD, 4), " (", benford_result$MAD.conformity, ")"), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(3), label = paste0("SSD: ",sprintf("%0.2f", sum((benford_result$bfd$data.dist*100 - benford_result$bfd$benford.dist*100)^2)), " (", SSD.conformity(sum((benford_result$bfd$data.dist*100 - benford_result$bfd$benford.dist*100)^2), test.type), ")"), size = 3, hjust = 1) +
    annotate("text", x = max_x-offset, y =max_y *lines(4), label = paste0("Chi-Squared: ", round(benford_result$stats$chisq$statistic, 1), ", p: ", as.character(signif(benford_result$stats$chisq$p.value,2)), ", cv: ", round(qchisq(p=0.05, df = benford_result$stats$chisq$parameter, lower.tail = F),1)),size = 3, hjust = 1) +
    geom_line(aes(x=as.factor(x_digits), y = expected_dist, group = 1, col = "red"), show.legend = FALSE) +
    xlab("First digit")+ ylab("Proportion") + ggtitle(title, subtitle = subtitle)
  
  if (test.type == "First-Two Digits") {
    p <- p +
      scale_x_discrete(breaks = seq(from = min_x, to = max_x, by = min_x), labels = seq(from = min_x, to = max_x, by = min_x)) +
      xlab("First-two digits")
  }
  
  line_count <- 4
  if (test_overBenford == T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("OverBenford: ", overBenford.result, " (", overBenford.p, ")"), size = 3, hjust = 1)
  }
  if (test_BREG == T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("BREG: ", BREG$verdict), size = 3, hjust = 1)
  }
  if (test_Kuiper ==T) {
    line_count <- line_count + 1
    p <- p + annotate("text", x = max_x-offset, y = max_y *lines(line_count), label = paste0("Kuiper: ", kuiper.result, " (", kuiper.p, ")"), size = 3, hjust = 1)
  }  
  
  if (test_overBenford == F) {
    overBenford.result = NA_character_
    overBenford.p = NA}
  
  if (test_BREG == F) {
    BREG <- data.frame(verdict = NA_character_, output_parameters = NA, ideal_parameters = NA, bootstrap_CIs=NA_character_, benford_CIs=NA_character_)
  }
  
  if (test_Kuiper ==F) {
    kuiper.result <- NA_character_
    kuiper.p <- NA}
  
  stats <- data.table(
    region = subtitle, 
    dataset = variable,
    test = test.type,
    total = sum(benford_result$data$data.used),
    range = max(fifelse(benford_result$data$data.used > 0, ceiling(log10(benford_result$data$data.used)), 0)),
    mantissa_difference = mantissa_distance(vector),
    oom_10 = oom_10,
    oom_5 = oom_5,
    oom_naive= oom_naive,
    ES12 = ES12,
    Nigrini = nigrini, 
    MAD = mad,
    Kossovsky = kossovsky,
    SSD = round(ssd, 3),
    Chi.squared = round(chisq, 3),
    n = N,
    df = df,
    p = signif(p_value,3),
    kuiper.result = kuiper.result,
    kuiper.p = kuiper.p,
    overBenford.result = overBenford.result,
    overBenford.p = overBenford.p,
    BREG_verdict = BREG$verdict,
    BREG_output = BREG$output_parameters[1], 
    BREG_ideal = BREG$ideal_parameters[1], 
    BREG_bootstrap_CI = BREG$bootstrap_CIs[1], 
    BREG_benford_CI = BREG$benford_CIs[1]
  )
  
  return(list(plot = p, stats = stats))
}
