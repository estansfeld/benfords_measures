
create_benfords_distribution <- function(n, digits = 1) {
  if (digits == 1) {
    return(sample(x = 1:9, size = n, prob = first_digit_frequency(), replace = T))
  } else {
    return(sample(x = 10:99, size = n, prob = first_two_digit_frequency(), replace = T))
  }
}


create_extreme_distribution <- function(n, digits = 1) {
  if (digits == 1) {
    return(sample(x = 1:9, size = n, prob = first_extreme_frequency(), replace = T))
  } else {
    return(sample(x = 10:99, size = n, prob = first_two_extreme_frequency(), replace = T))
  }
}

first_two_extreme_frequency <- function() {
  # v <- seq(10, 99, 1)
  # p <- log10(1 + 1 / v)
  # return(p)
}

first_extreme_frequency <- function() {
  p <- c(0.511, 0.215, 0.115, 0.068, 0.041, 0.025, 0.015, 0.007, 0.002)
  return(p)
}

first_two_digit_frequency <- function() {
  v <- seq(10, 99, 1)
  p <- log10(1 + 1 / v)
  return(p)
}

first_digit_frequency <- function() {
  v <- seq(1, 9, 1)
  p <- log10(1 + 1 / v)
  return(p)
}

second_digit_frequency <- function() {
  first <- rep(1:9, each = 10)
  second <- seq(0, 9, 1)
  p <- first_two_digit_frequency()
  DT <- data.table(first, second, p)
  DT <- dcast(DT, first~second, value.var = "p")[, -"first"]
  return(as.numeric(unlist(DT[j = lapply(.SD, sum)])))
}


get_leading_digits <- function(number, digits = 1) {
  number <- abs(number)
  # from Carlos Cinelli (benford.analysis package)
  if (digits == 1) {
    return(as.double(trunc((10^((floor(log10(number))*-1)))*number)))
  } else {
    return(as.double(trunc((10^((floor(log10(number))*-1) + 2 - 1))*number))) 
  }
}

get_expected_frequency <- function(digits = 1) {
  # returns frequency table for first and first two digits
  if (digits == 1) {
    return(first_digit_frequency())
  } else {
    return(first_two_digit_frequency()) 
  } 
}

get_extreme_frequency <- function(digits = 1) {
  # returns frequency table for first and first two digits
  if (digits == 1) {
    return(first_extreme_frequency())
  } else {
    return(first_two_extreme_frequency()) 
  } 
}

distribution_table <- function (vector, digits) {
  
  # look for gaps in the actual digit frequency table and fill them with zeros
  if (digits == 1) {
    template <- data.table(digit = as.numeric(seq(from = 1, to = 9)), N = 0) 
    vector <- data.table(get_leading_digits(vector, 1))    
  } else {
    template <- data.table(digit = as.numeric(seq(from = 10, to = 99)), N = 0)  
    vector <- vector[vector >= 10]
    vector <- data.table(get_leading_digits(vector, 2))
  }
  
  setnames(vector, c("digit"))
  missing <- template[!vector, on  = .(digit)][j = digit:=as.character(digit)]
  
  # find distribution of actuals (including the missing zeros)
  actual_dist <- data.table(table(vector))
  setnames(actual_dist, c("digit", "N"))
  
  # bind the actuals and the missing and sort
  combined <- rbindlist(list(actual_dist, missing), use.names = FALSE)[order(digit)]
  return(combined)
}

frequency_table <- function (vector_of_leading_digits, digits) {
  # look for gaps in the actual digit frequency table and fill them with zeros
  if (digits == 1) {
    template <- data.table(digit = as.numeric(seq(from = 1, to = 9)), frequency = 0) 
  } else {
    template <- data.table(digit = as.numeric(seq(from = 10, to = 99)), frequency = 0)  
  }
  vector <- data.table(vector_of_leading_digits)
  setnames(vector, c("digit"))
  missing <- template[!vector, on  = .(digit)][j = digit:=as.character(digit)]
  
  # find frequencies of actuals (including the missing zeros)
  actual_freq <- data.table(prop.table(table(vector)))
  setnames(actual_freq, c("digit", "frequency"))
  
  # bind the actuals and the missing,  sort and calculate cumulative frequency (cumulative distribution function)
  combined <- rbindlist(list(actual_freq, missing), use.names = FALSE)
  combined <- combined[order(digit)][j = cum_freq:=cumsum(frequency)]
  return(combined)
}

robust_magnitude_10 <- function(vector){
  #Find the quantiles (10th and 90th percentiles) of the vector
  quantiles <- quantile(vector[vector!=0], probs = c(0.1, 0.9))
  order_of_magnitude <- log10(quantiles[2] / quantiles[1])
  return(order_of_magnitude)
}

robust_magnitude_5 <- function(vector){
  #Find the quantiles (5th and 95th percentiles) of the vector
  quantiles <- quantile(vector[vector!=0], probs = c(0.05, 0.95))
  order_of_magnitude <- log10(quantiles[2] / quantiles[1])
  return(order_of_magnitude)
}

naive_magnitude <- function(vector){
  # no removal of outliers
  order_of_magnitude <- log10(max(vector) / max(min(vector),1))
  return(order_of_magnitude)
}

ES12 <- function(vector) {
  # actual frequency of digits one and two
  actual <- sum(frequency_table(get_leading_digits(vector,1),1)$frequency[1:2])
  # expected frequency of digits one and two
  expected <- sum(get_expected_frequency(1)[1:2])
  return((actual - expected)*100)
}
