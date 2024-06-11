# FUNCTIONS

# Bayesian
fit_model <- function(formula){
  model <- brm(formula, 
               data = Da2,
               data2 = list(phylogeny = phylogeny),
               control = list(adapt_delta = 0.95),
               chains = 4, cores = 4)  
  return(model)
}

# Folded normal distribution
fold <- function(data, nlevs, colu){
  nl <- 4+nlevs
  ma <- data.frame(sa[, 4:nl])
  seds <- apply(ma, 2, sd)
  mean_mpg <- rep(0, colu)
  lower <- mean_mpg <- rep(0, colu)
  upper <- mean_mpg <- rep(0, colu)
  
  for(i in 1:ncol(ma)){
    ma2 <- as.numeric(ma[, i])
    postfnorm <- stats::dnorm(ma2, 0, seds[i])*2*(seds[i]^2) + ma2*(2*stats::pnorm(ma2, 0, seds[i]) -1)
    mean_mpg[i] <- mean(postfnorm)
    lower[i] <- median_qi(postfnorm, .width = 0.95)[,2]
    upper[i] <- median_qi(postfnorm, .width = 0.95)[,3]
  }
  
  res <- data.frame(mean_mpg, lower, upper)
  colnames(res) <- c("mean", "lower", "upper")
  return(res)
}
