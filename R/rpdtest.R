

library(parallel)
library(doParallel)
library(foreach)




#' Randomized phi-divergence test: statistic part
#'
#' This is one of the auxiliary functions used to execute the rpdTest function.
#' This function calculates the statistic for a single  Randomized
#' phi-divergence test. Users generally do not need to call this function
#' except for testing purposes.
#'
#' @param data the same data structure that provided in \link{rpdTest}.
#' @param probability the same numeric vector that provided in \link{rpdTest}.
#' @param lambda the same parameter that provided in \link{rpdTest}.
#' @param nDim an integer indicating the dimension of the uniformly
#' distributed vectors generated during the computation of the statistic.
#' It is equal to the number of experiments for the multinomial distribution.
#' @param r an integer indicating the dimension of the data parameter.
#' It is equal to the number of possible outcomes of the multinomial distribution.
#' @param random.state a numeric that controls the randomness of the samples used
#' when generating uniformly distributed random vector on the n-sphere.
#'
#' @return a numeric value that reflects the statistic obtained after
#' an execution of rpdTest at that time.
#' @export
#'
#' @examples
#' d <- c(20,40)
#' #The next line is equivalent to rpdTest(d)$statistic
#'
#' rpdStat(d, c(1/2,1/2), nDim = sum(d), r = length(d))
#'
#' @importFrom stats D rnorm
rpdStat <- function(data, probability, lambda = 1, nDim, r,random.state = NULL) {
  if (is.data.frame(data))
    data <- as.matrix(data)
  if (is.matrix(data)) {
    if (min(dim(data)) == 1L)
      data <- as.vector(data)
    else
      stop("wrong data shape")
  }
  set.seed(random.state)
  Dnormal <- rnorm(nDim)
  Duniform <- Dnormal / sqrt(sum(Dnormal ^ 2))
  set.seed(NULL)
  phiF <- function(u, l = lambda) {
    if (l == 0) {
      return(u * log(u) - u + 1)
    }
    else if (l == -1) {
      return(-log(u) + u - 1)
    }
    else{
      return((u ^ (l + 1) - (l + 1) * (u - 1) - 1) / (l * (l + 1)))
    }
  }

  startIndex <- 1
  eta <- matrix(0, nrow = nDim, ncol = r)
  for (i in 1:r) {
    if(data[i] == 0){
      startIndex <- startIndex + data[i]
      next
    }
    # Create one-hot encoded vector for each element in data
    eta[startIndex:(startIndex + data[i] - 1), i] <- 1
    startIndex <- startIndex + data[i]
  }
  XTheta <- c(Duniform %*% (eta - probability))
  if (lambda == 0) {
    DDPhiF <-
      D(D(expression(u * log(u) - u + 1) ,   "u"),   "u")
  } else if (lambda == -1) {
    DDPhiF <-
      D(D(expression(-log(u) + u - 1) ,   "u"),   "u")
  } else{
    DDPhiF <-
      D(D(expression((
        u ^ (l + 1) - (l + 1) * (u - 1) - 1
      ) / (l * (
        l + 1
      ))) ,   "u"),   "u")
  }

  u <- 1
  l <- lambda
  DDPhiF1 <- eval(DDPhiF)
  #Use vectorized operations to compute rpdStat
  sum((2 * nDim / DDPhiF1) * probability * phiF(1 + XTheta / (sqrt(nDim) * probability)))
}


#' Randomized phi-divergence test: simulated p-value part
#'
#' This is one of the auxiliary functions used to execute the rpdTest function.
#' This function can be used to calculate p-values based on Monte Carlo simulation.
#' Users generally do not need to call this function except for testing purposes.
#' For more detailed description one can find in\link{rpdTest}.
#'
#' @param x the obtained multinomial distribution data.Same data structure
#' as the data parameter in \link{rpdTest}.
#' @param p the probability vector in the null hypothesis. It is necessary to
#' ensure beforehand that the vectors are valid.
#' @param lambda  a control parameter of the statistic calculation,
#' adjusting it will significantly change the final obtained statistic.
#' @param B an integer specifying the number of simulation data on the expected
#' null distribution (p) of the Monte Carlo simulation.
#' @param z an integer specifying the number by which to divide
#' the observation data group in a Monte Carlo simulation.
#' @param rs an integer that adjusts the number of statistics calculated in simulation.
#' @param n.cores  an integer used to specify the number of cores used
#' to perform parallel operations. The default is to use the maximum number
#' of cores available to the computer minus one.
#' @param nDim an integer indicating the dimension of the uniformly
#' distributed vectors generated during the computation of the statistic.
#' It is equal to the number of experiments for the multinomial distribution.
#' @param r an integer indicating the dimension of the data parameter.
#' It is equal to the number of possible outcomes of the multinomial distribution.
#'
#' @return an numeric value indicating simulated p-value.
#' @export
#'
#' @examples
#' d <- c(20,40)
#' #The next line is equivalent to rpdTest(d,sim.pValue = TRUE,n.cores = 2)$p.value
#' #It usually takes 1-2 minutes to perform this calculation process
#' \donttest{
#' pVals(d, c(1/2,1/2), B = 200, z = 40, rs = 1250, n.cores = 2, nDim = sum(d), r = length(d))
#' }
#' @importFrom stats ecdf pchisq rmultinom ks.test
pVals <- function(x, p, lambda = 1, B = 200, z = 40, rs = 1250,
                  n.cores, nDim, r) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1L)
      data <- as.vector(x)
    else
      stop("wrong data shape")
  }

  multd0 <- rmultinom(1,nDim,p)
  multd <- rmultinom(B,nDim,p)

  # set core number
  if (is.null(n.cores)) {
    n.cores <- detectCores()
    cluster <- makeCluster(n.cores - 1)
  }
  else{
    cluster <- makeCluster(n.cores)
  }
  registerDoParallel(cluster)
  getDoParRegistered()
  # initialize vector
  pv <- rep(NA, B)
  pv2 <- rep(NA, B)
  kl <- 0
  results <- foreach (kl = 1:200,.packages = c("RPDTest")) %dopar% {
    stat <- sapply(1:rs, function(i) rpdTest(multd0, p, lambda, random.state = i)$statistic)
    stat2 <- sapply(1:rs, function(i) rpdTest(multd[, kl], p, lambda, random.state = i)$statistic)
    stat3 <- sapply(1:rs, function(i) rpdTest(x, p, lambda, random.state = i)$statistic)
    list(pv = ks.test(stat, stat2)$statistic, pv2 = ks.test(stat3, stat2)$statistic)
  }

  # summarize result
  pv <- sapply(results, function(x) x$pv)
  pv2 <- sapply(results, function(x) x$pv2)

  stopCluster(cluster)
  stopImplicitCluster()

  splitVec <- split(pv, cut(seq_len(length(pv)), z, equal = FALSE))

  # Calculate the mean of each part
  mean0 <- lapply(splitVec, mean)

  benchmark <- mean(unlist(mean0))

  splitVec2 <- split(pv2, cut(seq_len(length(pv2)), z, equal = FALSE))

  # Calculate the mean of each part
  meanValues <- lapply(splitVec2, mean)

  (sum(meanValues <= benchmark) + 1)/ (z + 1)
}


#' Randomized phi-divergence test
#'
#' The most important part of the package:
#' a function for performing hypothesis testing ----
#' An analogue of Chi-square Goodness-of-Fit Test.
#' Accept a vector, matrix or a \link{data.frame} as observed data.
#' Then obtain a specific Randomized phi-divergence statistic,
#' which is computed based on a uniformly distributed random vector
#' on the n-sphere. This random vector is uniquely generated at runtime.
#' However, a p-values in Monte Carlo simulation is available as an option. It
#' executes in parallel way, comparing the empirical distribution function. In specific,
#' it simulates data under the null hypothesis and compares it to the observed data.
#' It generates B datasets based on the expected null distribution (p) and
#' the observed control data (v0). For each simulated dataset and the observed
#' data and v0, rs statistics are computed using different random seeds.
#' The Kolmogorov-Smirnov statistic is used to compare the distributions of the simulated and
#' observed data and the simulated and control data. We get B K-S statistics in both
#' observed data group and control data group.
#' The function then calculates a p-value based on how often the within-group mean of
#' the Kolmogorov-Smirnov statistic after dividing the observed data group into z groups
#' is more extreme than the mean of the statistic observed for the control vector group.
#' In the current version (0.0.2), this feature is still being debugged and improved,
#' so this option is not enabled by default.
#'
#' @param data a one-dimensional vector or matrix of this shape (data.frame)
#' in which observation data for some multinomial distribution are stored.
#' @param p the probability vector in the null hypothesis. Will check the
#' validity of this vector.
#' @param lambda a control parameter of the statistic calculation,
#' adjusting it will significantly change the final obtained statistic.
#' @param sim.pValue a logical variable. It decides whether to compute p-values
#' in Monte Carlo simulation.
#' @param B an integer specifying the number of simulation data on the expected
#' null distribution (p) of the Monte Carlo simulation.
#' @param z an integer specifying the number by which to divide
#' the observation data group in a Monte Carlo simulation.
#' @param rs an integer that adjusts the number of statistics calculated in simulation.
#' @param n.cores an integer used to specify the number of cores used
#' to perform parallel operations. The default is to use the maximum number
#' of cores available to the computer minus one.
#' @param random.state a numeric that controls the randomness of the samples used
#' when generating uniformly distributed random vector on the n-sphere.
#' @return standard list object with class "htest".
#' @export
#'
#' @examples
#' d <- rmultinom(1, 120, c(1/4,3/4))
#' #following will only obtain statistic
#' rpdTest(d)
#' #following will obtain sim.p.value either. You can also specify the number of
#' #cores to use. For example, two:
#' #It usually takes 1-2 minutes to perform this calculation process
#' \donttest{
#' rpdTest(d,sim.pValue = TRUE,n.cores = 2)
#' }
#' @import parallel doParallel foreach
#'
rpdTest <-
  function(data,
           p = rep(1 / length(data), length(data)),
           lambda = 1,
           sim.pValue = FALSE,
           B = 200, z = 40, rs = 1250, n.cores = NULL, random.state = NULL) {

    if (0 %in% p)
      warning("invalid 'p(probability)', statistic = NaN")
    if (length(data) == 1L)
      stop("'data' must at least have 2 elements")
    if (length(data) != length(p))
      stop("'data' and 'p(probability)' must have the same number of elements")
    if (any(p < 0))
      stop("probabilities must be non-negative.")
    if (abs(sum(p) - 1) > sqrt(.Machine$double.eps))
      stop("probabilities must sum to 1")

    if (is.data.frame(data))
      data <- as.matrix(data)
    if (is.matrix(data)) {
      if (min(dim(data)) == 1L)
        data <- as.vector(data)
      else
        stop("wrong data shape")
    }

    METHOD <- paste("Goodness-of-fit test using Randomized phi-divergence test statistics with parameter lambda = ",lambda)
    nDim <- sum(data)
    r <- length(data)
    PARAMETER <- r - 1
    STAT <- rpdStat(data, p, lambda, nDim, r,random.state)
    names(STAT) <- "Random phi-divergence"
    names(PARAMETER) <- "chisq df"
    PVAL <- NULL
    if (sim.pValue) {
      METHOD <- paste(METHOD," with simulated p-value (based on ",B," simulations and ",rs ," replicates)")
      PVAL <- pVals(data, p, lambda, B, z, rs, n.cores, nDim, r)
    }
    result <- structure(
      list(
        statistic = STAT,
        parameter = PARAMETER,
        p.value = PVAL,
        data.name = deparse(substitute(data)),
        method = METHOD,
        class = "htest"
      )
    )
    class(result) <- "htest"
    result
  }
