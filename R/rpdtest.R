

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
rpdStat <- function(data, probability, lambda = 1, nDim, r) {
  if (is.data.frame(data))
    data <- as.matrix(data)
  if (is.matrix(data)) {
    if (min(dim(data)) == 1L)
      data <- as.vector(data)
    else
      stop("wrong data shape")
  }
  Dnormal <- rnorm(nDim)
  Duniform <- Dnormal / sqrt(sum(Dnormal ^ 2))
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

  #Extensive use of apply family functions to cancel for loops
  eta <- diag(r)
  rStart <- cumsum(c(1, data[-length(data)]))
  rEnd <- cumsum(data)
  # Convert Duniform[rStart:rEnd] to a list
  Duniform_list <-
    lapply(1:r, function(i)
      Duniform[rStart[i]:rEnd[i]])
  # Convert eta[i,] - probability to a list
  eta_prob_list <- lapply(1:r, function(i)
    eta[i,] - probability)
  componentFun <- function(x, y)
    colSums(outer(x, y))
  # Use lapply function to calculate XTheta separately and add them up
  XTheta <-
    Reduce(`+`, lapply(1:r, function(i)
      componentFun(Duniform_list[[i]], eta_prob_list[[i]])))


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
#'
#' @param x the obtained multinomial distribution data.Same data structure
#' as the data parameter in \link{rpdTest}.
#' @param p the probability vector in the null hypothesis. It is necessary to
#' ensure beforehand that the vectors are valid.
#' @param lambda  a control parameter of the statistic calculation,
#' adjusting it will significantly change the final obtained statistic.
#' @param ll  an integer specifying the number of outer loops of the
#' Monte Carlo simulation.
#' @param simNum  an integer specifying the number of inner loops of the
#' Monte Carlo simulation.
#' @param edfLen  an integer that adjusts the number of points used to generate
#' the empirical distribution function used to perform the simulation.
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
#' pVals(d, c(1/2,1/2), ll = 5, simNum =  30, edfLen = 2500, n.cores = 2, nDim = sum(d), r = length(d))
#' }
#' @importFrom stats ecdf pchisq rmultinom
pVals <- function(x,
                  p,
                  lambda = 1,
                  ll,
                  simNum,
                  edfLen,
                  n.cores,
                  nDim,
                  r) {
  if (is.data.frame(x))
    x <- as.matrix(x)
  if (is.matrix(x)) {
    if (min(dim(x)) == 1L)
      data <- as.vector(x)
    else
      stop("wrong data shape")
  }
  if (is.null(n.cores)) {
    n.cores <- detectCores()
    cluster <- makeCluster(n.cores - 1)
  }
  else{
    cluster <- makeCluster(n.cores)
  }

  registerDoParallel(cl = cluster)
  getDoParRegistered()
  #aveP <- vector("numeric", ll*simNum)
  intv <- seq(0, 10, length = 5000)
  pdf <- pchisq(intv, df = r - 1)
  bV <-  rmultinom(ll,nDim, p)
  l <- numeric()
  aveP <- foreach(l = 1:ll, .combine = "cbind")%:%
    foreach(
      k = 1:simNum,
      .combine = "c",
      .packages = c("RPDTest"),
      .errorhandling = "pass"
    )%dopar%{
      rpdStat1 <- sapply(1:edfLen,function(i) rpdStat(bV[,l],p,lambda,nDim,r))
      rpdStat2 <- sapply(1:edfLen,function(i) rpdStat(x,p,lambda,nDim,r))
      edf1 <- ecdf(rpdStat1)
      edf2 <- ecdf(rpdStat2)
      d.1 <- abs(edf1(intv)-pdf)
      d.2 <- abs(edf2(intv)-pdf)
      (sum(abs(d.1 - d.2) < 0.001) + 1) / (edfLen + 1) > 0.5
    }
  if(ll == 1){
    return ((sum(aveP)+1)/(length(aveP) + 1))
  }
  else{
    aveP <- (colSums(aveP)+1)/(nrow(aveP) + 1)
  }
  stopCluster(cl = cluster)
  stopImplicitCluster()
  if(ll>=5){
  (sum(aveP)-max(aveP)-min(aveP))/(ll-2)
  }
  else{
    sum(aveP)/ll
  }
  #median?
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
#' No definite p-value is provided at current stage.
#' However, a p-values in Monte Carlo simulation is available as an option. It
#' executes in parallel within a nested for loop to reduce randomness.
#' In the current version (0.0.1), this feature is still being debugged and improved,
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
#' @param ll an integer specifying the number of outer loops of the
#' Monte Carlo simulation.
#' @param simNum an integer specifying the number of inner loops of the
#' Monte Carlo simulation.
#' @param edfLen an integer that adjusts the number of points used to generate
#' the empirical distribution function used to perform the simulation.
#' @param n.cores an integer used to specify the number of cores used
#' to perform parallel operations. The default is to use the maximum number
#' of cores available to the computer minus one.
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
           ll = 5,
           simNum = 30,
           edfLen = 2500,
           n.cores = NULL) {

    if (0 %in% p)
      stop("invalid 'p(probability)'")
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

    METHOD <- paste("Goodness-of-fit test using Randomized phi-divergence test statistics with parameter lambda = ",lambda,":")
    nDim <- sum(data)
    r <- length(data)
    PARAMETER <- r - 1
    STAT <- rpdStat(data, p, lambda, nDim, r)
    names(STAT) <- "Random phi-divergence"
    names(PARAMETER) <- "corresponding chisq df"
    PVAL <- NULL
    if (sim.pValue) {
      METHOD <- paste(METHOD,"total",ll,"*",simNum,"simulations")
      PVAL <- pVals(data, p, lambda, ll, simNum, edfLen, n.cores, nDim, r)
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
