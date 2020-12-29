# Load libraries
if (!require(lme4)) {install.packages('lme4')}
if (!require(ggplot2)) {install.packages('ggplot2')}
if (!require(knitr)) {install.packages('knitr')}
if (!require(psych)) {install.packages('psych')}
if (!require(puniform)) {install.packages('puniform')}
if (!require(reshape2)) {install.packages('reshape2')}
if (!require(kableExtra)) {install.packages('kableExtra')}
if (!require(lmerTest)) {install.packages('lmerTest')}
if (!require(pwr)) {install.packages('pwr')}

# Meta-analysis -----------------------------------------------------------

# Meta analysis run on a filtered dataset

# Custom robust multivariate RE meta-analytic model
# Needs specific naming of ES, variances and data on clustering; yi = yi, vi = vi, study, result
rmaCustom <- function(data = NA){
  rmaObject <- robust.rma.mv(rma.mv(yi = yi, V = vi, data = data, method = "REML", random = ~ 1|study/result), cluster = data$study)
  rmaObject
}

# 95% prediction interval -------------------------------------------------
pi95 <- function(rmaObject = NA){
  pi95Out <- c("95% PI LB" = round(predict.rma(rmaObject)$cr.lb, 3), "95% PI UB" = round(predict.rma(rmaObject)$cr.ub, 3))
  pi95Out
}

# Heterogeneity -----------------------------------------------------------

heterogeneity <- function(rmaObject = NA){
  
  # Total heterogeneity - tau
  tau <- sqrt(sum(rmaObject$sigma2))
  
  # I^2
  W <- diag(1/rmaObject$vi)
  X <- model.matrix(rmaObject)
  P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
  I2<- 100 * sum(rmaObject$sigma2) / (sum(rmaObject$sigma2) + (rmaObject$k-rmaObject$p)/sum(diag(P)))
  
  # Separate estimates of between- and within-cluster heterogeneity
  BW.hetero <- round(100 * rmaObject$sigma2 / (sum(rmaObject$sigma2) + (rmaObject$k-rmaObject$p)/sum(diag(P))), 2)
  
  studyID <- rmaObject$mf.r[[1]]$study
  resultID <- rmaObject$mf.r[[1]]$result
  resR <- rma.mv(yi = rmaObject$yi, V = rmaObject$vi, struct="UN", random = ~ 1|studyID/resultID)
  resF <- rma.mv(yi = rmaObject$yi, V = rmaObject$vi)
  
  # Jackson's approach to I^2
  JI2 <- round(c(100 * (vcov(resR)[1,1] - vcov(resF)[1,1]) / vcov(resR)[1,1]), 2)
  
  # Intra-class correlation of underlying true effects
  icc <- round(rmaObject$sigma2[1] / sum(rmaObject$sigma2), 2)
  
  c("Tau" = tau,
    "I^2" = I2,
    "Jackson's I^2" = JI2,
    "Between-cluster heterogeneity" = BW.hetero[1],
    "Within-cluster heterogeneity" = BW.hetero[2],
    "ICC" = icc)
}

# Proportion of significant results ---------------------------------------

propSig <- function(p.values = NA){
  as.integer(table(p.values < .05)[2])/length(p.values < .05)
}

# t-test for summary statistics -------------------------------------------

tTestSummary <- function(mean1, mean2, sd1, sd2, n1, n2)
{
  # pooled standard deviation, scaled by the sample sizes
  se <- sqrt((1/n1 + 1/n2) * ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2)/(n1 + n2 - 2)) 
  df <- n1 + n2 - 2
  t <- (mean1 - mean2)/se 
  dat <- c(mean1 - mean2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference in means", "SE", "t-statistic", "p-value")
  return(dat) 
}

# Random selection of effects ---------------------------------------------

# Choose effects from a single study by random (for the purpose of permutation p-curve)
duplicated.random = function(x, incomparables = FALSE, ...)
{
  if ( is.vector(x) )
  {
    permutation = sample(length(x))
    x.perm      = x[permutation]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
  else ( is.matrix(x) )
  {
    permutation = sample(nrow(x))
    x.perm      = x[permutation,]
    result.perm = duplicated(x.perm, incomparables, ...)
    result      = result.perm[order(permutation)]
    return(result)
  }
}

# Summary results ---------------------------------------------------------

maResults <- function(rmaObject = NA, data = NA, alpha = .05, briefBias = F, pcurve = T){
  list(
    "RMA results" = rmaObject,
    "Prediction interval" = pi95(rmaObject),
    "Heterogeneity" = heterogeneity(rmaObject),
    "p-curve" = ifelse(pcurve == TRUE, pcurve(data), paste("p-curve analysis not carried out")),
    "Proportion of significant results" = propSig(data$p),
    "Publication bias" = bias(data, rmaObject, briefBias = briefBias),
    "Power based on PEESE and 4PSM parameter estimates" = powerEst(data))
}

# General Grim Test -------------------------------------------------------

# Code adapted from https://osf.io/scpbz/ , by Nick Brown and 
# https://aurelienallard.netlify.com/post/anaytic-grimmer-possibility-standard-deviations/, by Aurélien Allard

grimTest <- function (n, mean, items = 1, decimals = 2) {
  # if(n>10^decimals){
  #   print("The sample size is too big compared to the precision of the reported mean, it is not possible to apply GRIM.")
  # } else {
  if(items == 0 | is.na(items)){
    return(NA)} else {
    N <- n*items
    dust <- 1e-12
    gMean <- mean
    int <- round(mean * N) # nearest integer; doesn't matter if this rounds up or down
    frac <- int / N
    dif <- abs(mean - frac)
    gran <- ((0.1 ^ decimals) / 2) + dust # allow for rounding errors
    gMean <- round(int / N, decimals)
    consistent <- ifelse(gMean == mean, TRUE, FALSE)
    return(consistent)
    #}
  }
}

# General Grimmer Test ----------------------------------------------------

# Result: -1 = GRIM inconsistent, 0 = GRIMMER inconsistent, 1 = mean & sd consistent
grimmerTest <- function(n, mean, SD, items = 1, decimals_mean = 2, decimals_SD = 2){
  # 
  # if(n>10^decimals_mean){
  #   print("Reported precision of mean too low, given N")
  # } else {
  # 
  #Applies the GRIM test, and computes the possible mean.
  if(items == 0 | is.na(items)){
    return(NA)} else {
  N <- n*items
  sum <- mean*N
  realsum <- round(sum)
  realmean <- realsum/N
  
  # Creates functions to round a number consistently up or down, when the last digit is 5
  
  round_down <- function(number, decimals=2){
    is_five <- number*10^(decimals+1)-floor(number*10^(decimals))*10
    number_rounded <- ifelse(is_five==5, floor(number*10^decimals)/10^decimals, round(number, digits = decimals))
    return(number_rounded)
  }
  
  round_up <- function(number, decimals=2){
    is_five <- number*10^(decimals+1)-floor(number*10^(decimals))*10
    number_rounded <- ifelse(is_five==5, ceiling(number*10^decimals)/10^decimals, round(number, digits = decimals))
    return(number_rounded)
  }
  
  # Applies the GRIM test, to see whether the reconstituted mean is the same as the reported mean (with both down and up rounding)
  
  consistency_down <- round_down(number = realmean, decimals = decimals_mean)==mean
  consistency_up <- round_up(number = realmean, decimals = decimals_mean)==mean
  
  if(consistency_down+consistency_up==0){
    return(-1)
  }
  
  #Computes the lower and upper bounds for the sd.
  
  Lsigma <- ifelse(SD<5/(10^decimals_SD), 0, SD-5/(10^decimals_SD))
  Usigma <- SD+5/(10^decimals_SD)
  
  #Computes the lower and upper bounds for the sum of squares of items.
  
  Lowerbound <- (N-1)*Lsigma^2+N*realmean^2
  Upperbound <- (N-1)*Usigma^2+N*realmean^2
  
  #Checks that there is at least an integer between the lower and upperbound
  
  FirstTest<- ifelse(ceiling(Lowerbound)>floor(Upperbound), FALSE, TRUE)
  
  if(FirstTest==FALSE){
    return(0)
  }
  
  #Takes a vector of all the integers between the lowerbound and upperbound
  
  Possible_Integers <- ceiling(Lowerbound):floor(Upperbound)
  
  #Creates the predicted variance and sd
  
  Predicted_Variance <- (Possible_Integers-N*realmean^2)/(N-1)
  Predicted_SD <- sqrt(Predicted_Variance)
  
  #Computes whether one Predicted_SD matches the SD (trying to round both down and up)
  
  Rounded_SD_down <- round_down(Predicted_SD, decimals_SD)
  Rounded_SD_up <- round_up(Predicted_SD, decimals_SD)
  
  Matches_SD <- Rounded_SD_down==SD | Rounded_SD_up==SD
  
  if(sum(Matches_SD)==0){
    return(0)
  }
  
  #Computes first whether there is any integer between lower and upper bound, and then whether there is 
  #an integer of the correct oddness between the lower and upper bounds.
  oddness <- realsum%%2
  Matches_Oddness <- Possible_Integers%%2==oddness
  Third_Test <- Matches_SD&Matches_Oddness
  return(ifelse(
    sum(Third_Test)==0, 0, 1)
  )
  }
}

# Plots -------------------------------------------------------------------

# Forest plot
# forest(x = data$yi, vi = data$vi,
#        xlim=c(-2.5,3.5),        ### adjust horizontal plot region limits
#        subset=order(data$vi),        ### order by size of yi
#        slab=NA, annotate=FALSE, ### remove study labels and annotations
#        efac=0,                  ### remove vertical bars at end of CIs
#        pch=19,                  ### changing point symbol to filled circle
#        col="gray40",            ### change color of points/CIs
#        psize=3,                 ### increase point size
#        cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
#        lty=c("solid","blank"))  ### remove horizontal line at top of plot
# title("Overall effect")
# addpoly(robma, row = 0, mlab = "", cex = 1, annotate = F)
# forest <- recordPlot()

# Contour enhanced funnel plot
#funnel(robma, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei")

# PET-PEESE plot
# if(fourPSM$value[4] < .05 & fourPSM$value[1] > 0)
# {plot(data$vi, data$yi, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27),xaxs="i")} else {plot(sqrt(data$vi), data$yi, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
# abline((if(fourPSM$value[4] < .05 & fourPSM$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
# PP.plot <- recordPlot()

# Publication bias --------------------------------------------------------

bias <- function(data = NA, rmaObject = NA, alpha = .05, briefBias = TRUE){
  # Correlation between the ES and precision (SE)
  esPrec <- cor(rmaObject$yi, sqrt(rmaObject$vi), method = "kendall")
  # Small-study effects correction
  # 4-parameter selection model
  result <- matrix(ncol = 4, nrow = nsim)
  model <- matrix(ncol = 4, nrow = 24)
  for(i in 1:nsim){
    model <- dat[!duplicated.random(dat$study),] %$% fourPSM.est(yi, vi)
    result[i,] <- c(model$value[1], "ciLB" = model$value[5], "ciUB" = model$value[6], "p-value" = model$value[4])
  }
  colnames(result) <- c("estimate", "ciLB", "ciUB", "p-value")
  fourPSM <- describe(result)[,3, drop = FALSE]
  
  # PET-PEESE
  pp <- pet.peese(rmaObject$yi, rmaObject$vi, rmaObject$mf.r[[1]]$study, rmaObject$mf.r[[1]]$result)
  
  # p-uniform
  puniform.out <- puniform(yi = rmaObject$yi, vi = rmaObject$vi, alpha = alpha, side = "right", method = "P")
  
  if(briefBias == TRUE){
    return(list("4PSM ES estimate" = fourPSM[1, 1],
                "4PSM confidence interval" = c(fourPSM[2, 1], fourPSM[3, 1]),
                "4PSM p-value" = fourPSM[4, 1],
                "Whether PET or PEESE was used" = ifelse(fourPSM[4, 1] < .05 & fourPSM[1, 1] > 0, "PEESE", "PET"),
                "PET-PEESE ES estimate" = as.numeric(pp[1]),
                "PET-PEESE confidence interval" = as.numeric(c(pp[5], pp[6])),
                "PET-PEESE p-value" = pp[4]))}
  else{
    return(list("4PSM" = fourPSM, 
                "PET-PEESE" = pp,
                "p-uniform" = puniform.out))
  }
}

# Power based on PEESE and 4PSM parameter estimates -----------------------

powerEst <- function(data = NA){
  powerPEESE <- NA
  power4PSM <- NA
  peeseEst <- with(data, pet.peese(yi, vi, study, result))[1]
  result <- matrix(ncol = 4, nrow = nsim)
  model <- matrix(ncol = 4, nrow = 24)
  for(i in 1:nsim){
    model <- data[!duplicated.random(data$study),] %$% fourPSM.est(yi, vi)
    result[i,] <- c(model$value[1], "ciLB" = model$value[5], "ciUB" = model$value[6], "p-value" = model$value[4])
  }
  colnames(result) <- c("estimate", "ciLB", "ciUB", "p-value")
  fourPSM <- describe(result)[,3, drop = FALSE]
  fpsmEst <- fourPSM[1, 1]
  
  powerPEESEresult <- median(pwr::pwr.t.test(n = data[!is.na(data$N),]$N, d = peeseEst)$power)
  power4PSMresult <- median(pwr::pwr.t.test(n = data[!is.na(data$N),]$N, d = fpsmEst)$power)
  c("Median power for detecting PET-PEESE estimate" = powerPEESEresult, 
    "Median power for detecting 4PSM estimate" = power4PSMresult)
}

# Multiple-parameter selection models -------------------------------------

# Code adapted from Carter, E. C., Schönbrodt, F. D., Hilgard, J., & Gervais, W. (2018). Correcting for bias in psychology: A comparison of meta-analytic methods. Retrieved from https://osf.io/rf3ys/.
# https://github.com/nicebread/meta-showdown/blob/master/MA-methods/7-Selection%20Models.R
# Return a result data frame either in wide or long format (for the 3PSM output)
returnRes <- function(res, long=TRUE, reduce=TRUE) {
  if (is.null(res)) return(NULL)
  
  # convert all factor columns to characters
  res %>% mutate_if(is.factor, as.character) -> res
  
  if (long==FALSE) {
    # return wide format
    return(res)
  } else {
    # transform to long format
    longRes <- melt(res, id.vars=c("method", "term"))
    if (reduce==TRUE & nrow(res) > 1) {longRes <- longRes %>% filter(!is.na(value)) %>% arrange(method, term, variable)}
    return(longRes)
  }
}

# 3-parameter selection model (3PSM)
# p-value intervals may be re-specified if they contain too few values
if (!require(weightr)) {
  install.packages('weightr')
}

threePSM.est <- function(d, v, min.pvalues=1, long=FALSE) {
  
  w1 <- tryCatch(
    weightfunct(d, v, steps = c(0.025, 1), mods = NULL, weights = NULL, fe = FALSE, table = TRUE),
    error = function(e) NULL
  )
  
  res.NA <- data.frame(
    method = "3PSM",
    term = c("tau2", "b0", "pr.nonsig"),
    estimate = NA,
    std.error = NA,
    statistic = NA,
    p.value = NA,
    conf.low = NA,
    conf.high = NA
  )
  
  if (is.null(w1)) return(returnRes(res.NA))
  
  # If <= 3 p-values in an interval: return NA
  p.table <- table(cut(w1$p, breaks=c(0, .025, 1)))
  if (any(p.table < min.pvalues)) {
    return(returnRes(res.NA))
  } else {
    est <- w1[[2]]$par
    
    # Compute standard errors from hessian
    std.err <- sqrt(abs(diag(solve(w1[[2]]$hessian))))
    
    
    res.wide <- data.frame(
      method = "3PSM",
      term = c("tau2", "b0", "pr.nonsig"),
      estimate = round(est, 4),
      std.error = round(std.err, 4),
      statistic = round(est/std.err, 4),
      p.value = round((pnorm(abs(est/std.err), lower.tail=FALSE)*2), 4),
      conf.low = round((est + qnorm(.025)*std.err), 4),
      conf.high = round((est + qnorm(1-.025)*std.err), 4)
    )
  }
  
  return(returnRes(res.wide))
}

fourPSM.est <- function(d, v, min.pvalues=0, long=TRUE, fallback = FALSE) {	
  w1 <- tryCatch(
    weightfunct(d, v, steps = c(0.025, 0.5, 1), mods = NULL, weights = NULL, fe = FALSE, table = TRUE),
    error = function(e) NULL
  )
  
  res.NA <- data.frame(
    method = "4PSM",
    term = c("tau2", "b0", "pr.nonsig", "pr.opposite"),
    estimate = NA,
    std.error = NA,
    statistic = NA,
    p.value = NA,
    conf.low = NA,
    conf.high = NA
  )
  
  if (is.null(w1)) return(returnRes(res.NA))
  
  # if <= min.pvalues p-values in an interval: return NA
  p.table <- table(cut(w1$p, breaks=c(0, .025, 0.5, 1)))
  if (any(p.table < min.pvalues)) {
    if (fallback==TRUE) {
      return(threePSM.est(d, v, min.pvalues=min.pvalues, long=long))
    } else {
      return(returnRes(res.NA))
    }	  
  } else {
    est <- w1[[2]]$par
    
    # compute standard errors from hessian
    std.err <- sqrt(abs(diag(solve(w1[[2]]$hessian))))
    
    res.wide <- data.frame(
      method = "4PSM",
      term = c("tau2", "b0", "pr.nonsig", "pr.opposite"),
      estimate = round(est, 4),
      std.error = round(std.err, 4),
      statistic = round(est/std.err, 4),
      p.value = round(pnorm(est/std.err, lower.tail=FALSE)*2, 4),
      conf.low = round(est + qnorm(.025)*std.err, 4),
      conf.high = round(est + qnorm(1-.025)*std.err, 4)
    )
  }
  
  return(returnRes(res.wide))
}

# PET-PEESE ---------------------------------------------------------------

#PET-PEESE with 3PSM as the conditional estimator instead of PET

pet.peese <- function(d, v, study, result){
  pet <<- robust.rma.mv(rma.mv(yi = d ~ sqrt(v), V = v, random = ~ 1|study/result, method="REML"), cluster = study)
  pet.out <- round(c(pet$b[1], pet$se[1], pet$zval[1], pet$pval[1], pet$ci.lb[1], pet$ci.ub[1]), 3)
  names(pet.out) <- c("PET estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
  pet.out
  
  peese <<- rma.mv(yi = d ~ v, V = v, random = ~ 1|study/result, method="REML")
  peese.out <- round(c(peese$b[1], peese$se[1], peese$zval[1], peese$pval[1], peese$ci.lb[1], peese$ci.ub[1]), 3)
  names(peese.out) <- c("PEESE estimate", "se", "zval", "pval", "ci.lb", "ci.ub")
  
  fourPSM <- fourPSM.est(d, v) # This is assuming independent effects
  
  ifelse(fourPSM$value[4] < .05 & fourPSM$value[1] > 0,
         return(peese.out),  return(pet.out))
}

# Permutation p-curve -----------------------------------------------------

# Evidential value; Permutation p-curve
# Subseting only the effects that are focal for the published study

pcurve <- function(data = NA){
  #p-curve data export
  set.seed(123)
  pcurveData <- na.omit(gsub("^.*?: ",": ", x = data[data$focalVariable == 1,]$label[!duplicated.random(data[data$focalVariable == 1,]$study)], replacement = ""))
  #write(pcurveData, "pcurve.data.txt")
  
  #Permutation p-curve
  cols <- 15
  pCurve <- matrix(ncol=cols, nrow=nsim)
  for(i in 1:nsim){
    pCurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data[data$focalVariable == 1,]$label[!duplicated.random(data[data$focalVariable == 1,]$study)], replacement = "")))
  }
  pCurveOut <- data.frame(pCurve)
  as.data.frame(describe(pCurveOut[c(4, 6, 8, 10)]))[,3, drop = FALSE]
}

# p-curve app code --------------------------------------------------------

# Adapted p-curve code from p-curve.com
#This is the R Code behind the p-curve app 4.052
#Last updated: 2017 03 17
#Written by Uri Simonsohn (urisohn@gmail.com)

#The app is a single function, pcurve_app() .
#To run it you need to store statistical tests in a textfile with the format from the example, see http://p-curve.com/app4/example.txt
#Then, you just run the function on that file. For exmaple, if you save the file "example2.txt" in the folder "c:\data" then you would run:
#pcurve_app("example2.txt","c://data"). The app generates various text files with the results in table format and saves the figures
# as .PNG files on the same folder as where you put the .txt file
#
#This R Code was written specifically to run in the back-end of the online app and thus it may do things in a way that's not the most intuitive or most efficient
#for a user actually interacting with the R Code. The goal is to exactly replicate what is done on the website, thus the R Code is left intact. The website runs the exact same code
#below.

if (!require(stringr)) {
  install.packages('stringr') #Library to process string variables (text of the entered tests)
}
if (!require(poibin)) {
  install.packages('poibin') #This library has the poisson binomial, the distribution of the sum of binomial with different underlying probabilities
}

pcurve_app=function(x)
{
  #0.1 Set up parameters
  #used to compute the binomial test given that each test has a (slightly) different probability of p<.025 depending on its own noncentral parameter
  #See Hong (2013) - "On computing the distribution function for the Poisson binomial distribution" Computational Statistics and Data Analysis, V59, p.41-51 - http://dx.doi.org/10.1016/j.csda.2012.10.006
  #setwd(dir1)                               #Set as Working Directory the folder on server where temporary files are saved
  #filek=substr(file1,1,nchar(file1)-4)      #filek is the name of the file entered by the user/server, dropping the extension
  x
  #Disable scientific notation
  options(scipen=999)
  
  ##############################################
  #(0) CREATE A FEW FUNCTIONS
  ##############################################
  #Function 1 - functions that find non-centrality parameter for f,chi distributions that gives some level of power
  
  #F-test
  #Note: starting with app 4.0, t() are converted to F() and Z to chisq() to avoid unnecessary repetition of code
  #So we only need code to find ncp for F() and chisq()
  
  getncp.f =function(df1,df2, power)   {
    error = function(ncp_est, power, x, df1,df2) pf(x, df1 = df1, df2=df2, ncp = ncp_est) - (1-power)
    xc=qf(p=.95, df1=df1,df2=df2)
    return(uniroot(error, c(0, 1000), x = xc, df1 = df1,df2=df2, power=power)$root)  }
  
  
  #chisq-test
  getncp.c =function(df, power)   {
    xc=qchisq(p=.95, df=df)
    error = function(ncp_est, power, x, df)      pchisq(x, df = df, ncp = ncp_est) - (1-power)
    return(uniroot(error, c(0, 1000), x = xc, df = df, power=power)$root)   }
  
  #Combine both in single function
  getncp=function(family,df1,df2,power) {
    if (family=="f") ncp=getncp.f(df1=df1,df2=df2,power=power)
    if (family=="c") ncp=getncp.c(df=df1,power=power)
    return(ncp)  }
  
  ###############################################################################
  #Function 2 - percent() : makes a number look like a percentage
  percent <- function(x, digits = 0, format = "f", ...)   {
    paste(formatC(100 * x, format = format, digits = digits, ...), "%", sep = "")
  }
  ###############################################################################
  
  
  ###############################################################################
  #Function 3 - pbound: bound p-values and pp-values by precision of measurement to avoid errors
  pbound=function(p) pmin(pmax(p,2.2e-16),1-2.2e-16)
  
  
  #Function 4 - prop33(pc) - Computes % of p-values that are expected to be smaller than pc,
  #for the tests submitted to p-curve, if power is 33%
  prop33=function(pc)
  {
    #pc: critical  p-value
    
    #Overview:
    #Creates a vector of the same length as the number of tests submitted to p-curve, significant and not,
    #    and computes the proportion of p-values expected to be smaller than {pc} given the d.f.
    #    and outputs the entire vector, with NA values where needed
    
    #F-tests (& thus  t-tests)
    prop=ifelse(family=="f" & p<.05,1-pf(qf(1-pc,df1=df1, df2=df2),df1=df1, df2=df2, ncp=ncp33),NA)
    #Chi2 (& thus Normal)
    prop=ifelse(family=="c" & p<.05,1-pchisq(qchisq(1-pc,df=df1),  df=df1, ncp=ncp33),prop)
    #output it
    prop
  }
  
  #Function 5 Stouffer test for a vector of pp-values
  stouffer=function(pp) sum(qnorm(pp),na.rm=TRUE)/sqrt(sum(!is.na(pp)))
  
  
  
  ###############################################################################
  #(1) PROCESS USER INPUT AND CONVERT IT INTO USABLE R DATA
  ###############################################################################
  
  
  #(1.1) Load data
  
  
  #From file;
  raw = x      #read file on server
  #  raw = scan(file=file1,what="")      #read file on server
  raw=tolower(raw)                                      #lower case
  ktot=length(raw)                                      #count studies
  
  #Create vector that numbers studies 1 to N,includes n.s. studies
  k=seq(from=1,to=length(raw))
  
  #1.2 Parse the entered text into usable statistical results
  #1.3 Create test type indicator
  stat=substring(raw,1,1)          #stat:   t,f,z,c,r
  test=ifelse(stat=="r","t",stat)  #test:   t,f,z,c      (r-->t)
  
  #1.4 Create family to turn t-->F and z-->chi2
  family=test
  family=ifelse(test=="t","f",family)
  family=ifelse(test=="z","c",family)
  
  #Note on terminology:
  #Stat:   t,f,c,z,r  is what the user entered, t,f,c,z,r
  #test:   t,f,c,z    is the test statistic, same as stat but with r-->t
  #family: f,c        converting t-->f and z-->c
  
  #1.5 Find comma,parentheses,equal sign
  par1 =str_locate(raw,"\\(")[,1]         #(  First  parenthesis
  par2 =str_locate(raw,"\\)")[,1]         #)  Second parenthesis
  comma=str_locate(raw,",")[,1]           #,  comma
  eq   =str_locate(raw,"=")[,1]           #=  equal
  
  #1.6 DF for t-tests
  df=as.numeric(ifelse(test=="t",substring(raw,par1+1,par2 -1),NA))             #t(df) later assigned to df2 in  F test with df1=1
  
  #1.7 DF1 for all tests
  #   recall, normal=sqrt(chi(1)) so df1=1 for Normal, same f(1,df)<-t(df)
  df1=as.numeric(ifelse(test=="f",substring(raw,par1+1,comma-1),NA))            #If F test, take value after comma, NA otherwise
  df1=as.numeric(ifelse(test=="z",1,df1))                                       #If Z replace missing value with a 1
  df1=as.numeric(ifelse(test=="t",1,df1))                                       #If t, replace missing value with a 1
  df1=as.numeric(ifelse(test=="c",substring(raw,par1+1,par2 -1),df1))           #If c, replace missing value with value in ()
  
  #1.8 DF2 for F(df1,df2) tests
  df2=as.numeric(ifelse(test=="f",substring(raw,comma+1,par2-1),NA))            #F-test
  df2=as.numeric(ifelse(test=="t",df,df2))                                      #t(df) into the df2 F test
  
  #1.9 Take value after equal sign, the value of the test-statistic, and put it in vector "equal"
  equal=abs(as.numeric(substring(raw,eq+1)))  #if not a r(), take the value after the ="
  
  #1.10  Go from "equal" (the value after = sign) to F or Chi2 value,
  value=ifelse((stat=="f" | stat=="c"),equal,NA)                      #For F and Chi2 test, equal=value
  value=ifelse(stat=="r", (equal/(sqrt((1-equal**2)/df2)))**2,value)  #For correlation, first turn value (r) to t, then square t. (using t=r/sqrt(1-r**2)/DF)
  value=ifelse(stat=="t", equal**2 ,value)                            #For t and Z, square it since t(df)**2=f(1,df) and z**2=chi(1)
  value=ifelse(stat=="z", equal**2 ,value)
  
  
  
  #1.11 Compute p-values
  p=ifelse(family=="f",1-pf(value,df1=df1,df2=df2),NA)
  p=ifelse(family=="c",1-pchisq(value,df=df1),p)
  p=pbound(p)  #Bound it to level of precision, see function 3 above
  
  #1.12 Count  studies
  #ktot is all studies
  ksig= sum(p<.05,na.rm=TRUE)     #significant studies
  khalf=sum(p<.025,na.rm=TRUE)    #half p-curve studies
  
  ###############################################################################
  #(2) COMPUTE PP-VALUES
  ##############################################################################
  
  #2.1 Right Skew, Full p-curve
  ppr=as.numeric(ifelse(p<.05,20*p,NA))            #If p<.05, ppr is 1/alpha*p-value, so 20*pvalue, otherwise missing.
  ppr=pbound(ppr)                                  #apply pbound function to avoid 0
  
  
  #2.2 Right Skew, half p-curve
  ppr.half=as.numeric(ifelse(p<.025,40*p,NA))    #If p<.05, ppr is 40*pvalue, otherwise missing.
  ppr.half=pbound(ppr.half)
  
  #2.3 Power of 33%
  #2.3.1 NCP for  f,c distributions
  # NCP33 (noncentrality parameter giving each test in p-curve 33% power given the d.f. of the test)
  ncp33=mapply(getncp,df1=df1,df2=df2,power=1/3,family=family)  #See function 1 above
  
  #2.3.2 Full-p-curve
  #Using the ncp33 compute pp33
  pp33=ifelse(family=="f" & p<.05,3*(pf(value, df1=df1, df2=df2, ncp=ncp33)-2/3),NA)
  pp33=ifelse(family=="c" & p<.05,3*(pchisq(value, df=df1, ncp=ncp33)-2/3),pp33)
  pp33=pbound(pp33)
  
  #2.3.3 HALF-p-curve
  #Share of p-values expected to be p<.025 if 33% power (using Function 4 from above, prop33() )
  prop25=3*prop33(.025)
  prop25.sig=prop25[p<.05]
  
  
  #Compute pp-values for the half
  pp33.half=ifelse(family=="f" & p<.025, (1/prop25)*(    pf(value,df1=df1,df2=df2,ncp=ncp33)-(1-prop25)),NA)
  pp33.half=ifelse(family=="c" & p<.025, (1/prop25)*(pchisq(value,df=df1,         ncp=ncp33)-(1-prop25)),pp33.half)
  pp33.half=pbound(pp33.half)
  
  
  ###############################################################################
  #(3) INFERENCE - STOUFFER & BINOMIAL
  ##############################################################################
  
  #3.1 Convert pp-values to Z scores, using Stouffer function above
  Zppr =     stouffer(ppr)            #right skew  - this is a Z value from Stouffer's test
  Zpp33=     stouffer(pp33)           #33% - idem
  Zppr.half= stouffer(ppr.half)       #right skew, half p-curve - idem
  Zpp33.half=stouffer(pp33.half)      #33% skew, half p-curve - idem
  
  #3.2 Overall p-values from Stouffer test
  p.Zppr =pnorm(Zppr)
  p.Zpp33=pnorm(Zpp33)
  p.Zppr.half =pnorm(Zppr.half)
  p.Zpp33.half=pnorm(Zpp33.half)
  
  #3.3 Save results to file
  main.results=as.numeric(c(ksig, khalf, Zppr, p.Zppr, Zpp33, p.Zpp33, Zppr.half, p.Zppr.half, Zpp33.half, p.Zpp33.half))
  
  #3.4 BINOMIAL
  #Observed share of p<.025
  prop25.obs=sum(p<.025)/sum(p<.05)
  #3.4.1 Flat null
  binom.r=1-pbinom(q=prop25.obs*ksig- 1, p=.5, size=ksig)     #The binomial in R computes the probability of x<=xo. We want prob(x>=x0) so we subtract one from x, and 1-prob()
  #3.4.2 Power of 33% null
  binom.33=ppoibin(kk=prop25.obs*ksig,pp=prop25[p<.05])
  
  #syntax for ppoibin():
  #   kk: is the proportion of succeses, a scalar, in this case, the share of p<.025
  #   pp: is the probabilty of success for each attempt, the number of attempts is determined
  #    by the length of the vector. For example ppoibin(kk=0,pp=c(.5,.5,.5)) is .125,
  #    if there are three attempts, each with probability .5, the odds of getting 0 succeses is .125
  #     ppoibin(kk=1,pp=c(1,.75)), in turn computes the probability of getting 1 success
  #     when one has a 100% of success, and the other 75%, and the solution is .25, since
  #     the first one succeeds for sure and the second would need to fail, with 25% chance.
  
  
  #3.4.3  Save binomial results
  binomial=c(binom.r, binom.33)
  
  ################################################
  #(4) POWER ESTIMATE
  ################################################
  
  #4.1 Function powerfit(power_est) - Returns the Stouffer Z of observing at least as right skewed a p-curve if  power=power_est
  #if Z=0, power_est is the best fit (p=50%).
  #if Z<0 the best fit is <power_est,
  #if Z>0 the best fit is >power_est
  powerfit=function(power_est)
  {
    #4.1.1 Get the implied non-centrality parameters (ncps) that give power_est to each test submitted to p-curve
    ncp_est=mapply(getncp,df1=df1,df2=df2,power=power_est,family=family)
    #4.1.2 Compute implied pp-values from those ncps_est,
    pp_est=ifelse(family=="f" & p<.05,(pf(value,df1=df1,df2=df2,ncp=ncp_est)-(1-power_est))/power_est,NA)
    pp_est=ifelse(family=="c" & p<.05,(pchisq(value,df=df1,ncp=ncp_est)-(1-power_est))/power_est,pp_est)
    pp_est=pbound(pp_est)
    #4.1.3 Aggregate pp-values for null that power=power_est via Stouffer
    return(stouffer(pp_est))   #This is a z score, so powerfit is expressed as the resulting Z score.
  }
  
  
  #4.2 COMPUTE FIT FOR EACH POWER for 5.1%, AND THEN 6-99%, AND PLOT IT. With power=5% boundary condition lead to errors
  #This becomes the diagnostic plot and gives us the best estimate, within 1%, of power.
  
  # Fit will be evaluated at every possible value of power between 5.1% and 99% in steps of 1%, stored in fit()
  fit=c()                                          #Create empty vector
  fit=abs(powerfit(.051))                      #Start it eavaluting fit of 5.1% power
  for (i in 6:99)   fit=c(fit,abs(powerfit(i/100))) #Now do 6% to 99%
  
  # Find the minimum
  #which ith power level considered leads to best estimate
  mini=match(min(fit,na.rm=TRUE),fit)
  #convert that into the power level, the ith value considered is (5+ith)/100
  hat=(mini+4)/100
  
  #4.3 Confidence interval for power estimate
  #4.3.1 Function get.power_pct(pct)
  get.power_pct =function(pct)   {
    #Function that finds power that gives p-value=pct for the Stouffer test
    #for example, get.power_pct(.5) returns the level of power that leads to p=.5  for the stouffer test.
    #half the time we would see p-curves more right skewed than the one we see, and half the time
    #less right-skewed, if the true power were that get.power_pct(.5). So it is the median estimate of power
    #similarliy, get.power_pct(.1) gives the 10th percentile estimate of power...
    #Obtain the normalized equivalent of pct, e.g., for 5% it is -1.64, for 95% it is 1.64
    z=qnorm(pct)  #convert to z because powerfit() outputs a z-score.
    #Quantify gap between computed p-value and desired pct
    error = function(power_est, z)  powerfit(power_est) - z
    #Find the value of power that makes that gap zero, (root)
    return(uniroot(error, c(.0501, .99),z)$root)   }
  
  #4.3.2 Boundary conditions (when the end of the ci=5% or 99% we cannot use root to find it,
  #use boundary value instead)
  
  #Boundary conditions
  p.power.05=pnorm(powerfit(.051)) #Proability p-curve would be at least at right-skewed if power=.051
  p.power.99=pnorm(powerfit(.99))  #Proability p-curve would be at least at right-skewed if power=.99
  
  #4.3.3 Find lower end of ci
  #Low boundary condition? If cannot reject 5% power, don't look for lower levels, use 5% as the end
  if (p.power.05<=.95) power.ci.lb=.05
  #High boundary condition? If we reject 99%, from below dont look for higher power, use 99% as the low end
  if (p.power.99>=.95) power.ci.lb=.99
  #If low bound is higher than 5.1% power and lower than 99% power, estimate it, find interior solution
  if (p.power.05>.95 && p.power.99<.95)  power.ci.lb=get.power_pct(.95)
  
  
  #4.3.4 Higher end of CI
  #If we reject 5% power from below, 5% is above the confidence interval, use 5% as the upper end of the confidence interval
  if (p.power.05<=.05) power.ci.ub=.05
  #If we do not reject that 99% power, don't look higher, use 99% as the higher end
  if (p.power.99>=.05) power.ci.ub=.99
  #If the the upper bound is between 5% and 99%, find it
  if (p.power.05>.05 && p.power.99<.05) power.ci.ub=get.power_pct(.05)
  
  
  #4.4 Save power fit estiate and ci
  power_results=c(power.ci.lb,hat,power.ci.ub)
  c(main.results, binomial, power_results)
  
  #Note, I use hat as the estimate of power, with powerfit(.5) we could get a more precise best fitting
  #level of power than the minimum in the figure above between .051 and .99, hat, but more precision than 1% in power is not informative.
}
