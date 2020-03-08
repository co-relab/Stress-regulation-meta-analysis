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
# Needs specific naming of ES, variances and data on clustering; yi = g.calc, vi = g.var.calc, study, result
rmaCustom <- function(data = NA){
  rmaObject <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result), cluster = data$study)
  rmaObject
}


# Permutation p-curve -----------------------------------------------------

# Evidential value; Permutation p-curve
# Subseting only the effects that are focal for the published study

# No of simulations for the permutation p-curve
nsim <- 5 # Set to 5 just to make code checking/running fast. For the final paper, it needs to be set to at least 1000 and run overnight.

pcurve <- function(data = NA){
  #p-curve data export
  set.seed(123)
  pcurveData <- na.omit(gsub("^.*?: ",": ", x = data[data$focal_variable == 1,]$label[!duplicated.random(data[data$focal_variable == 1,]$study)], replacement = ""))
  #write(pcurveData, "pcurve.data.txt")
  
  #Permutation p-curve
  cols <- 15
  pCurve <- matrix(ncol=cols, nrow=nsim)
  for(i in 1:nsim){
    pCurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data[data$focal_variable == 1,]$label[!duplicated.random(data[data$focal_variable == 1,]$study)], replacement = "")))
  }
  pCurveOut <- data.frame(pCurve)
  as.data.frame(describe(pCurveOut[c(4, 6, 8, 10)]))[,3, drop = FALSE]
}

# 95% prediction interval
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


# Publication bias --------------------------------------------------------

bias <- function(rmaObject = NA, alpha = .05, briefBias = TRUE){
  # Correlation between the ES and precision (SE)
  esPrec <- cor(rmaObject$yi, sqrt(rmaObject$vi), method = "kendall")
  # Small-study effects correction
  # 3-parameter selection model
  fourPSM <- fourPSM.est(rmaObject$yi, rmaObject$vi)
  
  # PET-PEESE
  pp <- pet.peese(rmaObject$yi, rmaObject$vi, rmaObject$mf.r[[1]]$study, rmaObject$mf.r[[1]]$result)
  
  # p-uniform
  puniform.out <- puniform(yi = rmaObject$yi, vi = rmaObject$vi, alpha = alpha, side = "right", method = "P")
  
  if(briefBias == TRUE){
         return(list("4PSM ES estimate" = fourPSM[1, 4],
              "4PSM confidence interval" = c(fourPSM[5, 4], fourPSM[6, 4]),
              "4PSM p-value" = fourPSM[4, 4],
              "Whether PET or PEESE was used" = ifelse(fourPSM$value[4] < .05 & fourPSM$value[1] > 0, "PEESE", "PET"),
              "PET-PEESE ES estimate" = as.numeric(pp[1]),
              "PET-PEESE confidence interval" = as.numeric(c(pp[5], pp[6])),
              "PET-PEESE p-value" = pp[4]))}
  else{
         return(list("4PSM" = fourPSM, 
               "PET-PEESE" = pp,
               "p-uniform" = puniform.out))
  }
}

# Power based on PEESE and 4PSM parameter estimates
powerEst <- function(rmaObject = NA, ni = NA){
  powerPEESE <- NA
  peeseEst <- pet.peese(rmaObject$yi, rmaObject$vi, rmaObject$mf.r[[1]]$study, rmaObject$mf.r[[1]]$result)[1]
  power4PSM <- NA
  fpsmEst <- fourPSM.est(rmaObject$yi, rmaObject$vi)$value[1]

  for(i in 1:rmaObject$k.all){
    powerPEESE[i] <- pwr::pwr.t.test(n = ni, d = peeseEst)$power[i]
    power4PSM[i] <- pwr::pwr.t.test(n = ni, d = fpsmEst)$power[i]
  }
  powerPEESEresult <- median(powerPEESE)*100
  power4PSMresult <- median(power4PSM)*100
  c("Median power for detecting PET-PEESE estimate" = powerPEESEresult, 
    "Median power for detecting 4PSM estimate" = power4PSMresult)
}


# Summary results ---------------------------------------------------------



maResults <- function(rmaObject = NA, data = NA, alpha = .05, briefBias = F){
  list(
    "RMA results" = rmaObject,
    "Prediction interval" = pi95(rmaObject),
    "Heterogeneity" = heterogeneity(rmaObject),
    "p-curve" = pcurve(data),
    "Proportion of significant results" = propSig(data$p),
    "Publication bias" = bias(rmaObject, briefBias = briefBias),
    "Power based on PEESE and 4PSM parameter estimates" = powerEst(rmaObject, ni = data$N))
}


# General Grim Test -------------------------------------------------------

# Code adapted from https://osf.io/scpbz/ , by Nick Brown and 
# https://aurelienallard.netlify.com/post/anaytic-grimmer-possibility-standard-deviations/, by Aurélien Allard

grimTest <- function (n, mean, items = 1, decimals = 2) {
  if(n>10^decimals){
    print("The sample size is too big compared to the precision of the reported mean, it is not possible to apply GRIM.")
  } else {
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
  }
}


# General Grimmer Test ----------------------------------------------------

# Result: -1 = GRIM inconsistent, 0 = GRIMMER inconsistent, 1 = mean & sd consistent
grimmerTest <- function(n, mean, SD, items = 1, decimals_mean = 2, decimals_SD = 2){
  
  if(n>10^decimals_mean){
    print("Reported precision of mean too low, given N")
  }
  
  #Applies the GRIM test, and computes the possible mean.
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

################
# Plots


# Forest plot
# forest(x = data$g.calc, vi = data$g.var.calc,
#        xlim=c(-2.5,3.5),        ### adjust horizontal plot region limits
#        subset=order(data$g.var.calc),        ### order by size of yi
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
# {plot(data$g.var.calc, data$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27),xaxs="i")} else {plot(sqrt(data$g.var.calc), data$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
# abline((if(fourPSM$value[4] < .05 & fourPSM$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
# PP.plot <- recordPlot()
