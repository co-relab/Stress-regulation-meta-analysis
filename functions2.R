
# Filtering options -------------------------------------------------------

# Not wishing to filter by a given criterion is represented by an NA
# mood.opt <- NA
# comp.prim.opt <- NA # 0 == compensatory, 1 = priming
# method.opt <- NA # 1 to 7
# methods.opt.label <- ifelse(is.na(method.opt), NA, c("Physical.Temperature.Manipulation.", "Visual.Verbal.Temperature.Prime.", "Outside.Temperature.", "Temperature.Estimate.", "Subjective.Warmth.Judgment", "Core.Temperature.Measurement.", "Skin.Temperature.Measurement.")[method.opt])
# category.opt <- NA # 1 to 10

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

# No of simulations for the permutation p-curve
nsim <- 5 # Set to 5 just to make code checking/running fast. For the final paper, it needs to be set to at least 1000 and run overnight.



# Meta-analysis -----------------------------------------------------------

# Meta analysis run on a filtered dataset

# Custom robust multivariate RE meta-analytic model
# Needs specific naming of ES, variances and data on clustering; yi = g.calc, vi = g.var.calc, study, result
rmaCustom <- function(data = NA){
  rmaObject <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result), cluster = data$study)
  rmaObject
}

# Evidential value
# Permutation p-curve
# Subseting only the effects that are focal for the published study
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


# Heterogeneity

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

# Proportion of significant results
propSig <- function(p.values = NA){
  as.integer(table(p.values < .05)[2])/length(p.values < .05)
}

# Publication bias

bias <- function(rmaObject = NA, alpha = .05, briefBias = TRUE){
  # Correlation between the ES and precision (SE)
  esPrec <- cor(rmaObject$yi, sqrt(rmaObject$vi), method = "kendall")
  # Small-study effects correction
  # 3-parameter selection model
  threePSM <- threePSM.est(rmaObject$yi, rmaObject$vi)
  
  # PET-PEESE
  pp <- pet.peese(rmaObject$yi, rmaObject$vi, rmaObject$mf.r[[1]]$study, rmaObject$mf.r[[1]]$result)
  
  # p-uniform
  puniform.out <- puniform(yi = rmaObject$yi, vi = rmaObject$vi, alpha = alpha, side = "right", method = "P")
  
  if(briefBias == TRUE){
         return(list("3PSM ES estimate" = threePSM[1, 4],
              "3PSM confidence interval" = c(threePSM[5, 4], threePSM[6, 4]),
              "3PSM p-value" = threePSM[4, 4],
              "Whether PET or PEESE was used" = ifelse(threePSM$value[4] < .05 & threePSM$value[1] > 0, "PEESE", "PET"),
              "PET-PEESE ES estimate" = as.numeric(pp[1]),
              "PET-PEESE confidence interval" = as.numeric(c(pp[5], pp[6])),
              "PET-PEESE p-value" = pp[4]))}
  else{
         return(list("Three PSM" = threePSM, 
               "PET-PEESE" = pp,
               "p-uniform" = puniform.out))
  }
}

# Power based on PEESE and 3PSM parameter estimates
powerEst <- function(rmaObject = NA, ni = NA){
  powerPEESE <- NA
  peeseEst <- pet.peese(rmaObject$yi, rmaObject$vi, rmaObject$mf.r[[1]]$study, rmaObject$mf.r[[1]]$result)[1]
  power3PSM <- NA
  tpsmEst <- threePSM.est(rmaObject$yi, rmaObject$vi)$value[1]

  for(i in 1:rmaObject$k.all){
    powerPEESE[i] <- pwr::pwr.t.test(n = ni, d = peeseEst)$power[i]
    power3PSM[i] <- pwr::pwr.t.test(n = ni, d = tpsmEst)$power[i]
  }
  powerPEESEresult <- median(powerPEESE)*100
  power3PSMresult <- median(power3PSM)*100
  c("Median power for detecting PET-PEESE estimate" = powerPEESEresult, 
    "Median power for detecting 3PSM estimate" = power3PSMresult)
}

maResults <- function(rmaObject = NA, data = NA, alpha = .05, briefBias = F){
  list(
    "RMA results" = rmaObject,
    "Prediction interval" = pi95(rmaObject),
    "Heterogeneity" = heterogeneity(rmaObject),
    "p-curve" = pcurve(data),
    "Proportion of significant results" = propSig(data$p),
    "Publication bias" = bias(rmaObject, briefBias = briefBias),
    "Power based on PEESE and 3PSM parameter estimates" = powerEst(rmaObject, ni = data$N))
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
# if(threePSM$value[4] < .05 & threePSM$value[1] > 0)
# {plot(data$g.var.calc, data$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27),xaxs="i")} else {plot(sqrt(data$g.var.calc), data$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
# abline((if(threePSM$value[4] < .05 & threePSM$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
# PP.plot <- recordPlot()
