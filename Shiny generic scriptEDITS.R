
# Filtering options -------------------------------------------------------

# Not wishing to filter by a given criterion is represented by an NA
mood.opt <- NA
comp.prim.opt <- NA # 0 == compensatory, 1 = priming
method.opt <- NA # 1 to 7
methods.opt.label <- ifelse(is.na(method.opt), NA, c("Physical.Temperature.Manipulation.", "Visual.Verbal.Temperature.Prime.", "Outside.Temperature.", "Temperature.Estimate.", "Subjective.Warmth.Judgment", "Core.Temperature.Measurement.", "Skin.Temperature.Measurement.")[method.opt])
category.opt <- NA # 1 to 10

# Data preparation --------------------------------------------------------

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
cor <- 0.5
corr <- 0.5

# Source scripts
source("functions.R")
source("ES_conversion.R")

# No of simulations for the permutation p-curve
nsim <- 5 # Set to 5 just to make code checking/running fast. For the final paper, it needs to be set to at least 1000 and run overnight.

# Outliers
# Excluding 3 outliers based on Baujat analysis, first 2 due to improbably big effect sizes, the last because it also excerts a big influence on the MA model due to combination of huge ES and small variance  (d = 3.23 with n=60, 2.55 with n=63, 1.72 with n=224, respectively).
dat[c(193, 213, 343),] <- NA

dat$et.temp <- dat$Effect.Type..Neither..0..Compensatory..1..Priming..2..Both..3.
dat$effect.type <- ifelse((dat$et.temp == 0) |  (dat$et.temp == 3), yes = NA, no = dat$et.temp)

# Recode "Y" and "N" dummies into a single "method" variable
d <- dat
dxx <- d[,44:50]
dxx <- apply(dxx, 2, as.factor)
dxx[dxx != "Y" & dxx != "N"] <- NA
dxx[dxx == "Y"] <- 1
dxx[dxx == "N"] <- 0
d$method <- colnames(dxx)[apply(dxx, 1, match, x = "1")]

# Recode "Y" and "N" dummies into 1 and 0 for Category variables (one result may be associated with several Categories)
dyy <- d[,52:61]
dyy <- apply(dyy, 2, as.factor)
dyy[dyy != "Y" & dyy != "N"] <- NA
dyy[dyy == "Y"] <- 1
dyy[dyy == "N"] <- 0
d[,52:61] <- dyy


# Filtering engine --------------------------------------------------------

data <-
  if(is.na(mood.opt)){
  if(is.na(comp.prim.opt)){d[d$PA.NA. != "Y",]} else(
    if(comp.prim.opt == 0){d[d$PA.NA. != "Y",] %>% filter(effect.type == 1)} else {d[d$PA.NA. != "Y",] %>% filter(effect.type == 2)})
  } else(d %>% filter(PA.NA. == "Y"))

data <-
  if(is.na(method.opt)){data} else(data %>% filter(method == methods.opt.label))

data <-
  if(is.na(category.opt)){data} else(data[unlist(lapply(data[,52:61], function(x)x == 1)[category.opt], use.names = F),])


# Meta-analysis -----------------------------------------------------------

# Meta analysis run on a filtered dataset

# Evidential value
# Permutation p-curve

pcurve <- function(data = NA){
  #p-curve data export
  set.seed(123)
  pcurveData <- na.omit(gsub("^.*?: ",": ", x = data$label[!duplicated.random(data$study)], replacement = ""))
  #write(pcurveData, "pcurve.data.txt")
  
  #Permutation p-curve
  cols <- 15
  pCurve <- matrix(ncol=cols, nrow=nsim)
  for(i in 1:nsim){
    pCurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data$label[!duplicated.random(data$study)], replacement = "")))
  }
  pCurveOut <- data.frame(pCurve)
  as.data.frame(describe(pCurveOut[c(4, 6, 8, 10)]))[,3, drop = FALSE]
}

pcurve(data = data)

# Now filter out results that don't go into the MA
data <- data[!is.na(data$g.calc),]

# Meta-analysis
ma <- rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result)
robma <- robust.rma.mv(x = ma, cluster = data$study)

# 95% prediction interval
pi95 <- function(rmaObject = NA){
  pi95Out <- c("95% PI LB" = round(predict.rma(robma)$cr.lb, 3), "95% PI UB" = round(predict.rma(robma)$cr.ub, 3))
  pi95Out
}

pi95(robma)


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

heterogeneity(robma)


# Proportion of significant results
sig.prop <- as.integer(table(data$p < .05)[2])/length(data$p < .05)

# Publication bias

bias <- function(rmaObject = NA, alpha = .05, brief = F){
  # Correlation between the ES and precision (SE)
  esPrec <- cor(rmaObject$yi, sqrt(rmaObject$vi), method = "kendall")
  # Small-study effects correction
  # 3-parameter selection model
  threePSM <- threePSM.est(rmaObject$yi, rmaObject$vi)
  
  # PET-PEESE
  pp <- pet.peese(rmaObject$yi, rmaObject$vi, rmaObject$mf.r[[1]]$study, rmaObject$mf.r[[1]]$result)
  
  # p-uniform
  puniform.out <- puniform(yi = rmaObject$yi, vi = rmaObject$vi, alpha = alpha, side = "right", method = "P")
  
  if(brief == TRUE){
         return(list("3PSM ES estimate" = threePSM[1, 4],
              "3PSM confidence interval" = c(threePSM[5, 4], threePSM[6, 4]),
              "3PSM p-value" = threePSM[4, 4],
              "Whether PET or PEESE was used" = ifelse(ThreePSM$value[4] < .05 & ThreePSM$value[1] > 0, "PEESE", "PET"),
              "PET-PEESE ES estimate" = as.numeric(pp[1]),
              "PET-PEESE confidence interval" = as.numeric(c(pp[5], pp[6])),
              "PET-PEESE p-value" = pp[4]))}
  else{
         return(list("Three PSM" = threePSM, 
               "PET-PEESE" = pp,
               "p-uniform" = puniform.out))
  }
}

bias(robma, brief = F)

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

powerEst(robma, ni = data$N)



maResults <- function(rmaObject = NA, data = data, alpha = .05, briefBias = T){
  list(
    "Prediction interval" = pi95(rmaObject),
    "Heterogeneity" = heterogeneity(rmaObject),
    "p-curve" = pcurve(data),
    "Proportion of significant results" = as.integer(table(data$p < .05)[2])/length(data$p < .05),
    "Publication bias" = bias(rmaObject, brief = F),
    "Power based on PEESE and 3PSM parameter estimates" = powerEst(rmaObject, ni = data$N))}

maResults(rmaObject = robma, data = data)


################
# Plots


# Forest plot
forest(x = data$g.calc, vi = data$g.var.calc,
       xlim=c(-2.5,3.5),        ### adjust horizontal plot region limits
       subset=order(data$g.var.calc),        ### order by size of yi
       slab=NA, annotate=FALSE, ### remove study labels and annotations
       efac=0,                  ### remove vertical bars at end of CIs
       pch=19,                  ### changing point symbol to filled circle
       col="gray40",            ### change color of points/CIs
       psize=3,                 ### increase point size
       cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
       lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Overall effect")
addpoly(robma, row = 0, mlab = "", cex = 1, annotate = F)
forest <- recordPlot()

# Contour enhanced funnel plot
funnel(robma, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei")



# PET-PEESE plot
if(threePSM$value[4] < .05 & threePSM$value[1] > 0)
{plot(data$g.var.calc, data$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27),xaxs="i")} else {plot(sqrt(data$g.var.calc), data$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(threePSM$value[4] < .05 & threePSM$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
PP.plot <- recordPlot()
