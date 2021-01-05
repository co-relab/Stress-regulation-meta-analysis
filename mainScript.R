#' ---
#' title: "Mindfulness & Biofeedback meta-analysis"
#' author: "Ivan Ropovik"
#' date: "`r Sys.Date()`"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       code_folding: show
#'       fig_retina: 2
#' always_allow_html: yes
#' ---

#+ setup, include = FALSE
knitr::opts_chunk$set(echo=FALSE, warning = FALSE)
# check with Alessandro if overall RoB was scored correctly

rm(list = ls())

# Settings ----------------------------------------------------------------

# Assumed default pre-post correlation for within-subjects design, .50.
# Here you can perform the sensitivity analysis to determine the impact of the assumed correlation on the overall effect size estimate.
# E.g., for corr = c(.10, .30, .50, .70, 90).
corr <- 0.5

# Assumed constant sampling correlation
rho <- 0.5

# Side argument for the p-uniform* and conditional estimator of PET-PEESE. If the target effect should be in negative values, set to "left", otherwise "right".
side <- "left"

# Define whether to use one-tailed or two-tailed test for PET-PEESE, 3PSM, and p-uniform*.
# Recommended by Stanley (2016) for literature where small sample-size studies are rather the norm.
# Assuming alpha level of .05 for the two-tailed test
test <- "one-tailed"

# No of simulations for the permutation-based bias correction models and p-curve specifically
nIterations <- 5 # Set to 5 just to make code checking/running fast. For the final paper, it needs to be set to at least 1000 and run overnight.
nIterationsPcurve <- 5

# Exclude studies having an overall Risk of Bias score of at least x.
acceptableRiskOfBias <- 2

# Sourcing and data -----------------------------------------------------------------
source("functions.R")
source("pcurvePlotOption.R")
source("esConversion.R")

# GRIM & GRIMMER Test -----------------------------------------------------
grimAndGrimmer(dat)

# Meta-analysis -----------------------------------------------------------

# Subset
dataMind <- dat[dat$strategy == 1 & !is.na(dat$yi),]
dataBio <- dat[dat$strategy == 2 & !is.na(dat$yi),]
dataMind <- dataMind %>% filter_all(any_vars(!is.na(.)))
dataBio <- dataBio %>% filter_all(any_vars(!is.na(.)))

#'# Meta-analysis results
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
namesObjects <- c("Self-administered mindfulness", "Biofeedback")
levels(dat$strategy) <- namesObjects
dataObjects <- list("Mind" = dataMind, "Bio" = dataBio)

rmaObjects <- setNames(lapply(dataObjects, function(x){rmaCustom(x)}), nm = namesObjects)

# Further results
results <- list(NA)
metaResultsPcurve <- list(NA)
for(i in 1:length(rmaObjects)){
  results[[i]] <- maResults(data = dataObjects[[i]], rmaObject = rmaObjects[[i]])
  metaResultsPcurve[[i]] <- metaResultPcurve
}

results <- setNames(results, nm = namesObjects)
metaResultsPcurve <- setNames(metaResultsPcurve, nm = namesObjects)

#+ include = TRUE
#'## Self-administered mindfulness
results$`Self-administered mindfulness`

#'## Biofeedback
results$Biofeedback


# Plots -------------------------------------------------------------------

#+ include = TRUE
#'# Plots
#'

#'## Contour enhanced funnel plot
#'### Mindfulness
dataMind %$% metafor::funnel.default(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei")

#'### Biofeedback
dataBio %$% metafor::funnel.default(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei")


#'## Forest plots
#'### Mindfulness
dataMind %$% forest(x = yi, vi = vi,
                    xlim=c(-2,2),            ### adjust horizontal plot region limits
                    subset=order(vi),        ### order by size of yi
                    slab=NA, annotate=FALSE, ### remove study labels and annotations
                    efac=0,                  ### remove vertical bars at end of CIs
                    pch=19,                  ### changing point symbol to filled circle
                    col="gray40",            ### change color of points/CIs
                    psize=2,                 ### increase point size
                    cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
                    lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Mindfulness")

#'### Biofeedback
dataBio %$% forest(yi, vi, subset=order(vi))

#'## p-curve plots
#'### Mindfulness
quiet(pcurveMod(metaResultsPcurve$`Self-administered mindfulness`, effect.estimation = FALSE, plot = TRUE))

#'### Biofeedback
quiet(pcurveMod(metaResultsPcurve$Biofeedback, effect.estimation = FALSE, plot = TRUE))

#'## PET-PEESE plots
#'### Mindfulness
quiet(petPeese(dataMind))
if(results[[1]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[1]]$`Publication bias`$`4/3PSM`["est"] > 0){
  dataObjects[[1]] %$% plot(vi, yi, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 2, cex = .30, xlim = c(0, .65), xaxs = "i")
} else {
  dataObjects[[1]] %$% plot(sqrt(vi), yi, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 2, cex = .3, xlim = c(0, .65), ylim = c(-2, 1), xaxs = "i")}
abline((if(results[[1]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[1]]$`Publication bias`$`4/3PSM`["est"] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#'### Biofeedback
quiet(petPeese(dataBio))
if(results[[2]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[2]]$`Publication bias`$`4/3PSM`["est"] > 0){
  dataObjects[[2]] %$% plot(vi, yi, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 2, cex = .30, xlim = c(0, .65), xaxs = "i")
} else {
  dataObjects[[2]] %$% plot(sqrt(vi), yi, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 2, cex = .3, xlim = c(0, .65), ylim = c(-2, 1), xaxs = "i")}
abline((if(results[[2]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[2]]$`Publication bias`$`4/3PSM`["est"] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")


# Moderator/sensitivity analyses ------------------------------------------

#'# Moderator/sensitivity analyses
#'## Published status

pubUnpub <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho, smooth_vi = TRUE)
  rmaObject <- rma.mv(yi ~ 0 + factor(published), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  pubUnpub[[i]] <- list("Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
pubUnpub <- setNames(pubUnpub, nm = namesObjects)
pubUnpub

#'## Excluding effects from non-randomized designs
rndNonrnd <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho, smooth_vi = TRUE)
  rmaObject <- rma.mv(yi ~ 0 + factor(researchDesign == 1), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  rndNonrnd[[i]] <- list("Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
rndNonrnd <- setNames(rndNonrnd, nm = namesObjects)
rndNonrnd

#'## Excluding effects due to inconsistent means or SDs
consIncons <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho, smooth_vi = TRUE)
  rmaObject <- rma.mv(yi ~ 0 + factor(as.logical(inconsistenciesCountGRIMMER)), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  consIncons[[i]] <- list("Count of GRIM/GRIMMER inconsistencies" = table(as.logical(dataObjects[[i]]$inconsistenciesCountGRIMMER)), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
consIncons <- setNames(consIncons, nm = namesObjects)
consIncons

#'## Excluding effects due to a high risk of bias

# # Probably need to edit to comply with that is given in the ms: "Following RoB 2 recommendations a study was categorized overall as a high risk of bias if one of two conditions are met: 
# # A) The study scores a  high risk of bias in at least one domain or B) the study is evaluated as having some concerns for more than one domain. 
# # A study was judged as having “some concern” whether it raised some concerns in at least one domain. 
# # Finally a study was assessed as having a low risk of bias if it was judged as having a low risk of bias in all of the five domains. 

highRoB <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho, smooth_vi = TRUE)
  rmaObject <- rma.mv(yi ~ 0 + factor(robOverall > acceptableRiskOfBias), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  highRoB[[i]] <- list("Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
highRoB <- setNames(highRoB, nm = namesObjects)
highRoB

#'## Comparison of strategies

#'### Model without covariates
viMatrixStratComp <- impute_covariance_matrix(dat$vi, cluster = dat$study, r = rho, smooth_vi = TRUE)
rmaObjectStratComp <- rma.mv(yi ~ 0 + factor(strategy), V = viMatrixStratComp, data = dat, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelStratComp <- conf_int(rmaObjectStratComp, vcov = "CR2", test = "z", cluster = dat$study)
list("Model results" = RVEmodelStratComp, "RVE Wald test" = Wald_test(rmaObjectStratComp, constraints = constrain_equal(1:2), vcov = "CR2"))

#'### Model with covariates
#' Controlling for design-related factors that are prognostic w.r.t. the effect sizes (i.e., might vary across moderator categories)
viMatrixStratComp <- impute_covariance_matrix(dat$vi, cluster = dat$study, r = rho, smooth_vi = TRUE)
rmaObjectStratComp <- rma.mv(yi ~ 0 + factor(strategy) + researchDesign + populationType + comparisonGroupType + published + robOverall, V = viMatrixStratComp, data = dat, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelStratComp <- conf_int(rmaObjectStratComp, vcov = "CR2", test = "z", cluster = dat$study)
list("Model results" = RVEmodelStratComp, "RVE Wald test" = Wald_test(rmaObjectStratComp, constraints = constrain_equal(1:2), vcov = "CR2"))
