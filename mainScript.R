#' ---
#' title: "Mindfulness & Biofeedback meta-analysis"
#' author: "Ivan Ropovik & Alessandro Sparacio"
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
# subgroup analysis A: if !is.na(affect) | !is.na(stressComponentType)) = 1; if !is.na(affectiveConsequencesStress) = 2
# subgroup analysis B: if (stressComponentType %in% c(1:4) | !is.na(affect)) ~ 1; stressComponentType = 5 ~ 2; stressComponentType = 6 ~ 3

rm(list = ls())

# Settings ----------------------------------------------------------------

# Assumed default pre-post correlation for within-subjects design, .50.
# Here you can perform the sensitivity analysis to determine the impact of the assumed correlation on the overall effect size estimate.
# E.g., for corr = c(.10, .30, .50, .70, 90).
rmCor <- 0.5

# Assumed constant sampling correlation
rho <- 0.5

# Side argument for the p-uniform* and conditional estimator of PET-PEESE. If the target effect should be in negative values, set to "left", otherwise "right".
side <- "right"

# Define whether to use one-tailed or two-tailed test for PET-PEESE, 3PSM, and p-uniform*.
# Recommended by Stanley (2016) for literature where small sample-size studies are rather the norm.
# Assuming alpha level of .05 for the two-tailed test
test <- "one-tailed"

# No of simulations for the permutation-based bias correction models and p-curve specifically
nIterations <- 5000 # Set to 5 just to make code checking/running fast. For the final analysis, it should be set to 5000.
nIterationsPcurve <- 2
nIterationVWsensitivity <- 2 # Number of iterations for the Vevea & Woods (2005) step function model sensitivity analysis 

# Controls for the multiple-parameter selection models 

# Whether to apply a 4- or 3-parameter selection model. If fallback == TRUE, the procedure falls back to the 3-parameter selection model. 
# This should be selected when too few effects in the opposite side make the estimate unstable.
fallback <- TRUE

# Even when fallback == FALSE, the 4-parameter selection model still falls back to 3 parameters for the given iteration if,
# (1) it fails to converge or (2) the number of p-values in each of the step intervals gets smaller than minPvalues.
minPvalues <- 4

# Steps and delta parameters for Vevea & Woods selection models 
# Can be adjusted if a different selection process is assumed. 
# Please note that steps vector represents one-tailed p-values.
stepsDelta <- data.frame(
  steps =     c(.0025, .005, .0125, .025, .05, .10, .25, .50, 1),
  moderateSelection = c(1, 0.99, 0.97, 0.95, 0.80, 0.60, 0.50, 0.50, 0.50),
  severeSelection =   c(1, 0.99, 0.97, 0.95, 0.65, 0.40, 0.25, 0.25, 0.25),
  extremeSelection =  c(1, 0.98, 0.95, 0.90, 0.50, 0.20, 0.10, 0.10, 0.10))

# Controls for PET-PEESE
condEst <- TRUE

# Exclude studies having an overall Risk of Bias score of at least x.
acceptableRiskOfBias <- 2

# Sourcing and data -----------------------------------------------------------------
source("functions.R")
source("pcurvePlotOption.R")
source("esConversion.R")
funnel <- metafor::funnel

# Meta-analysis -----------------------------------------------------------

#'# Meta-analysis results

#' **RMA results with model-based SEs**
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
#'
#' **RVE SEs with Satterthwaite small-sample correction**
#' Estimate based on a multilevel RE model with constant sampling correlation model (CHE - correlated hierarchical effects - working model) (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/). 
#' Interpretation of naive-meta-analysis should be based on these estimates.
#'
#' **Prediction interval**
#' Shows the expected range of true effects in similar studies.
#' As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between PI LB and PI UB.
#' Note that these are non-adjusted estimates. An unbiased newly conducted study will more likely fall in an interval centered around bias-adjusted estimate with a wider CI width.
#'
#' **Heterogeneity**
#' Tau can be interpreted as the total amount of heterogeneity in the true effects. 
#' I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates. Estimates calculated by 2 approaches are reported.
#' This is followed by separate estimates of between- and within-cluster heterogeneity and estimated intra-class correlation of underlying true effects.
#' 
#' **Proportion of significant results**
#' What proportion of effects were statistically at the alpha level of .05.
#' 
#' **ES-precision correlation**
#' Kendalls's correlation between the ES and precision
#' 
#' **4/3PSM**
#' Applies a permutation-based, step-function 4-parameter selection model (one-tailed p-value steps = c(.025, .5, 1)). 
#' Falls back to 3-parameter selection model if at least one of the three p-value intervals contains less than 4 p-values.
#' 
#' pvalue = p-value testing H0 that the effect is zero. ciLB and ciUB are lower and upper bound of the CI. k = number of studies. steps = 3 means that the 4PSM was applied, 2 means that the 3PSM was applied.
#' 
#' **PET-PEESE**
#' Estimated effect size of an infinitely precise study. Using 4/3PSM as the conditional estimator instead of PET (can be changed to PET). If the PET-PEESE estimate is in the opposite direction, the effect can be regarded nil. 
#' By default (can be changed to PET), the function employs a modified sample-size based estimator (see https://www.jepusto.com/pet-peese-performance/). 
#' It also uses the same RVE sandwich-type based estimator in a CHE (correlated hierarchical effects) working model with the identical random effects structure as the primary (naive) meta-analytic model.
#' 
#' Name of the estimate parameter denotes whether PET or PEESE was applied.
#' 
#' **WAAP-WLS**
#' Combined WAAP-WLS estimator (weighted average of the adequately powered - weighted least squares). The method tries to identify studies that are adequately powered to detect the meta-analytic effect.
#' If there's none or only one such study, the methods falls back to WLS estimator (Stanley & Doucouliagos, 2015).
#' If there are at least two, WAAP returns a WLS estimate based on only effects from those studies
#' 
#' type = 1: WAAP estimate, 2: WLS estimate. kAdequate = number of adequately powered studies
#' 
#' **p-uniform**
#' Permutation-based new version of p-uniform method, the so-called p-uniform* (van Aert, van Assen, 2021).
#' 
#' **p-curve**
#' Permutation-based p-curve method. Output should be pretty self-explanatory.
#' 
#' **Power for detecting SESOI and bias-corrected parameter estimates**
#' Estimates of the statistical power for detecting a smallest effect sizes of interest equal to .20, .50, and .70 in SD units (Cohen's d). 
#' A sort of a thought experiment, we also assumed that population true values equal the bias-corrected estimates (4/3PSM or PET-PEESE) and computed power for those.
#' 
#' **Handling of dependencies in bias-correction methods**
#' To handle dependencies among the effects, the 4PSM, p-curve, p-uniform are implemented using a permutation-based procedure, randomly selecting only one focal effect (i.e., excluding those which were not coded as being focal) from a single study and iterating nIterations times.
#' Lastly, the procedure selects the result with the median value of the ES estimate (4PSM, p-uniform) or median z-score of the full p-curve (p-curve).
#' 

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
dataMind %$% metafor::funnel.default(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline = 0, pch = 20, yaxis = "sei", digits = c(1, 2), xlab = expression(paste("Hedges' ",italic("g"))))

#'### Biofeedback
dataBio %$% metafor::funnel.default(yi, vi, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline = 0, pch = 20, yaxis = "sei", digits = c(1, 2), xlab = expression(paste("Hedges' ",italic("g"))))


#'## Forest plots
#'### Mindfulness
dataMind %$% forest(rmaObjects$`Self-administered mindfulness`$`RMA.MV object with RVE SEs with n/(n-p) small-sample correction`, order = "prec", slab = label, efac = 1, cex = .9, col="gray40", psize=1, cex.lab=1, cex.axis=1,
                    xlab = expression(paste("Hedges' ",italic("g"))), header="PaperID/StudyID/EffectID", mlab="", addpred = T)

#'### Biofeedback
dataBio %$% forest(rmaObjects$Biofeedback$`RMA.MV object with RVE SEs with n/(n-p) small-sample correction`, order = "prec", slab = label, efac = 1, cex = 1, col="gray40", psize=1, cex.lab=1, cex.axis=1,
                    xlab = expression(paste("Hedges' ",italic("g"))), header="PaperID/StudyID/EffectID", mlab="", addpred = T)

#'## p-curve plots
#'### Mindfulness
quiet(pcurveMod(metaResultsPcurve$`Self-administered mindfulness`, effect.estimation = FALSE, plot = TRUE))

#'### Biofeedback
quiet(pcurveMod(metaResultsPcurve$Biofeedback, effect.estimation = FALSE, plot = TRUE))

#'## PET-PEESE plots
#' Using the sqrt(2/n) and 2/n terms instead of SE and var for PET and PEESE, respectively since modified sample-size based estimator was implemented (see https://www.jepusto.com/pet-peese-performance/).
#' 

#'### Mindfulness
quiet(petPeese(dataMind))
if(results[[1]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[1]]$`Publication bias`$`4/3PSM`["est"] > 0){
  dataObjects[[1]] %$% plot(nTerm, yi, main="PEESE", xlab = "2/N", ylab = "Effect size", pch = 19, cex.main = 2, cex = .30, xlim = c(0, .4), xaxs = "i")
} else {
  dataObjects[[1]] %$% plot(sqrt(nTerm), yi, main="PET", xlab = "sqrt(2/n)", ylab = "Effect size", pch = 19, cex.main = 2, cex = .3, xlim = c(0, .4), ylim = c(-1.5, 2), xaxs = "i")}
abline((if(results[[1]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[1]]$`Publication bias`$`4/3PSM`["est"] > 0) {peese} else {pet}), lwd = 3, lty = 2, col = "red")

#'### Biofeedback
quiet(petPeese(dataBio))
if(results[[2]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[2]]$`Publication bias`$`4/3PSM`["est"] > 0){
  dataObjects[[2]] %$% plot(nTerm, yi, main="PEESE", xlab = "2/n", ylab = "Effect size", pch = 19, cex.main = 2, cex = .30, xlim = c(0, .2), ylim = c(-1, 3), xaxs = "i")
} else {
  dataObjects[[2]] %$% plot(sqrt(nTerm), yi, main="PET", xlab = "sqrt(2/n)", ylab = "Effect size", pch = 19, cex.main = 2, cex = .3, xlim = c(0, .2), ylim = c(-1.5, 3), xaxs = "i")}
abline((if(results[[2]]$`Publication bias`$`4/3PSM`["pvalue"] < alpha & ifelse(exists("side") & side == "left", -1, 1) * results[[2]]$`Publication bias`$`4/3PSM`["est"] > 0) {peese} else {pet}), lwd = 3, lty = 2, col = "red")

# Moderator/sensitivity analyses ------------------------------------------

#'## Moderator analysis for the number of intervention sessions

numberIntervention <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + interventionNumber, V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- coef_test(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  numberIntervention[[i]] <- list("Number of intervention sessions" = table(dataObjects[[i]]$interventionNumber), "Model results" = RVEmodel)
}
numberIntervention <- setNames(numberIntervention, nm = namesObjects)
numberIntervention

#'## Moderator analysis for the frequency of intervention

freqIntervention <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + interventionFrequency, V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- coef_test(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  freqIntervention[[i]] <- list("Frequency of intervention sessions" = table(dataObjects[[i]]$interventionFrequency), "Model results" = RVEmodel)
}
freqIntervention <- setNames(freqIntervention, nm = namesObjects)
freqIntervention

#'## Moderator analysis for the duration of intervention

durationIntervention <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + I(interventionDuration/60), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- coef_test(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  durationIntervention[[i]] <- list("Duration of intervention sessions" = table(dataObjects[[i]]$interventionDuration/60), "Model results" = RVEmodel)
}
durationIntervention <- setNames(durationIntervention, nm = namesObjects)
durationIntervention

#'## Moderator analysis for the proportion of females

femaleProportion <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + percFemale, V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- coef_test(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  femaleProportion[[i]] <- list("Proportion of females" = median(dataObjects[[i]]$percFemale, na.rm = TRUE), "Model results" = RVEmodel)
}
femaleProportion <- setNames(femaleProportion, nm = namesObjects)
femaleProportion

#'## Subgroup analysis for the type of the comparison group

compGroup <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + factor(comparisonGroupType), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- coef_test(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  compGroup[[i]] <- list(table(dataObjects[[i]]$comparisonGroupType), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(c(1:2)), vcov = "CR2"))
}
compGroup <- setNames(compGroup, nm = namesObjects)
compGroup

#'## Subgroup analysis for the population type

popType <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + factor(populationType), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- coef_test(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  popType[[i]] <- list(table(dataObjects[[i]]$populationType), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(c(1:3)), vcov = "CR2"))
}
popType <- setNames(popType, nm = namesObjects)
popType

#'## Subgroup analysis for timing of the effect
dataObjectsPost <- list("Mind" = dataMindPost, "Bio" = dataBioPost)

effTiming <- list(NA)
for(i in 1:length(dataObjectsPost)){
  viMatrix <- impute_covariance_matrix(dataObjectsPost[[i]]$vi, cluster = dataObjectsPost[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + factor(timingEffect), V = viMatrix, data = dataObjectsPost[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- coef_test(rmaObject, vcov = "CR2", test = "z", cluster = dataObjectsPost[[i]]$study)
  effTiming[[i]] <- list(table(dataObjectsPost[[i]]$timingEffect), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(c(1:2)), vcov = "CR2"))
}
effTiming <- setNames(effTiming, nm = namesObjects)
effTiming

#'## Subgroup analysis for type of source

viMatrixSourceType <- impute_covariance_matrix(dataMind$vi, cluster = dataMind$study, r = rho)
rmaObjectSourceType <- rma.mv(yi ~ 0 + factor(sourceType), V = viMatrixSourceType, data = dataMind, method = "REML", random = ~ 1|study/result, sparse = TRUE)
effTimingSourceType <- list(table(dataMind$sourceType), "Model results" = rmaObjectSourceType, "RVE Wald test" = Wald_test(rmaObjectSourceType, constraints = constrain_equal(c(1:5)), vcov = "CR2"))
effTimingSourceType


#'# Moderator/sensitivity analyses
#' The below reported meta-regressions are all implemented as a multivariate RVE-based models using the CHE working model (Pustejovsky & Tipton, 2020; https://osf.io/preprints/metaarxiv/vyfcj/).
#' Testing of contrasts is carried out using a robust Wald-type test testing the equality of estimates across levels of the moderator.
#' 

#'## Published status

pubUnpub <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + factor(published), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  pubUnpub[[i]] <- list(table(dataObjects[[i]]$published), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
pubUnpub <- setNames(pubUnpub, nm = namesObjects)
pubUnpub

#'## Excluding effects from non-randomized designs
rndNonrnd <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + factor(researchDesign == 1), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  rndNonrnd[[i]] <- list(table(dataObjects[[i]]$researchDesign), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
rndNonrnd <- setNames(rndNonrnd, nm = namesObjects)
rndNonrnd

#'## Excluding effects due to inconsistent means or SDs
dataObjectsIncon <- list("Mind" = dataMindIncon, "Bio" = dataBioIncon)

consIncons <- list(NA)
i <- 2 # Only for biofeedback, since there were 0 inconsistent means or SDs for mindfulness studies.
 #for(i in 1:length(dataObjectsIncon)){
  viMatrix <- impute_covariance_matrix(dataObjectsIncon[[i]]$vi, cluster = dataObjectsIncon[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + factor(as.logical(inconsistenciesCountGRIMMER)), V = viMatrix, data = dataObjectsIncon[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjectsIncon[[i]]$study)
  consIncons[[i]] <- list("Count of GRIM/GRIMMER inconsistencies" = table(as.logical(dataObjectsIncon[[i]]$inconsistenciesCountGRIMMER)), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
 #}
consIncons <- setNames(consIncons, nm = namesObjects)
consIncons

#'## Excluding effects due to a high risk of bias

# # Probably need to edit to comply with that is given in the ms: "Following RoB 2 recommendations a study was categorized overall as a high risk of bias if one of two conditions are met: 
# # A) The study scores a  high risk of bias in at least one domain or B) the study is evaluated as having some concerns for more than one domain. 
# # A study was judged as having “some concern” whether it raised some concerns in at least one domain. 
# # Finally a study was assessed as having a low risk of bias if it was judged as having a low risk of bias in all of the five domains. 
dataObjectsRob <- list("Mind" = dataMindRob, "Bio" = dataBioRob)

highRoB <- list(NA)
for(i in 1:length(dataObjectsRob)){
  viMatrix <- impute_covariance_matrix(dataObjectsRob[[i]]$vi, cluster = dataObjectsRob[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + factor(robOverall > acceptableRiskOfBias), V = viMatrix, data = dataObjectsRob[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjectsRob[[i]]$study)
  highRoB[[i]] <- list("RoB" = table(as.logical(dataObjectsRob[[i]]$robOverall2)), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
highRoB <- setNames(highRoB, nm = namesObjects)
highRoB

#'## Comparison of strategies

#'### Model without covariates
viMatrixStratComp <- impute_covariance_matrix(datIncluded$vi, cluster = datIncluded$study, r = rho)
rmaObjectStratComp <- rma.mv(yi ~ 0 + factor(strategy), V = viMatrixStratComp, data = datIncluded, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelStratComp <- conf_int(rmaObjectStratComp, vcov = "CR2", test = "z", cluster = datIncluded$study)
list("Model results" = RVEmodelStratComp, "RVE Wald test" = Wald_test(rmaObjectStratComp, constraints = constrain_equal(1:2), vcov = "CR2"))

#'### Model with covariates
#' Controlling for design-related factors that are prognostic w.r.t. the effect sizes (i.e., might vary across moderator categories)
viMatrixStratComp <- impute_covariance_matrix(datIncluded$vi, cluster = datIncluded$study, r = rho)
rmaObjectStratComp <- rma.mv(yi ~ 0 + factor(strategy) + researchDesign + populationType + comparisonGroupType + published + robOverall, V = viMatrixStratComp, data = datIncluded, method = "REML", random = ~ 1|study/result, sparse = TRUE)
RVEmodelStratComp <- conf_int(rmaObjectStratComp, vcov = "CR2", test = "z", cluster = datIncluded$study)
list("Model results" = RVEmodelStratComp, "RVE Wald test" = Wald_test(rmaObjectStratComp, constraints = constrain_equal(1:2), vcov = "CR2"))

#'## Subgroup analysis for stress vs affective consequences

stressAffectConseq <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + as.factor(stressAffective), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  stressAffectConseq[[i]] <- list("Number of included effects per category" = table(dataObjects[[i]]$stressAffective), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:2), vcov = "CR2"))
}
stressAffectConseq <- setNames(stressAffectConseq, nm = namesObjects)
stressAffectConseq

#'### Forest plots
#'#### Mindfulness
dataMind %>% filter(stressAffective == 1) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Mindfulness (StressResponse = 1)")
dataMind %>% filter(stressAffective == 2) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Mindfulness (AffectiveConsequences = 2)")
#'#### Biofeedback
dataBio %>% filter(stressAffective == 1) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Biofeedback (StressResponse = 1)")
dataBio %>% filter(stressAffective == 2) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Biofeedback (AffectiveConsequences = 2)")

#'## Subgroup analysis for stress vs affective consequences
stressComponentClusters <- list(NA)
for(i in 1:length(dataObjects)){
  viMatrix <- impute_covariance_matrix(dataObjects[[i]]$vi, cluster = dataObjects[[i]]$study, r = rho)
  rmaObject <- rma.mv(yi ~ 0 + as.factor(stressCompRecoded), V = viMatrix, data = dataObjects[[i]], method = "REML", random = ~ 1|study/result, sparse = TRUE)
  RVEmodel <- conf_int(rmaObject, vcov = "CR2", test = "z", cluster = dataObjects[[i]]$study)
  stressComponentClusters[[i]] <- list("Number of included effects per category" = table(dataObjects[[i]]$stressCompRecoded), "Model results" = RVEmodel, "RVE Wald test" = Wald_test(rmaObject, constraints = constrain_equal(1:3), vcov = "CR2"))
}
stressComponentClusters <- setNames(stressComponentClusters, nm = namesObjects)
stressComponentClusters

#'### Forest plots
#'#### Mindfulness
dataMind %>% filter(stressCompRecoded == 1) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Mindfulness (AffectiveComponent)")
dataMind %>% filter(stressCompRecoded == 2) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Mindfulness (CognititveComponent)")
dataMind %>% filter(stressCompRecoded == 3) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Mindfulness (PhysiologicalComponent)")
#'#### Biofeedback
dataBio %>% filter(stressCompRecoded == 1) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Biofeedback (AffectiveComponent)")
dataBio %>% filter(stressCompRecoded == 2) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Biofeedback (CognititveComponent)")
dataBio %>% filter(stressCompRecoded == 3) %$% forest(yi, vi, subset=order(vi), slab = label)
title("Biofeedback (PhysiologicalComponent)")

