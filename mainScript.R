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

# No of simulations for the permutation p-curve and 3PSM model
nsim <- 10 # Set to 5 just to make code checking/running fast. For the final paper, it needs to be set to at least 1000 and run overnight.

# Sourcing and data -----------------------------------------------------------------
source("functions.R")
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
briefBias <- FALSE # For a more elaborate output from the pub bias tests, set to FALSE
results <- list(NA)
for(i in 1:length(rmaObjects)){
  results[[i]] <- maResults(data = dataObjects[[i]], rmaObject = rmaObjects[[i]], briefBias = briefBias, pcurveOut = T)
}

results <- setNames(results, nm = namesObjects)

#+ include = TRUE
results

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


# Published vs unpublished studies ----------------------------------------
# 
# publishedRMA <- rmaCustom(dat[!is.na(dat$yi) & !is.na(dat$published) & dat$published == 1,])
# unpublishedRMA <- rmaCustom(dat[!is.na(dat$yi) & !is.na(dat$published) & dat$published == 0,])
# pubResults <- maResults(rmaObject = publishedRMA, data = dat[!is.na(dat$yi) & !is.na(dat$published) & dat$published == 1,], briefBias = T, pcurve = F)
# unpubResults <- maResults(rmaObject = unpublishedRMA, data = dat[!is.na(dat$yi) & !is.na(dat$published) & dat$published == 0,], briefBias = T, pcurve = F)
# pubResults
# unpubResults
# 
# # Sensitivity analysis excluding effects from non-randomized designs -------
# rmaRnd <- setNames(lapply(dataObjects, function(x){rmaCustom(x[x$researchDesign == 1,])}), nm = namesObjects)
# rndResults <- list(NA)
# for(i in 1:length(rmaRnd)){
#   rndResults[[i]] <- maResults(rmaObject = rmaRnd[[i]], data = dataObjects[[i]][dataObjects[[i]]$researchDesign == 1,], briefBias = T, pcurve = F)
# }
# rndResults <- setNames(rndResults, nm = namesObjects)
# rndResults
# 
# # Sensitivity analysis excluding effects based on inconsistent means or SDs -------
# 
# rmaRnd <- setNames(lapply(dataObjects, function(x){rmaCustom(x[x$inconsistenciesCount == 0,])}), nm = namesObjects)
# rndResults <- list(NA)
# for(i in 1:length(rmaRnd)){
#   rndResults[[i]] <- maResults(rmaObject = rmaRnd[[i]], data = dataObjects[[i]][dataObjects[[i]]$inconsistenciesCount == 0,], briefBias = T, pcurve = F)
# }
# rndResults <- setNames(rndResults, nm = namesObjects)
# rndResults
# 
# # Sensitivity analysis excluding effects based on a high risk of bias -------------
# 
# # Probably need to edit to comply with that is given in the ms: "Following RoB 2 recommendations a study was categorized overall as a high risk of bias if one of two conditions are met: 
# # A) The study scores a  high risk of bias in at least one domain or B) the study is evaluated as having some concerns for more than one domain. 
# # A study was judged as having “some concern” whether it raised some concerns in at least one domain. 
# # Finally a study was assessed as having a low risk of bias if it was judged as having a low risk of bias in all of the five domains. 
# 
# excludeRisk <- 3 # Exclude studies having a risk of bias of at least x
# rmaRoB <- setNames(lapply(dataObjects, function(x){rmaCustom(x[!is.na(x$robOverall) & x$robOverall < excludeRisk,])}), nm = namesObjects)
# RoBResults <- list(NA)
# for(i in 1:length(rmaRnd)){
#   RoBResults[[i]] <- maResults(rmaObject = rmaRnd[[i]], data = dataObjects[[i]][!is.na(dataObjects[[i]]$robOverall) & dataObjects[[i]]$robOverall < excludeRisk,], briefBias = T, pcurve = F)
# }
# RoBResults <- setNames(RoBResults, nm = namesObjects)
# RoBResults
# 
# # Moderator analysis for strategies ---------------------------------------
# # The other moderator analyses will follow the same analytic pipeline
# 
# # Comparison of categories after controlling for prognostic factors w.r.t. the effect sizes
# rmaCompare <- robust.rma.mv(rma.mv(yi = yi, V = vi, mods = ~factor(strategy), data = dat[!is.na(dat$yi),], method = "REML", random = ~ 1|study/result), cluster = dat[!is.na(dat$yi),]$study)
# rmaCompare
# 
# # Defining the null model for moderator analyses
# rmaNull <- robust.rma.mv(rma.mv(yi = yi, V = vi, mods = researchDesign + populationType + comparisonGroupType + published + robOverall - 1, struct="DIAG", data = dat[!is.na(dat$yi),], method = "ML", random = ~ factor(strategy) | result), cluster = dat[!is.na(dat$yi),]$study)
# 
# # Strategies
# # Comparison of categories of strategies after controlling for prognostic factors w.r.t. the effect sizes
# # What moderator/meta-regression analyses shall we conduct is a substantial question to discuss.
# rmaCat <- robust.rma.mv(rma.mv(yi = yi, V = vi, mods = ~factor(strategy) + researchDesign + populationType + comparisonGroupType + published + robOverall - 1,struct="DIAG", data = dat[!is.na(dat$yi),], method = "ML", random = ~ factor(strategy) | result), cluster = dat[!is.na(dat$yi),]$study)
# rmaCat
# 
# # Likelihood ratio test for the differences between categories
# # Omnibus test
# anova(rmaNull, rmaCat)
# 
# # Contrasts 
# # p-values adjusted using Holm's method
# summary(glht(rmaCat, linfct=cbind(contrMat(c("Self-administered mindfulness" = 1, "Biofeedback" = 1), type="Tukey"), 0, 0, 0, 0, 0)), test=adjusted("holm"))
# 
# # Components
# # Comparison of components after controlling for prognostic factors w.r.t. the effect sizes
# rmaComp <- robust.rma.mv(rma.mv(yi = yi, V = vi, mods = ~factor(stressComponentType) + researchDesign + populationType + comparisonGroupType + published + robOverall - 1, struct="DIAG", data = dat[!is.na(dat$yi),], method = "ML", random = ~ factor(strategy) | result), cluster = dat[!is.na(dat$yi),]$study)
# rmaComp
# 
# # Likelihood ratio test for the differences between categories
# # Omnibus test
# anova(rmaNull, rmaComp)
# 
# # Wald's robust F test
# # Wald_dv_multilevel <- Wald_test(rmaComp,
# #                                 constraints = constrain_equal(1:2), 
# #                                 vcov = "CR2")
# 
# # Contrasts 
# # p-values adjusted using Holm's method
# summary(glht(rmaComp, linfct=cbind(contrMat(c("AFloAneV" = 1, "AFhiAneV" = 1, "AFloApoV" = 1, "AFhiApoV" = 1, "cognitiveComp" = 1, "physiologicalComp" = 1), type="Tukey"), 0, 0, 0, 0, 0)), test=adjusted("holm"))
