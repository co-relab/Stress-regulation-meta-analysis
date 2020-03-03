

knitr::opts_chunk$set(echo=FALSE, warning = FALSE)
rm(list = ls())

# Assumed default pre-post correlation for within-subjects design, .50.
# Here you can perform the sensitivity analysis to determine the impact of the assumed correlation on the overall effect size estimate.
# E.g., for corr = c(.10, .30, .50, .70, 90).
corr <- 0.5

# Install required R libraries if not installed already
list.of.packages <- c("metafor", "lme4", "ggplot2", "knitr", "psych", "puniform", "reshape2", "kableExtra", "lmerTest", "pwr", "Amelia", "multcomp")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load required libraries
#+ include = FALSE
lapply(list.of.packages, require, quietly = TRUE, warn.conflicts = FALSE, character.only = TRUE)

#' Statistical analysis was carried out in R, version 3.4.3, using packages "metafor", "lme4", "ggplot2", "knitr", "psych", "puniform", "reshape2", "kableExtra", "lmerTest", "pwr", "Amelia".
#'


# Results -----------------------------------------------------------------
source("functions.R")
source("functions2.R")
source("ES_conversion.R")
vars <- subset(dat, select = c(g.calc, g.var.calc, p, N, study, result, label)) # For script testing purposes only
source("SimulateData.R")
dat <- cbind(dat, vars)

# Subset
dataMind <- dat[dat$category == 1 & !is.na(dat$g.calc),]
dataBio <- dat[dat$category == 2 & !is.na(dat$g.calc),]
dataNat <- dat[dat$category == 3 & !is.na(dat$g.calc),]
dataSoc <- dat[dat$category == 4 & !is.na(dat$g.calc),]

#'## Meta-analysis
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
namesObjects <- c("Self-administered mindfulness", "Biofeedback", "Being in nature", "Social support")
levels(dat$category) <- namesObjects
dataObjects <- list("Mind" = dataMind, "Bio" = dataBio, "Nat" = dataNat, "Soc" = dataSoc)

rmaObjects <- setNames(lapply(dataObjects, function(x){rmaCustom(x)}), nm = namesObjects)

# Further results
briefBias <- TRUE # For a more elaborate output from the pub bias tests, set to FALSE
results <- list(NA)
for(i in 1:length(rmaObjects)){
  results[[i]] <- maResults(rmaObject = rmaObjects[[i]], data = dataObjects[[i]], alpha = .05, briefBias = T)
}

results <- setNames(results, nm = namesObjects)
results

# Published vs unpublished studies 
publishedRMA <- rmaCustom(dat[!is.na(dat$g.calc) & dat$published == 1,])
unpublishedRMA <- rmaCustom(dat[!is.na(dat$g.calc) & dat$published == 2,])
pubResults <- maResults(rmaObject = publishedRMA, data = dat[!is.na(dat$g.calc) & dat$published == 1,], briefBias = T)
unpubResults <- maResults(rmaObject = unpublishedRMA, data = dat[!is.na(dat$g.calc) & dat$published == 2,], briefBias = T)
pubResults
unpubResults

# Sensitivity analysis excluding effects from non-randomized designs CHANGE THE CODING FROM 1,2 TO 0,1
rmaRnd <- setNames(lapply(dataObjects, function(x){rmaCustom(x[x$research_design == 1,])}), nm = namesObjects)
rndResults <- list(NA)
for(i in 1:length(rmaRnd)){
  rndResults[[i]] <- maResults(rmaObject = rmaRnd[[i]], data = dataObjects[[i]][dataObjects[[i]]$research_design == 1,], briefBias = T)
}
rndResults <- setNames(rndResults, nm = namesObjects)
rndResults

# Sensitivity analysis excluding effects reported in studies having a high risk of bias
# Probably need to edit to comply with that is given in the ms: "Following RoB 2 recommendations a study was categorized overall as a high risk of bias if one of two conditions are met: 
# A) The study scores a  high risk of bias in at least one domain or B) the study is evaluated as having some concerns for more than one domain. 
# A study was judged as having “some concern” whether it raised some concerns in at least one domain. 
# Finally a study was assessed as having a low risk of bias if it was judged as having a low risk of bias in all of the five domains. 

excludeRisk <- 3 # Exclude studies having a risk of bias of at least x
rmaRoB <- setNames(lapply(dataObjects, function(x){rmaCustom(x[x$Overall.risk.of.bias < excludeRisk,])}), nm = namesObjects)
RoBResults <- list(NA)
for(i in 1:length(rmaRnd)){
  RoBResults[[i]] <- maResults(rmaObject = rmaRnd[[i]], data = dataObjects[[i]][dataObjects[[i]]$Overall.risk.of.bias < excludeRisk,], briefBias = T)
}
RoBResults <- setNames(RoBResults, nm = namesObjects)
RoBResults

# Example: Comparison of categories after controlling for prognostic factors w.r.t. the effect sizes
# What moderator/meta-regression analyses shall we conduct is a substantial question to discuss.
rmaCompareNull <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, mods = research_design + type_of_population + type_of_comparison_group + published + Overall.risk.of.bias - 1, struct="DIAG", data = dat[!is.na(dat$g.calc),], method = "ML", random = ~ factor(category) | result), cluster = dat[!is.na(dat$g.calc),]$study)
rmaCompare <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, mods = ~factor(category) + research_design + type_of_population + type_of_comparison_group + published + Overall.risk.of.bias - 1,struct="DIAG", data = dat[!is.na(dat$g.calc),], method = "ML", random = ~ factor(category) | result), cluster = dat[!is.na(dat$g.calc),]$study)
rmaCompare


# Likelihood ratio test for the differences between categories
# Omnibus test
anova(rmaCompareNull, rmaCompare)

# Contrasts 
# p-values adjusted using Holm's method
summary(glht(rmaCompare, linfct=cbind(contrMat(c("Self-administered mindfulness" = 1, "Biofeedback" = 1, "Being in nature" = 1, "Social support" = 1), type="Tukey"), 0, 0, 0, 0, 0)), test=adjusted("holm"))


