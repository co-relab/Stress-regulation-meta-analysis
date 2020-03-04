

knitr::opts_chunk$set(echo=FALSE, warning = FALSE)
rm(list = ls())

# Assumed default pre-post correlation for within-subjects design, .50.
# Here you can perform the sensitivity analysis to determine the impact of the assumed correlation on the overall effect size estimate.
# E.g., for corr = c(.10, .30, .50, .70, 90).
corr <- 0.5

# Install required R libraries if not installed already
list.of.packages <- c("metafor", "lme4", "ggplot2", "knitr", "psych", "puniform", "reshape2", "kableExtra", "lmerTest", "pwr", "Amelia")
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
vars <- subset(dat, select = c(g.calc, g.var.calc, p, N, study, result, label))
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
# rmaObjects <- list("Mind" = rmaMind, "Bio" = rmaBio, "Nat" = rmaNat, "Soc" = rmaSoc)
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
excludeRisk <- 3 # Exclude studies having a risk of bias of at least x
rmaRoB <- setNames(lapply(dataObjects, function(x){rmaCustom(x[x$Overall.risk.of.bias < excludeRisk,])}), nm = namesObjects)
RoBResults <- list(NA)
for(i in 1:length(rmaRnd)){
  RoBResults[[i]] <- maResults(rmaObject = rmaRnd[[i]], data = dataObjects[[i]][dataObjects[[i]]$Overall.risk.of.bias < excludeRisk,], briefBias = T)
}
RoBResults <- setNames(RoBResults, nm = namesObjects)
RoBResults

# Comparison of categories after controlling for prognostic factors w.r.t. the effect sizes
rmaCompare <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, mods = ~factor(category), data = dat[!is.na(dat$g.calc),], method = "REML", random = ~ 1|study/result), cluster = dat[!is.na(dat$g.calc),]$study)
rmaCompare





# Diag NOT EVEN "UNDER CONSTRUCTION" YET

#+eval = FALSE
# Initial outlier diagnostics
# Univariate MA
ma.uni <- rma(yi = g.calc, vi = g.var.calc, data = dat, method = "REML", slab = result)

#+eval = FALSE
# MA diagnostics
baujat(ma.uni)

#+eval = FALSE
#fit FE model to all possible subsets
gosh.plot <- gosh(ma.uni, progbar = TRUE, subsets = 1000, parallel = "multicore")
# plot(gosh.plot, out = , breaks=50) # Testing the influence of single outliers

#+eval = FALSE
# Influence diagnostics
inf <- influence(ma.uni, progbar = T)

#+eval = FALSE
### Plot the influence diagnostics
plot(inf)

#+eval = TRUE
# Outlier removal in case of a need
# Excluding improbably big effect sizes or ES with improbably small SE, i.e. excerting a big influence on the MA model due to combination of huge ES and small variance.
# Sensitivity analysis with the outlying ESs included will be reported as well.
# dat[c(),] <- NA

#'##### Missing data
table(dat$Use.for.Meta == "Yes" & is.na(dat$g.calc))
table(is.na(dat$g.calc))

#'### Percentage of missing data
#'There is very little missing data. Regardless of what imputation procedure is applied, it won't have much effect.
#paste(round(sum(is.na(dat[,1:34]))/prod(dim(dat[,1:34]))*100, 3), "%", sep = "") # insert collumn numbers
#missmap(dat, rank.order = TRUE, margins = c(5, 0), legend = F)    # insert collumn numbers
e

