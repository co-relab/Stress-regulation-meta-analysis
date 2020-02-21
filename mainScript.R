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

# Source scripts
source("functions.R")
source("ES_conversion.R")

# No of simulations for the permutation p-curve
nsim <- 5 # Set to 5 just to make code checking/running fast. For the final paper, it needs to be set to at least 1000 and run overnight.

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


# Results -----------------------------------------------------------------

# Subset
dataMind <- dat[dat$category == 1 & !is.na(dat$g.calc),]
dataBio <- dat[dat$category == 2 & !is.na(dat$g.calc),]
dataNat <- dat[dat$category == 3 & !is.na(dat$g.calc),]
dataSoc <- dat[dat$category == 4 & !is.na(dat$g.calc),]

#'## Meta-analysis
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
maMind <- rma.mv(yi = g.calc, V = g.var.calc, data = dataMind, method = "REML", random = ~ 1|study/result)
rmaMind <- robust.rma.mv(x = maMind, cluster = dataMind$study)

maBio <- rma.mv(yi = g.calc, V = g.var.calc, data = dataBio, method = "REML", random = ~ 1|study/result)
rmaBio <- robust.rma.mv(x = maBio, cluster = dataBio$study)

maNat <- rma.mv(yi = g.calc, V = g.var.calc, data = dataNat, method = "REML", random = ~ 1|study/result)
rmaNat <- robust.rma.mv(x = maNat, cluster = dataNat$study)

maSoc <- rma.mv(yi = g.calc, V = g.var.calc, data = dataSoc, method = "REML", random = ~ 1|study/result)
rmaSoc <- robust.rma.mv(x = maSoc, cluster = dataSoc$study)





