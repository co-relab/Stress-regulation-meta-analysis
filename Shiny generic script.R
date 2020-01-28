
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

#p-curve data export
set.seed(123)
pcurve.data <- na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y",]$study)], replacement = ""))
#write(pcurve.data, "pcurve.data.txt")

#Permutation p-curve
variables <- 15
pcurve <- matrix(ncol=variables, nrow=nsim)
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y",]$study)], replacement = "")))
}
pcurve.out <- data.frame(pcurve)

# Now filter out results that don't go into the MA
data <- data[!is.na(data$g.calc),]

# Meta-analysis
ma <- rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result)
robma <- robust.rma.mv(x = ma, cluster = data$study)

# 95% prediction interval
pi.lb <- round(predict.rma(robma)$cr.lb, 3)
pi.ub <- round(predict.rma(robma)$cr.ub, 3)

# Heterogeneity
# Total heterogeneity - tau
tau <- sqrt(sum(robma$sigma2))


# I^2
W <- diag(1/data$g.var.calc)
X <- model.matrix(robma)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2<- 100 * sum(robma$sigma2) / (sum(robma$sigma2) + (robma$k-robma$p)/sum(diag(P)))

# Separate estimates of between- and within-cluster heterogeneity
BW.hetero <- round(100 * robma$sigma2 / (sum(robma$sigma2) + (robma$k-robma$p)/sum(diag(P))), 2)

res.R <- rma.mv(yi = g.calc, V = g.var.calc, data = data, struct="UN", random = ~ 1|study/result)
res.F <- rma.mv(yi = g.calc, V = g.var.calc, data = data)

# Jackson's approach to I^2
JI2 <- round(c(100 * (vcov(res.R)[1,1] - vcov(res.F)[1,1]) / vcov(res.R)[1,1]), 2)

# Proportion of significant results
sig.prop <- as.integer(table(data$p < .05)[2])/length(data$p < .05)

# Intra-class correlation of underlying true effects
icc <- round(robma$sigma2[1] / sum(robma$sigma2), 2)

# Contour enhanced funnel plot
funnel(robma, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei")
funnel <- recordPlot()

# Correlation between the ES and precision
es.prec <- cor(data$g.calc, sqrt(data$g.var.calc), method = "kendall")

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

# Small-study effects correction
# 3-parameter selection model
ThreePSM <- threePSM.est(data$g.calc, data$g.var.calc)

# PET-PEESE
pp <- with(data, pet.peese(g.calc, g.var.calc, study, result))

# PET-PEESE plot
if(ThreePSM$value[4] < .05 & ThreePSM$value[1] > 0)
{plot(data$g.var.calc, data$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27),xaxs="i")} else {plot(sqrt(data$g.var.calc), data$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM$value[4] < .05 & ThreePSM$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
PP.plot <- recordPlot()

# p-uniform
puniform.out <- puniform(yi = data[!(data$result %in% c(211)),]$g.calc, vi = data[!(data$result %in% c(211)),]$g.var.calc, alpha = .05, side = "right", method = "P")

# Power based on PEESE and 3PSM parameter estimates
power.peese <- NA
peese.est <- with(data, pet.peese(g.calc, g.var.calc, study, result))[1]
power.3PSM <- NA
TPSM.est <- ThreePSM$value[1]

for(i in 1:length(data$N)){
  power.peese[i] <- pwr::pwr.t.test(n = data$N, d = peese.est)$power[i]
  power.3PSM[i] <- pwr::pwr.t.test(n = data$N, d = TPSM.est)$power[i]
}
pwr.peese.result <- median(power.peese)*100
pwr.3PSM.result <- median(power.3PSM)*100


# Nested list -------------------------------------------------------------

Results.ES <- list("Overall effect size estimate" = as.numeric(robma$beta),
                           "Confidence interval" = c(robma$ci.lb, robma$ci.ub),
                           "p-value" = robma$pval,
                           "Tau" = tau,
                           "I^2" = I2,
                           "Proportion of significant results" = sig.prop,
                           "Prediction Interval" = c(pi.lb, pi.ub),
                           "Forest plot" = forest)

Publication.bias <- list("p-curve mean p-values" = as.data.frame(describe(pcurve.out[c(4, 6, 8, 10)]))[,3, drop = FALSE],
                                 "Funnel plot" = funnel,
                                 "3PSM ES estimate" = ThreePSM[1, 4],
                                 "3PSM confidence interval" = c(ThreePSM[5, 4], ThreePSM[6, 4]),
                                 "3PSM p-value" = ThreePSM[4, 4],
                                 "Whether PET or PEESE was used" = ifelse(ThreePSM$value[4] < .05 & ThreePSM$value[1] > 0, "PEESE", "PET"),
                                 "PET-PEESE ES estimate" = as.numeric(pp[1]),
                                 "PET-PEESE confidence interval" = as.numeric(c(pp[5], pp[6])),
                                 "PET-PEESE p-value" = pp[4],
                                 "Median power for detecting PEESE estimate" = pwr.peese.result,
                                 "Median power for detecting 3PSM estimate" = pwr.3PSM.result)
