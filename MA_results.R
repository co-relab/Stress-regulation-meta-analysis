#' ---
#' title: "Social Thermoregulation: A Meta-Analysis"
#' author: "Ivan Ropovik"
#' output:
#'    html_document:
#'       toc: true
#'       toc_float: true
#'       number_sections: true
#'       code_folding: hide
#'       fig_retina: 2
#'       theme: paper
#' always_allow_html: yes
#' ---
#+ setup, include=FALSE
knitr::opts_chunk$set(echo=FALSE, warning = FALSE)

# further code edits: addcred = T

rm(list = ls())

# Assumed default pre-post correlation for within-subjects design, .50.
# Here you can perform the sensitivity analysis to determine the impact of the assumed correlation on the overall effect size estimate.
# E.g., for corr = c(.10, .30, .50, .70, 90).
corr <- 0.5

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

# Source scripts
source("functions.R")
source("ES_conversion.R")

# No of simulations for the permutation p-curve
nsim <- 5 # Set to 5 just to make code checking/running fast. For the final paper, it needs to be set to at least 1000 and run overnight.

#+eval = FALSE
# Initial outlier diagnostics
# Univariate MA
#ma.uni <- rma(yi = g.calc, vi = g.var.calc, data = dat, method = "REML", slab=result)

#+eval = FALSE
# MA diagnostics
#baujat(ma.uni)

#+eval = FALSE
#fit FE model to all possible subsets
#gosh.plot <- gosh(ma.uni, progbar = TRUE, subsets = 1000, parallel = "multicore")

#+eval = FALSE
# GOSH plot
# MA-fitted outliers for gosh plot
# huge yi 193=120, 213=134, 343=199
# huge vi 211=132
#plot(gosh.plot, out = 120, breaks=50)
#plot(gosh.plot, out = 134, breaks=50)
#plot(gosh.plot, out = 199, breaks=50)
#plot(gosh.plot, out = 132, breaks=50)
# Study 211 causes large shift (∆I^2 = 11) in heterogeneity.

#+eval = FALSE
# Influence diagnostics
#inf <- influence(ma.uni, progbar = T)

#+eval = FALSE
### Plot the influence diagnostics
#plot(inf)

#+eval = FALSE
# View the outliers
#View(dat[c(193, 213, 343, 211),])

#+eval = TRUE
# Outliers
# Excluding 3 outliers based on Baujat analysis, first 2 due to improbably big effect sizes, the last because it also excerts a big influence on the MA model due to combination of huge ES and small variance  (d = 3.23 with n=60, 2.55 with n=63, 1.72 with n=224, respectively).
dat[c(193, 213, 343),] <- NA

#'##### Missing data
table(dat$Use.for.Meta == "Yes" & is.na(dat$g.calc))
table(is.na(dat$g.calc))

###############
#RESULTS
###############

#'# Mood
data.mood <- dat[dat$PA.NA. == "Y" & !is.na(dat$g.calc),]

#'## Meta-analysis
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
ma.mood <- rma.mv(yi = g.calc, V = g.var.calc, data = data.mood, method = "REML", random = ~ 1|study/result)
robma.mood <- robust.rma.mv(x = ma.mood, cluster = data.mood$study)
robma.mood

#'#### 95% prediction interval
pi.lb.mood <- round(predict.rma(robma.mood)$cr.lb, 3)
pi.ub.mood <- round(predict.rma(robma.mood)$cr.ub, 3)

#'As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between `r pi.lb.mood` and `r pi.ub.mood`.
#'Note that these are non-adjusted estimates. An unbiased newly conducted study will likely fall in an interval centered around PET-PEESE estimate with a similar CI width of `r pi.ub.mood - pi.lb.mood`.
#'

#'### Heterogeneity
#'#### Total heterogeneity - tau
#'
tau.mood <- {{round(sqrt(sum(robma.mood$sigma2)), 3)}}
#'The sum of the two variance components is equal to
{{round(tau.mood, 3)}}
#'. That can be interpreted as the total amount of heterogeneity in the true effects.
#'

#'#### $I^2$
W <- diag(1/data.mood$g.var.calc)
X <- model.matrix(robma.mood)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'$I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates.
#'
#'Here, total relative heterogeneity was
I2.mood <- {{round(100 * sum(robma.mood$sigma2) / (sum(robma.mood$sigma2) + (robma.mood$k-robma.mood$p)/sum(diag(P))), 2)}}
{{round(I2.mood, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * robma.mood$sigma2 / (sum(robma.mood$sigma2) + (robma.mood$k-robma.mood$p)/sum(diag(P))), 2)}}
#'%, respectively.

res.R <- rma.mv(yi = g.calc, V = g.var.calc, data = data.mood, struct="UN", random = ~ 1|study/result)
res.F <- rma.mv(yi = g.calc, V = g.var.calc, data = data.mood)
#'Jackson's approach to $I^2$ yields a relative heterogeneity estimate of
{{round(c(100 * (vcov(res.R)[1,1] - vcov(res.F)[1,1]) / vcov(res.R)[1,1]), 2)}}
#'%.

#Confidence intervals and profile plots
#+eval = FALSE
#confint(robma.mood, digits=3)
#par(mfrow=c(2,1))
#profile(robma.mood, sigma2=1, progbar = F)
#profile(robma.mood, sigma2=2, progbar = F)
#Both profile likelihood plots are peaked at the respective parameter estimates (as indicated by the vertical dotted lines) and the log likelihoods quickly decrease (i.e., become more negative) as the values of the components are moved away from the actual REML estimates. Hence, we can be fairly confident that both variance components are identifiable.

#+eval = FALSE
#Testing significance of variance components
#ma0 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.mood, method = "ML", random = ~ 1|study/result)
#ma1 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.mood, method = "ML", random = ~ 1|study/result, sigma2 = c(NA,0))
#ma2 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.mood, method = "ML", random = ~ 1|study/result, sigma2 = c(0,NA))
#anova(ma1, ma0)
#anova(ma2, ma0)

#'#### Proportion of significant results
sig.prop.mood <- as.integer(table(data.mood$p < .05)[2])/length(data.mood$p < .05)
round(sig.prop.mood, 2)

#'#### Intra-class correlation of underlying true effects
round(robma.mood$sigma2[1] / sum(robma.mood$sigma2), 3)

#'### Contour enhanced funnel plot
funnel(robma.mood, level=c(90, 95, 99), shade=c("white", "lightgray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g", steps = 7, digits = c(1,2))
title("Funnel plot for Mood", cex.main = 1)
funnel.mood <- recordPlot()

#'Correlation between the ES and precision
cor(data.mood$g.calc, sqrt(data.mood$g.var.calc), method = "kendall")

#'### Forest plot
forest(x = data.mood$g.calc, vi = data.mood$g.var.calc, subset=order(data.mood$g.calc), slab = data.mood$result,
       xlab = "Hedges' g", xlim = c(-1.6, 3.5), at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
       cex = .8)
title("Forest plot for Mood", cex.main = 1)
addpoly(robma.mood, row = 0, mlab = "", cex = .5, annotate = F)
forest.mood <- recordPlot()


#'## Small-study effects correction
#'
#'### 3-parameter selection model
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.mood <- threePSM.est(data.mood$g.calc, data.mood$g.var.calc)
kable(ThreePSM.mood[c(1, 4, 5, 6),], "html", digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'### PET-PEESE
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.mood <- with(data.mood, pet.peese(g.calc, g.var.calc, study, result))
kable(pp.mood, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### PET-PEESE plot
#'
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.mood$value[4] < .05 & ThreePSM.mood$value[1] > 0)
{plot(data.mood$g.var.calc, data.mood$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .5), ylim = c(-.2, .8), xaxs="i")} else {plot(sqrt(data.mood$g.var.calc), data.mood$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .5), ylim = c(-.2, .8), xaxs="i")}
abline((if(ThreePSM.mood$value[4] < .05 & ThreePSM.mood$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#+include = FALSE
#'#### p-uniform
#'Additional bias-corrected estimate. Because it's far less precise than PET-PEESE, when the n of studies is small, look just at the CI width and p-value.
puniform(yi = data.mood$g.calc, vi = data.mood$g.var.calc, alpha = .05, side = "right", method = "P")

#########################################################################
#########################################################################
#########################################################################

#'# Overall effect (excluding mood)
data <- dat[dat$PA.NA. != "Y" & !is.na(dat$g.calc),]

#'## Evidential value
#'#### Permutation p-curve

#p-curve data export
set.seed(123)
pcurve.data <- na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y",]$study)], replacement = ""))
write(pcurve.data, "pcurve.data.txt")

#Permutation p-curve
variables <- 15
pcurve <- matrix(ncol=variables, nrow=nsim)

#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y",]$study)], replacement = "")))
}

pcurve.out.overall <- data.frame(pcurve)
colnames(pcurve.out.overall) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
#'
#'NOTE: For now, all the p-curve permutations are based just on 100 sets of draws, because they are quite computationally intensive.

kable(describe(pcurve.out.overall, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#hist(pcurve.out$fullz, breaks = nrow(pcurve.out))
#hist(pcurve.out$halfz, breaks = nrow(pcurve.out))
#hist(pcurve.out$fullz33, breaks = nrow(pcurve.out))
#hist(pcurve.out$halfz33, breaks = nrow(pcurve.out))
#hist(pcurve.out$fullp33, breaks = nrow(pcurve.out))
#hist(pcurve.out$halfp33, breaks = nrow(pcurve.out))

#'## Meta-analysis
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
ma.overall <- rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result)
robma.overall <- robust.rma.mv(x = ma.overall, cluster = data$study)
robma.overall

#'#### 95% prediction interval
pi.lb.overall <- round(predict.rma(robma.overall)$cr.lb, 3)
pi.ub.overall <- round(predict.rma(robma.overall)$cr.ub, 3)

#'As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between `r pi.lb.overall` and `r pi.ub.overall`.
#'Note that these are non-adjusted estimates. An unbiased newly conducted study will likely fall in an interval centered around PET-PEESE estimate with a similar CI width of `r pi.ub.overall - pi.lb.overall`.
#'

#'### Heterogeneity
#'#### Total heterogeneity - tau
#'
tau.overall <- sqrt(sum(robma.overall$sigma2))
#'The sum of the two variance components is equal to
{{round(tau.overall, 3)}}
#'. That can be interpreted as the total amount of heterogeneity in the true effects.
#'

#'#### $I^2$
W <- diag(1/data$g.var.calc)
X <- model.matrix(robma.overall)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'$I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates.
#'
I2.overall<- 100 * sum(robma.overall$sigma2) / (sum(robma.overall$sigma2) + (robma.overall$k-robma.overall$p)/sum(diag(P)))
#'Here, total relative heterogeneity was
{{round(I2.overall, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * robma.overall$sigma2 / (sum(robma.overall$sigma2) + (robma.overall$k-robma.overall$p)/sum(diag(P))), 2)}}
#'%, respectively.

res.R.overall <- rma.mv(yi = g.calc, V = g.var.calc, data = data, struct="UN", random = ~ 1|study/result)
res.F.overall <- rma.mv(yi = g.calc, V = g.var.calc, data = data)
#'Jackson's approach to $I^2$ yields a relative heterogeneity estimate of
{{round(c(100 * (vcov(res.R.overall)[1,1] - vcov(res.F.overall)[1,1]) / vcov(res.R.overall)[1,1]), 2)}}
#'%.

# Confidence intervals and profile plots
#+eval = FALSE
#confint(robma.overall, digits=3)
#par(mfrow=c(2,1))
#profile(robma.overall, sigma2=1, progbar = F)
#profile(robma.overall, sigma2=2, progbar = F)
#Both profile likelihood plots are peaked at the respective parameter estimates (as indicated by the vertical dotted lines) and the log likelihoods quickly decrease (i.e., become more negative) as the values of the components are moved away from the actual REML estimates. Hence, we can be fairly confident that both variance components are identifiable.

#+eval = FALSE
# Testing significance of variance components
#ma0.overall <- rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "ML", random = ~ 1|study/result)
#ma1.overall <- rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "ML", random = ~ 1|study/result, sigma2 = c(NA,0))
#ma2.overall <- rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "ML", random = ~ 1|study/result, sigma2 = c(0,NA))
#anova(ma1.overall, ma0.overall)
#anova(ma2.overall, ma0.overall)

#+eval = FALSE
# Univariate MA, qqnorm plot, permutation test
#ma.uni.overall <- rma(yi = g.calc, sei = sqrt(g.var.calc), data = data, method = "REML", slab=result)
#metafor::qqnorm.rma.uni(ma.uni.overall)
#permutest(ma.uni.overall)

#'#### Proportion of significant results
sig.prop.overall <- as.integer(table(data$p < .05)[2])/length(data$p < .05)
round(sig.prop.overall, 2)

#'#### Intra-class correlation of underlying true effects
round(robma.overall$sigma2[1] / sum(robma.overall$sigma2), 2)

#'### Contour enhanced funnel plot
funnel(robma.overall, level=c(90, 95, 99), shade=c("white", "lightgray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g", xlim = c(-1.4, 1.7), steps = 7, digits = c(1,2))
title("Funnel plot", cex.main = 1)
funnel.overall <- recordPlot()

#'Correlation between the ES and precision
cor(data$g.calc, sqrt(data$g.var.calc), method = "kendall")

#'### Forest plot
forest(x = data$g.calc, vi = data$g.var.calc, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
                  xlim = c(-3.5,4.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
                  subset = order(data$g.var.calc),        ### order by size of yi
                  slab = NA, annotate = FALSE, ### remove study labels and annotations
                  efac = 0,                  ### remove vertical bars at end of CIs
                  pch = 19,                  ### changing point symbol to filled circle
                  col = "gray40",            ### change color of points/CIs
                  psize = 5,                 ### increase point size
                  cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
                  lty=c("solid","blank"))  ### remove horizontal line at top of plot

                  title("Forest plot", cex.main = 1)
                  addpoly(robma.overall, row = -3, mlab = "", cex = .5, annotate = FALSE)
                  forest.overall <- recordPlot()

#'## Small-study effects correction
#'
#'### 3-parameter selection model
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.overall <- threePSM.est(data$g.calc, data$g.var.calc)
kable(ThreePSM.overall[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'### PET-PEESE
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.overall <- with(data, pet.peese(g.calc, g.var.calc, study, result))
kable(pp.overall, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### PET-PEESE plot
#'
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.overall$value[4] < .05 & ThreePSM.overall$value[1] > 0)
{plot(data$g.var.calc, data$g.calc, main="PEESE", xlab = "Variance", ylab = "Hedges' g", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27),xaxs="i")} else {plot(sqrt(data$g.var.calc), data$g.calc, main="PET", xlab = "Standard error", ylab = "Hedges' g", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.overall$value[4] < .05 & ThreePSM.overall$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#+include = FALSE
#'#### p-uniform
#'Additional bias-corrected estimate. Because it's far less precise than PET-PEESE, when the n of studies is small, look just at the CI width and p-value. Leaving out study id# 211 because p-uniform won't converge due to huge variance.
puniform(yi = data[!(data$result %in% c(211)),]$g.calc, vi = data[!(data$result %in% c(211)),]$g.var.calc, alpha = .05, side = "right", method = "P")

#'#### Power based on PEESE and 3PSM parameter estimates
power.peese <- NA
peese.est.overall <- with(data, pet.peese(g.calc, g.var.calc, study, result))[1]
power.3PSM <- NA
TPSM.est.overall <- ThreePSM.overall$value[1]

for(i in 1:length(data$N)){
  power.peese[i] <- pwr::pwr.t.test(n = data$N, d = peese.est.overall)$power[i]
  power.3PSM[i] <- pwr::pwr.t.test(n = data$N, d = TPSM.est.overall)$power[i]
}
pwr.peese.result <- median(power.peese)*100
pwr.3PSM.result <- median(power.3PSM)*100

paste("Power to detect PEESE estimate = ", round(pwr.peese.result, 2), "%", sep = "")
paste("Power to detect 3PSM estimate = ", round(pwr.3PSM.result, 2), "%", sep = "")

#' How many studies had more than 50% power to detect the overall PEESE estimate?
table(power.peese>.50)

#########################################################################

# Recode the effect type into compensatory and priming
dat$et.temp <- dat$Effect.Type..Neither..0..Compensatory..1..Priming..2..Both..3.
dat$effect.type <- ifelse((dat$et.temp == 0) |  (dat$et.temp == 3), yes = NA, no = dat$et.temp)
data <- dat[dat$PA.NA. != "Y" & !is.na(dat$g.calc),]

data.comp <- dat[dat$PA.NA. != "Y" & dat$effect.type == 1,]
data.prim <- dat[dat$PA.NA. != "Y" & dat$effect.type == 2,]
data$effect.type <- relevel(factor(data$effect.type), ref="1")

#'# Effect type: compensatory vs priming

#'## Evidential value for effect type = compensatory
#'#### Permutation p-curve

#Permutation p-curve
variables <- 15
pcurve <- matrix(ncol=variables, nrow=nsim)

#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data.comp$label[!duplicated.random(data.comp$study)], replacement = "")))
}
pcurve.out.comp <- data.frame(pcurve)
colnames(pcurve.out.comp) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)

kable(describe(pcurve.out.comp, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#'## Evidential value for effect type = priming
#'#### Permutation p-curve
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data.prim$label[!duplicated.random(data.prim$study)], replacement = "")))
}
pcurve.out.prim <- data.frame(pcurve)
colnames(pcurve.out.prim) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.prim, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")


#'## Meta-analysis
#'#### Compensatory vs. priming effect type
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
ma <- rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], method = "REML", random = ~ 1|study/result, mods = ~ effect.type)
robma <- robust.rma.mv(x = ma, cluster = data$study[!is.na(data$effect.type)])

ma.comp <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], method = "ML", random = ~ 1|study/result, subset = effect.type == 1), cluster = data$study[!is.na(data$effect.type)])
ma.prim <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], method = "ML", random = ~ 1|study/result, subset = effect.type == 2), cluster = data$study[!is.na(data$effect.type)])

tau.comp <- sqrt(sum(ma.comp$sigma2))
tau.prim <- sqrt(sum(ma.prim$sigma2))

result <- data.frame(meta = c("Compensatory","Priming"),
                     k = c(ma.comp$k, ma.prim$k),
                     estimate = round(c(coef(ma.comp), coef(ma.prim)), 3), stderror = round(c(ma.comp$se, ma.prim$se), 3),
                     tau = round(c(sqrt(sum(ma.comp$sigma2)), sqrt(sum(ma.prim$sigma2))) ,3))
kable(result, "html", digits = 3) %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### Wald test p-value
#'Testing the difference in the uncorrected MA estimates between effect types
comp.vs.prim.z <- with(result, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))
comp.vs.prim.p <- 2*pnorm(abs(as.numeric(comp.vs.prim.z)), lower.tail = F)
comp.vs.prim.p

#'#### 95% prediction interval for compensatory effect
pi.lb.comp <- round(predict.rma(ma.comp)$cr.lb, 3)
pi.ub.comp <- round(predict.rma(ma.comp)$cr.ub, 3)

#'As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between `r pi.lb.comp` and `r pi.ub.comp`.
#'Note that these are non-adjusted estimates. An unbiased newly conducted study will likely fall in an interval centered around PET-PEESE estimate with a similar CI width of `r pi.ub.comp - pi.lb.comp`.
#'

#'#### 95% prediction interval for priming effect
pi.lb.prim <- round(predict.rma(ma.prim)$cr.lb, 3)
pi.ub.prim <- round(predict.rma(ma.prim)$cr.ub, 3)

#'As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between `r pi.lb.prim` and `r pi.ub.prim`.
#'An unbiased newly conducted study will likely fall in an interval centered around PET-PEESE estimate with a similar CI width of `r pi.ub.prim - pi.lb.prim`.
#'

#'#### Relative heterogeneity
#'#### $I^2$ for Compensatory
W <- diag(1/data.comp$g.var.calc[!is.na(data.comp$g.var.calc)])
X <- model.matrix(ma.comp)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'$I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates.
#'
#'Here, total relative heterogeneity was
I2.comp<- {{round(100 * sum(ma.comp$sigma2) / (sum(ma.comp$sigma2) + (ma.comp$k-ma.comp$p)/sum(diag(P))), 2)}}
{{round(I2.comp, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * ma.comp$sigma2 / (sum(ma.comp$sigma2) + (ma.comp$k-ma.comp$p)/sum(diag(P))), 2)}}
#'%, respectively.

res.R <- rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], struct="UN", random = ~ 1|study/result, subset = effect.type == 1)
res.F <- rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], subset = effect.type == 1)
#'Jackson's approach to $I^2$ yields a relative heterogeneity estimate of
{{round(c(100 * (vcov(res.R)[1,1] - vcov(res.F)[1,1]) / vcov(res.R)[1,1]), 2)}}
#'%.
#'
#'
#'#### $I^2$ for Priming
W <- diag(1/data.prim$g.var.calc[!is.na(data.prim$g.var.calc)])
X <- model.matrix(ma.prim)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'Here, total relative heterogeneity was
I2.prim <- {{round(100 * sum(ma.prim$sigma2) / (sum(ma.prim$sigma2) + (ma.prim$k-ma.prim$p)/sum(diag(P))), 2)}}
{{round(I2.prim, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * ma.prim$sigma2 / (sum(ma.prim$sigma2) + (ma.prim$k-ma.prim$p)/sum(diag(P))), 2)}}
#'%, respectively.

res.R <- rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], struct="UN", random = ~ 1|study/result, subset = effect.type == 2)
res.F <- rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], subset = effect.type == 2)
#'Jackson's approach to $I^2$ yields a relative heterogeneity estimate of
{{round(c(100 * (vcov(res.R)[1,1] - vcov(res.F)[1,1]) / vcov(res.R)[1,1]), 2)}}
#'%.

#+eval = FALSE
# Confidence intervals and profile plots
confint(robma, digits=3)
par(mfrow=c(2,1))
profile(robma, sigma2=1, progbar = F)
profile(robma, sigma2=2, progbar = F)
# Both profile likelihood plots are peaked at the respective parameter estimates (as indicated by the vertical dotted lines) and the log likelihoods quickly decrease (i.e., become more negative) as the values of the components are moved away from the actual REML estimates. Hence, we can be fairly confident that both variance components are identifiable.

#+eval = FALSE
# Testing significance of variance components for Compensatory
ma0 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.comp[!is.na(data.comp$effect.type),], method = "ML", random = ~ 1|study/result, subset = effect.type == 1)
ma1 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.comp[!is.na(data.comp$effect.type),], method = "ML", random = ~ 1|study/result, sigma2 = c(NA,0), subset = effect.type == 1)
ma2 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.comp[!is.na(data.comp$effect.type),], method = "ML", random = ~ 1|study/result, sigma2 = c(0,NA), subset = effect.type == 1)
anova(ma1, ma0)
anova(ma2, ma0)

#+eval = FALSE
# Testing significance of variance components for Priming
ma0 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.prim[!is.na(data.prim$effect.type),], method = "ML", random = ~ 1|study/result, subset = effect.type == 2)
ma1 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.prim[!is.na(data.prim$effect.type),], method = "ML", random = ~ 1|study/result, sigma2 = c(NA,0), subset = effect.type == 2)
ma2 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.prim[!is.na(data.prim$effect.type),], method = "ML", random = ~ 1|study/result, sigma2 = c(0,NA), subset = effect.type == 2)
anova(ma1, ma0)
anova(ma2, ma0)

#'#### Intra-class correlation of underlying true effects
#'For compensatory
round(ma.comp$sigma2[1] / sum(ma.comp$sigma2), 2)
#'For priming
round(ma.prim$sigma2[1] / sum(ma.prim$sigma2), 2)

#'### Contour enhanced funnel plot
par(mfrow=c(1,2))
funnel(ma.comp, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g", ylim = c(0, 0.6), xlim = c(-1, 1.7), steps = 7, digits = c(1,2))
title("Compensatory", cex.main = 1)
funnel(ma.prim, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Hedges' g", ylim = c(0, 0.6), xlim = c(-1, 1.7), steps = 7, digits = c(1,2))
title("Priming", cex.main = 1)
funnel.comp.prim <- recordPlot()

#'**Correlation between the ES and precision**
#'
#'For compensatory
with(data.comp[!is.na(data.comp$g.calc),], cor(g.calc, sqrt(g.var.calc), method = "kendall"))
#'For priming
with(data.prim[!is.na(data.prim$g.calc),], cor(g.calc, sqrt(g.var.calc), method = "kendall"))

#'### Forest plot
par(mar=c(4,4,1,2), mfrow=c(1,2))
forest(x = data.comp$g.calc, vi = data.comp$g.var.calc, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
       xlim = c(-1.5,2.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
       subset = order(data.comp$g.var.calc),        ### order by size of yi
       slab = NA, annotate = FALSE, ### remove study labels and annotations
       efac = 0,                  ### remove vertical bars at end of CIs
       pch = 19,                  ### changing point symbol to filled circle
       col = "gray40",            ### change color of points/CIs
       psize = 1.5,                 ### increase point size
       cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
       lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Compensatory", cex.main = 1)
addpoly(ma.comp, row = 0, mlab = "", cex = 1, annotate = F)

forest(x = data.prim$g.calc, vi = data.prim$g.var.calc, at = c(-1, -.5, 0, .5, 1, 1.5, 2, 2.5),
       xlim = c(-1.5,2.5), alim = c(-1, 2.5), xlab = "Hedges' g",        ### adjust horizontal plot region limits
       subset = order(data.prim$g.var.calc),        ### order by size of yi
       slab = NA, annotate = FALSE, ### remove study labels and annotations
       efac = 0,                  ### remove vertical bars at end of CIs
       pch = 19,                  ### changing point symbol to filled circle
       col = "gray40",            ### change color of points/CIs
       psize = 3,                 ### increase point size
       cex.lab=.7, cex.axis=.7,   ### increase size of x-axis title/labels
       lty=c("solid","blank"))  ### remove horizontal line at top of plot
title("Priming", cex.main = 1)
addpoly(ma.prim, row = -1, mlab = "", cex = 1, annotate = F)


forest.overall <- recordPlot()

################

#'## Small-study effects correction for compensatory effect
#'
#'### 3-parameter selection model for compensatory effect
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.comp <- threePSM.est(data.comp[!is.na(data.comp$g.calc),]$g.calc, data.comp[!is.na(data.comp$g.calc),]$g.var.calc)
kable(ThreePSM.comp[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'### PET-PEESE for compensatory effect
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.comp <- with(data.comp[!is.na(data.comp$g.calc),], pet.peese(g.calc, g.var.calc, study, result))
kable(pp.comp, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'**PET-PEESE plot for compensatory effect**
#'
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.comp$value[4] < .05 & ThreePSM.comp$value[1] > 0)
{plot(data.comp$g.var.calc, data.comp$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .12), xaxs="i")} else {plot(sqrt(data.comp$g.var.calc), data.comp$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.comp$value[4] < .05 & ThreePSM.comp$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#+include = FALSE
#'#### p-uniform for compensatory effect
#'Additional bias-corrected estimate. Because it's far less precise than PET-PEESE, when the n of studies is small, look just at the CI width and p-value.
puniform(yi = data.comp[!is.na(data.comp$g.calc),]$g.calc, vi = data.comp[!is.na(data.comp$g.calc),]$g.var.calc, alpha = .05, side = "right", method = "P")

#'## Small-study effects correction for priming effect
#'
#'### 3-parameter selection model for priming effect
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.prim <- threePSM.est(data.prim$g.calc, data.prim$g.var.calc)
kable(ThreePSM.prim[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'### PET-PEESE for priming effect
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.prim <- with(data.prim[!is.na(data.prim$g.calc),], pet.peese(g.calc, g.var.calc, study, result))
kable(pp.prim, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'**PET-PEESE plot for priming effect**
#'
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.prim$value[4] < .05 & ThreePSM.prim$value[1] > 0)
{plot(data.prim$g.var.calc, data.prim$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27), xaxs="i")} else {plot(sqrt(data.prim$g.var.calc), data.prim$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.prim$value[4] < .05 & ThreePSM.prim$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#+include = FALSE
#'#### p-uniform for priming effect
#'Additional bias-corrected estimate. Because it's far less precise than PET-PEESE, when the n of studies is small, look just at the CI width and p-value. Leaving out study id# 211 because p-uniform won't converge due to huge variance.
puniform(yi = data.prim$g.calc[!(data.prim$result %in% c(211)) & !is.na(data.prim$g.calc)], vi = data.prim$g.var.calc[!(data.prim$result %in% c(211)) & !is.na(data.prim$g.calc)], alpha = .05, side = "right", method = "P")

#########################################################################
#########################################################################
#########################################################################

# Recode the temperature manipulation into Y/N
dat$Physical.Temperature.Manipulation. <- relevel(factor(dat$Physical.Temperature.Manipulation.), ref="N")
data <- dat[dat$PA.NA. != "Y" & !is.na(dat$g.calc),]
data.physN <- data[data$Physical.Temperature.Manipulation. == "N",]
data.physY <- data[data$Physical.Temperature.Manipulation. == "Y",]

#'# Physical manipulation
#'
#'## Evidential value for non-physical manipulation
#'#### Permutation p-curve

#Permutation p-curve
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Physical.Temperature.Manipulation. == "N",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Physical.Temperature.Manipulation. == "N",]$study)], replacement = "")))
}
pcurve.out.physN <- data.frame(pcurve)
colnames(pcurve.out.physN) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
kable(describe(pcurve.out.physN, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#'## Evidential value for physical manipulation
#'#### Permutation p-curve
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Physical.Temperature.Manipulation. == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Physical.Temperature.Manipulation. == "Y",]$study)], replacement = "")))
}
pcurve.out.physY <- data.frame(pcurve)
colnames(pcurve.out.physY) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
kable(describe(pcurve.out.physY, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#'## Meta-analysis
#'#### Physical vs. non-physical manipulation
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate

ma.physN <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Physical.Temperature.Manipulation.),], method = "ML", random = ~ 1|study/result, subset = Physical.Temperature.Manipulation. == "N"), cluster = data$study[!is.na(data$Physical.Temperature.Manipulation.)])
ma.physY <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Physical.Temperature.Manipulation.),], method = "ML", random = ~ 1|study/result, subset = Physical.Temperature.Manipulation. == "Y"), cluster = data$study[!is.na(data$Physical.Temperature.Manipulation.)])

tau.physN <- sqrt(sum(ma.physN$sigma2))
tau.physY <- sqrt(sum(ma.physY$sigma2))

result <- data.frame(meta = c("(Physical manipulation = N)","(Physical manipulation = Y)"),
                     k = c(ma.physN$k, ma.physY$k),
                     estimate = round(c(coef(ma.physN), coef(ma.physY)), 3),
                     stderror = round(c(ma.physN$se, ma.physY$se), 3),
                     tau = round(c(sqrt(sum(ma.physN$sigma2)), sqrt(sum(ma.physY$sigma2))) ,3))

kable(result, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### 95% prediction interval for non-physical manipulation
pi.lb.physN <- round(predict.rma(ma.physN)$cr.lb, 3)
pi.ub.physN <- round(predict.rma(ma.physN)$cr.ub, 3)

#'As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between `r pi.lb.physN` and `r pi.ub.physN`.
#'Note that these are non-adjusted estimates. An unbiased newly conducted study will likely fall in an interval centered around PET-PEESE estimate with a similar CI width of `r pi.ub.physN - pi.lb.physN`.
#'

#'#### 95% prediction interval for physical manipulation
pi.lb.physY <- round(predict.rma(ma.physY)$cr.lb, 3)
pi.ub.physY <- round(predict.rma(ma.physY)$cr.ub, 3)

#'As an approximation, in 95% of cases the true effect in a new *published* study can be expected to fall between `r pi.lb.physY` and `r pi.ub.physY`.
#'An unbiased newly conducted study will likely fall in an interval centered around PET-PEESE estimate with a similar CI width of `r pi.ub.physY - pi.lb.physY`.
#'

#'#### Wald test p-value
#'Testing the difference in the uncorrected MA estimates between physical and non-physical manipulation
physN.vs.physY.z <- with(result, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))
physN.vs.physY.p <- 2*pnorm(abs(as.numeric(physN.vs.physY.z)), lower.tail = F)
physN.vs.physY.p

#'#### Relative heterogeneity
#'#### $I^2$ for non-physical manipulation
W <- diag(1/data.physN$g.var.calc[!is.na(data.physN$g.var.calc)])
X <- model.matrix(ma.physN)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'$I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates.
#'
#'Here, total relative heterogeneity was
I2.physN <- {{round(100 * sum(ma.physN$sigma2) / (sum(ma.physN$sigma2) + (ma.physN$k-ma.physN$p)/sum(diag(P))), 2)}}
{{round(I2.physN, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * ma.physN$sigma2 / (sum(ma.physN$sigma2) + (ma.physN$k-ma.physN$p)/sum(diag(P))), 2)}}
#'%, respectively.

res.R <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physN, struct="UN", random = ~ 1|study/result)
res.F <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physN)
#'Jackson's approach to $I^2$ yields a relative heterogeneity estimate of
{{round(c(100 * (vcov(res.R)[1,1] - vcov(res.F)[1,1]) / vcov(res.R)[1,1]), 2)}}
#'%.
#'

#'#### $I^2$ for physical manipulation
W <- diag(1/data.physY$g.var.calc[!is.na(data.physY$g.var.calc)])
X <- model.matrix(ma.physY)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'$I^2$ represents the ratio of true heterogeneity to total variance across the observed effect estimates.
#'
#'Here, total relative heterogeneity was
I2.physY <- {{round(100 * sum(ma.physY$sigma2) / (sum(ma.physY$sigma2) + (ma.physY$k-ma.physY$p)/sum(diag(P))), 2)}}
{{round(I2.physY, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * ma.physY$sigma2 / (sum(ma.physY$sigma2) + (ma.physY$k-ma.physY$p)/sum(diag(P))), 2)}}
#'%, respectively.

res.R <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physY, struct="UN", random = ~ 1|study/result)
res.F <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physY)
#'Jackson's approach to $I^2$ yields a relative heterogeneity estimate of
{{round(c(100 * (vcov(res.R)[1,1] - vcov(res.F)[1,1]) / vcov(res.R)[1,1]), 2)}}
#'%.

#+eval = FALSE
#Testing significance of variance components for physical manipulation == N
ma0 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physN, method = "ML", random = ~ 1|study/result)
ma1 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physN, method = "ML", random = ~ 1|study/result, sigma2 = c(NA,0))
ma2 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physN, method = "ML", random = ~ 1|study/result, sigma2 = c(0,NA))
anova(ma1, ma0)
anova(ma2, ma0)

#+eval = FALSE
#Testing significance of variance components for physical manipulation == Y
ma0 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physY, method = "ML", random = ~ 1|study/result)
ma1 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physY, method = "ML", random = ~ 1|study/result, sigma2 = c(NA,0))
ma2 <- rma.mv(yi = g.calc, V = g.var.calc, data = data.physY, method = "ML", random = ~ 1|study/result, sigma2 = c(0,NA))
anova(ma1, ma0)
anova(ma2, ma0)

#'#### Intra-class correlation of underlying true effects
#'For physical manipulation == N
round(ma.physN$sigma2[1] / sum(ma.physN$sigma2), 2)
#'For physical manipulation == Y
round(ma.physY$sigma2[1] / sum(ma.physY$sigma2), 2)

#'### Contour enhanced funnel plot
par(mfrow=c(1,2), mar=c(6,2,6,2))
funnel(ma.physN, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor physical manipulation == N", ylim = c(0, 0.6), xlim = c(-1, 1.8))
funnel(ma.physY, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor physical manipulation == Y", ylim = c(0, 0.6), xlim = c(-1, 1.8))
funnel.phys <- recordPlot()

#'Correlation between the ES and precision
#'
#'For physical manipulation == N
with(data.physN[!is.na(data.physN$g.calc),], cor(g.calc, sqrt(g.var.calc), method = "kendall"))
#'For physical manipulation == Y
with(data.physY[!is.na(data.physY$g.calc),], cor(g.calc, sqrt(g.var.calc), method = "kendall"))

#'### Forest plot
par(mar=c(4,4,1,2), mfrow=c(1,2))
forest(data.physN$g.calc, vi = data.physN$g.var.calc, subset=order(data.physN$g.calc), slab = data.physN$result, addcred = T)
title("Non-physical manipulation")
addpoly(ma.physN, row = 0, mlab = "", cex = 1, annotate = F)
forest(data.physY$g.calc, vi = data.physY$g.var.calc, subset=order(data.physY$g.calc), slab = data.physY$result, addcred = T)
title("Physical manipulation")
addpoly(ma.physY, row = 0, mlab = "", cex = 1, annotate = F)

#'## Small-study effects correction for non-physical manipulation
#'
#'### 3-parameter selection model for non-physical manipulation
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.physN <- threePSM.est(data.physN[!is.na(data.physN$g.calc),]$g.calc, data.physN[!is.na(data.physN$g.calc),]$g.var.calc)
kable(ThreePSM.physN[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'### PET-PEESE for non-physical manipulation
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.physN <- with(data.physN[!is.na(data.physN$g.calc),], pet.peese(g.calc, g.var.calc, study, result))
kable(pp.physN, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'**PET-PEESE plot for non-physical manipulation**
#'
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.physN$value[4] < .05 & ThreePSM.physN$value[1] > 0)
{plot(data.physN$g.var.calc, data.physN$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .23), xaxs="i")} else {plot(sqrt(data.physN$g.var.calc), data.physN$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.physN$value[4] < .05 & ThreePSM.physN$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#+include = FALSE
#'#### p-uniform for non-physical manipulation
#'Additional bias-corrected estimate. Because it's far less precise than PET-PEESE, when the n of studies is small, look just at the CI width and p-value. Leaving out study id# 211 because p-uniform won't converge due to huge variance.
puniform(yi = data.physN[!(data.physN$result %in% c(211)) & !is.na(data.physN$g.calc),]$g.calc, vi = data.physN[!(data.physN$result %in% c(211)) & !is.na(data.physN$g.calc),]$g.var.calc, alpha = .05, side = "right", method = "P")

#'## Small-study effects correction for physical manipulation
#'
#'### 3-parameter selection model for physical manipulation
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.physY <- threePSM.est(data.physY$g.calc, data.physY$g.var.calc)
kable(ThreePSM.physY[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'### PET-PEESE for physical manipulation
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.physY <- with(data.physY[!is.na(data.physY$g.calc),], pet.peese(g.calc, g.var.calc, study, result))
kable(pp.physY, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'**PET-PEESE plot for physical manipulation**
#'
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.physY$value[4] < .05 & ThreePSM.physY$value[1] > 0)
{plot(data.physY$g.var.calc, data.physY$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27), xaxs="i")} else {plot(sqrt(data.physY$g.var.calc), data.physY$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.physY$value[4] < .05 & ThreePSM.physY$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#+include = FALSE
#'#### p-uniform for physical manipulation
#'Additional bias-corrected estimate. Because it's far less precise than PET-PEESE, when the n of studies is small, look just at the CI width and p-value.
puniform(yi = data.physY$g.calc[!is.na(data.physY$g.calc)], vi = data.physY$g.var.calc[!is.na(data.physY$g.calc)], alpha = .05, side = "right", method = "P")

#########################################################################
#########################################################################

#'# Method effects
data <- dat[dat$PA.NA. != "Y" & !is.na(dat$g.calc),]

#'## Evidential value
#'### Permutation p-curve for Visual.Verbal.Temperature.Prime.
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Visual.Verbal.Temperature.Prime. == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Visual.Verbal.Temperature.Prime. == "Y",]$study)], replacement = "")))
}
pcurve.out.VVTP <- data.frame(pcurve)
colnames(pcurve.out.VVTP) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
#'
kable(describe(pcurve.out.VVTP, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Outside.Temperature.
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Outside.Temperature. == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Outside.Temperature. == "Y",]$study)], replacement = "")))
}
pcurve.out.OT <- data.frame(pcurve)
colnames(pcurve.out.OT) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.OT, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Temperature.Estimate.
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Temperature.Estimate. == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Temperature.Estimate. == "Y",]$study)], replacement = "")))
}
pcurve.out.TE <- data.frame(pcurve)
colnames(pcurve.out.TE) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.TE, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Subjective.Warmth.Judgment
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Subjective.Warmth.Judgment == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Subjective.Warmth.Judgment == "Y",]$study)], replacement = "")))
}
pcurve.out.SWJ <- data.frame(pcurve)
colnames(pcurve.out.SWJ) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.SWJ, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#'## Meta-analysis for different methods
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
#'
#For Visual.Verbal.Temperature.Prime.
ma.VVTP <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Visual.Verbal.Temperature.Prime.),], method = "REML", random = ~ 1|study/result, subset = Visual.Verbal.Temperature.Prime. == "Y"), cluster = data$study[!is.na(data$Visual.Verbal.Temperature.Prime.)])
#$I^2$
W <- diag(1/data[data$Visual.Verbal.Temperature.Prime. == "Y",]$g.var.calc)
X <- model.matrix(ma.VVTP)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.VVTP <- {{round(100 * sum(ma.VVTP$sigma2) / (sum(ma.VVTP$sigma2) + (ma.VVTP$k-ma.VVTP$p)/sum(diag(P))), 2)}}

#For Outside.Temperature.
ma.OT <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Outside.Temperature.),], method = "REML", random = ~ 1|study/result, subset = Outside.Temperature. == "Y"), cluster = data$study[!is.na(data$Outside.Temperature.)])
#$I^2$
W <- diag(1/data[data$Outside.Temperature. == "Y",]$g.var.calc)
X <- model.matrix(ma.OT)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.OT <- {{round(100 * sum(ma.OT$sigma2) / (sum(ma.OT$sigma2) + (ma.OT$k-ma.OT$p)/sum(diag(P))), 2)}}

#For Temperature.Estimate.
ma.TE <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Temperature.Estimate.),], method = "REML", random = ~ 1|study/result, subset = Temperature.Estimate. == "Y"), cluster = data$study[!is.na(data$Temperature.Estimate.)])
#$I^2$
W <- diag(1/data[data$Temperature.Estimate. == "Y",]$g.var.calc)
X <- model.matrix(ma.TE)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.TE <- {{round(100 * sum(ma.TE$sigma2) / (sum(ma.TE$sigma2) + (ma.TE$k-ma.TE$p)/sum(diag(P))), 2)}}

#For Subjective.Warmth.Judgment
ma.SWJ <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Subjective.Warmth.Judgment),], method = "REML", random = ~ 1|study/result, subset = Subjective.Warmth.Judgment == "Y"), cluster = data$study[!is.na(data$Subjective.Warmth.Judgment)])
#$I^2$
W <- diag(1/data[data$Subjective.Warmth.Judgment == "Y",]$g.var.calc)
X <- model.matrix(ma.SWJ)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.SWJ <- {{round(100 * sum(ma.SWJ$sigma2) / (sum(ma.SWJ$sigma2) + (ma.SWJ$k-ma.SWJ$p)/sum(diag(P))), 2)}}

######################################
method.result <- data.frame(meta = c("Visual.Verbal.Temperature.Prime.","Outside.Temperature.", "Temperature.Estimate.", "Subjective.Warmth.Judgment"),
                     k = c(ma.VVTP$k, ma.OT$k, ma.TE$k, ma.SWJ$k),
                     estimate = round(c(coef(ma.VVTP), coef(ma.OT), coef(ma.TE), coef(ma.SWJ)), 3),
                     stderror = round(c(ma.VVTP$se, ma.OT$se, ma.TE$se, ma.SWJ$se), 3),
                     tau = round(c(sqrt(sum(ma.VVTP$sigma2)), sqrt(sum(ma.OT$sigma2)), sqrt(sum(ma.TE$sigma2)), sqrt(sum(ma.SWJ$sigma2))) ,3),
                     I2 = c(I2.VVTP, I2.OT, I2.TE, I2.SWJ)
                     )

kable(method.result, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

######################################

#'### Contour enhanced funnel plots
par(mfrow=c(2,2), mar=c(4,3,4,3))
funnel(ma.VVTP, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor Visual.Verbal.Temperature.Prime.", ylim = c(0, 0.43), xlim = c(-.5, 1.5))
funnel(ma.OT, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor Outside.Temperature.", ylim = c(0, 0.43), xlim = c(-.5, 1.5))
funnel(ma.TE, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor Temperature.Estimate.", ylim = c(0, 0.43), xlim = c(-.5, 1.5))
funnel(ma.SWJ, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor Subjective.Warmth.Judgment", ylim = c(0, 0.43), xlim = c(-.5, 1.5))
funnel.method <- recordPlot()

#'### Forest plots
par(mar=c(4,4,1,2), mfrow=c(2,2))
forest(data[data$Visual.Verbal.Temperature.Prime. == "Y", ]$g.calc, vi = data[data$Visual.Verbal.Temperature.Prime. == "Y",]$g.var.calc, subset=order(data[data$Visual.Verbal.Temperature.Prime. == "Y",]$g.calc), slab = data[data$Visual.Verbal.Temperature.Prime. == "Y",]$result, addcred = T)
title("Visual.Verbal.Temperature.Prime.")
addpoly(ma.VVTP, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Outside.Temperature. == "Y",]$g.calc, vi = data[data$Outside.Temperature. == "Y",]$g.var.calc, subset=order(data[data$Outside.Temperature. == "Y",]$g.calc), slab = data[data$Outside.Temperature. == "Y",]$result, addcred = T)
title("Outside.Temperature.")
addpoly(ma.OT, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Temperature.Estimate. == "Y",]$g.calc, vi = data[data$Temperature.Estimate. == "Y",]$g.var.calc, subset=order(data[data$Temperature.Estimate. == "Y",]$g.calc), slab = data[data$Temperature.Estimate. == "Y",]$result, addcred = T)
title("Temperature.Estimate.")
addpoly(ma.TE, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Subjective.Warmth.Judgment == "Y",]$g.calc, vi = data[data$Subjective.Warmth.Judgment == "Y",]$g.var.calc, subset=order(data[data$Subjective.Warmth.Judgment == "Y",]$g.calc), slab = data[data$Subjective.Warmth.Judgment == "Y",]$result, addcred = T)
title("Subjective.Warmth.Judgment")
addpoly(ma.SWJ, row = 0, mlab = "", cex = 1, annotate = F)
######################################

#'## Small-study effects correction for different methods
#'
#### 3-parameter selection model
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
#
#For Visual.Verbal.Temperature.Prime.
ThreePSM.VVTP <- threePSM.est(data[data$Visual.Verbal.Temperature.Prime. == "Y", ]$g.calc, data[data$Visual.Verbal.Temperature.Prime. == "Y",]$g.var.calc)
#For Outside.Temperature.
ThreePSM.OT <- threePSM.est(data[data$Outside.Temperature. == "Y",]$g.calc, data[data$Outside.Temperature. == "Y",]$g.var.calc)
#For Temperature.Estimate.
ThreePSM.TE <- threePSM.est(data[data$Temperature.Estimate. == "Y",]$g.calc, data[data$Temperature.Estimate. == "Y",]$g.var.calc)
#For Subjective.Warmth.Judgment
ThreePSM.SWJ <- threePSM.est(data[data$Subjective.Warmth.Judgment == "Y",]$g.calc, data[data$Subjective.Warmth.Judgment == "Y",]$g.var.calc)

result <- rbind(Visual.Verbal.Temperature.Prime. = ThreePSM.VVTP[c(1, 4, 5, 6),], Outside.Temperature. = ThreePSM.OT[c(1, 4, 5, 6),],
                Temperature.Estimate. = ThreePSM.TE[c(1, 4, 5, 6),], Subjective.Warmth.Judgment = ThreePSM.SWJ[c(1, 4, 5, 6),])

kable(result, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")


#'### PET-PEESE
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
#'
#'y-axis intercept represents the estimated bias-corrected ES.
#'
#For Visual.Verbal.Temperature.Prime.
par(mar=c(4,4,1,2), mfrow=c(2,2))
PP.VVTP <- with(data[data$Visual.Verbal.Temperature.Prime. == "Y", ], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.VVTP$value[4] < .05 & ThreePSM.VVTP$value[1] > 0)
{plot(data[data$Visual.Verbal.Temperature.Prime. == "Y", ]$g.var.calc, data[data$Visual.Verbal.Temperature.Prime. == "Y", ]$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")} else {plot(sqrt(data[data$Visual.Verbal.Temperature.Prime. == "Y", ]$g.var.calc), data[data$Visual.Verbal.Temperature.Prime. == "Y", ]$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .45), ylim = c(0, 1.25), xaxs="i")}
abline((if(ThreePSM.VVTP$value[4] < .05 & ThreePSM.VVTP$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#For Outside.Temperature.
PP.OT <- with(data[data$Outside.Temperature. == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.OT$value[4] < .05 & ThreePSM.OT$value[1] > 0)
{plot(data[data$Outside.Temperature. == "Y",]$g.var.calc, data[data$Outside.Temperature. == "Y",]$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .08), xaxs="i")} else {plot(sqrt(data[data$Outside.Temperature. == "Y",]$g.var.calc), data[data$Outside.Temperature. == "Y",]$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.OT$value[4] < .05 & ThreePSM.OT$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#For Temperature.Estimate.
PP.TE <- with(data[data$Temperature.Estimate. == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.TE$value[4] < .05 & ThreePSM.TE$value[1] > 0)
{plot(data[data$Temperature.Estimate. == "Y",]$g.var.calc, data[data$Temperature.Estimate. == "Y",]$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")} else {plot(sqrt(data[data$Temperature.Estimate. == "Y",]$g.var.calc), data[data$Temperature.Estimate. == "Y",]$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .45), ylim = c(-0.1, 1.5), xaxs="i")}
abline((if(ThreePSM.TE$value[4] < .05 & ThreePSM.TE$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#For Subjective.Warmth.Judgment
PP.SWJ <- with(data[data$Subjective.Warmth.Judgment == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.SWJ$value[4] < .05 & ThreePSM.SWJ$value[1] > 0)
{plot(data[data$Subjective.Warmth.Judgment == "Y",]$g.var.calc, data[data$Subjective.Warmth.Judgment == "Y",]$g.calc, main="PEESE", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")} else {plot(sqrt(data[data$Subjective.Warmth.Judgment == "Y",]$g.var.calc), data[data$Subjective.Warmth.Judgment == "Y",]$g.calc, main="PET", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, 0.2), xaxs="i")}
abline((if(ThreePSM.SWJ$value[4] < .05 & ThreePSM.SWJ$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

PP.categ.results <- rbind(Visual.Verbal.Temperature.Prime. = PP.VVTP, Outside.Temperature. = PP.OT,
                Temperature.Estimate. = PP.TE, Subjective.Warmth.Judgment = PP.SWJ)
colnames(PP.categ.results)[1] <- "PET-PEESE estimate"

kable(PP.categ.results, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#########################################################################
#########################################################################
#########################################################################

#'# Category effects
#'
#'## Evidential value
#'### Permutation p-curve for Category..Emotion
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Category..Emotion == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Category..Emotion == "Y",]$study)], replacement = "")))
}
pcurve.out.emot <- data.frame(pcurve)
colnames(pcurve.out.emot) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
#'
kable(describe(pcurve.out.emot, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Category...Interpersonal
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Category...Interpersonal == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Category...Interpersonal == "Y",]$study)], replacement = "")))
}
pcurve.out.interp <- data.frame(pcurve)
colnames(pcurve.out.interp) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.interp, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Category..Person.Perception
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Category..Person.Perception == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Category..Person.Perception == "Y",]$study)], replacement = "")))
}
pcurve.out.pers.perc <- data.frame(pcurve)
colnames(pcurve.out.pers.perc) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.pers.perc, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Category..Group.Processes
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Category..Group.Processes == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Category..Group.Processes == "Y",]$study)], replacement = "")))
}
pcurve.out.group.process <- data.frame(pcurve)
colnames(pcurve.out.group.process) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.group.process, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Category..Self.Regulation
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Category..Self.Regulation == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Category..Self.Regulation == "Y",]$study)], replacement = "")))
}
pcurve.out.self.reg <- data.frame(pcurve)
colnames(pcurve.out.self.reg) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.self.reg, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Category..Cognitive.Processes
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Category..Cognitive.Processes == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Category..Cognitive.Processes == "Y",]$study)], replacement = "")))
}
pcurve.out.cogn.proc <- data.frame(pcurve)
colnames(pcurve.out.cogn.proc) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.cogn.proc, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#############

#'### Permutation p-curve for Category..Economic.Decision.Making
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Category..Economic.Decision.Making == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Category..Economic.Decision.Making == "Y",]$study)], replacement = "")))
}
pcurve.out.econ.decis.making <- data.frame(pcurve)
colnames(pcurve.out.econ.decis.making) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'
kable(describe(pcurve.out.econ.decis.making, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

############

#'## Meta-analysis for effect categories
#'
#'k = number of studies; sqrt in "Variance components" = tau, the standard deviation of true effects; estimate in "Model results" = naive MA estimate
#
#For Category..Emotion
ma.emot <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Category..Emotion),], method = "REML", random = ~ 1|study/result, subset = Category..Emotion == "Y"), cluster = data$study[!is.na(data$Category..Emotion)])
#$I^2$
W <- diag(1/data[data$Category..Emotion == "Y",]$g.var.calc)
X <- model.matrix(ma.emot)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.emot <- {{round(100 * sum(ma.emot$sigma2) / (sum(ma.emot$sigma2) + (ma.emot$k-ma.emot$p)/sum(diag(P))), 2)}}

#For Category...Interpersonal
ma.interp <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Category...Interpersonal),], method = "REML", random = ~ 1|study/result, subset = Category...Interpersonal == "Y"), cluster = data$study[!is.na(data$Category...Interpersonal)])
#$I^2$
W <- diag(1/data[data$Category...Interpersonal == "Y",]$g.var.calc)
X <- model.matrix(ma.interp)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.interp <- {{round(100 * sum(ma.interp$sigma2) / (sum(ma.interp$sigma2) + (ma.interp$k-ma.interp$p)/sum(diag(P))), 2)}}

#For Category..Person.Perception
ma.pers.perc <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Category..Person.Perception),], method = "REML", random = ~ 1|study/result, subset = Category..Person.Perception == "Y"), cluster = data$study[!is.na(data$Category..Person.Perception)])
#$I^2$
W <- diag(1/data[data$Category..Person.Perception == "Y",]$g.var.calc)
X <- model.matrix(ma.pers.perc)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.pers.perc <- {{round(100 * sum(ma.pers.perc$sigma2) / (sum(ma.pers.perc$sigma2) + (ma.pers.perc$k-ma.pers.perc$p)/sum(diag(P))), 2)}}

#For Category..Group.Processes
ma.group.process <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Category..Group.Processes),], method = "REML", random = ~ 1|study/result, subset = Category..Group.Processes == "Y"), cluster = data$study[!is.na(data$Category..Group.Processes)])
#$I^2$
W <- diag(1/data[data$Category..Group.Processes == "Y",]$g.var.calc)
X <- model.matrix(ma.group.process)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.group.process <- {{round(100 * sum(ma.group.process$sigma2) / (sum(ma.group.process$sigma2) + (ma.group.process$k-ma.group.process$p)/sum(diag(P))), 2)}}

#For Category..Self.Regulation
ma.self.reg <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Category..Self.Regulation),], method = "REML", random = ~ 1|study/result, subset = Category..Self.Regulation == "Y"), cluster = data$study[!is.na(data$Category..Self.Regulation)])
#$I^2$
W <- diag(1/data[data$Category..Self.Regulation == "Y",]$g.var.calc)
X <- model.matrix(ma.self.reg)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.self.reg <- {{round(100 * sum(ma.self.reg$sigma2) / (sum(ma.self.reg$sigma2) + (ma.self.reg$k-ma.self.reg$p)/sum(diag(P))), 2)}}

#For Category..Cognitive.Processes
ma.cogn.proc <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Category..Cognitive.Processes),], method = "REML", random = ~ 1|study/result, subset = Category..Cognitive.Processes == "Y"), cluster = data$study[!is.na(data$Category..Cognitive.Processes)])
#$I^2$
W <- diag(1/data[data$Category..Cognitive.Processes == "Y",]$g.var.calc)
X <- model.matrix(ma.cogn.proc)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.cogn.proc <- {{round(100 * sum(ma.cogn.proc$sigma2) / (sum(ma.cogn.proc$sigma2) + (ma.cogn.proc$k-ma.cogn.proc$p)/sum(diag(P))), 2)}}

#For Category..Economic.Decision.Making
ma.econ.decis.making <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Category..Economic.Decision.Making),], method = "REML", random = ~ 1|study/result, subset = Category..Economic.Decision.Making == "Y"), cluster = data$study[!is.na(data$Category..Economic.Decision.Making)])
#$I^2$
W <- diag(1/data[data$Category..Economic.Decision.Making == "Y",]$g.var.calc)
X <- model.matrix(ma.econ.decis.making)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
I2.econ.decis.making <- {{round(100 * sum(ma.econ.decis.making$sigma2) / (sum(ma.econ.decis.making$sigma2) + (ma.econ.decis.making$k-ma.econ.decis.making$p)/sum(diag(P))), 2)}}

######################################
category.result <- data.frame(meta = c("Category..Emotion","Category...Interpersonal", "Category..Person.Perception", "Category..Group.Processes", "Category..Self.Regulation", "Category..Cognitive.Processes", "Category..Economic.Decision.Making"),
                     k = c(ma.emot$k, ma.interp$k, ma.pers.perc$k, ma.group.process$k, ma.self.reg$k, ma.cogn.proc$k, ma.econ.decis.making$k),
                     estimate = round(c(coef(ma.emot), coef(ma.interp), coef(ma.pers.perc), coef(ma.group.process), coef(ma.self.reg), coef(ma.cogn.proc), coef(ma.econ.decis.making)), 3),
                     stderror = round(c(ma.emot$se, ma.interp$se, ma.pers.perc$se, ma.group.process$se, ma.self.reg$se, ma.cogn.proc$se, ma.econ.decis.making$se), 3),
                     tau = round(c(sqrt(sum(ma.emot$sigma2)), sqrt(sum(ma.interp$sigma2)), sqrt(sum(ma.pers.perc$sigma2)), sqrt(sum(ma.group.process$sigma2)),
                                   sqrt(sum(ma.self.reg$sigma2)), sqrt(sum(ma.cogn.proc$sigma2)), sqrt(sum(ma.econ.decis.making$sigma2))) ,3),
                     I2 = c(I2.emot, I2.interp, I2.pers.perc, I2.group.process, I2.self.reg, I2.cogn.proc, I2.econ.decis.making)
)

kable(category.result, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")
######################################

#'### Contour enhanced funnel plots
par(mfrow=c(4,2), mar=c(4,4,0,4))
funnel(ma.emot, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20,yaxis = "sei", xlab = "Emotion", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(ma.interp, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Interpersonal", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(ma.pers.perc, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Person.Perception", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(ma.group.process, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Group.Processes", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(ma.self.reg, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Self.Regulation", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(ma.cogn.proc, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Cognitive.Processes", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel(ma.econ.decis.making, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "Economic.Decision.Making", ylim = c(0, 0.51), xlim = c(-1, 1.7))
funnel.category <- recordPlot()

#'#### Forest plots
par(mar=c(4,4,1,2))
forest(data[data$Category..Emotion == "Y", ]$g.calc, vi = data[data$Category..Emotion == "Y",]$g.var.calc, subset=order(data[data$Category..Emotion == "Y",]$g.calc), slab = data[data$Category..Emotion == "Y",]$result, addcred = T)
title("Category..Emotion")
addpoly(ma.emot, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Category...Interpersonal == "Y",]$g.calc, vi = data[data$Category...Interpersonal == "Y",]$g.var.calc, subset=order(data[data$Category...Interpersonal == "Y",]$g.calc), slab = data[data$Category...Interpersonal == "Y",]$result, addcred = T)
title("Category...Interpersonal")
addpoly(ma.interp, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Category..Person.Perception == "Y",]$g.calc, vi = data[data$Category..Person.Perception == "Y",]$g.var.calc, subset=order(data[data$Category..Person.Perception == "Y",]$g.calc), slab = data[data$Category..Person.Perception == "Y",]$result, addcred = T)
title("Category..Person.Perception")
addpoly(ma.pers.perc, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Category..Group.Processes == "Y",]$g.calc, vi = data[data$Category..Group.Processes == "Y",]$g.var.calc, subset=order(data[data$Category..Group.Processes == "Y",]$g.calc), slab = data[data$Category..Group.Processes == "Y",]$result, addcred = T)
title("Category..Group.Processes")
addpoly(ma.group.process, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Category..Self.Regulation == "Y",]$g.calc, vi = data[data$Category..Self.Regulation == "Y",]$g.var.calc, subset=order(data[data$Category..Self.Regulation == "Y",]$g.calc), slab = data[data$Category..Self.Regulation == "Y",]$result, addcred = T)
title("Category..Self.Regulation")
addpoly(ma.self.reg, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Category..Cognitive.Processes == "Y",]$g.calc, vi = data[data$Category..Cognitive.Processes == "Y",]$g.var.calc, subset=order(data[data$Category..Cognitive.Processes == "Y",]$g.calc), slab = data[data$Category..Cognitive.Processes == "Y",]$result, addcred = T)
title("Category..Cognitive.Processes")
addpoly(ma.cogn.proc, row = 0, mlab = "", cex = 1, annotate = F)
forest(data[data$Category..Economic.Decision.Making == "Y",]$g.calc, vi = data[data$Category..Economic.Decision.Making == "Y",]$g.var.calc, subset=order(data[data$Category..Economic.Decision.Making == "Y",]$g.calc), slab = data[data$Category..Economic.Decision.Making == "Y",]$result, addcred = T)
title("Category..Economic.Decision.Making")
addpoly(ma.econ.decis.making, row = 0, mlab = "", cex = 1, annotate = F)
######################################

#'## Small-study effects correction for effect category
#'
#'### 3-parameter selection model
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
#
#For Category..Emotion
ThreePSM.emot <- threePSM.est(data[data$Category..Emotion == "Y", ]$g.calc, data[data$Category..Emotion == "Y",]$g.var.calc)
#For Category...Interpersonal
ThreePSM.interp <- threePSM.est(data[data$Category...Interpersonal == "Y",]$g.calc, data[data$Category...Interpersonal == "Y",]$g.var.calc)
#For Category..Person.Perception
ThreePSM.pers.perc <- threePSM.est(data[data$Category..Person.Perception == "Y",]$g.calc, data[data$Category..Person.Perception == "Y",]$g.var.calc)
#For Category..Group.Processes
ThreePSM.group.process <- threePSM.est(data[data$Category..Group.Processes == "Y",]$g.calc, data[data$Category..Group.Processes == "Y",]$g.var.calc)
#For Category..Self.Regulation
ThreePSM.self.reg <- threePSM.est(data[data$Category..Self.Regulation == "Y",]$g.calc, data[data$Category..Self.Regulation == "Y",]$g.var.calc)
#For Category..Cognitive.Processes
ThreePSM.cogn.proc <- threePSM.est(data[data$Category..Cognitive.Processes == "Y",]$g.calc, data[data$Category..Cognitive.Processes == "Y",]$g.var.calc)
#For Category..Economic.Decision.Making
ThreePSM.econ.decis.making <- threePSM.est(data[data$Category..Economic.Decision.Making == "Y",]$g.calc, data[data$Category..Economic.Decision.Making == "Y",]$g.var.calc)


result <- rbind(Category..Emotion = ThreePSM.emot[c(1, 4, 5, 6),], Category...Interpersonal = ThreePSM.interp[c(1, 4, 5, 6),],
                Category..Person.Perception = ThreePSM.pers.perc[c(1, 4, 5, 6),], Category..Group.Processes = ThreePSM.group.process[c(1, 4, 5, 6),],
                Category..Self.Regulation = ThreePSM.self.reg[c(1, 4, 5, 6),], Category..Cognitive.Processes = ThreePSM.cogn.proc[c(1, 4, 5, 6),],
                Category..Economic.Decision.Making = ThreePSM.econ.decis.making[c(1, 4, 5, 6),])
kable(result, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")


#'### PET-PEESE
#'
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
#'
#'y-axis intercept represents the estimated bias-corrected ES.
par(mar=c(4,4,1,2), mfrow=c(3,3))
#For Category..Emotion
PP.emot <- with(data[data$Category..Emotion == "Y", ], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.emot$value[4] < .05 & ThreePSM.emot$value[1] > 0)
{plot(data[data$Category..Emotion == "Y", ]$g.var.calc, data[data$Category..Emotion == "Y", ]$g.calc, main="PEESE for Emotion", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")} else {plot(sqrt(data[data$Category..Emotion == "Y", ]$g.var.calc), data[data$Category..Emotion == "Y", ]$g.calc, main="PET for Emotion", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .4), ylim = c(-.1, .9), xaxs="i")}
abline((if(ThreePSM.emot$value[4] < .05 & ThreePSM.emot$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
#For Category...Interpersonal
PP.interp <- with(data[data$Category...Interpersonal == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.interp$value[4] < .05 & ThreePSM.interp$value[1] > 0)
{plot(data[data$Category...Interpersonal == "Y",]$g.var.calc, data[data$Category...Interpersonal == "Y",]$g.calc, main="PEESE for Interpersonal", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .25), xaxs="i")} else {plot(sqrt(data[data$Category...Interpersonal == "Y",]$g.var.calc), data[data$Category...Interpersonal == "Y",]$g.calc, main="PET for Interpersonal", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.interp$value[4] < .05 & ThreePSM.interp$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
#For Category..Person.Perception
PP.pers.perc <- with(data[data$Category..Person.Perception == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.pers.perc$value[4] < .05 & ThreePSM.pers.perc$value[1] > 0)
{plot(data[data$Category..Person.Perception == "Y",]$g.var.calc, data[data$Category..Person.Perception == "Y",]$g.calc, main="PEESE for Person Perception", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .27), xaxs="i")} else {plot(sqrt(data[data$Category..Person.Perception == "Y",]$g.var.calc), data[data$Category..Person.Perception == "Y",]$g.calc, main="PET for Person Perception", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.pers.perc$value[4] < .05 & ThreePSM.pers.perc$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
#For Category..Group.Processes
PP.group.process <- with(data[data$Category..Group.Processes == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.group.process$value[4] < .05 & ThreePSM.group.process$value[1] > 0)
{plot(data[data$Category..Group.Processes == "Y",]$g.var.calc, data[data$Category..Group.Processes == "Y",]$g.calc, main="PEESE for Group Processes", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .25),xaxs="i")} else {plot(sqrt(data[data$Category..Group.Processes == "Y",]$g.var.calc), data[data$Category..Group.Processes == "Y",]$g.calc, main="PET for Group Processes", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.group.process$value[4] < .05 & ThreePSM.group.process$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
#For Category..Self.Regulation
PP.self.reg  <- with(data[data$Category..Self.Regulation == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.self.reg$value[4] < .05 & ThreePSM.self.reg$value[1] > 0)
{plot(data[data$Category..Self.Regulation == "Y",]$g.var.calc, data[data$Category..Self.Regulation == "Y",]$g.calc, main="PEESE for Self-regulation", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .2), xaxs="i")} else {plot(sqrt(data[data$Category..Self.Regulation == "Y",]$g.var.calc), data[data$Category..Self.Regulation == "Y",]$g.calc, main="PET for Self-regulation", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.self.reg$value[4] < .05 & ThreePSM.self.reg$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
#For Category..Cognitive.Processes
PP.cogn.proc <- with(data[data$Category..Cognitive.Processes == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.cogn.proc$value[4] < .05 & ThreePSM.cogn.proc$value[1] > 0)
{plot(data[data$Category..Cognitive.Processes == "Y",]$g.var.calc, data[data$Category..Cognitive.Processes == "Y",]$g.calc, main="PEESE for Cognitive Processes", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .2), xaxs="i")} else {plot(sqrt(data[data$Category..Cognitive.Processes == "Y",]$g.var.calc), data[data$Category..Cognitive.Processes == "Y",]$g.calc, main="PET for Cognitive Processes", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.cogn.proc$value[4] < .05 & ThreePSM.cogn.proc$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
#For Category..Economic.Decision.Making
PP.econ.decis.making <- with(data[data$Category..Economic.Decision.Making == "Y",], pet.peese(g.calc, g.var.calc, study, result))
if(ThreePSM.econ.decis.making$value[4] < .05 & ThreePSM.econ.decis.making$value[1] > 0)
{plot(data[data$Category..Economic.Decision.Making == "Y",]$g.var.calc, data[data$Category..Economic.Decision.Making == "Y",]$g.calc, main="PEESE for Economic Decision-making", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .2), xaxs="i")} else {plot(sqrt(data[data$Category..Economic.Decision.Making == "Y",]$g.var.calc), data[data$Category..Economic.Decision.Making == "Y",]$g.calc, main="PET for Economic Decision-making", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.econ.decis.making$value[4] < .05 & ThreePSM.econ.decis.making$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#'
PP.categ.results <- rbind(Category..Emotion = PP.emot, Category...Interpersonal = PP.interp,
                          Category..Person.Perception = PP.pers.perc, Category..Group.Processes = PP.group.process,
                          Category..Self.Regulation = PP.self.reg, Category..Cognitive.Processes = PP.cogn.proc,
                          Category..Economic.Decision.Making = PP.econ.decis.making)
colnames(PP.categ.results)[1] <- "PET-PEESE estimate"

kable(PP.categ.results, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

##################
##################
##################

#'# Meta-regression
#'
#'## Moderation by citations and IF
#'#### Overall effect moderated by citations and IF
ma.cit.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result, mods = ~ scale(Publication.Year) + scale(Citations.March.1.2016..GS.) + scale(H5.Index.GS.Journal.March.2016)), cluster = data$study)
ma.cit.mod

#'## Moderation by lattitude
#'#### Overall effect moderated by lattitude
ma.lat.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result, mods = ~ scale(Latitude.University..proxy.for.climate.)), cluster = data$study)
ma.lat.mod

#'#### Priming effects moderated by lattitude
ma.lat.prim.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], method = "REML", random = ~ 1|study/result, subset = effect.type == 2, mods = ~ scale(Latitude.University..proxy.for.climate.)), cluster = data$study[!is.na(data$effect.type)])
ma.lat.prim.mod

#'#### Compensatory effects moderated by lattitude
ma.lat.comp.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], method = "REML", random = ~ 1|study/result, subset = effect.type == 1, mods = ~ scale(Latitude.University..proxy.for.climate.)), cluster = data$study[!is.na(data$effect.type)])
ma.lat.comp.mod

#'#### Mood effects moderated by lattitude
ma.lat.mood.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data.mood, method = "REML", random = ~ 1|study/result, mods = ~ scale(Latitude.University..proxy.for.climate.)), cluster = data.mood$study)
ma.lat.mood.mod

#####################################
#'## Moderation by gender
#'#### Overall effect moderated by gender
ma.gender.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result, mods = ~ scale(gender.ratio)), cluster = data$study)
ma.gender.mod

#'#### Priming effects moderated by gender
ma.prim.gender.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], method = "REML", random = ~ 1|study/result, subset = effect.type == 2, mods = ~ scale(gender.ratio)), cluster = data$study[!is.na(data$effect.type)])
ma.prim.gender.mod

#'#### Compensatory effects moderated by gender
ma.comp.gender.mod <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$effect.type),], method = "REML", random = ~ 1|study/result, subset = effect.type == 1, mods = ~ scale(gender.ratio)), cluster = data$study[!is.na(data$effect.type)])
ma.comp.gender.mod

#####################################

#'## Overall effect moderated by published
#'

#'### Evidential value
#'#### Permutation p-curve for Unpublished
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data[data$published == "unpublished",]$label[!duplicated.random(data[data$published == "unpublished",]$study)], replacement = "")))
}

pcurve.out.unpub <- data.frame(pcurve)
colnames(pcurve.out.unpub) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
#'
kable(describe(pcurve.out.unpub, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#'### Evidential value
#'#### Permutation p-curve for Published
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = data[data$published == "published",]$label[!duplicated.random(data[data$published == "published",]$study)], replacement = "")))
}

pcurve.out.pub <- data.frame(pcurve)
colnames(pcurve.out.pub) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'P-Curve analysis combines the half and full p-curve to make inferences about evidential value. In particular, if the half p-curve test is right-skewed (halfp) with p<.05 or both the half and full test (fullp) are right-skewed with p < .1, then p-curve analysis indicates the presence of evidential value.
#'Similarly, p-curve analysis indicates that evidential value is inadequate or absent if the 33% power test is p < .05 for the full p-curve (fullp33) or both the half p-curve (halfp33) and binomial 33% power test (binomp33) are p < .1.
#'
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
#'
kable(describe(pcurve.out.pub, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

##############

data <- dat[dat$PA.NA. != "Y" & !is.na(dat$g.calc),]

ma.unpub <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$published) & !is.na(data$g.var.calc),], method = "ML", random = ~ 1|study/result, subset = published == "unpublished"), cluster = data$study[!is.na(data$published)])
ma.pub <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$published) & !is.na(data$g.var.calc),], method = "ML", random = ~ 1|study/result, subset = published == "published"), cluster = data$study[!is.na(data$published)])

result.pub.vs.unpub <- data.frame(meta = c("Unpublished","Published"), k = c(ma.unpub$k, ma.pub$k), estimate = round(c(coef(ma.unpub), coef(ma.pub)), 3), stderror = round(c(ma.unpub$se, ma.pub$se), 3),
                       tau = round(c(sqrt(sum(ma.unpub$sigma2)), sqrt(sum(ma.pub$sigma2))) ,3))
kable(result.pub.vs.unpub, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### $I^2$ for unpublished
W <- diag(1/data[data$published == "unpublished",]$g.var.calc)
X <- model.matrix(ma.unpub)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'
I2.unpub <- {{round(100 * sum(ma.unpub$sigma2) / (sum(ma.unpub$sigma2) + (ma.unpub$k-ma.unpub$p)/sum(diag(P))), 2)}}
#'Here, total relative heterogeneity was
{{round(I2.unpub, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * ma.unpub$sigma2 / (sum(ma.unpub$sigma2) + (ma.unpub$k-ma.unpub$p)/sum(diag(P))), 2)}}
#'%, respectively.
#'

#'#### $I^2$ for published studies
W <- diag(1/data[data$published == "published",]$g.var.calc)
X <- model.matrix(ma.pub)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'
I2.pub <- {{round(100 * sum(ma.pub$sigma2) / (sum(ma.pub$sigma2) + (ma.pub$k-ma.pub$p)/sum(diag(P))), 2)}}
#'Here, total relative heterogeneity was
{{round(I2.pub, 2)}}
#'%.

#'Separate estimates of between- and within-cluster heterogeneity were
{{round(100 * ma.pub$sigma2 / (sum(ma.pub$sigma2) + (ma.pub$k-ma.pub$p)/sum(diag(P))), 2)}}
#'%, respectively.
#'

#'#### Wald test p-value
#'Testing the difference in the uncorrected MA estimates between published and unpublished
unpub.vs.pub.z <- with(result.pub.vs.unpub, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))
unpub.vs.pub.p <- 2*pnorm(abs(as.numeric(unpub.vs.pub.z)), lower.tail = F)
unpub.vs.pub.p

#'### Contour enhanced funnel plots
par(mfrow=c(1,2))
funnel(ma.unpub, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor unpublished studies", xlim = c(-1, 1.8), ylim = c(0, 0.6))
funnel(ma.pub, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor published studies", xlim = c(-1, 1.8), ylim = c(0, 0.6))
funnel.pub.unpub <- recordPlot()

#'### Small-study effects correction for unpublished effects

#'#### 3-parameter selection model for unpublished effects
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.unpub <- threePSM.est(data[!is.na(data$g.calc) & data$published == "unpublished",]$g.calc, data[!is.na(data$g.calc) & data$published == "unpublished",]$g.var.calc)
kable(ThreePSM.unpub[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### PET-PEESE for unpublished effects
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.unpub <- with(data[data$published == "unpublished",], pet.peese(g.calc, g.var.calc, study, result))
kable(pp.unpub, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#**PET-PEESE plot for unpublished effects**
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.unpub$value[4] < .05 & ThreePSM.unpub$value[1] > 0)
{plot(data[data$published == "unpublished",]$g.var.calc, data[data$published == "unpublished",]$g.calc, main="PEESE for unpublished studies", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .17) , xaxs="i")} else {plot(sqrt(data[data$published == "unpublished",]$g.var.calc), data[data$published == "unpublished",]$g.calc, main="PET for unpublished studies", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.unpub$value[4] < .05 & ThreePSM.unpub$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#'### Small-study effects correction for published effects

#'#### 3-parameter selection model for published effects
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.published <- threePSM.est(data[data$published == "published",,]$g.calc, data[data$published == "published",]$g.var.calc)
kable(ThreePSM.published[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### PET-PEESE for published effects
#Estimated effect size of an infinitely precise study. If the estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
pp.pub <- with(data[data$published == "published",], pet.peese(g.calc, g.var.calc, study, result))
kable(pp.pub, "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#**PET-PEESE plot for published effects**
#y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.published$value[4] < .05 & ThreePSM.published$value[1] > 0)
{plot(data[data$published == "published",]$g.var.calc, data[data$published == "published",]$g.calc, main="PEESE for published studies", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .28) , xaxs="i")} else {plot(sqrt(data[data$published == "published",]$g.var.calc), data[data$published == "published",]$g.calc, main="PET for published studies", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.published$value[4] < .05 & ThreePSM.published$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")
########################

#'## Randomization
#'
#'Overall effect moderated by the presence of randomization.
data <- dat[dat$PA.NA. != "Y" & !is.na(dat$g.calc),]

ma.obs <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Randomized.Controlled.Study...Y.vs..N.) & !is.na(data$g.var.calc),], method = "ML", random = ~ 1|study/result, subset = Randomized.Controlled.Study...Y.vs..N. == -.5), cluster = data$study[!is.na(data$Randomized.Controlled.Study...Y.vs..N.)])
ma.rand <- robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data[!is.na(data$Randomized.Controlled.Study...Y.vs..N.) & !is.na(data$g.var.calc),], method = "ML", random = ~ 1|study/result, subset = Randomized.Controlled.Study...Y.vs..N. == .5), cluster = data$study[!is.na(data$Randomized.Controlled.Study...Y.vs..N.)])

result <- data.frame(meta = c("Observational","Randomized"), k = c( ma.obs$k, ma.rand$k), estimate = round(c(coef(ma.obs), coef(ma.rand)), 3), stderror = round(c(ma.obs$se, ma.rand$se), 3),
                     tau = round(c(sqrt(sum(ma.obs$sigma2)), sqrt(sum(ma.rand$sigma2))) ,3))
kable(result, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'COMMENT: Observational studies seem to give somewhat higher ESs but not significantly so. Note especially huge heterogeneity of observational studies. That causes the huge SE of MA estimate, primarily producing ns difference btw those two sets.
#'

#'#### Wald test p-value
#'Testing the difference in the uncorrected MA estimates between effects from observational and randomized studies.
z.obs.vs.rand <- with(result, round(c(zval = (estimate[1] - estimate[2])/sqrt(stderror[1]^2 + stderror[2]^2)), 3))
2*pnorm(abs(as.numeric(z.obs.vs.rand)), lower.tail = F)

#'#### $I^2$ for observed
W <- diag(1/data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,]$g.var.calc)
X <- model.matrix(ma.obs)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'
I2.obs <- {{round(100 * sum(ma.obs$sigma2) / (sum(ma.obs$sigma2) + (ma.obs$k-ma.obs$p)/sum(diag(P))), 2)}}
#'Here, total relative heterogeneity was
{{round(I2.obs, 2)}}
#'%.

#'#### $I^2$ for randomized
W <- diag(1/data[data$Randomized.Controlled.Study...Y.vs..N. == .5,]$g.var.calc)
X <- model.matrix(ma.rand)
P <- W - W %*% X %*% solve(t(X) %*% W %*% X) %*% t(X) %*% W
#'
I2.rand <- {{round(100 * sum(ma.rand$sigma2) / (sum(ma.rand$sigma2) + (ma.rand$k-ma.rand$p)/sum(diag(P))), 2)}}
#'Here, total relative heterogeneity was
{{round(I2.rand, 2)}}
#'%.

#'#### Contour enhanced funnel plots
par(mfrow=c(1,2))
funnel(ma.obs, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor observational studies", ylim = c(0, 0.5), xlim = c(-1, 1.8))
funnel(ma.rand, level=c(90, 95, 99), shade=c("white", "gray", "darkgray"), refline=0, pch = 20, yaxis = "sei", xlab = "ES\nFor randomized studies", ylim = c(0, 0.5), xlim = c(-1, 1.8))

#'### Small-study effects correction for observational studies

#'#### 3-parameter selection model for observational studies
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.unpub <- threePSM.est(data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,]$g.calc, data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,]$g.var.calc)
kable(ThreePSM.unpub[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### PET-PEESE for observational studies
#'Estimated effect size of an infinitely precise study. Using 3PSM as the conditional estimator instead of PET. If the PET-PEESE estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
kable(with(data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,], pet.peese(g.calc, g.var.calc, study, result)), "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#**PET-PEESE plot for observational studies**
#'y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.unpub$value[4] < .05 & ThreePSM.unpub$value[1] > 0)
{plot(data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,]$g.var.calc, data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,]$g.calc, main="PEESE for observational studies", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .17) , xaxs="i")} else {plot(sqrt(data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,]$g.var.calc), data[data$Randomized.Controlled.Study...Y.vs..N. == -.5,]$g.calc, main="PET for observational studies", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.unpub$value[4] < .05 & ThreePSM.unpub$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

#'### Small-study effects correction for randomized studies

#'#### 3-parameter selection model for randomized studies
#'Bias-corrected estimate, note especially the CI (conf.low, conf.high).
ThreePSM.published <- threePSM.est(data[data$Randomized.Controlled.Study...Y.vs..N. == .5,,]$g.calc, data[data$Randomized.Controlled.Study...Y.vs..N. == .5,]$g.var.calc)
kable(ThreePSM.published[c(1, 4, 5, 6),], "html", digits = 3) %>%   kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'#### PET-PEESE for randomized studies
#Estimated effect size of an infinitely precise study. If the estimate is negative, the effect can be regarded 0. pval = p-value testing H0 that the effect is zero. cil.lb and ci.ub are upper and lower bound of the CI.
kable(with(data[data$Randomized.Controlled.Study...Y.vs..N. == .5,], pet.peese(g.calc, g.var.calc, study, result)), "html", digits = 3, col.names = "") %>% kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#**PET-PEESE plot for randomized studies**
#y-axis intercept represents the estimated bias-corrected ES.
if(ThreePSM.published$value[4] < .05 & ThreePSM.published$value[1] > 0)
{plot(data[data$Randomized.Controlled.Study...Y.vs..N. == .5,]$g.var.calc, data[data$Randomized.Controlled.Study...Y.vs..N. == .5,]$g.calc, main="PEESE for randomized studies", xlab = "Variance", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xlim = c(0, .28) , xaxs="i")} else {plot(sqrt(data[data$Randomized.Controlled.Study...Y.vs..N. == .5,]$g.var.calc), data[data$Randomized.Controlled.Study...Y.vs..N. == .5,]$g.calc, main="PET for randomized studies", xlab = "Standard error", ylab = "Effect size", pch = 19, cex.main = 1.3, cex = .6, xaxs="i")}
abline((if(ThreePSM.published$value[4] < .05 & ThreePSM.published$value[1] > 0) {peese} else {pet}), lwd=3, lty = 2, col = "red")

###################
#'## Year of Publication
#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
year.pub <- summary(lmer(scale(sqrt(g.var.calc)) ~ scale(H5.Index.GS.Journal.March.2016) + scale(Publication.Year) + (1|study), data = data, REML = T))$coefficients
kable(year.pub, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")
#'Comment: all the variables were centered for easier interpretation of model coefficients. See the negative beta for Publication Year. The higher the publication year, the lower the variance (better precision), controlling for H5.
#'
#'Thus, practices regarding the precision of studies (mainly due to N) seem to have improved throughout last years.
#'

#'#### Scatterplot year <-> precision
#'
#'Size of the points indicate the H5 index (the bigger the higher) of the journal that the ES is published in.
s <- ggplot(data, aes(Publication.Year, sqrt(g.var.calc)))
s + geom_point(aes(size = H5.Index.GS.Journal.March.2016), alpha=.70, colour="#80afce") +
  geom_smooth(method=lm) +
  scale_x_continuous(breaks=2005:2017) +
  xlab("Year of publication") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")

#'## Citations
#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
citations <- summary(lmer(scale(sqrt(g.var.calc)) ~ scale(Publication.Year) + scale(H5.Index.GS.Journal.March.2016) + scale(Citations.March.1.2016..GS.) + (1|study), data = data, REML = T))$coefficients
kable(citations, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")


#'#### Scatterplot precision <-> citations
#'
#'The relationship between precision (sqrt of variance) and number of citations.

s <- ggplot(data, aes(log(Citations.March.1.2016..GS. + 1), sqrt(g.var.calc)))
s + geom_point(alpha=.70, colour="#80afce") +
  geom_smooth(method=lm) +
  xlim(0, 8) +
  xlab("Citations") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")

#'#### Scatterplot precision <-> journal H5
#'

#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study. Using square root of variance to make the distribution normal.
H5.index <- summary(lmer(scale(sqrt(g.var.calc)) ~ scale(H5.Index.GS.Journal.March.2016) + (1|study), data = data, REML = T))$coefficients
kable(H5.index, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'The relationship between precision (sqrt of variance) and H5 index of the journal.

s <- ggplot(data, aes(H5.Index.GS.Journal.March.2016, sqrt(g.var.calc)))
s + geom_point(alpha=.70, colour="#80afce") +
  geom_smooth(method=lm) +
  xlim(20,165) +
  xlab("H5 index") +
  ylab("Precision") +
  theme_bw() +
  theme(legend.position = "none")

#'## Decline effect
#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study.
MA.estimate <- as.numeric(robust.rma.mv(rma.mv(yi = g.calc, V = g.var.calc, data = data, method = "REML", random = ~ 1|study/result), cluster = data$study)[1])
decline.eff <- summary(lmer(scale(abs(MA.estimate - g.calc)) ~ scale(sqrt(g.var.calc)) + scale(Publication.Year) + (1|study), data = data, REML = T))$coefficients
kable(decline.eff, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'## Citation bias
#'Linear mixed-effects model. Taking into effect clustering of ESs due to originating from the same study.
cit.bias <- summary(lmer(scale(abs(MA.estimate - g.calc)) ~ scale(sqrt(g.var.calc)) + scale(Citations.March.1.2016..GS.) + (1|study), data = data, REML = T))$coefficients
kable(cit.bias, "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left")

#'## P-curve for interaction effects
#+include = FALSE
for(i in 1:nsim){
  pcurve[i,] <- pcurve_app(na.omit(gsub("^.*?: ",": ", x = dat[dat$PA.NA. != "Y" & dat$Moderated.Effect. == "Y",]$label[!duplicated.random(dat[dat$PA.NA. != "Y" & dat$Moderated.Effect. == "Y",]$study)], replacement = "")))
}
pcurve.out.mod <- data.frame(pcurve)
colnames(pcurve.out.mod) <- c("ksig", "khalf", "fullz", "fullp", "fullz33", "fullp33", "halfz", "halfp", "halfz33", "halfp33", "binomp", "binomp33", "power.ci.lb", "power.est", "power.ci.up")
#'ksig = average number of effects associated with p < .05; khalf = average number of effects associated with p < .025;...z = average z-values; power.est = average estimated statistical power of the studies (with lower bound and upper bound)
#'
kable(describe(pcurve.out.mod, skew = FALSE, ranges = FALSE), "html", digits = 3) %>%
  kable_styling(bootstrap_options = "striped", full_width = F, font_size = 12, position = "left") %>%
  row_spec(c(4, 6, 8, 10, 14), bold = T, color = "white", background = "#bbd4f0")

#+include = FALSE
hist(dat[dat$Moderated.Effect. == "Y",]$p, xlim = c(0, .0481), breaks = 200, main = "p-values for moderated effects", xlab = "p-value")

