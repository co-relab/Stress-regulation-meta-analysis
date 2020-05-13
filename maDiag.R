# These procedures will be used to do the diagnostics

#+eval = FALSE
# Initial outlier diagnostics
# Univariate MA
ma.uni <- rma(yi = yi, vi = vi, data = dat, method = "REML", slab = result)

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
#paste(round(sum(is.na(dat[,1:34]))/prod(dim(dat[,1:34]))*100, 3), "%", sep = "") # insert collumn numbers
#missmap(dat, rank.order = TRUE, margins = c(5, 0), legend = F)    # insert collumn numbers


