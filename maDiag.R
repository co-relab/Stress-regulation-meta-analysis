# Outlier diagnostics -----------------------------------------------------

#'# For mindfulness studies

#+eval = FALSE
# Initial outlier diagnostics
# Univariate MA
ma.uni <- dataMind %$% rma(yi = yi, vi = vi, method = "REML", slab = result)

#+eval = FALSE
# MA diagnostics
baujat(ma.uni)

#+eval = FALSE
#fit FE model to all possible subsets
gosh.plot <- gosh(ma.uni, progbar = TRUE, subsets = 1000, parallel = "multicore")
plot(gosh.plot, out = 45, breaks = 50) # Testing the influence of single outliers

#+eval = FALSE
# Influence diagnostics
inf <- influence(ma.uni, progbar = T)

#+eval = FALSE
### Plot the influence diagnostics
plot(inf)

#'# For biofeedback studies

#+eval = FALSE
# Initial outlier diagnostics
# Univariate MA
ma.uni <- dataBio %$% rma(yi = yi, vi = vi, method = "REML", slab = result)

#+eval = FALSE
# MA diagnostics
baujat(ma.uni)

#+eval = FALSE
#fit FE model to all possible subsets
gosh.plot <- gosh(ma.uni, progbar = TRUE, subsets = 1000, parallel = "multicore")
plot(gosh.plot, out = 12, breaks = 50) # Testing the influence of single outliers

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

# Missing data ------------------------------------------------------------

# For mindfulness studies
dat %>% filter(strategy == 1) %$% table(is.na(.$yi))

# For biofeedback studies
dat %>% filter(strategy == 2) %$% table(is.na(.$yi))

#'### Percentage of missing data overall
# dat %>% filter(strategy == 2) %$% paste(round(sum(is.na(.))/prod(dim(.))*100, 3), "%", sep = "") # insert collumn numbers
dat %>% filter(strategy == 2) %>% missmap(rank.order = TRUE, margins = c(10, 0), legend = F)
