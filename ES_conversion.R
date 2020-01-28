# Load required libraries
if (!require(metafor)) {
  install.packages('metafor')
}
if (!require(dplyr)) {
  install.packages('dplyr')
}
if (!require(esc)) {
  install.packages('esc')
}
if (!require(magrittr)) {
  install.packages('magrittr')
}

# Read in the data
dat <- read.csv("MA data oct 2018 IR.csv", sep = ";")

# Some data wrangling to get the right type of data (formatting of the raw dataset in Excel introduces a lot of junk otherwise)
dat$F <- as.numeric(as.character(dat$F))
dat$t <- as.numeric(as.character(dat$t))
dat$df1 <- as.numeric(as.character(dat$df1))
dat$df2 <- as.numeric(as.character(dat$df2))
dat$N.warm <- as.numeric(as.character(dat$Nwarm))
dat$N.cold <- as.numeric(as.character(dat$Ncold))
dat$N.control <- as.numeric(as.character(dat$Ncontrol))
dat$Design <- as.factor(as.character(dat$Design))
dat$published <- as.factor(as.character(dat$published))
dat$Randomized.Controlled.Study...Y.vs..N. <- as.factor(as.character(dat$Randomized.Controlled.Study...Y.vs..N.))
dat$Publication.Year <- as.numeric(as.character(dat$Publication.Year))
dat$n.rep <- as.numeric(as.character(dat$Reported.Overall.N))
dat$d.reported <- as.numeric(as.character(dat$d.reported))
dat$p.reported <- as.numeric(as.character(gsub("[^0-9.]", "", dat$p.reported)))
dat$PA.NA. <- as.factor(as.character(dat$PA.NA.))
dat$Latitude.University..proxy.for.climate. <- as.numeric(as.character(dat$Latitude.University..proxy.for.climate.))
dat$H5.Index.GS.Journal.March.2016 <- as.numeric(as.character(dat$H5.Index.GS.Journal.March.2016))

# Transform moderators into factors and specify contrast coding
dat$published <- relevel(factor(ifelse(dat$published == "Y", "published", "unpublished")), ref= "unpublished")
dat$Randomized.Controlled.Study...Y.vs..N. <- relevel(factor(ifelse(dat$Randomized.Controlled.Study...Y.vs..N. == "Y", .5, -.5)), ref= "-0.5")

# Which designs are present?
table(dat$Design, useNA="ifany")

# Initialize new variables
dat$final.design <- "NA"
dat$bias.test <- "NA"
dat$N <- NA
dat$t.from.r <- NA
dat$t.from.beta <- NA
dat$r.var <- NA
dat$beta.var <- NA
dat$d.calc <- NA
dat$g.calc <- NA
dat$g.var.calc <- NA
dat$p <- NA
dat$label <- NA

# Create result and study ID
dat$result <- 1:nrow(dat)
dat$study <- paste(dat$paperID, "/", dat$Study.Indicator, sep = "")

# For 2-cell designs, establish cell sizes
dat$cell1n <- dat$N.warm
dat$cell2n <- ifelse(!is.na(dat$N.cold), dat$N.cold, dat$N.control)

# Compute gender ratio (% of female) ??added code
dat$N.male <- as.numeric(as.character(dat$N.male))
dat$N.female <- as.numeric(as.character(dat$N.female))
dat$gender.ratio <- dat$N.female/(dat$N.female + dat$N.male)

####
# F-Test between with df1 == 1 ---------------------------------------------------------------------
####

# Specify the design, compute N and p
dat[dat$Use.for.Meta == "Yes" & !is.na(dat$F) & !is.na(dat$df1) & !is.na(dat$df2) & dat$df1 == 1, "final.design"] <- "F1"
dat[dat$final.design == "F1", ]$N <- dat[dat$final.design == "F1", ]$df2 + 2
dat[dat$final.design == "F1", ]$p <- with(dat[dat$final.design == "F1", ], 1-pf(F, df1, df2))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat$UseCellN <- ifelse((dat$cell1n + dat$cell2n) >= (dat$N - 2) & (dat$cell1n + dat$cell2n) <= (dat$N + 2), 1, 0)
dat$UseCellN[is.na(dat$UseCellN)] <- 0

# Compute ES and var based on total N
dat[dat$final.design == "F1", ]$g.calc <- esc_f(f = dat[dat$final.design == "F1", ]$F, totaln = dat[dat$final.design == "F1", ]$df2 + 2, es.type = "g")$es
dat[dat$final.design == "F1", ]$g.var.calc <- esc_f(f = dat[dat$final.design == "F1", ]$F, totaln = dat[dat$final.design == "F1", ]$df2 + 2, es.type = "g")$var

# Compute ES and var based on n1 and n2 if available
dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$g.calc <- esc_f(f = dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$F, grp1n = dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$cell1n, grp2n = dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$cell2n, es.type = "g")$es
dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$g.var.calc <- esc_f(f = dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$F, grp1n = dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$cell1n, grp2n = dat[dat$final.design == "F1" & dat$UseCellN == 1, ]$cell2n, es.type = "g")$var

# Show the converted ESs
dat %>% filter(final.design == "F1") %>% select(g.calc, g.var.calc, Design, d.reported, df2, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & !is.na(dat$F) & !is.na(dat$df1) & !is.na(dat$df2), "bias.test"] <- "F"
dat$label[dat$bias.test == "F"] <- paste(dat[dat$bias.test == "F",]$paperID, "/", dat[dat$bias.test == "F",]$Study.Indicator, "/", dat[dat$bias.test == "F",]$Variable.Indicator, ": ",
                                                    "F(", dat[dat$bias.test == "F",]$df1, ",",dat[dat$bias.test == "F",]$df2, ")=", dat[dat$bias.test == "F",]$F, sep = "")

####
# t-tests between ---------------------------------------------------------
####

# Specify the design, compute N and p
dat[dat$Use.for.Meta == "Yes" & dat$Design == "Between" & !is.na(dat$t) & !is.na(dat$df2), "final.design"] <- "between.t"
dat[dat$final.design == "between.t", ]$N <- dat[dat$final.design == "between.t", ]$df2 + 2
dat[dat$final.design == "between.t", ]$p <- with(dat[dat$final.design == "between.t", ], 2*pt(abs(t), df2, lower.tail = FALSE))

#Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat$UseCellN <- ifelse((dat$cell1n + dat$cell2n) >= (dat$N - 2) & (dat$cell1n + dat$cell2n) <= (dat$N + 2), 1, 0)
dat$UseCellN[is.na(dat$UseCellN)] <- 0

#Compute ES and var based on total N
dat[dat$final.design == "between.t", ]$g.calc <- esc_t(t = abs(dat[dat$final.design == "between.t", ]$t), totaln = dat[dat$final.design == "between.t", ]$df2 + 2, es.type = "g")$es
dat[dat$final.design == "between.t", ]$g.var.calc <- esc_t(t = dat[dat$final.design == "between.t", ]$t, totaln = dat[dat$final.design == "between.t", ]$df2 + 2, es.type = "g")$var

#Compute ES and var based on n1 and n2 if available
dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$g.calc <- esc_t(t = abs(dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$t), grp1n = dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$cell1n, grp2n = dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$cell2n, es.type = "g")$es
dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$g.var.calc <- esc_t(t = dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$t, grp1n = dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$cell1n, grp2n = dat[dat$final.design == "between.t" & dat$UseCellN == 1, ]$cell2n, es.type = "g")$var

# Show the converted ESs
dat %>% filter(final.design == "between.t") %>% select(g.calc, g.var.calc, Design, d.reported, df2, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & dat$Design == "Between" & !is.na(dat$t) & !is.na(dat$df2), "bias.test"] <- "between.t"
dat$label[dat$bias.test == "between.t"] <- paste(dat[dat$bias.test == "between.t",]$paperID, "/", dat[dat$bias.test == "between.t",]$Study.Indicator, "/", dat[dat$bias.test == "between.t",]$Variable.Indicator, ": ",
                                                    "t(", dat[dat$bias.test == "between.t",]$df2, ")=", dat[dat$bias.test == "between.t",]$t, sep = "")

####
# Correlation -------------------------------------------------------------
####

# Specify the design
dat[dat$Use.for.Meta == "Yes" & !is.na(dat$r) & !is.na(dat$df2), "final.design"] <- "correlation"

# Compute ES, var, N, and p
dat[dat$final.design == "correlation", ] %<>% mutate(
  t.from.r = abs(r)*sqrt(df2 / (1 - r^2)),
  g.calc = ((2*abs(r))/sqrt(1-r^2))*(1 - (3/(4*df2 - 1))),
  N = df2 + 2,
  r.var = escalc(measure = "COR", ri = r, ni = df2 + 2, data = dat[dat$final.design == "correlation", ])$vi,
  g.var.calc = (1 - (3/(4*df2 - 1))) * (4 * r.var/(1 - r^2)^3),
  p = 2*pt(abs(t.from.r), df2, lower.tail=FALSE)
)

# Show the converted ESs
dat %>% filter(final.design == "correlation") %>% select(g.calc, g.var.calc, Design, d.reported, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & !is.na(dat$r) & !is.na(dat$df2), "bias.test"] <- "correlation"
dat$label[dat$bias.test == "correlation"] <- paste(dat[dat$bias.test == "correlation",]$paperID, "/", dat[dat$bias.test == "correlation",]$Study.Indicator, "/", dat[dat$bias.test == "correlation",]$Variable.Indicator, ": ",
                                                    "r(", dat[dat$bias.test == "correlation",]$df2, ")=", dat[dat$bias.test == "correlation",]$r, sep = "")

####
# Within-subjects design, ES based on t-distribution ----------------------
####

# Identify within-subjects design reporting F and compute t
dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep), "final.design"] <- "within.t"
dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep),]$t <- sqrt(dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep),]$F)

# Compute ES using t
dat[dat$Use.for.Meta == "Yes" & dat$Design == "Within" & !is.na(dat$t) & !is.na(dat$n.rep), "final.design"] <- "within.t"
dat[dat$final.design == "within.t", ] %<>% mutate(
  d.calc = abs(t)*sqrt((2 * (1 - corr)) / n.rep),
  g.calc = (1 - (3/(4*n.rep - 3))) * d.calc,
  g.var.calc = (1 - (3/(4*n.rep - 3)))^2 * ((1 / n.rep) + ((d.calc^2) / (2 * n.rep))) * 2 * (1 - corr),
  N = n.rep,
  p = 2*pt(abs(t), n.rep - 1, lower.tail = FALSE)
)

# Show the converted ESs
dat %>% filter(final.design == "within.t") %>% select(g.calc, g.var.calc, Design, d.reported, N, df2)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & dat$Design == "Within" & !is.na(dat$t) & !is.na(dat$df2), "bias.test"] <- "within.t"
dat$label[dat$bias.test == "within.t"] <- paste(dat[dat$bias.test == "within.t",]$paperID, "/", dat[dat$bias.test == "within.t",]$Study.Indicator, "/", dat[dat$bias.test == "within.t",]$Variable.Indicator, ": ",
                                                    "t(", dat[dat$bias.test == "within.t",]$df2, ")=", dat[dat$bias.test == "within.t",]$t, sep = "")

# PaperID 71 used mixed-effects models, couldn't convert, so using the reported d (converted to g)
dat$g.calc[266:270] <- (1 - (3/(4*dat$n.rep[266:270] - 3))) * dat$d.reported[266:270]
dat$g.var.calc[266:270] = (1 - (3/(4*dat$n.rep[266:270] - 3))) * ((dat$n.rep[266:270])/(dat$n.rep[266:270]/2 * dat$n.rep[266:270]/2) + (dat$d.reported[266:270]^2)/(2 * (dat$n.rep[266:270])))

####
# Betas in between-subjects designs ---------------------------------------
####

#Betas in between-subjects designs (in ML, beta considered as covariate/confounding adjusted r, then using r to d conversion)
dat[dat$Use.for.Meta == "Yes" & (dat$Design == "Between" | dat$Design == "Continuous") & !is.na(dat$beta) & !is.na(dat$df2), "final.design"] <- "Beta.between"

dat[dat$final.design == "Beta.between", ] %<>% mutate(
  t.from.beta = abs(beta)*sqrt(df2 / (1 - beta^2)),
  g.calc = ((2*abs(beta))/sqrt(1-beta^2))*(1 - (3/(4*df2 - 1))),
  N = df2 + 2,
  beta.var = escalc(measure = "COR", ri = beta, ni = df2 + 2, data = dat[dat$final.design == "Beta.between", ])$vi,
  g.var.calc = (1 - (3/(4*df2 - 1))) * (4 * beta.var)/((1 - beta^2)^3),
  p = 2*pt(abs(t.from.beta), df2, lower.tail=FALSE)
)

# Show the converted ESs
dat %>% filter(final.design == "Beta.between") %>% select(g.calc, beta, g.var.calc, final.design, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & (dat$Design == "Between" | dat$Design == "Continuous") & !is.na(dat$beta) & !is.na(dat$p.reported), "bias.test"] <- "Beta.between"
dat$label[dat$bias.test == "Beta.between"] <- paste(dat[dat$bias.test == "Beta.between",]$paperID, "/", dat[dat$bias.test == "Beta.between",]$Study.Indicator, "/", dat[dat$bias.test == "Beta.between",]$Variable.Indicator, ": ",
                                                      "Z=", qnorm(1-dat[dat$bias.test == "Beta.between",]$p.reported/2), sep = "")

####
# chi^2 -------------------------------------------------------------------
####

# Specify the design, compute ES, var, N, and p
dat[dat$Use.for.Meta == "Yes" & dat$Design == "Between" & !is.na(dat$Chisq) & !is.na(dat$n.rep), "final.design"] <- "between.chisq"

dat[dat$final.design == "between.chisq", ]$g.calc <- esc_chisq(chisq = dat[dat$final.design == "between.chisq", ]$Chisq, totaln = dat[dat$final.design == "between.chisq", ]$n.rep, es.type = "g")$es
dat[dat$final.design == "between.chisq", ]$g.var.calc <- esc_chisq(chisq = dat[dat$final.design == "between.chisq", ]$Chisq, totaln = dat[dat$final.design == "between.chisq", ]$n.rep, es.type = "g")$var
dat[dat$final.design == "between.chisq", ]$N <- dat[dat$final.design == "between.chisq", ]$n.rep
dat[dat$final.design == "between.chisq", ]$p <- with(dat[dat$final.design == "between.chisq", ], 1 - pchisq(Chisq, 1))

# Show the converted ESs
dat %>% filter(final.design == "between.chisq") %>% select(g.calc, g.var.calc, Chisq, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & dat$Design == "Between" & !is.na(dat$Chisq) & !is.na(dat$n.rep), "bias.test"] <- "between.chisq"
dat$label[dat$bias.test == "between.chisq"] <- paste(dat[dat$bias.test == "between.chisq",]$paperID, "/", dat[dat$bias.test == "between.chisq",]$Study.Indicator, "/", dat[dat$bias.test == "between.chisq",]$Variable.Indicator, ": ",
                                                        "chi2(", dat[dat$bias.test == "between.chisq",]$df1, ")=", dat[dat$bias.test == "between.chisq",]$Chisq, sep = "")

####
# t for B in continuous designs -------------------------------------------
####

# Specify the design, compute ES, var, N, and p
dat[dat$Use.for.Meta == "Yes" & (dat$Design == "Continuous" |  dat$Design == "Correlation") & !is.na(dat$B) & !is.na(dat$t) & !is.na(dat$df2) & is.na(dat$beta), "final.design"] <- "t.from.B"

dat[dat$final.design == "t.from.B", ]$g.calc <- esc_t(t = abs(dat[dat$final.design == "t.from.B", ]$t), totaln = dat[dat$final.design == "t.from.B", ]$df2 + 2, es.type = "g")$es
dat[dat$final.design == "t.from.B", ]$g.var.calc <- esc_t(t = dat[dat$final.design == "t.from.B", ]$t, totaln = dat[dat$final.design == "t.from.B", ]$df2 + 2, es.type = "g")$var
dat[dat$final.design == "t.from.B", ]$N <- dat[dat$final.design == "t.from.B", ]$df2 + 2
dat[dat$final.design == "t.from.B", ]$p <- with(dat[dat$final.design == "t.from.B", ], 2*pt(abs(t), N - 1, lower.tail = FALSE))

# Show the converted ESs
dat %>% filter(final.design == "t.from.B") %>% select(g.calc, g.var.calc, B, t, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$Use.for.Bias.Test == "Yes" & (dat$Design == "Continuous" |  dat$Design == "Correlation") & !is.na(dat$B) & !is.na(dat$t) & !is.na(dat$df2) & is.na(dat$beta), "bias.test"] <- "t.from.B"
dat$label[dat$bias.test == "t.from.B"] <- paste(dat[dat$bias.test == "t.from.B",]$paperID, "/", dat[dat$bias.test == "t.from.B",]$Study.Indicator, "/", dat[dat$bias.test == "t.from.B",]$Variable.Indicator, ": ",
                                                    "t(", dat[dat$bias.test == "t.from.B",]$df2, ")=", dat[dat$bias.test == "t.from.B",]$t, sep = "")

# Multiply the ES by -1 if not in the predicted direction
dat$g.calc <- ifelse(test = (dat$Predicted.Direction == -1 & dat$g.calc > 0), yes = (dat$g.calc*-1), no = dat$g.calc)


####
# Create a results label for the rest of bias=yes and meta=no, based on p-values, converted to z-score (p-curve app works with z but not with p).
dat[dat$Use.for.Bias.Test == "Yes" & is.na(dat$label) & !is.na(dat$paperID), "bias.test"] <- "z.from.p"
dat$p <- ifelse(is.na(dat$p), dat$p.reported, dat$p)
dat$label[dat$bias.test == "z.from.p"] <- paste(dat[dat$bias.test == "z.from.p",]$paperID, "/", dat[dat$bias.test == "z.from.p",]$Study.Indicator, "/", dat[dat$bias.test == "z.from.p",]$Variable.Indicator, ": ",
                                              "Z=", qnorm(1-(dat[dat$bias.test == "z.from.p",]$p)/2), sep = "")

# Remove unused variables
delvars <- names(dat) %in% c(
  "t.from.r", "t.from.beta", "r.var", "beta.var", "d.calc")
dat <- dat[!delvars]
rm(delvars)

