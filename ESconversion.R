# Need paperID studyID nMale nfemale useMeta useBias predictedDirection journalH5
# Need to simulate: mean, sd, n, items


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


##### under construction
dat$Design <- as.factor(as.character(dat$Design))
##### continue here



dat$published <- as.factor(as.character(dat$published))
dat$design <- as.factor(as.character(dat$design))
dat$Publication.Year <- as.numeric(as.character(dat$Publication.Year))
dat$n.rep <- as.numeric(as.character(dat$Reported.Overall.N))
dat$d.reported <- as.numeric(as.character(dat$d.reported))
dat$p.reported <- as.numeric(as.character(gsub("[^0-9.]", "", dat$p.reported)))
dat$journalH5 <- as.numeric(as.character(dat$journalH5))

# Transform moderators into factors and specify contrast coding
dat$published <- relevel(factor(ifelse(dat$published == "1", "published", "unpublished")), ref= "unpublished")
dat$design <- relevel(factor(ifelse(dat$design == 1, -.5, 5)), ref= "-0.5")

# Which designs are present?
table(dat$Design, useNA="ifany")

# Initialize new variables
dat$finalDesign <- "NA"
dat$biasTest <- "NA"
dat$N <- NA
dat$t.from.r <- NA
dat$t.from.beta <- NA
dat$r.var <- NA
dat$beta.var <- NA
dat$d.calc <- NA
dat$gCalc <- NA
dat$gVarCalc <- NA
dat$p <- NA
dat$label <- NA

# Create result and study ID
dat$result <- 1:nrow(dat)
dat$study <- paste(dat$paperID, "/", dat$studyID, sep = "")

# For 2-cell designs, establish cell sizes
dat$cell1n <- dat$N.warm
dat$cell2n <- ifelse(!is.na(dat$N.cold), dat$N.cold, dat$N.control)

# Compute gender ratio (% of female) ??added code
dat$nMale <- as.numeric(as.character(dat$nMale))
dat$nFemale <- as.numeric(as.character(dat$nFemale))
dat$genderRatio <- dat$nFemale/(dat$nFemale + dat$nMale)

####
# F-Test between with df1 == 1 ---------------------------------------------------------------------
####

# Specify the design, compute N and p
dat[dat$useMeta == 1 & !is.na(dat$F) & !is.na(dat$df1) & !is.na(dat$df2) & dat$df1 == 1, "finalDesign"] <- "F1"
dat[dat$finalDesign == "F1", ]$N <- dat[dat$finalDesign == "F1", ]$df2 + 2
dat[dat$finalDesign == "F1", ]$p <- with(dat[dat$finalDesign == "F1", ], 1-pf(F, df1, df2))

# Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat$useCellN <- ifelse((dat$cell1n + dat$cell2n) >= (dat$N - 2) & (dat$cell1n + dat$cell2n) <= (dat$N + 2), 1, 0)
dat$useCellN[is.na(dat$useCellN)] <- 0

# Compute ES and var based on total N
dat[dat$finalDesign == "F1", ]$gCalc <- esc_f(f = dat[dat$finalDesign == "F1", ]$F, totaln = dat[dat$finalDesign == "F1", ]$df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "F1", ]$gVarCalc <- esc_f(f = dat[dat$finalDesign == "F1", ]$F, totaln = dat[dat$finalDesign == "F1", ]$df2 + 2, es.type = "g")$var

# Compute ES and var based on n1 and n2 if available
dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$gCalc <- esc_f(f = dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$F, grp1n = dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$cell1n, grp2n = dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$cell2n, es.type = "g")$es
dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$gVarCalc <- esc_f(f = dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$F, grp1n = dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$cell1n, grp2n = dat[dat$finalDesign == "F1" & dat$useCellN == 1, ]$cell2n, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "F1") %>% select(gCalc, gVarCalc, Design, d.reported, df2, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$useBias == 1 & !is.na(dat$F) & !is.na(dat$df1) & !is.na(dat$df2), "biasTest"] <- "F"
dat$label[dat$biasTest == "F"] <- paste(dat[dat$biasTest == "F",]$paperID, "/", dat[dat$biasTest == "F",]$studyID, "/", dat[dat$biasTest == "F",]$Variable.Indicator, ": ",
                                                    "F(", dat[dat$biasTest == "F",]$df1, ",",dat[dat$biasTest == "F",]$df2, ")=", dat[dat$biasTest == "F",]$F, sep = "")

####
# t-tests between ---------------------------------------------------------
####

# Specify the design, compute N and p
dat[dat$useMeta == 1 & dat$Design == "Between" & !is.na(dat$t) & !is.na(dat$df2), "finalDesign"] <- "between.t"
dat[dat$finalDesign == "between.t", ]$N <- dat[dat$finalDesign == "between.t", ]$df2 + 2
dat[dat$finalDesign == "between.t", ]$p <- with(dat[dat$finalDesign == "between.t", ], 2*pt(abs(t), df2, lower.tail = FALSE))

#Decide whether n1+n2 approximately corresponds to the reported df (in order for n1 and n2 be used in equations).
dat$useCellN <- ifelse((dat$cell1n + dat$cell2n) >= (dat$N - 2) & (dat$cell1n + dat$cell2n) <= (dat$N + 2), 1, 0)
dat$useCellN[is.na(dat$useCellN)] <- 0

#Compute ES and var based on total N
dat[dat$finalDesign == "between.t", ]$gCalc <- esc_t(t = abs(dat[dat$finalDesign == "between.t", ]$t), totaln = dat[dat$finalDesign == "between.t", ]$df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "between.t", ]$gVarCalc <- esc_t(t = dat[dat$finalDesign == "between.t", ]$t, totaln = dat[dat$finalDesign == "between.t", ]$df2 + 2, es.type = "g")$var

#Compute ES and var based on n1 and n2 if available
dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$gCalc <- esc_t(t = abs(dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$t), grp1n = dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$cell1n, grp2n = dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$cell2n, es.type = "g")$es
dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$gVarCalc <- esc_t(t = dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$t, grp1n = dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$cell1n, grp2n = dat[dat$finalDesign == "between.t" & dat$useCellN == 1, ]$cell2n, es.type = "g")$var

# Show the converted ESs
dat %>% filter(finalDesign == "between.t") %>% select(gCalc, gVarCalc, Design, d.reported, df2, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$useBias == 1 & dat$Design == "Between" & !is.na(dat$t) & !is.na(dat$df2), "biasTest"] <- "between.t"
dat$label[dat$biasTest == "between.t"] <- paste(dat[dat$biasTest == "between.t",]$paperID, "/", dat[dat$biasTest == "between.t",]$studyID, "/", dat[dat$biasTest == "between.t",]$Variable.Indicator, ": ",
                                                    "t(", dat[dat$biasTest == "between.t",]$df2, ")=", dat[dat$biasTest == "between.t",]$t, sep = "")

####
# Correlation -------------------------------------------------------------
####

# Specify the design
dat[dat$useMeta == 1 & !is.na(dat$r) & !is.na(dat$df2), "finalDesign"] <- "correlation"

# Compute ES, var, N, and p
dat[dat$finalDesign == "correlation", ] %<>% mutate(
  t.from.r = abs(r)*sqrt(df2 / (1 - r^2)),
  gCalc = ((2*abs(r))/sqrt(1-r^2))*(1 - (3/(4*df2 - 1))),
  N = df2 + 2,
  r.var = escalc(measure = "COR", ri = r, ni = df2 + 2, data = dat[dat$finalDesign == "correlation", ])$vi,
  gVarCalc = (1 - (3/(4*df2 - 1))) * (4 * r.var/(1 - r^2)^3),
  p = 2*pt(abs(t.from.r), df2, lower.tail=FALSE)
)

# Show the converted ESs
dat %>% filter(finalDesign == "correlation") %>% select(gCalc, gVarCalc, Design, d.reported, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$useBias == 1 & !is.na(dat$r) & !is.na(dat$df2), "biasTest"] <- "correlation"
dat$label[dat$biasTest == "correlation"] <- paste(dat[dat$biasTest == "correlation",]$paperID, "/", dat[dat$biasTest == "correlation",]$studyID, "/", dat[dat$biasTest == "correlation",]$Variable.Indicator, ": ",
                                                    "r(", dat[dat$biasTest == "correlation",]$df2, ")=", dat[dat$biasTest == "correlation",]$r, sep = "")

####
# Within-subjects design, ES based on t-distribution ----------------------
####

# Identify within-subjects design reporting F and compute t
dat[dat$useMeta == 1 & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep), "finalDesign"] <- "within.t"
dat[dat$useMeta == 1 & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep),]$t <- sqrt(dat[dat$useMeta == 1 & dat$Design == "Within" & !is.na(dat$F) & !is.na(dat$n.rep),]$F)

# Compute ES using t
dat[dat$useMeta == 1 & dat$Design == "Within" & !is.na(dat$t) & !is.na(dat$n.rep), "finalDesign"] <- "within.t"
dat[dat$finalDesign == "within.t", ] %<>% mutate(
  d.calc = abs(t)*sqrt((2 * (1 - corr)) / n.rep),
  gCalc = (1 - (3/(4*n.rep - 3))) * d.calc,
  gVarCalc = (1 - (3/(4*n.rep - 3)))^2 * ((1 / n.rep) + ((d.calc^2) / (2 * n.rep))) * 2 * (1 - corr),
  N = n.rep,
  p = 2*pt(abs(t), n.rep - 1, lower.tail = FALSE)
)

# Show the converted ESs
dat %>% filter(finalDesign == "within.t") %>% select(gCalc, gVarCalc, Design, d.reported, N, df2)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$useBias == 1 & dat$Design == "Within" & !is.na(dat$t) & !is.na(dat$df2), "biasTest"] <- "within.t"
dat$label[dat$biasTest == "within.t"] <- paste(dat[dat$biasTest == "within.t",]$paperID, "/", dat[dat$biasTest == "within.t",]$studyID, "/", dat[dat$biasTest == "within.t",]$Variable.Indicator, ": ",
                                                    "t(", dat[dat$biasTest == "within.t",]$df2, ")=", dat[dat$biasTest == "within.t",]$t, sep = "")

# PaperID 71 used mixed-effects models, couldn't convert, so using the reported d (converted to g)
dat$gCalc[266:270] <- (1 - (3/(4*dat$n.rep[266:270] - 3))) * dat$d.reported[266:270]
dat$gVarCalc[266:270] = (1 - (3/(4*dat$n.rep[266:270] - 3))) * ((dat$n.rep[266:270])/(dat$n.rep[266:270]/2 * dat$n.rep[266:270]/2) + (dat$d.reported[266:270]^2)/(2 * (dat$n.rep[266:270])))

####
# Betas in between-subjects designs ---------------------------------------
####

#Betas in between-subjects designs (in ML, beta considered as covariate/confounding adjusted r, then using r to d conversion)
dat[dat$useMeta == 1 & (dat$Design == "Between" | dat$Design == "Continuous") & !is.na(dat$beta) & !is.na(dat$df2), "finalDesign"] <- "Beta.between"

dat[dat$finalDesign == "Beta.between", ] %<>% mutate(
  t.from.beta = abs(beta)*sqrt(df2 / (1 - beta^2)),
  gCalc = ((2*abs(beta))/sqrt(1-beta^2))*(1 - (3/(4*df2 - 1))),
  N = df2 + 2,
  beta.var = escalc(measure = "COR", ri = beta, ni = df2 + 2, data = dat[dat$finalDesign == "Beta.between", ])$vi,
  gVarCalc = (1 - (3/(4*df2 - 1))) * (4 * beta.var)/((1 - beta^2)^3),
  p = 2*pt(abs(t.from.beta), df2, lower.tail=FALSE)
)

# Show the converted ESs
dat %>% filter(finalDesign == "Beta.between") %>% select(gCalc, beta, gVarCalc, finalDesign, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$useBias == 1 & (dat$Design == "Between" | dat$Design == "Continuous") & !is.na(dat$beta) & !is.na(dat$p.reported), "biasTest"] <- "Beta.between"
dat$label[dat$biasTest == "Beta.between"] <- paste(dat[dat$biasTest == "Beta.between",]$paperID, "/", dat[dat$biasTest == "Beta.between",]$studyID, "/", dat[dat$biasTest == "Beta.between",]$Variable.Indicator, ": ",
                                                      "Z=", qnorm(1-dat[dat$biasTest == "Beta.between",]$p.reported/2), sep = "")

####
# chi^2 -------------------------------------------------------------------
####

# Specify the design, compute ES, var, N, and p
dat[dat$useMeta == 1 & dat$Design == "Between" & !is.na(dat$Chisq) & !is.na(dat$n.rep), "finalDesign"] <- "between.chisq"

dat[dat$finalDesign == "between.chisq", ]$gCalc <- esc_chisq(chisq = dat[dat$finalDesign == "between.chisq", ]$Chisq, totaln = dat[dat$finalDesign == "between.chisq", ]$n.rep, es.type = "g")$es
dat[dat$finalDesign == "between.chisq", ]$gVarCalc <- esc_chisq(chisq = dat[dat$finalDesign == "between.chisq", ]$Chisq, totaln = dat[dat$finalDesign == "between.chisq", ]$n.rep, es.type = "g")$var
dat[dat$finalDesign == "between.chisq", ]$N <- dat[dat$finalDesign == "between.chisq", ]$n.rep
dat[dat$finalDesign == "between.chisq", ]$p <- with(dat[dat$finalDesign == "between.chisq", ], 1 - pchisq(Chisq, 1))

# Show the converted ESs
dat %>% filter(finalDesign == "between.chisq") %>% select(gCalc, gVarCalc, Chisq, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$useBias == 1 & dat$Design == "Between" & !is.na(dat$Chisq) & !is.na(dat$n.rep), "biasTest"] <- "between.chisq"
dat$label[dat$biasTest == "between.chisq"] <- paste(dat[dat$biasTest == "between.chisq",]$paperID, "/", dat[dat$biasTest == "between.chisq",]$studyID, "/", dat[dat$biasTest == "between.chisq",]$Variable.Indicator, ": ",
                                                        "chi2(", dat[dat$biasTest == "between.chisq",]$df1, ")=", dat[dat$biasTest == "between.chisq",]$Chisq, sep = "")

####
# t for B in continuous designs -------------------------------------------
####

# Specify the design, compute ES, var, N, and p
dat[dat$useMeta == 1 & (dat$Design == "Continuous" |  dat$Design == "Correlation") & !is.na(dat$B) & !is.na(dat$t) & !is.na(dat$df2) & is.na(dat$beta), "finalDesign"] <- "t.from.B"

dat[dat$finalDesign == "t.from.B", ]$gCalc <- esc_t(t = abs(dat[dat$finalDesign == "t.from.B", ]$t), totaln = dat[dat$finalDesign == "t.from.B", ]$df2 + 2, es.type = "g")$es
dat[dat$finalDesign == "t.from.B", ]$gVarCalc <- esc_t(t = dat[dat$finalDesign == "t.from.B", ]$t, totaln = dat[dat$finalDesign == "t.from.B", ]$df2 + 2, es.type = "g")$var
dat[dat$finalDesign == "t.from.B", ]$N <- dat[dat$finalDesign == "t.from.B", ]$df2 + 2
dat[dat$finalDesign == "t.from.B", ]$p <- with(dat[dat$finalDesign == "t.from.B", ], 2*pt(abs(t), N - 1, lower.tail = FALSE))

# Show the converted ESs
dat %>% filter(finalDesign == "t.from.B") %>% select(gCalc, gVarCalc, B, t, N)

# Create a "result label" to be used as an input for p-curve analysis
dat[dat$useBias == 1 & (dat$Design == "Continuous" |  dat$Design == "Correlation") & !is.na(dat$B) & !is.na(dat$t) & !is.na(dat$df2) & is.na(dat$beta), "biasTest"] <- "t.from.B"
dat$label[dat$biasTest == "t.from.B"] <- paste(dat[dat$biasTest == "t.from.B",]$paperID, "/", dat[dat$biasTest == "t.from.B",]$studyID, "/", dat[dat$biasTest == "t.from.B",]$Variable.Indicator, ": ",
                                                    "t(", dat[dat$biasTest == "t.from.B",]$df2, ")=", dat[dat$biasTest == "t.from.B",]$t, sep = "")

# Multiply the ES by -1 if not in the predicted direction
dat$gCalc <- ifelse(test = (dat$predictedDirection == -1 & dat$gCalc > 0), 1 = (dat$gCalc*-1), 0 = dat$gCalc)


####
# Create a results label for the rest of bias=1 and meta=0, based on p-values, converted to z-score (p-curve app works with z but not with p).
dat[dat$useBias == 1 & is.na(dat$label) & !is.na(dat$paperID), "biasTest"] <- "z.from.p"
dat$p <- ifelse(is.na(dat$p), dat$p.reported, dat$p)
dat$label[dat$biasTest == "z.from.p"] <- paste(dat[dat$biasTest == "z.from.p",]$paperID, "/", dat[dat$biasTest == "z.from.p",]$studyID, "/", dat[dat$biasTest == "z.from.p",]$Variable.Indicator, ": ",
                                              "Z=", qnorm(1-(dat[dat$biasTest == "z.from.p",]$p)/2), sep = "")

# Remove unused variables
delvars <- names(dat) %in% c(
  "t.from.r", "t.from.beta", "r.var", "beta.var", "d.calc")
dat <- dat[!delvars]
rm(delvars)

