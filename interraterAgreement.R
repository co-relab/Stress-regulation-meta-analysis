#+ setup, include=FALSE

library(tidyverse)
library(irr)

firstCoder <- read_delim("FIRST - SAM.csv", delim = ",")
secondCoder <- read_delim("SECOND - SAM.csv", delim = ",")

# Check of the first 10% (remove it from the comment to calculate the interrater agreement of the first 10%)
# ratersData <- list(firstCoder[1:11,], secondCoder[1:11,])
# 20% overall
ratersData <- list(firstCoder[1:25,], secondCoder[1:25,])

variableNames <- names(secondCoder)[c(5,6,10,15,17:21,23,25:29, 31, 34,37:40, 42:43, 45)]

# kappas <- list()
# for(i in variableNames){
#   kappas[[i]] <- kappa2(cbind(select(.data = ratersData[[1]], contains(variableNames[i])), select(.data = ratersData[[2]], contains(variableNames[i]))))
# }
# names(kappas) <- variableNames

agreement <- list()
for(i in variableNames){
  agreement[[i]] <- agree(cbind(ratersData[[1]][i], ratersData[[2]][i]))
}
names(agreement) <- variableNames


#+ include=TRUE
#'# Inter-rater agreement Mindfulness
agreement

rm(list = ls())
#+ include=FALSE
firstCoder <- read_delim("FIRST - BIO.csv", delim = ",")
secondCoder <- read_delim("SECOND - BIO.csv", delim = ",")

# Check of the first 10% (remove it from the comment to calculate the interrater agreement of the first 10%)
#ratersData <- list(firstCoder[1:12,], secondCoder[1:12,])
# Overall 20% 
ratersData <- list(firstCoder[1:27,], secondCoder[1:27,])

variableNames <- names(secondCoder)[c(5,6,10,15,17:22,24,26:29, 31, 34,38:41, 43:44, 46)]

# kappas <- list()
# for(i in variableNames){
#   kappas[[i]] <- kappa2(cbind(select(.data = ratersData[[1]], contains(variableNames[i])), select(.data = ratersData[[2]], contains(variableNames[i]))))
# }
# names(kappas) <- variableNames

agreement <- list()
for(i in variableNames){
  agreement[[i]] <- agree(cbind(ratersData[[1]][i], ratersData[[2]][i]))
}
names(agreement) <- variableNames


#+ include=TRUE
#'# Inter-rater agreement Biofeedback
agreement

