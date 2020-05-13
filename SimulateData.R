library(dplyr)
#simulate a dataset for the stress meta-analysis 

#General info of a study

set.seed(123)

research_design <-sample(c(1,2),346, replace=T)
# 1=RCT, 2=observational
journal <-sample(c("Psychological Bulletin","Trends in Cognitive Sciences","American Psychologist","Child Development"),346, replace=T)
database <-sample(c(1,2,3,4,5),346, replace=T)
# 1=psychinfo, 2=pubmed, 3=scopus, 4=google scholar, 5=hand search
country <-sample(c(NA,"italy","france","singapore","US","slovakia","netherlands"),346, replace=T)
nMale <-as.integer(runif(346,min=1, max=60))
nFemale <-as.integer(runif(346,min=1, max=60))
journalH5 <-as.integer(runif(346,min=0, max=10))
useMeta <-sample(c(0,1),346, replace=T)
useBias <-sample(c(0,1),346, replace=T)
predictedDirection <-sample(c(-1,1),346, replace=T)
items <-as.integer(runif(346,min=0, max=10))
# for items=  (0 if the DV is not a discrete/likert variable, 
#otherwise an integer count of, e.g., Likert scale items comprising the DV)

mean1 <- round(rnorm(346,5,1), 2)
mean2 <- round(rnorm(346,5,2), 2)
#mean of group1 (experimental) and of group2 (control)
sd1 <- round(rnorm(346,2,1), 2)
sd2 <- round(rnorm(346,2,1), 2)
#sd of group1 (experimental) and of group2 (control)
n1 <- as.integer(runif(346,min=1, max=100))
n2 <- as.integer(runif(346,min=1, max=100))
n3 <- as.integer(runif(346,min=1, max=100))
#sample size per cells, n1,n2,n3

############################################################################################

#Intervention characteristics 



presence_of_individual_differences <-sample(c(1,0),346, replace=T)
# 1=if there are individual differences in the study, 0= If there aren't individual differences in the study
focal_variable <-sample(c(0,1),346, replace=T)
# Added the focal variable; 0= not a focal variable, 1= it is a focal variable
timing_intervention <-sample(c(1,2,3),346, replace=T)
# 1= before stressor, 2= after the stressor, 3= both
presence_of_stressTest <-sample(c(1,2),346, replace=T)
# 1=present, 2=absent
MASdata<-as.data.frame(presence_of_stressTest)
Type_of_StressTest<-ifelse(MASdata$presence_of_stressTest==1, 
                    c(1,2,3,NA), NA)
MASdata$Type_of_StressTest <-Type_of_StressTest
#If the stressTest is present I code: 1= TSST, 2=stroop task, 3=others 999= no stress test
type_of_population <-sample(c(1,2),346, replace=T)
# 1=normal, 2= clinical
type_of_comparison_group <-sample(c(1,2),346, replace=T)
# 1=passive control group, 2= active control group
type_of_component <-sample(c(1,2,3,4,5,6),346, replace=T)
#1= Affect: low arousal, negative valence, 2= Affect: high arousal, negative valence, 3= Affect: low arousal, positive valence 4= Affect: high arousal, positive valence, 
#5= cognitive component 6=physiological component 
exact_type_of_population <- sample(c(NA,"nursersy students","employers","psychology students","patients with depression","COPD patients","managers"),346, replace=T)
frequency_of_intervention <- as.integer(runif(346,min=0, max=5))
#how many times each week the participants receive the intervention
duration_of_intervention <- as.integer(runif(346,min=0, max=80))
#total duration of the intervention in hours
number_of_intervention <- as.integer(runif(346,min=0, max=90))
#total number of SEM session, biofeedback trainings, exposure to nature, social support received
Instrument <-sample(c("PSS","STAI","cortisol","heart-rate","RSS","HADS"),346, replace=T)
#Type of scale used in the experiment (e.g. PSS)

nationality <-sample(c(NA,"italian","french","american","dutch","slovakian"),346, replace=T)

######################################################################

#Specific coding for categories

category <-sample(c(1,2,3,4),346, replace=T)
#Divide each article by categories for the meta-analysis: 1= self-administered mindfulness, 2=biofeedback 3= being in nature 4=social support
data1<-as.data.frame(category)

Type_of_Sam<-ifelse(data1$category==1, 
                    c(1,2,3,4,5), NA)
#If the article is on self-administered mindfulness, I code 1 = internet 2= smartphone app 3=book 4= audio 5= mixed
data1$Type_of_Sam<-Type_of_Sam
type_of_environmet <-ifelse(data1$category==2, 
                            c(1,2), NA)
#If the article is on being in nature, I code 1 = nature env, 2= built environment
data1$type_of_environmet <-type_of_environmet 

type_of_exposure <-ifelse(data1$category==2, 
                          c(1,2,3), NA)
#If the article is on being in nature, I code 1 = outdoor walk 2 = nature viewing 3= outdoor gardening 
data1$type_of_exposure<-type_of_exposure


type_of_SocialSupport <-ifelse(data1$category==4, 
                               c(1,2,3,4), NA) 
data1$type_of_SocialSupport<-type_of_SocialSupport
#If the article is on social support I code 1 =no support 2= physical contact 3= verbal social support 4= mixed 
source_of_SocialSupport <-ifelse(data1$category==4,                               
                                 c(1,2,3), NA) 
data1$source_of_SocialSupport <- source_of_SocialSupport
#If the article is on social support I code 1 =partner 2= friends 3=stranger                            
#Code the source of SocialSupport 

MetaData <-cbind(nMale,nFemale,journalH5,useMeta,useBias,predictedDirection,items,mean1,mean2,sd1,sd2,n1,n2,n3,research_design,focal_variable,journal,database,country,number_of_intervention,Instrument,presence_of_individual_differences,timing_intervention,MASdata,type_of_population,type_of_comparison_group,type_of_component,exact_type_of_population,frequency_of_intervention,duration_of_intervention,nationality,data1)
#Create the first simulated dataset with the info encoded until now
paperID <- 1:nrow(MetaData)
studyID <-1:nrow(MetaData)

dat <- read.csv("MA data oct 2018 IR.csv", sep = ";")
# Read in the social thermoregulation dataset in order to merge it with the one simulated for stress.

published <-sample(c(0,1),346, replace=T)
#1=published, 0=unpublished

StressData <- cbind(MetaData,studyID,paperID,dat$F, dat$log.reg.B,dat$B,dat$t,dat$r,dat$Chisq,dat$beta,dat$Waldchisq,dat$df1,dat$df2,dat$Design,published)
#Merging a simulated dataset on stress, with the real dataset of social thermoregulation (in order to take the effect sizes from that and other statistics)



rob2 <- read.csv("Rob_2.csv", sep = ";")
#bringing-in the Rob2 for each study, after having used the Rob2 excel sheet for each study



rob2 <- rob2[ , which(names(rob2) %in% c("Domain.1.risk.of.bias","Domain.2.risk.of.bias","Domain.3.risk.of.bias","Domain.4.risk.of.bias","Domain.5.risk.of.bias","Overall.risk.of.bias"))]
#taking from the raw Rob2 dataset just the columns we need (Rob for all domains and overall)

#rob2 <- rob2 %>%
#select("Domain.1.risk.of.bias","Domain.2.risk.of.bias","Domain.3.risk.of.bias","Domain.4.risk.of.bias","Domain.5.risk.of.bias","Overall.risk.of.bias")
#View(rob2)
#Here I tried to do the same with dplyr, it was just and exercise

StressData <-cbind(StressData, rob2)
#View(rob2)

#Merging the simulated dataset with the Rob2 simulated dataset 
gender <- c("MALE","FEMALE","FEMALE","UNKNOWN","MALE")
ifelse(gender == "MALE", 1, ifelse(gender == "FEMALE", 2, 3))

StressData$Domain.1.risk.of.bias<- ifelse(StressData$Domain.1.risk.of.bias == "Low", 1, ifelse(StressData$Domain.1.risk.of.bias == "High",3,2))
StressData$Domain.2.risk.of.bias<- ifelse(StressData$Domain.2.risk.of.bias == "Low", 1, ifelse(StressData$Domain.2.risk.of.bias == "High",3,2))
StressData$Domain.3.risk.of.bias<- ifelse(StressData$Domain.3.risk.of.bias == "Low", 1, ifelse(StressData$Domain.3.risk.of.bias == "High",3,2))
StressData$Domain.4.risk.of.bias<- ifelse(StressData$Domain.4.risk.of.bias == "Low", 1, ifelse(StressData$Domain.4.risk.of.bias == "High",3,2))
StressData$Domain.5.risk.of.bias<- ifelse(StressData$Domain.5.risk.of.bias == "Low", 1, ifelse(StressData$Domain.5.risk.of.bias == "High",3,2))
StressData$Overall.risk.of.bias<- ifelse(StressData$Overall.risk.of.bias == "Low", 1, ifelse(StressData$Overall.risk.of.bias == "High", 3, 2))
StressData$`dat$Design`<- ifelse(StressData$`dat$Design` == "Within", 1,2)

#recoded to numeric values items for RoB2 and type of design 
dat <- StressData
# View(dat)
# StressData <- StressData %>%
  #mutate (Domain.1.risk.of.bias = ifelse(Domain.1.risk.of.bias == "Low", 1,
                                       # ifelse(Domain.1.risk.of.bias == "High",3,2)) altro modo di fare ifelse
