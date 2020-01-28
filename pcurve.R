###################################################################################################################
# Adapted p-curve code from p-curve.com
#
#This is the R Code behind the p-curve app 4.052
#Last updated: 2017 03 17
#Written by Uri Simonsohn (urisohn@gmail.com)
#
###################################################################################################################

#The app is a single function, pcurve_app() .
#To run it you need to store statistical tests in a textfile with the format from the example, see http://p-curve.com/app4/example.txt
#Then, you just run the function on that file. For exmaple, if you save the file "example2.txt" in the folder "c:\data" then you would run:
#pcurve_app("example2.txt","c://data"). The app generates various text files with the results in table format and saves the figures
# as .PNG files on the same folder as where you put the .txt file
#
#This R Code was written specifically to run in the back-end of the online app and thus it may do things in a way that's not the most intuitive or most efficient
#for a user actually interacting with the R Code. The goal is to exactly replicate what is done on the website, thus the R Code is left intact. The website runs the exact same code
#below.

###################################################################################################################


library(stringr)  #Library to process string variables (text of the entered tests)
library(poibin)   #This library has the poisson binomial, the distribution of the sum of binomial with different underlying probabilities

pcurve_app=function(x)
{
  #0.1 Set up parameters
  #used to compute the binomial test given that each test has a (slightly) different probability of p<.025 depending on its own noncentral parameter
  #See Hong (2013) - "On computing the distribution function for the Poisson binomial distribution" Computational Statistics and Data Analysis, V59, p.41-51 - http://dx.doi.org/10.1016/j.csda.2012.10.006
  #setwd(dir1)                               #Set as Working Directory the folder on server where temporary files are saved
  #filek=substr(file1,1,nchar(file1)-4)      #filek is the name of the file entered by the user/server, dropping the extension
  x
  #Disable scientific notation
  options(scipen=999)

  ##############################################
  #(0) CREATE A FEW FUNCTIONS
  ##############################################
  #Function 1 - functions that find non-centrality parameter for f,chi distributions that gives some level of power

  #F-test
  #Note: starting with app 4.0, t() are converted to F() and Z to chisq() to avoid unnecessary repetition of code
  #So we only need code to find ncp for F() and chisq()

  getncp.f =function(df1,df2, power)   {
    error = function(ncp_est, power, x, df1,df2) pf(x, df1 = df1, df2=df2, ncp = ncp_est) - (1-power)
    xc=qf(p=.95, df1=df1,df2=df2)
    return(uniroot(error, c(0, 1000), x = xc, df1 = df1,df2=df2, power=power)$root)  }


  #chisq-test
  getncp.c =function(df, power)   {
    xc=qchisq(p=.95, df=df)
    error = function(ncp_est, power, x, df)      pchisq(x, df = df, ncp = ncp_est) - (1-power)
    return(uniroot(error, c(0, 1000), x = xc, df = df, power=power)$root)   }

  #Combine both in single function
  getncp=function(family,df1,df2,power) {
    if (family=="f") ncp=getncp.f(df1=df1,df2=df2,power=power)
    if (family=="c") ncp=getncp.c(df=df1,power=power)
    return(ncp)  }

  ###############################################################################
  #Function 2 - percent() : makes a number look like a percentage
  percent <- function(x, digits = 0, format = "f", ...)   {
    paste(formatC(100 * x, format = format, digits = digits, ...), "%", sep = "")
  }
  ###############################################################################


  ###############################################################################
  #Function 3 - pbound: bound p-values and pp-values by precision of measurement to avoid errors
  pbound=function(p) pmin(pmax(p,2.2e-16),1-2.2e-16)


  #Function 4 - prop33(pc) - Computes % of p-values that are expected to be smaller than pc,
  #for the tests submitted to p-curve, if power is 33%
  prop33=function(pc)
  {
    #pc: critical  p-value

    #Overview:
    #Creates a vector of the same length as the number of tests submitted to p-curve, significant and not,
    #    and computes the proportion of p-values expected to be smaller than {pc} given the d.f.
    #    and outputs the entire vector, with NA values where needed

    #F-tests (& thus  t-tests)
    prop=ifelse(family=="f" & p<.05,1-pf(qf(1-pc,df1=df1, df2=df2),df1=df1, df2=df2, ncp=ncp33),NA)
    #Chi2 (& thus Normal)
    prop=ifelse(family=="c" & p<.05,1-pchisq(qchisq(1-pc,df=df1),  df=df1, ncp=ncp33),prop)
    #output it
    prop
  }

  #Function 5 Stouffer test for a vector of pp-values
  stouffer=function(pp) sum(qnorm(pp),na.rm=TRUE)/sqrt(sum(!is.na(pp)))



  ###############################################################################
  #(1) PROCESS USER INPUT AND CONVERT IT INTO USABLE R DATA
  ###############################################################################


  #(1.1) Load data


  #From file;
  raw = x      #read file on server
  #  raw = scan(file=file1,what="")      #read file on server
  raw=tolower(raw)                                      #lower case
  ktot=length(raw)                                      #count studies

  #Create vector that numbers studies 1 to N,includes n.s. studies
  k=seq(from=1,to=length(raw))

  #1.2 Parse the entered text into usable statistical results
  #1.3 Create test type indicator
  stat=substring(raw,1,1)          #stat:   t,f,z,c,r
  test=ifelse(stat=="r","t",stat)  #test:   t,f,z,c      (r-->t)

  #1.4 Create family to turn t-->F and z-->chi2
  family=test
  family=ifelse(test=="t","f",family)
  family=ifelse(test=="z","c",family)

  #Note on terminology:
  #Stat:   t,f,c,z,r  is what the user entered, t,f,c,z,r
  #test:   t,f,c,z    is the test statistic, same as stat but with r-->t
  #family: f,c        converting t-->f and z-->c

  #1.5 Find comma,parentheses,equal sign
  par1 =str_locate(raw,"\\(")[,1]         #(  First  parenthesis
  par2 =str_locate(raw,"\\)")[,1]         #)  Second parenthesis
  comma=str_locate(raw,",")[,1]           #,  comma
  eq   =str_locate(raw,"=")[,1]           #=  equal

  #1.6 DF for t-tests
  df=as.numeric(ifelse(test=="t",substring(raw,par1+1,par2 -1),NA))             #t(df) later assigned to df2 in  F test with df1=1

  #1.7 DF1 for all tests
  #   recall, normal=sqrt(chi(1)) so df1=1 for Normal, same f(1,df)<-t(df)
  df1=as.numeric(ifelse(test=="f",substring(raw,par1+1,comma-1),NA))            #If F test, take value after comma, NA otherwise
  df1=as.numeric(ifelse(test=="z",1,df1))                                       #If Z replace missing value with a 1
  df1=as.numeric(ifelse(test=="t",1,df1))                                       #If t, replace missing value with a 1
  df1=as.numeric(ifelse(test=="c",substring(raw,par1+1,par2 -1),df1))           #If c, replace missing value with value in ()

  #1.8 DF2 for F(df1,df2) tests
  df2=as.numeric(ifelse(test=="f",substring(raw,comma+1,par2-1),NA))            #F-test
  df2=as.numeric(ifelse(test=="t",df,df2))                                      #t(df) into the df2 F test

  #1.9 Take value after equal sign, the value of the test-statistic, and put it in vector "equal"
  equal=abs(as.numeric(substring(raw,eq+1)))  #if not a r(), take the value after the ="

  #1.10  Go from "equal" (the value after = sign) to F or Chi2 value,
  value=ifelse((stat=="f" | stat=="c"),equal,NA)                      #For F and Chi2 test, equal=value
  value=ifelse(stat=="r", (equal/(sqrt((1-equal**2)/df2)))**2,value)  #For correlation, first turn value (r) to t, then square t. (using t=r/sqrt(1-r**2)/DF)
  value=ifelse(stat=="t", equal**2 ,value)                            #For t and Z, square it since t(df)**2=f(1,df) and z**2=chi(1)
  value=ifelse(stat=="z", equal**2 ,value)



  #1.11 Compute p-values
  p=ifelse(family=="f",1-pf(value,df1=df1,df2=df2),NA)
  p=ifelse(family=="c",1-pchisq(value,df=df1),p)
  p=pbound(p)  #Bound it to level of precision, see function 3 above

  #1.12 Count  studies
  #ktot is all studies
  ksig= sum(p<.05,na.rm=TRUE)     #significant studies
  khalf=sum(p<.025,na.rm=TRUE)    #half p-curve studies

  ###############################################################################
  #(2) COMPUTE PP-VALUES
  ##############################################################################

  #2.1 Right Skew, Full p-curve
  ppr=as.numeric(ifelse(p<.05,20*p,NA))            #If p<.05, ppr is 1/alpha*p-value, so 20*pvalue, otherwise missing.
  ppr=pbound(ppr)                                  #apply pbound function to avoid 0


  #2.2 Right Skew, half p-curve
  ppr.half=as.numeric(ifelse(p<.025,40*p,NA))    #If p<.05, ppr is 40*pvalue, otherwise missing.
  ppr.half=pbound(ppr.half)

  #2.3 Power of 33%
  #2.3.1 NCP for  f,c distributions
  # NCP33 (noncentrality parameter giving each test in p-curve 33% power given the d.f. of the test)
  ncp33=mapply(getncp,df1=df1,df2=df2,power=1/3,family=family)  #See function 1 above

  #2.3.2 Full-p-curve
  #Using the ncp33 compute pp33
  pp33=ifelse(family=="f" & p<.05,3*(pf(value, df1=df1, df2=df2, ncp=ncp33)-2/3),NA)
  pp33=ifelse(family=="c" & p<.05,3*(pchisq(value, df=df1, ncp=ncp33)-2/3),pp33)
  pp33=pbound(pp33)

  #2.3.3 HALF-p-curve
  #Share of p-values expected to be p<.025 if 33% power (using Function 4 from above, prop33() )
  prop25=3*prop33(.025)
  prop25.sig=prop25[p<.05]


  #Compute pp-values for the half
  pp33.half=ifelse(family=="f" & p<.025, (1/prop25)*(    pf(value,df1=df1,df2=df2,ncp=ncp33)-(1-prop25)),NA)
  pp33.half=ifelse(family=="c" & p<.025, (1/prop25)*(pchisq(value,df=df1,         ncp=ncp33)-(1-prop25)),pp33.half)
  pp33.half=pbound(pp33.half)


  ###############################################################################
  #(3) INFERENCE - STOUFFER & BINOMIAL
  ##############################################################################

  #3.1 Convert pp-values to Z scores, using Stouffer function above
  Zppr =     stouffer(ppr)            #right skew  - this is a Z value from Stouffer's test
  Zpp33=     stouffer(pp33)           #33% - idem
  Zppr.half= stouffer(ppr.half)       #right skew, half p-curve - idem
  Zpp33.half=stouffer(pp33.half)      #33% skew, half p-curve - idem

  #3.2 Overall p-values from Stouffer test
  p.Zppr =pnorm(Zppr)
  p.Zpp33=pnorm(Zpp33)
  p.Zppr.half =pnorm(Zppr.half)
  p.Zpp33.half=pnorm(Zpp33.half)

  #3.3 Save results to file
  main.results=as.numeric(c(ksig, khalf, Zppr, p.Zppr, Zpp33, p.Zpp33, Zppr.half, p.Zppr.half, Zpp33.half, p.Zpp33.half))

  #3.4 BINOMIAL
  #Observed share of p<.025
  prop25.obs=sum(p<.025)/sum(p<.05)
  #3.4.1 Flat null
  binom.r=1-pbinom(q=prop25.obs*ksig- 1, p=.5, size=ksig)     #The binomial in R computes the probability of x<=xo. We want prob(x>=x0) so we subtract one from x, and 1-prob()
  #3.4.2 Power of 33% null
  binom.33=ppoibin(kk=prop25.obs*ksig,pp=prop25[p<.05])

  #syntax for ppoibin():
  #   kk: is the proportion of succeses, a scalar, in this case, the share of p<.025
  #   pp: is the probabilty of success for each attempt, the number of attempts is determined
  #    by the length of the vector. For example ppoibin(kk=0,pp=c(.5,.5,.5)) is .125,
  #    if there are three attempts, each with probability .5, the odds of getting 0 succeses is .125
  #     ppoibin(kk=1,pp=c(1,.75)), in turn computes the probability of getting 1 success
  #     when one has a 100% of success, and the other 75%, and the solution is .25, since
  #     the first one succeeds for sure and the second would need to fail, with 25% chance.


  #3.4.3  Save binomial results
  binomial=c(binom.r, binom.33)

  ################################################
  #(4) POWER ESTIMATE
  ################################################

  #4.1 Function powerfit(power_est) - Returns the Stouffer Z of observing at least as right skewed a p-curve if  power=power_est
  #if Z=0, power_est is the best fit (p=50%).
  #if Z<0 the best fit is <power_est,
  #if Z>0 the best fit is >power_est
  powerfit=function(power_est)
  {
    #4.1.1 Get the implied non-centrality parameters (ncps) that give power_est to each test submitted to p-curve
    ncp_est=mapply(getncp,df1=df1,df2=df2,power=power_est,family=family)
    #4.1.2 Compute implied pp-values from those ncps_est,
    pp_est=ifelse(family=="f" & p<.05,(pf(value,df1=df1,df2=df2,ncp=ncp_est)-(1-power_est))/power_est,NA)
    pp_est=ifelse(family=="c" & p<.05,(pchisq(value,df=df1,ncp=ncp_est)-(1-power_est))/power_est,pp_est)
    pp_est=pbound(pp_est)
    #4.1.3 Aggregate pp-values for null that power=power_est via Stouffer
    return(stouffer(pp_est))   #This is a z score, so powerfit is expressed as the resulting Z score.
  }


  #4.2 COMPUTE FIT FOR EACH POWER for 5.1%, AND THEN 6-99%, AND PLOT IT. With power=5% boundary condition lead to errors
  #This becomes the diagnostic plot and gives us the best estimate, within 1%, of power.

  # Fit will be evaluated at every possible value of power between 5.1% and 99% in steps of 1%, stored in fit()
  fit=c()                                          #Create empty vector
  fit=abs(powerfit(.051))                      #Start it eavaluting fit of 5.1% power
  for (i in 6:99)   fit=c(fit,abs(powerfit(i/100))) #Now do 6% to 99%

  # Find the minimum
  #which ith power level considered leads to best estimate
  mini=match(min(fit,na.rm=TRUE),fit)
  #convert that into the power level, the ith value considered is (5+ith)/100
  hat=(mini+4)/100

  #4.3 Confidence interval for power estimate
  #4.3.1 Function get.power_pct(pct)
  get.power_pct =function(pct)   {
    #Function that finds power that gives p-value=pct for the Stouffer test
    #for example, get.power_pct(.5) returns the level of power that leads to p=.5  for the stouffer test.
    #half the time we would see p-curves more right skewed than the one we see, and half the time
    #less right-skewed, if the true power were that get.power_pct(.5). So it is the median estimate of power
    #similarliy, get.power_pct(.1) gives the 10th percentile estimate of power...
    #Obtain the normalized equivalent of pct, e.g., for 5% it is -1.64, for 95% it is 1.64
    z=qnorm(pct)  #convert to z because powerfit() outputs a z-score.
    #Quantify gap between computed p-value and desired pct
    error = function(power_est, z)  powerfit(power_est) - z
    #Find the value of power that makes that gap zero, (root)
    return(uniroot(error, c(.0501, .99),z)$root)   }

  #4.3.2 Boundary conditions (when the end of the ci=5% or 99% we cannot use root to find it,
  #use boundary value instead)

  #Boundary conditions
  p.power.05=pnorm(powerfit(.051)) #Proability p-curve would be at least at right-skewed if power=.051
  p.power.99=pnorm(powerfit(.99))  #Proability p-curve would be at least at right-skewed if power=.99

  #4.3.3 Find lower end of ci
  #Low boundary condition? If cannot reject 5% power, don't look for lower levels, use 5% as the end
  if (p.power.05<=.95) power.ci.lb=.05
  #High boundary condition? If we reject 99%, from below dont look for higher power, use 99% as the low end
  if (p.power.99>=.95) power.ci.lb=.99
  #If low bound is higher than 5.1% power and lower than 99% power, estimate it, find interior solution
  if (p.power.05>.95 && p.power.99<.95)  power.ci.lb=get.power_pct(.95)


  #4.3.4 Higher end of CI
  #If we reject 5% power from below, 5% is above the confidence interval, use 5% as the upper end of the confidence interval
  if (p.power.05<=.05) power.ci.ub=.05
  #If we do not reject that 99% power, don't look higher, use 99% as the higher end
  if (p.power.99>=.05) power.ci.ub=.99
  #If the the upper bound is between 5% and 99%, find it
  if (p.power.05>.05 && p.power.99<.05) power.ci.ub=get.power_pct(.05)


  #4.4 Save power fit estiate and ci
  power_results=c(power.ci.lb,hat,power.ci.ub)
  c(main.results, binomial, power_results)

  #Note, I use hat as the estimate of power, with powerfit(.5) we could get a more precise best fitting
  #level of power than the minimum in the figure above between .051 and .99, hat, but more precision than 1% in power is not informative.
}