
library(quantreg)
library(numDeriv)
library(rootSolve)
library(msm)

source("required_funcs.r") ## required functions for the proposed method
source("qor_bs.R")         ## qor_bs is the main function for the proposed method

write.table(Data,"sampledata.txt", col.names = T, row.names = F, sep=" ")
Data<-read.table("sampledata.txt",header = T)

Y1<-Data$Y1
Y2<-Data$Y2

delta1<-Data$delta1
delta2<-Data$delta2
delta.c<-1-delta1*delta2
sum(delta.c)/length(delta.c)      # overall censoring rate= 0.21

tauseq<-c( .2, .3, .4, .5, .6)    # the quantile of interest, total quantiles: m
Zmat<-cbind(1, Data$Z1, Data$Z2)  # the design matrix : nxp

#####################################################################
#######  conduct the proposed method using qor_bs function  #########
#####################################################################

output<-qor_bs(Y1, Y2, delta1, delta2, Zmat, tauseq) 

################
output$est.beta1 ## estimate of beta1 (pxm)
output$est.beta2 ## estimate of beta2 (pxm)
output$est.gamma ## estimate of gamma (pxm)

output$est.sd.beta1 ## estimate of sd.beta1 (pxm)
output$est.sd.beta2 ## estimate of sd.beta2 (pxm)
output$est.sd.gamma ## estimate of sd.gamma (pxm)


