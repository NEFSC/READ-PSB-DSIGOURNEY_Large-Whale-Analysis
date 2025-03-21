library(tidyverse)
library(furrr)
library(mgcv)
library(readr)
library(tictoc)


final.estimates.all<-data.frame(chain=NA,grid=NA, sp=NA, seasonum=NA,
                            density=NA,low=NA,high=NA,cv=NA)

Extraps.all<-data.frame(Year=NA,Layer=NA,Grid=NA,ExDet=NA,Estimate=NA,
                    Variance=NA,gIVH.Indicator=NA,keep=NA,obs=NA)

#Pick extrapolation cutoffs
U.ExDet=1
L.ExDet=0.00023
# gIVh.Off=1 #If 1 don't us gIVH if 0 use gIVH

tic()
for(z in 1:20){
  
  rm(list=setdiff(ls(),c("z","final.estimates.all","Extraps.all","U.ExDet",
                         "L.ExDet")))

#Choose number of cores
n.cores<-availableCores()-2
#set to parallel
plan("multisession", workers=n.cores)  #To run things in parallel

setwd("C:/Large Whale Analysis/BHDSM")
#get total number of available cores
# availableCores()

# n.cores<-availableCores()

#Pick species
cSpecies<-'SEWH'

#Pick season
cSeason<-'spring'
cSeason.num=1

#Pick Model
 cModel<-2

 
 #HUWH 
   # Mod 20
  # covs<-c("sstmur","dist2GSSw","dist125","depth","btemp")
  # covs<-c("sstmur","btemp","salinity","dist125","depth","pic")
 
   # Mod 17
     # covs<-c("sstmur","dist2GSSw","dist125", "depth","btemp","pp",
     #      "mlp", "dist1000")
 #MIWH 
  #Mod 1
 # covs<-c("salinity","btemp","sstmur","poc","depth","dist125",
 #         "chla","slope","pic")
  #Mod 10
  # covs<-c("salinity","btemp","sstmur","poc","dist2shore","dist125",
  #  "chla","slope","pic")
 
#SPWH Mod 18
# covs<-c("slope","depth","btemp","salinity","dist1000","sstmur",
#         "chlfma","chla","sstfma", "dist200","picma",
#         "sla")

#FIWH Mod 8
 # covs<-c("btemp","lat","dist125","dist1000")
 
 #SEWH 
 # Mod 2
 covs<-c("sstmur","lon","dist200","mlp","depth",
         "dist1000","pic","pp")
 # Mod 14
 # covs<-c("lat","dist200","mlp","slope")
 

# load(paste("/net/work7/dsigourney/AMAPPS Sims/New DSM Analyses/RData files/",cSeason,"_10-17.RData", sep=""))

#True Coastal
 # load(paste("/net/work7/dsigourney/AMAPPS Sims/New DSM Analyses/RData files/",cSeason,"_TRUE_coastal_10-17.RData", sep=""))
 load("../RData files/amapps1021.RData")


season.all<-amapps1021 %>% filter(seasonum==cSeason.num) %>% 
dplyr::rename(lat=latdd, lon=londd)
rm(amapps1021)


#First subset
n.grids1<-length(unique(season.all$grid))
subsets1<-20
n.subs1<-ceiling(n.grids1/subsets1)
ss.all1<-rep(seq(1,subsets1,1),n.subs1)
ss1<-ss.all1[1:n.grids1]
sub1<-data.frame(grid=unique(season.all$grid),subsets.all=ss1)
season<-merge(season.all,sub1,by='grid', all.x=TRUE) %>% filter(subsets.all==z) %>%
  select(-subsets.all)
rm(season.all)

#Make second subset for parallel
n.grids<-length(unique(season$grid))
n.subs<-ceiling(n.grids/n.cores)
ss.all<-rep(seq(1,n.cores,1),n.subs)
ss<-ss.all[1:n.grids]

sub<-data.frame(grid=unique(season$grid),subsets=ss)

allcovs<-merge(season,sub,by='grid', all.x=TRUE)


#Load data out put from Nimble
# load(paste("/net/work7/dsigourney/AMAPPS Sims/New DSM Analyses/Large Whale Analysis/BHDSM/",cSpecies,"/",cSpecies,"_Model_",cModel,"_Final.RData",sep=""))
 load(paste("./",cSpecies,"/",cSpecies,"_Model_",cModel,"_NewCovs.RData",sep=""))


n.chains<-length(output)


n.betas<-length(covs)*4+1
n.cols<-n.betas+3
df<-matrix(NA,nrow=1,ncol=n.cols)
mcmcsamps<-dim(cSamps<-output[[1]]$samples)[1]


for (i in 1:n.chains){
  
  cSamps<-output[[i]]$samples
  cols.betas<-grep(colnames(cSamps), pattern='beta')
  col.gs<-grep(colnames(cSamps), pattern='lambda.size')
  col.psi<-grep(colnames(cSamps), pattern='psi')
  cols.mu<-grep(colnames(cSamps), pattern='mu')
 

  cMat<-matrix(NA,nrow=mcmcsamps,ncol= n.cols)
  cMat[,1]=i
  cMat[,2:(n.betas+1)]=cSamps[,cols.betas]
  cMat[,n.betas+2]=cSamps[,col.gs]
  cMat[,n.betas+3]=cSamps[,col.psi]
  
  df<-rbind(df,cMat)
}

n.rows<-dim(df)[1] 
mcmcsamps<-n.rows-1
all.pars<-df[2:n.rows,]

# all.pars<-sample_n(as.data.frame(all.pars),mcmcsamps, replace=FALSE)

b.pars<-t(all.pars[,2:(n.betas+1)])
lambda.size<-all.pars[,n.betas+2]

getEst<-function(i){

  cSamps<-output[[i]]$samples
  cols.mu<-grep(colnames(cSamps), pattern='mu')
  cEstimates<-as.data.frame(cSamps[,cols.mu])

  return(cEstimates)

}

Estimates<-future_map_dfr(1:n.chains, ~getEst(.),.id='chain',options(future.globals.maxSize = 1000 * 1024^2))

Est_subset <- Estimates[,sapply(Estimates, is.numeric)]
varEst<-sapply(Est_subset, var)
  maxVar<-quantile(varEst, 0.999)
  maxVar<-max(varEst)
rm(Estimates)
rm(Est_subset)
maxVar<-1

#Read in All Covariate Data
load("../Final_Cov_Mat_new.rda")
Final_Cov_Mat_All<-Final_Cov_Mat_All.new
# #Subset species
Final_Cov_Mat<-subset(Final_Cov_Mat_All,SPECIES==cSpecies)
rm(Final_Cov_Mat_All)
rm(Final_Cov_Mat_All.new)

cov_vec<-c("lat","lon","depth","dist2shore","slope","dist200","dist125","dist1000",
           "sstmur","chla","salinity","btemp","mlp","pic",
           "poc","pp","dist2GSNw","dist2GSSw","chlaf","sstf")
#Subset covariates
# DS<-Final_Cov_Mat[,cov_vec] #Subset
# Final_Cov_Mat<-Final_Cov_Mat[complete.cases(DS),] %>% filter(EFFORT_KM >0.75)
# rm(DS)


source("Make_Predictions_func.R")


rm(season)
rm(output)
rm(all.pars)
rm(cMat)
# b.pars<-dplyr::sample_n(as.data.frame(t(b.pars)), 5000, replace=FALSE)
# b.pars<-t(b.pars)
# lambda.size<-lambda.size[1:5000]
# mcmcsamps<-5000

out_df<-future_map_dfr(1:n.cores,~makePreds(.,b.pars,lambda.size), .id='Rep',seed=TRUE,options(future.globals.maxSize = 10000 * 1024^2))
# save(out_df, file=paste("preds_",cSpecies,"_",cSeason,"_",z,".rda", sep=""))

source("summarize_estimates.R")

# save(final.estimates,file=
# paste0("./",cSpecies,"/","final_estimates_",cSpecies,"_",cSeason,"_Model_",cModel,"_subset_",z,".rda"))
# 
# save(Extraps,file=
# paste0("../Extrapolations/",cSpecies,"_Extrapolations/",cSpecies,"_",cSeason,"_Extraps_",z,".rda"))

final.estimates.all<-rbind(final.estimates.all,final.estimates)
Extraps.all<-rbind(Extraps.all,Extraps)

print(z)
z=z+1
print(z)


}
toc()

final.estimates.all<-final.estimates.all[!is.na(final.estimates$chain),]
Extraps.all<-Extraps.all[!is.na(Extraps.all$Year),]

save(final.estimates.all,file=
paste0("./",cSpecies,"/",cSpecies,"_Predictions/","final_estimates_",cSpecies,"_",cSeason,"_Model_",cModel,"_NewCovs_rep2.rda"))

save(Extraps.all,file=
paste0("../Extrapolations/",cSpecies,"_Extrapolations/",cSpecies,"_",cSeason,"_Model_",cModel,"_Extraps_NewCovs_rep2.rda"))





