library(tidyverse)
library(furrr)
library(mgcv)
library(readr)
library(tictoc)


# rm(list=ls())
#Choose number of cores
# n.cores<-14
#set to parallel
plan("multisession", workers=n.cores)  #To run things in parallel


#get total number of available cores
# availableCores()

#Subsets for processing final density estimates
n.subs.P=10


n.grids<-length(unique(out_df$Grid))
n.subs<-ceiling(n.grids/n.subs.P)
ss.all<-rep(seq(1,n.subs.P,1),n.subs)
ss<-ss.all[1:n.grids]
df<-data.frame(Grid=unique(out_df$Grid),subsets=ss)
est_all.F<-merge(out_df,df,by='Grid', all.x=TRUE)

est.func<-function(i){
  
  
  cSub<-est_all.2 %>% dplyr::filter(subsets==i & keep==1) %>%
    dplyr::select(-Year,-Layer,-ExDet,-keep,
                  -Rep,-Estimate,-Variance,-gIVH.Indicator,
                  -NA.Indicator,-subsets)
  
  
  test<- cSub%>% tidyr::gather("rep","estimates",2:5001)%>%
    group_by(Grid,rep)%>%summarise(M.estimates=mean(estimates)) %>% 
    summarise(Estimate=median(M.estimates), meanEst=mean(M.estimates), SD=sd(M.estimates),
              Upper=quantile(M.estimates,0.975), Lower=quantile(M.estimates,0.025))%>%
    dplyr::mutate(CV=SD/meanEst) 
  
  out<-data.frame(grid=test$Grid, sp=cSpecies, seasonum=cSeason.num,
                  density=test$Estimate,low=test$Lower,high=test$Upper,
                  cv=test$CV)
  
  return(out)
}

final.estimates<-data.frame(chain=NA,grid=NA, sp=NA, seasonum=NA,
                            density=NA,low=NA,high=NA,cv=NA)

Extraps<-data.frame(Year=NA,Layer=NA,Grid=NA,ExDet=NA,Estimate=NA,
                    Variance=NA,gIVH.Indicator=NA,keep=NA,obs=NA)

 for (j in 1:n.subs.P){
   
   est_all<-est_all.F %>% filter(subsets==j) %>% select(-subsets)
   
   source('../Extrapolations/Extrap.R')
   
   rm(est_all)
   
n.grids<-length(unique(extraps$Grid))
n.subs<-ceiling(n.grids/n.cores)
ss.all<-rep(seq(1,n.cores,1),n.subs)
ss<-ss.all[1:n.grids]
df<-data.frame(Grid=unique(extraps$Grid),subsets=ss)
est_all.2<-merge(extraps,df,by='Grid', all.x=TRUE)

rm(extraps)
rm(ss.all)
rm(ss)
rm(df)


est_df<-future_map_dfr(1:n.cores,~est.func(.), .id='chain', options(future.globals.maxSize = 10000 * 1024^2))

final.estimates<-rbind(final.estimates,est_df)
Extraps<-rbind(Extraps,sum.extraps)
rm(est_all.2)

 }

final.estimates<-final.estimates[!is.na(final.estimates$chain),]
Extraps<-Extraps[!is.na(Extraps$Year),]


