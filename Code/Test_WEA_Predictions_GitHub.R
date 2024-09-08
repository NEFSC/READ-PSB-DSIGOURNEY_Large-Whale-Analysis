
library(tidyverse)
library(MASS)
library(mgcv)
library(ggplot2)
library(patchwork)
library(matrixStats)
   
setwd("C:/Large Whale Analysis")

#For github

#Pick species
cSpecies<-'SPWH'

#Pick Model
cModel<-10

#Pick trheshold
 # threshold=11 #PBR fin whales=11
 # threshold=9.7 #PBR humpback whales=1
  # threshold=11 #PBR minke whales=170 (Canadian East Coastal stock 2021 PBR)
 #threshold=1 #PBR sei whales=6.2 (Nova Scotia stock 2021 PBR)
 threshold=3.2 #PBR sperm whales=3.2

 #FIWH Mod 8
 # covs<-c("btemp","lat","dist125","dist1000")
 
 #HUWH Model 9
 # covs<-c("sstmur","dist2GSSw","dist125",
 #         "depth","btemp")
 
 #MIWH Mod 1
 # covs<-c("salinity","btemp","sstmur","poc","depth","dist125",
 #         "chla","slope","pic")
 # 
 #MIWH Mod 9
 # covs<-c("salinity","btemp","sstmur","depth","dist125",
 #         "chla","slope","pic")
 

 #SEWH Mod 2
 # covs<-c("sstmur","lon","dist200","mlp","depth","dist1000",
 # "pic","pp")
 
 #SPWH Mod 10
   covs<-c("slope","depth","btemp","salinity","dist1000","dist200")
         
 
#Pick Min/Max ExDet values
 minExDet=0.00023
 maxExDet=1
 
#Load model
load(paste("./BHDSM/",cSpecies,"/",cSpecies,"_Model_",cModel,"_NewCovs.RData",sep=""))

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
psi<-all.pars[,n.betas+3]


#Subset WEA Covariates    
 # WEA<-read.csv("../WEA grid with buffer.csv")
 # WEA1<-subset(WEA,Block==2)
 # grids1<-WEA1$Grid

WEA<-read.csv("./Mapping code/WEA_Grids_1_for_Plotting.csv")
WEA1<-subset(WEA,NewID==2)
grids1<-WEA1$TenKomerc

    

load("./RData files/amapps1021.RData")
 all.grids<-amapps1021 %>% rename(lat=latdd, lon=londd)
# all.grids<-subset(amapps1021, year<2017)
rm(amapps1021)

prediction_grids<- all.grids[all.grids$grid%in%grids1,]  
PG.covs<-prediction_grids[,covs] #Subset
prediction_grids<-prediction_grids[complete.cases(PG.covs),] 
rm(all.grids)  

#Remove extrapolations
source("./Extrapolations/Remove_Extraps_WEAs.R")
prediction_grids<-prediction_grids %>% left_join(wea.extraps, by=c("year","layer","grid")) %>% 
   dplyr::filter(remove==0) %>% 
   dplyr::select(-ExDet,-mic,-prob,-keep,-obs,-remove,-seasonum.y,-x,-y) %>% rename(seasonum=seasonum.x)

#Read in All Covariate Data
load("./Final_Cov_Mat_new.rda")
Final_Cov_Mat_All<-Final_Cov_Mat_All.new

Final_Cov_Mat<-subset(Final_Cov_Mat_All,SPECIES==cSpecies)
rm(Final_Cov_Mat_All)
rm(Final_Cov_Mat_All.new)
 
 cov_vec<-c("lat","lon","depth","dist2shore","slope","dist200","dist125","dist1000",
             "sstmur","chla","salinity","btemp","mlp","pic",
             "poc","pp","dist2GSNw","dist2GSSw","chlaf","sstf")
 
 DS<-Final_Cov_Mat[,cov_vec] #Subset
 Final_Cov_Mat<-Final_Cov_Mat[complete.cases(DS),] %>% filter(EFFORT_KM >0.75)
 N_t<-dim(Final_Cov_Mat)[1]
 

 #FIWH Mod 8
 # gam.mod <- gam(Nsights~s(btemp,bs="ts", k=5)+s(lat,bs="ts", k=5)+
 #                  s(dist125,bs="ts", k=5)+ s(dist1000,bs="ts", k=5),
 #                data=Final_Cov_Mat,family = poisson(link=log),
 #                file = "FIWH_jagam")
 
 
#HUWH Mod 9
 # gam.mod <- gam(Nsights~s(sstmur,bs="ts", k=5)+s(dist2GSSw,bs="ts", k=5)+
 #                   s(dist125,bs="ts", k=5)+s(depth,bs="ts", k=5)+
 #                   s(btemp,bs="ts", k=5),data=Final_Cov_Mat,
 #                   family = poisson(link=log), file = "HUWH_jagam")
 
 #MIWH Mod 1
  # gam.mod <- gam(Nsights~ s(salinity,bs="ts", k=5)+ s(btemp,bs="ts", k=5)+
  #               s(sstmur,bs="ts", k=5)+s(poc,bs="ts", k=5)+
  #               s(depth,bs="ts", k=5)+  s(dist125,bs="ts", k=5)+
  #               s(chla,bs="ts", k=5)+s(slope,bs="ts", k=5)+
  #               s(pic,bs="ts", k=5),data=Final_Cov_Mat,
  #                   family = poisson(link=log), file = "MIWH_jagam")
  
  #MIWH Mod 9
  # gam.mod <- gam(Nsights~ s(salinity,bs="ts", k=5)+ s(btemp,bs="ts", k=5)+
  #                  s(sstmur,bs="ts", k=5)+s(depth,bs="ts", k=5)+
  #                  s(dist125,bs="ts", k=5)+s(chla,bs="ts", k=5)+
  #                  s(slope,bs="ts", k=5)+s(pic,bs="ts", k=5),data=Final_Cov_Mat,
  #                family = poisson(link=log), file = "MIWH_jagam")
 
 #SEWH Mod 2
  # gam.mod <- gam(Nsights~s(sstmur,bs="ts", k=5)+ s(lon,bs="ts", k=5)+
  #              s(dist200,bs="ts", k=5)+s(mlp,bs="ts", k=5)+
  #              s(depth,bs="ts", k=5)+s(dist1000,bs="ts", k=5)+
  #              s(pic,bs="ts", k=5)+ s(pp,bs="ts", k=5),
  #              data=Final_Cov_Mat,family = poisson(link=log),
  #              file = "SEWH_jagam")
 
 #SPWH Mod10
 gam.mod <- gam(Nsights~s(slope,bs="ts", k=5)+ s(depth,bs="ts", k=5)+
                  s(btemp,bs="ts", k=5)+  s(salinity,bs="ts", k=5)+
                  s(dist1000,bs="ts", k=5)+s(dist200,bs="ts", k=5),
                data=Final_Cov_Mat,family = poisson(link=log),
                file = "SPWH_jagam")
 
 

pd<-prediction_grids %>% dplyr::select(all_of(covs))
 
   Xp<- predict(gam.mod, type = "lpmatrix", newdata = pd)
   # load("Xp_test.rda")
   # load("Xp_fin.rda")
 # N.tot.grids<-dim(Xp)[1]
 
 # Z<-matrix(rep(z.pars,N.tot.grids),ncol=mcmcsamps,byrow=TRUE )
 # preds <- Xp %*% b.pars+Z
 
mu.all <- exp(Xp %*% b.pars)
 
Estvars<-rowVars(mu.all)  
keep.vars<-which(Estvars<0.05)
mu<-mu.all[keep.vars,]
N.tot.grids<-dim(mu)[1]

   GS_Mat<-matrix(NA, nrow=N.tot.grids, ncol=mcmcsamps) 
   
    mat<-matrix(NA,nrow=N.tot.grids,ncol=mcmcsamps)
   for (i in 1:mcmcsamps){
  # r<-r_mat[i]
   r<-psi[i]
   gs<-(lambda.size[i]+1)
   p<-r/(mu[,i]*100+r)
   Ng<-rnbinom(N.tot.grids,size=r,prob=p)
   N<-Ng*gs
   mat[,i]<-N
     # mat[,i]<-mu[,i]*100*gs
   } 
  gridCols<-prediction_grids[keep.vars,]%>%dplyr::select(grid,layer,year,seasonum)
 df.final<-cbind(gridCols,mat)  

 #Summarize number of grids removed by year-layer
 gridCount<-df.final%>%dplyr::select(grid,year,layer,seasonum)%>%group_by(year,layer,seasonum)%>%dplyr::count()
 
    ncols<-mcmcsamps+3
   Abun1<-df.final%>%dplyr::select(-grid)%>%group_by(year,layer,seasonum)%>%summarise_all(funs(sum))

  # Prob1<-as.numeric(Abun1[,4:ncols]<threshold)
    Prob1<-matrix(as.numeric(Abun1[,4:ncols]>threshold), nrow=dim(Abun1)[1], ncol=mcmcsamps)
   Prob_mat<-cbind(Abun1[,1:3], rowSums(Prob1, na.rm=TRUE)/mcmcsamps)
   Prob_mat1<-as.data.frame(Prob_mat)
   colnames(Prob_mat1)<-c("year", "layer", "seasonum", "prob")
   
   
#Year specific Abundance and Prob
   cYear=2021
   Abun.year<-subset(Abun1, year==cYear)
   Prob.year<-matrix(as.numeric(Abun.year[,4:ncols]>threshold), nrow=dim(Abun.year)[1], ncol=mcmcsamps)
   Prob_mat.year<-cbind(Abun.year[,1:3], rowSums(Prob.year, na.rm=TRUE)/mcmcsamps)
   Prob_mat.year<-as.data.frame(Prob_mat.year)
   colnames(Prob_mat.year)<-c("year", "layer", "seasonum", "prob.year")
 
 #Add mean probability among years  
   Prob_mat2<-Prob_mat1%>%dplyr::select(-year)%>%group_by(layer)%>%summarise_each(mean)
   Prob_mat2$year<-'Mean'
    Prob_mat3<-Prob_mat2[,c(4,1,2,3)]
   
   Prob_mat.Final<-rbind(Prob_mat1,Prob_mat3)
  
   ProbPlot1<-ggplot(Prob_mat.Final, aes(x = layer, y = prob, colour = as.factor(year), group=year)) +
   geom_line()
   
   
#Abundance Plots   
year.vec<-unique(Abun1$year)
n.years<-length(year.vec)

layer.vec<-unique(Abun1$layer)
n.layers<-length(layer.vec)

new.mat<-data.frame<-(Layer1=seq(1,46,1))
for(j in 1:n.years){
cAbun.mat<-subset(Abun1, year==year.vec[j])
new.mat<-cbind(new.mat, cAbun.mat[,4:ncols])
}


# plot total abundance by summing "fits" and "predictions" of all sites

#plot(c(0,46),c(0, 100), type = 'n', xaxt ='n', xlab = 'Layer',
#	ylab = 'Abundance', cex.lab = 1.5, main = '',
#	cex.axis = 1.5, cex.main = 1.5)
#axis(1, at = c(1,10,20,30,4,46), labels = c(1,10,20,30,4,46),
#     cex.axis = 1.5)
#
#for(i in 1:(dim(new.mat)[2]-1)) {
#
#lines(1:46, new.mat[,i+1],
#	col = rgb(0, 0, 0, alpha = .05))
#}
#
#lines(1:n.layers, apply(new.mat,1,quantile, 0.975, na.rm=TRUE), lty = 2, lwd = 3, col = 'coral1')
#lines(1:n.layers, apply(new.mat,1,quantile, 0.5, na.rm=TRUE), lty = 1, lwd = 4, col = 'red')
#lines(1:n.layers, apply(new.mat,1,quantile, 0.025, na.rm=TRUE), lty = 2, lwd = 3, col = 'coral1')
#
 new.abun_mat<-data.frame(Layer=seq(1,46,1),Abundance=apply(new.mat,1,quantile, 0.5, na.rm=TRUE),
                           High= apply(new.mat,1,quantile, 0.975, na.rm=TRUE),
                           Low= apply(new.mat,1,quantile, 0.025, na.rm=TRUE))
                           
 AbunPlot<-ggplot(new.abun_mat, aes(Layer, Abundance)) +
 geom_line(lwd=2)+
    geom_hline(yintercept=threshold,  color = "red", size=2)+
 geom_ribbon(aes(ymin = Low, ymax = High), alpha = 0.2)+theme(text=element_text(face="bold", color="black",size=14),
 axis.title.x = element_blank(),axis.text.x = element_text(face="bold", color="black", size=16),
 axis.text.y = element_text(face="bold", color="black", size=16),
 axis.title.y = element_text(face="bold", color="black", size=16)) +ylab("Abundance")+
    scale_x_continuous(breaks=c(1,13,25,37), labels=c("Jan","Apr","Jul","Oct"))  +
 theme(axis.text.x = element_text(face="bold", size=14, angle=45))+
    labs(tag = "A") +
    theme(plot.tag = element_text(),
          plot.tag.position = c(0.87, 0.97))
 
 Prob.Mat2<-Prob_mat.Final%>%group_by(layer)%>%
# summarize(Med.Prob=median(prob, na.rm=TRUE), SD=sd(prob, na.rm=TRUE))%>%
#mutate(Diff1=Med.Prob-SD, Diff2=Med.Prob+SD)#%>%mutate(Low=max(0,Diff1), High=min(1,Diff2))
summarize(Med.Prob=median(prob, na.rm=TRUE), High=quantile(prob, 0.975, na.rm=TRUE),
Low=quantile(prob, 0.025, na.rm=TRUE))
 
 
 #max.fun<-function(x){
#c<-max(0,x) 
# return(c)
# }
# 
# min.fun<-function(x){
#a<-min(1,x)
# return(a)
# }
#
# Prob.Mat2$Low<-sapply(Prob.Mat2$Diff1,max.fun)
# Prob.Mat2$High<-sapply(Prob.Mat2$Diff2,min.fun)
 
 # ggplot(economics, aes(x=date)) + 
 #   geom_line(aes(y = psavert), color = "darkred") + 
 #   geom_line(aes(y = uempmed), color="steelblue", linetype="twodash") 
 # 
#
 Prob.Mat.All<-cbind(Prob.Mat2,Prob_mat.year$prob.year)
 
 ProbPlot2<-ggplot(Prob.Mat.All, aes(x=layer)) +
 geom_line(aes(y = Med.Prob), color = "black", lwd=2)+
   geom_line(aes(y = Prob_mat.year$prob.year), color = "green", lwd=2)+
 geom_ribbon(aes(ymin =Low, ymax = High), alpha = 0.2)+theme(text=element_text(face="bold", color="black",size=14),
 axis.title.x = element_blank(),axis.text.x = element_text(face="bold", color="black", size=16),
 axis.text.y = element_text(face="bold", color="black", size=16),
 axis.title.y = element_text(face="bold", color="black", size=16)) +ylab("Probability")+
    scale_x_continuous(breaks=c(1,13,25,37), labels=c("Jan","Apr","Jul","Oct"))+
 theme(axis.text.x = element_text(face="bold", size=14, angle=45))+
    labs(tag = "B") +
    theme(plot.tag = element_text(),
          plot.tag.position = c(0.87, 0.97))
 
 #ProbPlot1 + AbunPlot

#fiwh<- AbunPlot + ProbPlot2
#huwh<- AbunPlot + ProbPlot2
miwh<- AbunPlot + ProbPlot2
#sewh<- AbunPlot + ProbPlot2
# spwh<- AbunPlot + ProbPlot2
# 
#  fiwh/huwh/miwh/sewh/spwh
 
 
 
 
 
 
 
####################
#2-step predictions
#####################
 
 #Old covariates
 covs_old<-c("lat","lon","depth","dist2shore","slope","dist200","dist125","dist1000",
             "naoim","sstmur","chla","sla","mld","salinity","btemp","mlp","picma",
             "pocma","pp","dist2GSNw","dist2GSSw","chlfma","sstfma","sstfmt")
 
 #New covariates
 covs_new<-c("lat","lon","depth","dist2shore","slope","dist200","dist125","dist1000",
             "sstmur","chla","salinity","btemp","mlp","pic",
             "poc","pp","dist2GSNw","dist2GSSw","chlaf","sstf")
 

 
 new.covs<-Final_Cov_Mat %>% dplyr::select(Grid,Year,Layer, all_of(covs_new)) %>% distinct()
 
   
 ######################
 #HUWH
 #######################
 
 HUWH.data<-read.csv("./2-STEP/HUWH_2_Step_input_data.csv")
 
 #Get rid of NAs to make model slection Consistent
 HUWH.data<-HUWH.data%>%rename(Effort=EFFORT_KM)%>% dplyr::select(-covs_old) %>%
    left_join(new.covs,by=c("Grid","Layer","Year")) %>%
    dplyr::select(Grid,Nsights, W, phat, Effort,all_of(covs_new))%>%
    filter(Effort >0.75)%>% na.omit()

 #Fit best model 
 HUWH.best <- gam(Nsights~s(sstmur,bs="ts", k=5)+ s(btemp,bs="ts", k=5)+
                     s(salinity,bs="ts", k=5)+s(dist125,bs="ts", k=5)+
                     s(depth,bs="ts", k=5)+ s(pic,bs="ts", k=5),
               family=nb(theta=NULL, link="log"), select=TRUE,method="REML",
               offset=log(phat*2*W*Effort), data= HUWH.data)
 
 summary(HUWH.best)
 
 pd.test<-prediction_grids %>% dplyr::filter(year==2020 & layer==22) %>% dplyr::select(all_of(covs))
 preds.test<-predict(HUWH.best,newdata= pd.test,se=TRUE, type="response")
 sum(preds.test$fit*100*1.5)
 
 
 
 