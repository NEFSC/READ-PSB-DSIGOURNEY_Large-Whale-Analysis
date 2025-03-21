
  path <- Sys.getenv('PATH')
newPath <- paste("C:\\Rtools\\bin;C:\\Rtools\\mingw_64\\bin;",
path, sep = "")
Sys.setenv(PATH = newPath)

#library(rjags)  
library(mgcv)
library(nimble)
library(dplyr)
library(furrr) # also attaches 'foreach' and 'parallel'

#Choose number of cores
n.cores<-5

#Set to parallel
plan("multisession", workers=n.cores)

#Run 1
setwd("/net/work7/dsigourney/AMAPPS Sims/New DSM Analyses")


#Run Models
cSpecies<-"FIWH"

#MRDS Groups
  
   # NES_Group<-10
#    SES_Group<-5
    # NEA_Group<-1
    # NEACB_Group<-1
    SEA_Group<-4
#Import MRDS DataFrame
    
 Dat_Sightings<-read.csv("MRDS_Data_All.csv")
 

 #NES Group
 nes.group <- c("FIWH", "FISE", "SEWH", "MIWH", "BLWH")  #Took out RIWH
 N.species_1_1 <- length(nes.group)
 SpeciesID_1_1<-seq(1:N.species_1_1)
 df1_1<-data.frame(Species=nes.group, SpeciesID=SpeciesID_1_1)
 Dat_Sightings_1_1<-Dat_Sightings%>%filter(Survey=="NES" & Species%in%nes.group)%>%
   full_join(df1_1, by="Species")
 
 W_1_1<-unique(Dat_Sightings_1_1$W.R)
 W.L_1_1<-unique(Dat_Sightings_1_1$W.L)
 
 #NEA Groups
 #Final Grouping
 nea.group_Final <- c("FIWH", "MIWH","HUWH","SEWH","SPWH","RIWH", "FISE", "BEWH", "GOBW", "SOBW","NBWH")
 
 nea.group <-  nea.group_Final

 N.species_2 <- length(nea.group)
 SpeciesID_2<-seq(1:N.species_2)
 df2<-data.frame(Species=nea.group, SpeciesID=SpeciesID_2)
 Dat_Sightings_2_1<-Dat_Sightings%>%filter((Survey=="NEA" | Survey=="CB") & Species%in%nea.group)%>%
   full_join(df2, by="Species")%>%mutate(SpeciesID=ifelse(Species=="FISE",4,SpeciesID))%>%
   mutate(SpeciesID=ifelse(Species=="RIWH",5,SpeciesID))%>%
   mutate(SpeciesID=ifelse(Species=="BEWH",2,SpeciesID))%>%
   mutate(SpeciesID=ifelse(Species=="GOBW",2,SpeciesID))%>%
   mutate(SpeciesID=ifelse(Species=="SOBW",2,SpeciesID))%>%
   mutate(SpeciesID=ifelse(Species=="NBWH",5,SpeciesID))
 
 N.species_2_1 <- length(unique(Dat_Sightings_2_1$SpeciesID))
 
 #Set new truncation distances
 W_2_1<-900
 W.L_2_1<-35
 
 Dat_Sightings_2_1<-Dat_Sightings_2_1%>%filter(Distance>W.L_2_1 & Distance<W_2_1)
 
   
 #SES Group
 ses.group <- c("FIWH", "HUWH", "RIWH", "SPWH",  "MIWH", "FISE", "SEWH")

 N.species_3 <- length(ses.group)
 SpeciesID_3<-seq(1:N.species_3)
 df3<-data.frame(Species=ses.group, SpeciesID=SpeciesID_3)
 Dat_Sightings_1_2<-Dat_Sightings%>%filter(Survey=="SES" & Species%in%ses.group)%>%
   full_join(df3, by="Species")%>%mutate(SpeciesID=ifelse(Species=="SEWH",5,SpeciesID))
 #merged MIWH and SEWH into one SpeciesID
 
 N.species_1_2 <- length(unique(Dat_Sightings_1_2$SpeciesID))
 
 W_1_2<-unique(Dat_Sightings_1_2$W.R)
 W.L_1_2<-unique(Dat_Sightings_1_2$W.L)
 
 #SE Aerial 
 Dat_Sightings_SE_Plane<-subset(Dat_Sightings, Survey=="SEA" & Group==SEA_Group)
 Dat_Sightings_SE_Plane<-subset(Dat_Sightings_SE_Plane, y01==0) #Get rid of back team observations
 
 Dat_Sightings_2_2.All<-Dat_Sightings_SE_Plane
 W_2_2<-Dat_Sightings_2_2.All$W.R[1]
 W.L_2_2<-Dat_Sightings_2_2.All$W.L[1]
 Dat_Sightings_2_2<- subset(Dat_Sightings_2_2.All, Distance>=W.L_2_2)  
 N_Sights_2_2<-dim(Dat_Sightings_2_2)[1]
 ones.dist_2_2<-array(NA,N_Sights_2_2)
 
 for (i in 1:N_Sights_2_2){
   ones.dist_2_2[i]<-1
 } 
 
          
Dat_Sightings_sp<-subset(Dat_Sightings, Species==cSpecies)
N.size<-dim(Dat_Sightings_sp)[1]
gs<-Dat_Sightings_sp$Group.size
gs[309]=1

#Read in All Covariate Data
load("Final_Cov_Mat_new.rda")
Final_Cov_Mat_All<-Final_Cov_Mat_All.new

#Subset species
source("Separate_Fin_Sei_AddFin.R")

# Final_Cov_Mat<-subset(Final_Cov_Mat_All,SPECIES==cSpecies)
rm(Final_Cov_Mat_All)
rm(Final_Cov_Mat_All.new)

#Subset only covariates in final model
cov_vec<-c("lat","lon","depth","dist2shore","slope","dist200","dist125","dist1000",
           "sstmur","chla","salinity","btemp","mlp","pic",
           "poc","pp","dist2GSNw","dist2GSSw","chlaf","sstf")

DS<-Final_Cov_Mat[,cov_vec] #Subset
Final_Cov_Mat<-Final_Cov_Mat[complete.cases(DS),] %>% filter(EFFORT_KM >0.75)
N_t<-dim(Final_Cov_Mat)[1]


Final_Cov_Mat$TOD_AVG[is.na(Final_Cov_Mat$TOD_AVG)]<-mean(Final_Cov_Mat$TOD_AVG, na.rm=TRUE)
Indicator_1_1<-as.numeric(Final_Cov_Mat$survey=="NES") #NE Ship Indicator
Indicator_1_2<-as.numeric(Final_Cov_Mat$survey=="SES") #SE Ship Indicator
Indicator_2_1<-as.numeric(Final_Cov_Mat$survey=="NEA" |Final_Cov_Mat$survey=="CB" ) #NE Aerial Indicator
Indicator_2_2<-as.numeric(Final_Cov_Mat$survey=="SEA") #SE Aerial Indicator


#Pick GAM 
#FIWH   
jd <-jagam( Nsights~s(btemp,bs="ts", k=5)+ s(lat,bs="ts", k=5)+
              s(dist125,bs="ts", k=5)+s(dist1000,bs="ts", k=5),
              data=Final_Cov_Mat, family = poisson(link=log),
             file = "FIWH_jagam")

#GAM Data  
X<-jd$jags.data$X
S1<-jd$jags.data$S1
S2<-jd$jags.data$S2
S3<-jd$jags.data$S3
S4<-jd$jags.data$S4


zero<-jd$jags.data$zero

  #GAM Inits
  b<-jd$jags.ini$b
  lambda<-jd$jags.ini$lambda
#  b[1]=-8.16
  
  
  Sightings=Final_Cov_Mat$Nsights
  NonZeroObs<-which(Sightings>0)
  Nnonzero<-length(NonZeroObs)

  ZeroObs<-which(Sightings==0)
  Nzero<-length(ZeroObs)
 
 #g(0)
     #NE Surveys
        g0_est_NE=0.67
        g0_CV_NE=0.16
         
           c_g0_NE=(g0_est_NE*(1-g0_est_NE)/(g0_est_NE*g0_CV_NE)^2)-1
           a_g0_NE=c_g0_NE*g0_est_NE  
           b_g0_NE=c_g0_NE*(1-g0_est_NE)  
      
        #SE Surveys
          g0_est_SE=0.86
           g0_CV_SE=0.179
      
           c_g0_SE=(g0_est_SE*(1-g0_est_SE)/(g0_est_SE*g0_CV_SE)^2)-1
           a_g0_SE=c_g0_SE*g0_est_SE  
           b_g0_SE=c_g0_SE*(1-g0_est_SE)  
      
      
      a_g0<-c(a_g0_NE, a_g0_SE)
      b_g0<-c(b_g0_NE, b_g0_SE)
     
  #Availability  bias
      Av_est=0.374
      Av_CV=0.336
  # Av_CV=0.01
     c_Av=((Av_est*(1-Av_est))/((Av_est*Av_CV)^2))-1
     a_Av=c_Av*Av_est
     b_Av=c_Av*(1-Av_est)  
     
      #Perception Bias Plane for FIWH
     P_est=0.752
     P_CV=0.171
     c_SEA=((P_est*(1-P_est))/((P_est*P_CV)^2))-1
     a_SEA=c_SEA*P_est
     b_SEA=c_SEA*(1-P_est)  
      
     zeros<-matrix(0,nrow=N_t, ncol=1)
         
       
 #For Hazard Rate Model
    steps<-50


#For temporal autocorrleation    
 # xcoord = unique(unique(Final_Cov_Mat$Year))       
 # year.mat=-as.matrix(dist(xcoord))   
 year<-Final_Cov_Mat$Year-min(Final_Cov_Mat$Year)+1
 n.years<-length(unique(Final_Cov_Mat$Year))      
 # mu_vec<-array(0,n.years)  
 # Z=runif(n.years,-0.1,0.1) 
    
    
 #Run Species-specific model    

#NB
    source("Large Whale Analysis/BHDSM/FIWH/FIWH_Model_8_NewCovs.R") 
# # 
   save(output, file="Large Whale Analysis/BHDSM/FIWH/FIWH_Model_8_NewCovs_with_FinSeiSep.RData")

#CPG
# source("Large Whale Analysis/FIWH/FIWH_Model_CPG.R")
# save(samples, file="Large Whale Analysis/FIWH/FIWH_Nimble_CPG_1.RData")

