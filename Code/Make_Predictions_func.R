
makePreds<-function(i,b.pars,lambda.size){
  
  # rm(All_Dat.rmNA)
  
  #Read in All Covariate Data
  # load("../Final_Cov_Mat_new.rda")
  load("C:/Large Whale Analysis/Final_Cov_Mat_new.rda")
  #
  # #Subset species
  Final_Cov_Mat_All<-Final_Cov_Mat_All.new
  Final_Cov_Mat<-subset(Final_Cov_Mat_All,SPECIES==cSpecies)
  rm(Final_Cov_Mat_All)
  rm(Final_Cov_Mat_All.new)

  #Subset covariates
  
  DS<-Final_Cov_Mat[,cov_vec] #Subset
  Final_Cov_Mat<-Final_Cov_Mat[complete.cases(DS),] %>% filter(EFFORT_KM >0.75)
  # 
  # Final_Cov_Mat$TOD_AVG[is.na(Final_Cov_Mat$TOD_AVG)]<-mean(Final_Cov_Mat$TOD_AVG, na.rm=TRUE)
  # Final_Cov_Mat$SWELL_AVG[is.na(Final_Cov_Mat$SWELL_AVG)]<-mean(Final_Cov_Mat$SWELL_AVG, na.rm=TRUE)
  # Final_Cov_Mat$SEASTATE_WAVG[is.na(Final_Cov_Mat$SEASTATE_WAVG)]<-mean(Final_Cov_Mat$SEASTATE_WAVG, na.rm=TRUE)
  # Final_Cov_Mat$GLAREMAG_WAVG[is.na(Final_Cov_Mat$GLAREMAG_WAVG)]<-mean(Final_Cov_Mat$GLAREMAG_WAVG, na.rm=TRUE)
  # Final_Cov_Mat$SUBJ_WAVG[is.na(Final_Cov_Mat$SUBJ_WAVG)]<-mean(Final_Cov_Mat$SUBJ_WAVG, na.rm=TRUE)
  
  #
  #HUWH 
   # Mod 17
     gam.mod <- gam(Nsights~s(sstmur,bs="ts", k=5)+ s(dist2GSSw,bs="ts", k=5)+
                      s(dist125,bs="ts", k=5)+ s(depth,bs="ts", k=5)+
                      s(btemp,bs="ts", k=5)+ s(pp,bs="ts", k=5)+ s(mlp,bs="ts", k=5)+
                      s(dist1000,bs="ts", k=5),data=Final_Cov_Mat, family = poisson(link=log))
                   
  
  # Mod 9
  # gam.mod <- gam(Nsights~s(sstmur,bs="ts", k=5)+s(dist2GSSw,bs="ts", k=5)+
  #                  s(dist125,bs="ts", k=5)+ s(depth,bs="ts", k=5)+
  #                  s(btemp,bs="ts", k=5),data=Final_Cov_Mat, family = poisson(link=log))
 
  # gam.mod <- gam(Nsights~s(sstmur,bs="ts", k=5)+ s(btemp,bs="ts", k=5)+
  #               s(salinity,bs="ts", k=5)+s(dist125,bs="ts", k=5)+
  #               s(depth,bs="ts", k=5)+ s(pic,bs="ts", k=5),
  #             data=Final_Cov_Mat, family = poisson(link=log),
  #             file = "HUWH_jagam")
   #MIWH Mod 1
  # gam.mod <- gam(Nsights~ s(salinity,bs="ts", k=5)+ s(btemp,bs="ts", k=5)+
  #                  s(sstmur,bs="ts", k=5)+s(poc,bs="ts", k=5)+
  #                  s(depth,bs="ts", k=5)+  s(dist125,bs="ts", k=5)+s(chla,bs="ts", k=5)+
  #                  s(slope,bs="ts", k=5)+ s(pic,bs="ts", k=5),
  #                data=Final_Cov_Mat, family = poisson(link=log))
  #MIWH Mod 10
  # gam.mod <- gam(Nsights~ s(salinity,bs="ts", k=5)+ s(btemp,bs="ts", k=5)+
  #                  s(sstmur,bs="ts", k=5)+s(poc,bs="ts", k=5)+
  #                  s(dist2shore,bs="ts", k=5)+  s(dist125,bs="ts", k=5)+s(chla,bs="ts", k=5)+
  #                  s(slope,bs="ts", k=5)+ s(pic,bs="ts", k=5),
  #                data=Final_Cov_Mat, family = poisson(link=log))
  
  #SPWH Mod 18
  # gam.mod <- gam( Nsights~s(slope,bs="ts", k=5)+ s(depth,bs="ts", k=5)+
  #                   s(btemp,bs="ts", k=5)+  s(salinity,bs="ts", k=5)+ s(dist1000,bs="ts", k=5)+
  #                   s(sstmur,bs="ts", k=5)+ s(chlfma,bs="ts", k=5)+
  #                   s(chla,bs="ts", k=5)+ s(sstfma,bs="ts", k=5)+s(dist200,bs="ts", k=5)+
  #                   s(picma,bs="ts", k=5)+ s(sla,bs="ts", k=5),
  #                   data=Final_Cov_Mat, family = poisson(link=log))
  
  #FIWH Model 19
  # gam.mod <- gam(Nsights~s(btemp,bs="ts", k=5)+s(lat,bs="ts", k=5)+
  #              s(dist125,bs="ts", k=5)+ s(dist1000,bs="ts", k=5),
  #              data=Final_Cov_Mat, family = poisson(link=log),
  #            file = "FIWH_jagam")
  
  
     #SEWH Mod 14
     gam.mod <- gam(Nsights~s(sstmur,bs="ts", k=5)+ s(lon,bs="ts", k=5)+
                      s(dist200,bs="ts", k=5)+s(mlp,bs="ts", k=5)+
                      s(depth,bs="ts", k=5)+s(dist1000,bs="ts", k=5)+
                      s(pic,bs="ts", k=5)+ s(pp,bs="ts", k=5))
     
  #SEWH Mod 14
  # gam.mod <- gam(Nsights~s(lat,bs="ts", k=5)+ s(dist200,bs="ts", k=5)+
  #                  s(mlp,bs="ts", k=5)+s(slope,bs="ts", k=5),
  #                data=Final_Cov_Mat, family = poisson(link=log))
     
     
     
  

 allcovs.sub<-allcovs %>% dplyr::filter(subsets==i) 
 pd<-allcovs.sub %>% dplyr::select(all_of(covs))
     
  Xp<- predict(gam.mod, type = "lpmatrix", newdata = pd)
  N.tot.grids<-dim(Xp)[1]

 # Z<-matrix(rep(z.pars,N.tot.grids),ncol=mcmcsamps,byrow=TRUE )
 # preds <- Xp %*% b.pars+Z
  preds <- Xp %*% b.pars


   
   GS_Mat<-matrix(NA, nrow=N.tot.grids, ncol=mcmcsamps)  
   
      for (j in 1:mcmcsamps){ 
     
      lam<-lambda.size[j]
      gs<-lam+1
     
      GS_Mat[,j]<-gs
      
      }
      
      pred_mat<-preds
     Lambda_Mat<-exp(pred_mat)
     Pred_Mat<-Lambda_Mat*GS_Mat
  
 
 
   Details_Mat<-matrix(NA, nrow=N.tot.grids, ncol=7)
   
     Grid<-allcovs.sub$grid
     Year<-allcovs.sub$year
     Layer<-allcovs.sub$layer
          
  for (s in 1:N.tot.grids){
      
       mean.grid<-mean(as.numeric(Lambda_Mat[s,]))
       var.grid<-var(as.numeric(Lambda_Mat[s,]))
       gIVH_Ind<-as.numeric(var.grid>maxVar)  
        NA_Ind<-as.numeric(is.na(mean.grid))
       
       
       Details_Mat[s,1]<-Grid[s]
       Details_Mat[s,2]<-Year[s]
       Details_Mat[s,3]<-Layer[s]
       Details_Mat[s,4]<-mean.grid
       Details_Mat[s,5]<-var.grid
       Details_Mat[s,6]<-gIVH_Ind
       Details_Mat[s,7]<-NA_Ind
       
     #  Details_Mat[i,8]<-pred.gam$fit
#       Details_Mat[i,9]<-pred.gam$se
       
         }
      
      
         Details_Mat<-as.data.frame(Details_Mat)
       #  colnames(Details_Mat)<-c("Grid","Year","Layer","Estimate","Variance","gIVH.Indicator", "NA.Indicator", "gam.Estimate", "gam.SE")
       colnames(Details_Mat)<-c("Grid","Year","Layer","Estimate","Variance","gIVH.Indicator", "NA.Indicator")
         
           All_Dat<-cbind(Details_Mat,Pred_Mat)
           All_Dat.rmNA<-subset(All_Dat,NA.Indicator==0)
            
            # rm(list=setdiff(ls(),c("All_Dat.rmNA")))
            
           return(All_Dat.rmNA)
          
}