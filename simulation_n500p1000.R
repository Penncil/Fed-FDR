########simulation GLM
source("Fed_FDR_functions.R", echo=TRUE)

########simulation GLM varying signal strength
set.seed(123)
K <- 5 # Suppose there are 10 sites (include the local)
n <- 500#sample size in each site
p<-1000
supp<-50
Corr<-0.3
signal<-7
theta<-c(rep(1,supp),rep(0,p-1-supp))


result1<-NULL
for (signal in c(7,8,9,10,11)) {
  re<-1
  while(re<=50) {
    tryCatch({
      data.tran<-data.generate(cov.type=2,
                               betashift =FALSE, 
                               K=K,
                               n=n,
                               p=p,
                               Gamma=rep(-0.5,K),
                               Corr=rep(Corr,K),
                               # Gamma=seq(-0.5,0,length.out=K),
                               # Corr=seq(0.2,0.8,length.out=K),
                               supp=supp,
                               signal=signal)
      
      fit0<-local.glm.lasso(Treatall=data.tran$Treatall,
                            Xall=data.tran$Xall,
                            site=data.tran$site)
      ####our method Fed-FDR
      fit<-MFDR(fit0$Theta,0.1)
      Fed_FDR<-Tpr(fit$FDR,theta)
      Fed_FDR
      
      fit0<-local.glm(Treatall=data.tran$Treatall,
                      Xall=data.tran$Xall,
                      site=data.tran$site)
      #####Fed-FDR-N
      fit<-MFDR(fit0$delasso,0.1)
      Fed_FDR_N<-Tpr(fit$FDR,theta)
      Fed_FDR_N
      
      ####metha+BHq method
      fit<-BH(fit0$pvalue,0.1)
      Meta_BHq<-Tpr(fit$FDR,theta)
      Meta_BHq
      
      
      ####Fed-FDR-Adapt
      fit1<-local.glm.scad(Treatall=data.tran$Treatall,
                           Xall=data.tran$Xall,
                           site=data.tran$site)
      
      fit<-MFDR(fit1$Theta,0.1)
      Fed_FDR_scad<-Tpr(fit$FDR,theta)
      Fed_FDR_scad
      
      
      ####BHq###
      fit0<-BHq(Treatall=data.tran$Treatall,
                Xall=data.tran$Xall)
      
      fit<-BH(fit0$pvalue,0.1)
      FDR_BHq<-Tpr(fit$FDR,theta)
      
      result<-rbind(Fed_FDR,Fed_FDR_N,Meta_BHq,Fed_FDR_scad,FDR_BHq)
      #method<-c("Fed","metaBH","FedBH","MDS","Knockoff","BHq" )
      result<-cbind(result,signal)
      result1<-rbind(result1,result)
      re<-re+1
    }, error=function(e){
      cat("At repnum:", re, "Error:", conditionMessage(e), "\n")
    })
    cat("Rep=",re,"\n")
  }
  
  cat("signal=",signal,"\n")
}

name<-paste("n",n,"p",p,"Corr",Corr,"sig.rds",sep = "")
#name<-paste("result-cov1-h5-k",k,"-sepF.rds",sep = "")
saveRDS(result1,file=name)



########simulation GLM varying correlation
set.seed(123)
K <- 5 # Suppose there are 10 sites (include the local)
n <- 500#sample size in each site
p<-1000
supp<-50
Corr<-0.7
signal<-9
theta<-c(rep(1,supp),rep(0,p-1-supp))


result2<-NULL
for (Corr in c(0.3,0.4,0.5,0.6,0.7)) {
  re<-1
  while(re<=50) {
    tryCatch({
      data.tran<-data.generate(cov.type=2,
                               betashift =FALSE, 
                               K=K,
                               n=n,
                               p=p,
                               Gamma=rep(-0.5,K),
                               Corr=rep(Corr,K),
                               # Gamma=seq(-0.5,0,length.out=K),
                               # Corr=seq(0.2,0.8,length.out=K),
                               supp=supp,
                               signal=signal)
      
      fit0<-local.glm.lasso(Treatall=data.tran$Treatall,
                            Xall=data.tran$Xall,
                            site=data.tran$site)
      ####our method Fed-FDR
      fit<-MFDR(fit0$Theta,0.1)
      Fed_FDR<-Tpr(fit$FDR,theta)
      Fed_FDR
      
      fit0<-local.glm(Treatall=data.tran$Treatall,
                      Xall=data.tran$Xall,
                      site=data.tran$site)
      #####Fed-FDR-N
      fit<-MFDR(fit0$delasso,0.1)
      Fed_FDR_N<-Tpr(fit$FDR,theta)
      Fed_FDR_N
      
      ####metha+BHq method
      fit<-BH(fit0$pvalue,0.1)
      Meta_BHq<-Tpr(fit$FDR,theta)
      Meta_BHq
      
      
      ####Fed-FDR-Adapt
      fit1<-local.glm.scad(Treatall=data.tran$Treatall,
                           Xall=data.tran$Xall,
                           site=data.tran$site)
      
      fit<-MFDR(fit1$Theta,0.1)
      Fed_FDR_scad<-Tpr(fit$FDR,theta)
      Fed_FDR_scad
      
      
      ####BHq###
      fit0<-BHq(Treatall=data.tran$Treatall,
                Xall=data.tran$Xall)
      
      fit<-BH(fit0$pvalue,0.1)
      FDR_BHq<-Tpr(fit$FDR,theta)
      
      result<-rbind(Fed_FDR,Fed_FDR_N,Meta_BHq,Fed_FDR_scad,FDR_BHq)
      #method<-c("Fed","metaBH","FedBH","MDS","Knockoff","BHq" )
      result<-cbind(result,Corr)
      
      result2<-rbind(result2,result)
      re<-re+1
    }, error=function(e){
      cat("At repnum:", re, "Error:", conditionMessage(e), "\n")
    })
    cat("Rep=",re,"\n")
  }
  cat("Corr=",Corr,"\n")
}

name<-paste("n",n,"p",p,"sig",signal,"corr.rds",sep = "")
#name<-paste("result-cov1-h5-k",k,"-sepF.rds",sep = "")
saveRDS(result2,file=name)




