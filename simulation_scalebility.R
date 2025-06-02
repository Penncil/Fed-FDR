########simulation GLM
source("Fed_FDR_functions.R", echo=TRUE)


###sparsity level 
set.seed(123)
K <- 5 # Suppose there are 10 sites (include the local)
n <- 500#sample size in each site
p<-500
supp<-50
Corr<-0.5
signal<-7
theta<-c(rep(1,supp),rep(0,p-1-supp))


result3<-NULL
for (supp in c(30,35,40,45,50)) {
  re<-1
  while(re<=50) {
    tryCatch({
      data.tran<-data.generate(cov.type=2,
                               betashift =FALSE,
                               K=K,
                               n=n,
                               p=p,
                               #Gamma=rep(-0.5,K),
                               #Corr=rep(Corr,K),
                               Gamma=seq(-0.5,0.5,length.out=K),
                               Corr=seq(0.4,0.6,length.out=K),
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
      result<-cbind(result,supp)

      result3<-rbind(result3,result)
      re<-re+1
    }, error=function(e){
      cat("At repnum:", re, "Error:", conditionMessage(e), "\n")
    })
    cat("Rep=",re,"\n")
  }
  cat("supp=",supp,"\n")
}

name<-paste("n",n,"p",p,"corr",Corr,"signal",signal,"supp.rds",sep = "")
#name<-paste("result-cov1-h5-k",k,"-sepF.rds",sep = "")
saveRDS(result3,file=name)



###local sample size 
set.seed(123)
K <- 5 # Suppose there are 10 sites (include the local)
n <- 500#sample size in each site
p<-500
supp<-50
Corr<-0.5
signal<-7
theta<-c(rep(1,supp),rep(0,p-1-supp))


result2<-NULL
for (n in c(500,750,1000,1250,1500)) {
  re<-1
  while(re<=50) {
    tryCatch({
      data.tran<-data.generate(cov.type=2,
                               betashift =FALSE,
                               K=K,
                               n=n,
                               p=p,
                               #Gamma=rep(-0.5,K),
                               #Corr=rep(Corr,K),
                               Gamma=seq(-0.5,0.5,length.out=K),
                               Corr=seq(0.4,0.6,length.out=K),
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
      result<-cbind(result,n)
      
      result2<-rbind(result2,result)
      re<-re+1
    }, error=function(e){
      cat("At repnum:", re, "Error:", conditionMessage(e), "\n")
    })
    cat("Rep=",re,"\n")
  }
  cat("supp=",supp,"\n")
}

name<-paste("p",p,"corr",Corr,"signal",signal,"samplesize.rds",sep = "")
#name<-paste("result-cov1-h5-k",k,"-sepF.rds",sep = "")
saveRDS(result2,file=name)



########simulation GLM varying site number K
set.seed(123)

K <- 5 # Suppose there are 10 sites (include the local)
n <- 500#sample size in each site
p<-500
supp<-50
Corr<-0.3
signal<-7
theta<-c(rep(1,supp),rep(0,p-1-supp))


result4<-NULL
for (K in c(4,6,8,10,12)) {
  re<-1
  while(re<=50) {
    tryCatch({
      data.tran<-data.generate(cov.type=2,
                               betashift =FALSE, 
                               K=K,
                               n=n,
                               p=p,
                               #Gamma=rep(-0.5,K),
                               #Corr=rep(Corr,K),
                               Gamma=seq(-0.5,0.5,length.out=K),
                               Corr=seq(0.3,0.5,length.out=K),
                               supp=supp,
                               signal=signal)
      
      fit0<-local.glm.lasso(Treatall=data.tran$Treatall,
                            Xall=data.tran$Xall,
                            site=data.tran$site,
                            weight = FALSE)
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
      result<-cbind(result,K)
      
      result4<-rbind(result4,result)
      re<-re+1
    }, error=function(e){
      cat("At repnum:", re, "Error:", conditionMessage(e), "\n")
    })
    cat("Rep=",re,"\n")
  }
  cat("K=",K,"\n")
}

name<-paste("n",n,"p",p,"corr",Corr,"signal",signal,"K.rds",sep = "")
#name<-paste("result-cov1-h5-k",k,"-sepF.rds",sep = "")
saveRDS(result4,file=name)





###lambda tuning 
set.seed(123)
K <- 5 # Suppose there are 10 sites (include the local)
n <- 500#sample size in each site
p<-500
supp<-50
Corr<-0.5
signal<-7
theta<-c(rep(1,supp),rep(0,p-1-supp))


result3<-NULL

re<-1
while(re<=50) {
  tryCatch({
    data.tran<-data.generate(cov.type=2,
                             betashift =FALSE,
                             K=K,
                             n=n,
                             p=p,
                             #Gamma=rep(-0.5,K),
                             #Corr=rep(Corr,K),
                             Gamma=seq(-0.5,0.5,length.out=K),
                             Corr=seq(0.4,0.6,length.out=K),
                             supp=supp,
                             signal=signal)
    for (lambda_scale in c(0.5,1,1.5,2,2.5)) {
      fit0<-local.glm.lasso(Treatall=data.tran$Treatall,
                            Xall=data.tran$Xall,
                            site=data.tran$site,
                            lambda_scale = lambda_scale)
      ####our method Fed-FDR
      fit<-MFDR(fit0$Theta,0.1)
      Fed_FDR<-Tpr(fit$FDR,theta)
      Fed_FDR
      
      # fit0<-local.glm(Treatall=data.tran$Treatall,
      #                 Xall=data.tran$Xall,
      #                 site=data.tran$site)
      # 
      # #####Fed-FDR-N
      # fit<-MFDR(fit0$delasso,0.1)
      # Fed_FDR_N<-Tpr(fit$FDR,theta)
      # Fed_FDR_N
      # 
      # ####metha+BHq method
      # fit<-BH(fit0$pvalue,0.1)
      # Meta_BHq<-Tpr(fit$FDR,theta)
      # Meta_BHq
      
      
      ####Fed-FDR-Adapt
      fit1<-local.glm.scad(Treatall=data.tran$Treatall,
                           Xall=data.tran$Xall,
                           site=data.tran$site,
                           lambda_scale = lambda_scale)
      
      fit<-MFDR(fit1$Theta,0.1)
      Fed_FDR_scad<-Tpr(fit$FDR,theta)
      Fed_FDR_scad
      
      
      ####BHq###
      fit0<-BHq(Treatall=data.tran$Treatall,
                Xall=data.tran$Xall)
      
      fit<-BH(fit0$pvalue,0.1)
      FDR_BHq<-Tpr(fit$FDR,theta)
      
      result<-rbind(Fed_FDR,Fed_FDR_scad)
      #method<-c("Fed","metaBH","FedBH","MDS","Knockoff","BHq" )
      result<-cbind(result,lambda_scale)
      
      result3<-rbind(result3,result)
    }
    re<-re+1
  }, error=function(e){
    cat("At repnum:", re, "Error:", conditionMessage(e), "\n")
  })
  cat("Rep=",re,"\n")
}


name<-paste("n",n,"p",p,"corr",Corr,"signal",signal,"lambda.rds",sep = "")
#name<-paste("result-cov1-h5-k",k,"-sepF.rds",sep = "")
saveRDS(result3,file=name)



