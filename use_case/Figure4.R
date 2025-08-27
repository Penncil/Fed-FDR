results <- readRDS("output_table1.rds")
fit_BHq = results$gold_standard
fit_local_wide_lasso = results$lasso
fit_SCAD_final = results$SCAD
fit_meta_BH = results$meta
fit_proposed = results$proposed


## read in data
sample_data = read.csv("sample_data_to_run.csv")[,-1]

## site
sites = sample_data$new_site

## covariates matrix
design_matrix = sample_data[,-c(1,2,ncol(sample_data))]
new_ind <- which(apply(design_matrix,2,mean) > 0.1)
X = design_matrix[,new_ind] %>% as.matrix()

## outcome vector
treat = sample_data[,"outcome"] %>% as.matrix()

##########################################################
########### ROC ------ Generate testing data
##########################################################

set.seed(1)

ind = sample(c(1:nrow(treat)), 0.8*nrow(treat))
y.train = treat[ind]
x.train = X[ind,]

y.test = y.train[-ind]
x.test = x.train[-ind,]

##########################################################
########### PLOT ROC curve
##########################################################
####function for obtain roc and auc under selected support
fit.glm.roc<-function(support,y.train,x.train,y.test,x.test){
  x.data<-data.frame(x.train[,support])  
  fit<-glm(y.train~.,data = x.data,family = "binomial")  
  
  dat<-data.frame(x.test[,support])
  prdict<-predict(fit,newdata = dat,type = "response")
  pred<-prediction(prdict,y.test)
  roc<-performance(pred,"tpr","fpr")
  auc<-performance(pred,measure = "auc")@y.values[[1]]
  return(list(roc=roc,auc=auc))
}


########obtain roc for each method

sup_BHq<- which(fit_BHq$FDR == 1)
fit.BHq<-fit.glm.roc(support = sup_BHq,y.train,x.train,y.test,x.test)
print(paste("Gold standard auc is", fit.BHq$auc))
roc.BHq<-fit.BHq$roc


sup_Fed_FDR_N<-which(fit_local_wide_lasso$FDR == 1)
fit.N<-fit.glm.roc(sup_Fed_FDR_N,y.train,x.train,y.test,x.test)
print(paste("Lasso auc is", fit.N$auc))
roc.N<-fit.N$roc

sup_meta_BHq<-which(fit_meta_BH$FDR == 1)
fit.meta_BHq<-fit.glm.roc(sup_meta_BHq,y.train,x.train,y.test,x.test)
print(paste("meta auc is", fit.meta_BHq$auc))
roc.meta_BHq<-fit.meta_BHq$roc


sup_Fed_SCAD<-which(fit_SCAD_final$FDR  == 1)
fit.SCAD<-fit.glm.roc(sup_Fed_SCAD,y.train,x.train,y.test,x.test)
print(paste("SCAD auc is", fit.SCAD$auc))
roc.SCAD<-fit.SCAD$roc


sup_Fed_FDR_C<-which(fit_proposed$FDR == 1)
fit.C<-fit.glm.roc(sup_Fed_FDR_C,y.train,x.train,y.test,x.test)
print(paste("our_1 auc is", fit.C$auc))
roc.C<-fit.C$roc


x.BHq<-roc.BHq@x.values[[1]]
y.BHq<-roc.BHq@y.values[[1]]
loess_fit <- loess(y.BHq ~ x.BHq, span = 0.5)  
y.BHq.smooth <-predict(loess_fit)


x.C<-roc.C@x.values[[1]]
y.C<-roc.C@y.values[[1]]
loess_fit <- loess(y.C ~ x.C, span = 0.5)  
y.C.smooth <-predict(loess_fit)

x.SCAD<-roc.SCAD@x.values[[1]]
y.SCAD<-roc.SCAD@y.values[[1]]
loess_fit <- loess(y.SCAD ~ x.SCAD, span = 0.5)  
y.SCAD.smooth <-predict(loess_fit)

x.N<-roc.N@x.values[[1]]
y.N<-roc.N@y.values[[1]]
loess_fit <- loess(y.N ~ x.N, span = 0.5)  
y.N.smooth <-predict(loess_fit)

x.meta_BHq<-roc.meta_BHq@x.values[[1]]
y.meta_BHq<-roc.meta_BHq@y.values[[1]]
loess_fit <- loess(y.meta_BHq ~ x.meta_BHq, span = 0.5)  
y.meta_BHq.smooth <-predict(loess_fit)


par(mgp=c(2,0.5,0),mar=c(3,3,1,1))
plot(x.BHq, y.BHq.smooth,
     xlab="FPR",ylab="TPR", xlim=c(0,1),ylim=c(0,1),
     col="#FFB000",type="l",lty=1,lwd=2)
par(new=TRUE)
plot(x.C, y.C.smooth,xlab="",ylab="", axes=FALSE, xlim=c(0,1),ylim=c(0,1), col="#8C3333",type="l",lty=1,lwd=2)
par(new=TRUE)
plot(x.SCAD, y.SCAD.smooth,xlab="",ylab="",axes=FALSE,xlim=c(0,1),ylim=c(0,1), col="#D2691E",type="l",lty=1,lwd=2)
par(new=TRUE)
plot(x.N, y.N.smooth,xlab="",ylab="",axes=FALSE,xlim=c(0,1),ylim=c(0,1), col="#016A70",type="l",lty=2,lwd=2)
par(new=TRUE)
plot(x.meta_BHq, y.meta_BHq.smooth,xlab="",ylab="",axes=FALSE,xlim=c(0,1),ylim=c(0,1), col="#7A9D54",type="l",lty=2,lwd=2)



labelname<-c(paste("BHq, AUC = ", round(fit.BHq$auc,3)),
             paste("Fed-FDR, AUC = ", round(fit.C$auc,3)),
             paste("Fed-FDR-SCAD, AUC = ",round(fit.SCAD$auc,3)),
             paste("Fed-FDR-N, AUC = ", round(fit.N$auc,3)),
             paste("Meta-BHq, AUC = ", round(fit.meta_BHq$auc,3)))
cols <- c("#FFB000","#8C3333","#D2691E","#016A70","#7A9D54")
legend(x=0.31,y=0.45, legend=labelname,bty="o", col=cols,lwd=2,lty=c(1,1,1,2,2))




