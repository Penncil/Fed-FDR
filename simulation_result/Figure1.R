
library(ggplot2)
library(ggh4x)


######Signal plot

#####n500 p500
result1<-readRDS("n500p500Corr0.5sig.rds")
result1[which(is.na(result1))]<-0
n<-dim(result1)[1]
value<-c(result1[,1:2])
type<-rep(c("FDR","Power"),each=n)
labelname<-c("Fed-FDR","Fed-FDR-N","Meta-BHq","Fed-FDR-SCAD","BHq")
method<-rep(labelname,n/length(labelname))
Correlation<-rep(result1[,3],2)
a<-data.frame(value,type,method,Correlation)
res1<-aggregate(a$value,by=list(a$method,a$type,a$Correlation),mean )
samp<-rep(c("n=500, p=500"),nrow(res1))
res1<-data.frame(res1,samp)
colnames(res1)<-c("method","type","Correlation","value","scale")

#####n500 p1000
result1<-readRDS("n500p1000Corr0.3sig.rds")
result1[which(is.na(result1))]<-0
n<-dim(result1)[1]
value<-c(result1[,1:2])
type<-rep(c("FDR","Power"),each=n)
labelname<-c("Fed-FDR","Fed-FDR-N","Meta-BHq","Fed-FDR-SCAD","BHq")
method<-rep(labelname,n/length(labelname))
Correlation<-rep(result1[,3],2)
a<-data.frame(value,type,method,Correlation)
res2<-aggregate(a$value,by=list(a$method,a$type,a$Correlation),mean )
samp<-rep(c("n=500, p=1000"),nrow(res2))
res2<-data.frame(res2,samp)
colnames(res2)<-c("method","type","Correlation","value","scale")

#####n1000 p500
result1<-readRDS("n1000p500Corr0.5sig.rds")
result1[which(is.na(result1))]<-0
n<-dim(result1)[1]
value<-c(result1[,1:2])
type<-rep(c("FDR","Power"),each=n)
labelname<-c("Fed-FDR","Fed-FDR-N","Meta-BHq","Fed-FDR-SCAD","BHq")
method<-rep(labelname,n/length(labelname))
Correlation<-rep(result1[,3],2)
a<-data.frame(value,type,method,Correlation)
res3<-aggregate(a$value,by=list(a$method,a$type,a$Correlation),mean )
samp<-rep(c("n=1000, p=500"),nrow(res3))
res3<-data.frame(res3,samp)
colnames(res3)<-c("method","type","Correlation","value","scale")

res<-rbind(res1,res2,res3)


#cols <- c("#FFB000","#8C3333","#D2691E","#FF7256","#016A70","#7A9D54","#000000","#F2EE9D")
cols <- c("#FFB000","#8C3333","#016A70","#D2691E","#7A9D54")
plot <- ggplot2::ggplot(res, ggplot2::aes(x = Correlation, y = value, color = method))
plot1 <- plot + 
  geom_point(aes(shape = method,size=method))+
  scale_shape_manual(values = c(10,16,3,17,7))+
  scale_size_manual(values =rep(3.5,7))+
  scale_color_manual(values = cols)+
  geom_line(aes(linetype = method),size = 0.8) + 
  #scale_linetype_manual(values = c("solid","solid","solid","twodash","solid","solid","twodash"))+
  scale_linetype_manual(values = c("solid","solid","twodash","solid","twodash"))+
  facet_grid(type~scale,scales = "free")+
  # facet_manual(vars(type,scale),design =matrix(c(1,2,1,2,1,2),2,3),
  #              strip.position = "left",
  #              scales = "free_y",
  #              heights =c(1,2))+
  labs(x = "Singnal strength", y="")+
  #ggtitle(c("n=500, p=500"))+
  #scale_y_continuous(breaks=seq(0,0.1,0.02),limits = c(0,0.1))
  #scale_colour_manual(values=cols)+
  #scale_x_continuous(breaks =res$Correlation,labels =rep(seq(2,6,1),2))+
  ggplot2::theme(legend.position = "top",
                 plot.title = element_text(hjust = 0.5),
                 legend.title = ggplot2::element_blank(),
                 legend.key = ggplot2::element_blank(),
                 panel.grid.major.y = ggplot2::element_line(color = "#DDDDDD"), 
                 panel.grid.major.x = ggplot2::element_blank(),
                 panel.background = ggplot2::element_blank(),
                 #axis.title.y = ggplot2::element_blank(),
                 panel.border = ggplot2::element_rect(fill = NA))+
  theme(text = element_text(size =12),
        strip.text = element_text(size = 12))
plot1








