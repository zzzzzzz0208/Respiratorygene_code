
###===Figure S3===###

rm (list=ls())
gc()

#====loading data====####

library(psych)
library(ggplot2)
library(dplyr)
library(forcats)

####====Fig. S3A====####

fit<-survfit(Surv(OS.time,OS)~group,data = MAOB_survival)
Fig.S3A <- ggsurvplot(fit,
                      pval = TRUE, conf.int = TRUE,
                      risk.table = FALSE, 
                      risk.table.col = "strata", 
                      linetype = "strata", 
                      surv.median.line = "hv", 
                      ggtheme = theme_bw(),  
                      palette = c("#f5b7b1", "#a2d9ce"))
Fig.S3A

####====Fig. S3B====####

dat<-read.table("D:/Code/Data/FigS3/GEO-VCAN-HR_value.csv",header=T,sep=",")
colnames(dat)[1]<-'cancer'

sample[,3]<-as.numeric(sample[,3])
sample[,4]<-as.numeric(sample[,4])
sample[,5]<-as.numeric(sample[,5])
sample[,6]<-as.numeric(sample[,6])

##定义因子顺序
group<-read.table("D:/Code/Data/FigS3/因子顺序txt",header=F,sep="\t")
group<-as.matrix(group)
sample$cancer<- factor(sample$cancer,levels=group) 

ggplot(data=sample,aes(x=HR,y=cancer,color=p))+
  geom_errorbarh(aes(xmax=High.95.CI,xmin=Low.95.CI),color="black",height=0.3,size=0.8)+
  geom_point(aes(x=HR,y=cancer),size=5.5,shape=18)+
  geom_vline(xintercept=1,linetype="dashed",size=0.5)+
  scale_x_continuous(breaks=c(1,2))+
  scale_y_discrete(limit=c(as.character(sample$cancer)))+
  scale_color_continuous(high="#FFEBEE",low="#B71C1C")+  
  ylab("")+xlab("Hazard ratios")+
  labs(color="P value",title="")+
  theme(panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(),axis.line = element_line(colour = "grey"))+  
  theme(axis.text=element_text(colour="black",size=12))+  
  theme(legend.text=element_text(size=12))+  
  theme(legend.title=element_text(size=10))+  
  theme(axis.title=element_text(colour="black",size=12))


