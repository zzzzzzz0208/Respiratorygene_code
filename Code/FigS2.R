
###===Figure S2===###

rm (list=ls())
gc()

#====loading data====####

library(psych)
library(ggplot2)
library(dplyr)
library(forcats)

####====Fig. S2A====####

gly<-read.table("D:/Code/Data/glycolysis.txt",header=F,sep="\t")
tca<-read.table("D:/Code/Data/TCA cycle.txt",header=F,sep="\t")
oxphos<-read.table("D:/Code/Data/oxphos.txt",header=F,sep="\t")

load('D:/Code/Data/53gene-TCGA_pancancer_expression.RData')
texp<-t(exp)

##两两计算相关性
texp[texp==0] = NA
newexp=na.omit(texp)

gly<-newexp[,c(1:32)]
tca<-newexp[,c(33:41)]
oxp<-newexp[,c(42:53)]

#TCA-OXPHOS
xiangguan1<-matrix(c("tca","oxp","pearson","p-value"),1,4)
for(i in 1:ncol(tca)){
  for(j in 1:ncol(oxp)){
    cor<-cor.test(tca[,i],oxp[,j],method = "pearson")
    xg=c(colnames(tca)[i],colnames(oxp)[j],cor$estimate,cor$p.value)
    xiangguan1<-rbind(xiangguan1,xg)
  }
}
colnames(xiangguan1)<-c("tca","oxp","cor","p-value")
xiangguan1<-as.data.frame(xiangguan1[-1,])

xiangguan1$`p-value`<-as.numeric(xiangguan1$`p-value`)
xiangguan1$cor<-as.numeric(xiangguan1$cor)
xiangguan1$cor<-round(xiangguan1$cor,4)

xiangguan1$pstar <- ifelse(xiangguan1$`p-value` < 0.05,
                           ifelse(xiangguan1$`p-value` < 0.01,
                                  ifelse(xiangguan1$`p-value` < 0.001,"***","**"),"*"),"")

#GLYCOLYSIS-TCA
xiangguan2<-matrix(c("gly","tca","pearson","p-value"),1,4)
for(i in 1:ncol(gly)){
  for(j in 1:ncol(tca)){
    cor<-cor.test(gly[,i],tca[,j],method = "pearson")
    xg=c(colnames(gly)[i],colnames(tca)[j],cor$estimate,cor$p.value)
    xiangguan2<-rbind(xiangguan2,xg)
  }
}
colnames(xiangguan2)<-c("gly","tca","cor","p-value")
xiangguan2<-as.data.frame(xiangguan2[-1,])

xiangguan2$`p-value`<-as.numeric(xiangguan2$`p-value`)
xiangguan2$cor<-as.numeric(xiangguan2$cor)
xiangguan2$cor<-round(xiangguan2$cor,4)

xiangguan2$pstar <- ifelse(xiangguan2$`p-value` < 0.05,
                           ifelse(xiangguan2$`p-value` < 0.01,
                                  ifelse(xiangguan2$`p-value` < 0.001,"***","**"),"*"),"")

#GLYCOLYSIS-OXPHOS
xiangguan3<-matrix(c("gly","oxp","pearson","p-value"),1,4)
for(i in 1:ncol(gly)){
  for(j in 1:ncol(oxp)){
    cor<-cor.test(gly[,i],oxp[,j],method = "pearson")
    xg=c(colnames(gly)[i],colnames(oxp)[j],cor$estimate,cor$p.value)
    xiangguan3<-rbind(xiangguan3,xg)
  }
}
colnames(xiangguan3)<-c("gly","oxp","cor","p-value")
xiangguan3<-as.data.frame(xiangguan3[-1,])

xiangguan3$`p-value`<-as.numeric(xiangguan3$`p-value`)
xiangguan3$cor<-as.numeric(xiangguan3$cor)
xiangguan3$cor<-round(xiangguan3$cor,4)

xiangguan3$pstar <- ifelse(xiangguan3$`p-value` < 0.05,
                           ifelse(xiangguan3$`p-value` < 0.01,
                                  ifelse(xiangguan3$`p-value` < 0.001,"***","**"),"*"),"")

Fig.S2A <- ggplot(xiangguan2, aes(gly,tca)) + 
  geom_tile(aes(fill = cor),colour = "white",size=1)+
  scale_fill_gradient2(low = "#2874C5",mid = "white",high = "#EE0000")+
  geom_text(aes(label=pstar),col ="#424949",size = 4)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(size = 12))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))
Fig.S2A

####====Fig. S2B====####

Fig.S2B <- ggplot(xiangguan1, aes(tca,oxp)) + 
  geom_tile(aes(fill = cor),colour = "white",size=1)+
  scale_fill_gradient2(low = "#2874C5",mid = "white",high = "#EE0000")+
  geom_text(aes(label=pstar),col ="#424949",size = 4)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1,size = 12),
        axis.text.y = element_text(size = 12))+
  labs(fill =paste0(" * p < 0.05","\n\n","** p < 0.01","\n\n","*** p < 0.001","\n\n","Correlation"))
Fig.S2B

