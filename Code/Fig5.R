
###===Figure 5===###

rm (list=ls())
gc()

#====loading data====####

library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)

####====Fig. 5B====####

drug_target <-  read.table("D:/Code/Data/Fig5/药靶基因.txt",header=F,sep = "\t")

dir1 = "D:/TCGA-glycolysis_gene_expression"
file_list1 = list.files(path =dir1, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
dir2 = "D:/TCGA-oxhpos_gene_expression"
file_list2 = list.files(path =dir2, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
dir3 = "D:/TCGA-drug_target_expression"
file_list3 = list.files(path =dir3, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)

for(f in 1:length(file_list1)){
  
  ##glycolysis
  gly_res<-read.table(file_list1[f],header = T,row.names=1,sep = ",")
 
  rowmean=apply(gly_res,1,mean)                   
  rowsd=apply(gly_res,1,sd)                               
  res1=sweep(gly_res,1,rowmean)              
  res2=sweep(res1,1,rowsd,'/') 

  gly_colmean<-as.data.frame(colMeans(res2))
  colnames(gly_colmean)[1]<-"score"
  
  ##oxphos
  oxp_res<-read.table(file_list2[f],header = T,row.names=1,sep = ",")

  rowmean=apply(oxp_res,1,mean)                   
  rowsd=apply(oxp_res,1,sd)                               
  res3=sweep(oxp_res,1,rowmean)              
  res4=sweep(res3,1,rowsd,'/') 

  oxp_colmean<-as.data.frame(colMeans(res4))
  colnames(oxp_colmean)[1]<-"score"
  
  setwd("D:/样本得分/GLYCOLYSIS")
  cname<-strsplit(file_list1[f],split="\\/")[[1]][2]
  write.csv(gly_colmean,cname)
  setwd("D:/样本得分/OXPHOS")
  write.csv(oxp_colmean,cname)
  
  ##计算样本得分与药靶基因表达相关性
  drug_exp <-read.table(file_list3[f],header = T,row.names=1,sep = ",")
  del <- which(apply(drug_exp==0, 1, sum) > 0.3*ncol(drug_exp))
  drug_exp <- drug_exp[-del,]
  
  #glycolysis
  xiangguan1<-matrix(c("drug_target","pearson","p-value"),1,3)
  for(j in 1:nrow(drug_exp)){
    cor<-cor.test(as.numeric(gly_colmean[,1]),as.numeric(drug_exp[j,]),method = "pearson")
    if(cor$p.value<0.05&cor$estimate>0.3){
      xg=c(rownames(drug_exp)[j],cor$estimate,cor$p.value)
      xiangguan1<-rbind(xiangguan1,xg)
    }
  }
  colnames(xiangguan1)<-c("drug_target","cor","p-value")
  xiangguan1<-xiangguan1[-1,]
  
  #oxphos
  xiangguan2<-matrix(c("drug_target","pearson","p-value"),1,3)
  for(j in 1:nrow(drug_exp)){
    cor<-cor.test(as.numeric(oxp_colmean[,1]),as.numeric(drug_exp[j,]),method = "pearson")
    if(cor$p.value<0.05&cor$estimate>0.3){
      xg=c(rownames(drug_exp)[j],cor$estimate,cor$p.value)
      xiangguan2<-rbind(xiangguan2,xg)
    }
  }
  colnames(xiangguan2)<-c("drug_target","cor","p-value")
  xiangguan2<-xiangguan2[-1,]
  
  setwd("D:/respiratory_signature-drug_target-correlation/glycolysis")
  write.csv(xiangguan1,cname,row.names=F)
  setwd("D:/respiratory_signature-drug_target-correlation/oxphos")
  write.csv(xiangguan2,cname,row.names=F)
}

signature<-read.table("D:/Code/Data/Fig5/kidney-respiratory_gene_signature.txt",header=F,sep="\t")

library(VennDiagram)
Fig.5B <- venn.diagram(list(target_gene=drug_target[,1], respiratory_gene_signature=signature[,1]), 
                       filename = NULL,
                       fill=c("#fbeee6","#cd6155"), 
                       alpha=c(0.5,0.5),cex=2, cat.fontface=4, cat.cex=1.4)
pdf(file="venn.pdf")
grid.draw(Fig.5B)
dev.off()

####====Fig. 5C====####



####====Fig. 5D====####

dat<-read.table("D:/KICH-VCAN_expression.csv",header=T,row.names=1,sep=",")

library(ggplot2)
Fig.5D <- ggplot(dat, aes(x=group, y=value,fill=group)) + 
  geom_violin(trim=F,color="white") + 
  geom_boxplot(width=0.2,position=position_dodge(0.8))+  
  scale_fill_manual(values = c("#F1948A", "#85C1E9"))+
  theme_bw()+ 
  theme(axis.text.x=element_text(colour="black",size=12),
        axis.text.y=element_text(size=12,face="plain"), 
        axis.title.y=element_text(size = 12,face="plain"), 
        panel.border = element_blank(),axis.line = element_line(colour = "black"), 
        legend.text=element_text(face="italic",colour="black", 
                                 size=12),
        legend.title=element_text(face="italic",colour="black", 
                                  size=12),
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank())+ 
  ylab("KICH")+xlab("") 
Fig.5D

####====Fig. 5E====####

fit<-survfit(Surv(OS.time,OS)~group,data= VCAN_survival)
Fig.5E <- ggsurvplot(fit, data = VCAN_survival, 
                     risk.table = F,
                     palette = c("#A93226","#2874A6"),
                     xlab="Time(days)",pval =T,
                     ggtheme =theme_light())
Fig.5E

####====Fig. 5F====####

dat<-read.table('D:/Code/Data/Fig5/KICH-VCAN_exp-glycolysis_score.csv',header=T,row.names=1,sep='')

Fig.5F <- ggplot(dat, aes(x = VCAN, y = gly_colmean))+
  theme(axis.title.x = element_text(size = 12),axis.title.y = element_text(size = 12))+
  geom_point(shape = 20,size=2,color="#5DADE2")+ 
  geom_smooth(method="lm", color = "#F1948A", fill = "lightgray",formula = y~x)+ 
  stat_cor(size=5)+   
  theme_gray() 
Fig.5F













