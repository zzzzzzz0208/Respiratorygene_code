
###===Figure 3===###

rm (list=ls())
gc()

#====loading data====####

library(survival)
library(survminer)

####====Fig. 3A====####

dir = "D:/TCGA-53gene_expression+survival" 
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)    
for(f in 1:length(file_list)){
  
  sur<-read.table(file_list[f],header=T,row.names=1,sep=",")
  sur_data<-subset(sur, select=c(OS,OS.time,X_PATIENT)) 
  exp<-sur[1:(ncol(sur)-3)] 
  
  cox_result<-matrix(c("gene","p","HR","Low.95.CI","High.95.CI"),1,5)
  for (i in 1:ncol(exp)) {
    df<-cbind(sur_data,exp[,i])
    colnames(df)[ncol(df)]<-"exp"
    res.cox<-coxph(Surv(OS.time,OS)~exp, data =df )
    a<-summary(res.cox)
    b<-exp(confint(res.cox))
    b<- as.vector(b)
    
    if(a$coefficients[5]<0.05){
    xg=c(colnames(exp)[i],a$coefficients[5],a$coefficients[2],b[1],b[2])
    cox_result<-rbind(cox_result,xg)
    }
  }
  
  colnames(cox_result)<-cox_result[1,]
  cox_result<-cox_result[drop=F,-1,]
  
  setwd("D:/单因素cox_result")
  cname1<-strsplit(file_list[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split=' \\.')[[1]][1]
  name<-paste(cname2,"_cox.csv",sep='')
  write.csv(cox_result,name,row.names=F)
}

##整理画图数据
geneset<-read.table("D:/Code/Data/RespiratorySig.txt",header=F,sep="\t")
geneset<-as.matrix(geneset)

dir = "D:/单因素cox_result"  
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)     

cox <- read.table(file_list[1],header = T,row.names=1,sep=",")

z<-matrix(nr = 53, nc = 4)
rownames(z)=geneset[,1]  
colnames(z)<-colnames(cox)
cox<-as.matrix(cox)

for (i in 1:nrow(z)) {
  a<-which(rownames(cox)==rownames(z)[i])

  if(length(a)!=0) {
    z[i,]=cox[which(rownames(cox)==rownames(z)[i]),]
  } 
}         
merge_data<-as.data.frame(z[,2])

cname1<-strsplit(file_list[1],split="/")[[1]][2]
cname2<-strsplit(cname1,split=' _')[[1]][1]
colnames(merge_data)<-c(cname2)

for(f in 2:length(file_list)){
  cox <- read.table(file_list[f],header = T,row.names=1,sep=",")
  
  z<-matrix(nr = 53, nc = 4)
  rownames(z)=geneset[,1]  
  colnames(z)<-colnames(cox)
  cox<-as.matrix(cox)
  
  for (i in 1:nrow(z)) {
    a<-which(rownames(cox)==rownames(z)[i])
    
    if(length(a)!=0) {
      z[i,]=cox[which(rownames(cox)==rownames(z)[i]),]
    } 
  }         
  data<-as.data.frame(z[,2])
  cname1<-strsplit(file_list[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split=' _')[[1]][1]
  colnames(data)<-c(cname2)
  
  merge_data<-cbind(merge_data,data)
}
setwd("D:/单因素cox_result/画图数据")
write.csv(merge_data,"53genes-HR_value.csv")

##画图
dat<-as.matrix(merge_data)
dat[dat >= 1]= 1
dat[dat>0 & dat<1] =(-1)
dat[is.na(dat)]<-0

GLY<-dat[c(1:32),]
TCA<-dat[c(33:41),]
OXPHOS<-dat[c(42:53),]

library(pheatmap)
Fig.3A<-pheatmap(t(GLY),
                 cellwidth = 5, cellheight = 10, 
                 legend = F, 
                 border_color="#FFFFFF",
                 color = colorRampPalette(c("#FFCDD2", "#FFFFFF", "#D32F2F"))(50),
                 cluster_rows = F,cluster_cols= F,
                 fontsize_row = 10,fontsize_col=7 ,angle_col=90
) 
Fig.3A<-pheatmap(t(TCA),
                 cellwidth = 5, cellheight = 10, 
                 legend = F, 
                 border_color="#FFFFFF",
                 color = colorRampPalette(c("#FFECB3", "#FFFFFF", "#FFB300"))(50),
                 cluster_rows = F,cluster_cols= F,
                 fontsize_row = 10,fontsize_col=7 ,angle_col=90
) 
Fig.3A<-pheatmap(t(OXPHOS),
                 cellwidth = 5, cellheight = 10, 
                 legend = F, 
                 border_color="#FFFFFF",
                 color = colorRampPalette(c("#C5CAE9", "#FFFFFF", "#5C6BC0"))(50),
                 cluster_rows = F,cluster_cols= F,
                 fontsize_row = 10,fontsize_col=7 ,angle_col=90
) 
Fig.3A

####====Fig. 3B====####

vacn_HR<-read.table("D:/Code/Data/Fig3/TCGA_pancancer-VCAN-HR_value.csv",header=T,sep=",")
colnames(vacn_HR)[1]<-"cancer"

Fig.3B <- ggplot(data=vacn_HR,aes(x=HR,y=cancer,color=p))+
  geom_errorbarh(aes(xmax=High.95.CI,xmin=Low.95.CI),color="black",height=0.3,size=0.8)+
  geom_point(aes(x=HR,y=cancer),size=4,shape=18)+
  geom_vline(xintercept=1,linetype="dashed",size=0.5)+
  scale_x_continuous(breaks=c(0.5,1,1.5))+
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
Fig.3B

####====Fig. 3C====####

load('D:/Code/Data/Fig3/kidney-53gene_expression.RData')
load('D:/Code/Data/Fig3/kidney_survival.RData')

##一致性聚类
texp<-t(exp)
z_normal = sweep(texp,1, apply(texp,1,median,na.rm=T))
z_normal<-as.matrix(z_normal)

library(ConsensusClusterPlus)
title="D:/肾癌样本一致性聚类" 
#title=tempdir()
results = ConsensusClusterPlus(z_normal,maxK=3,reps=50,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="hc",distance="pearson",
                               seed=1262118388.71279,plot="pdf") #Fig.3C

cluster<-as.data.frame(results[[2]]$consensusClass)
colnames(cluster)[1]<-c("group")
cluster$group[which(cluster$group ==2)] <- "High-Expression"
cluster$group[which(cluster$group ==1)] <- "Low-Expression"

setwd("D:/肾癌样本一致性聚类")
write.table(cluster,"一致性聚类分组.txt",row.names=T,col.names=T,sep="\t",quote=F)
dat<-cbind(exp,cluster)
write.csv(dat,"expression-consensus_group.csv")

datall<-cbind(dat,sur_dat)
write.csv(datall,"expression-consensus_group-survival.csv")

##survival analysis
fit <- survfit(Surv(OS.time,OS)~group,data = datall)
Fig.3D <- ggsurvplot(fit, data = datall, 
                     risk.table = F,
                     palette = c("#A93226","#2874A6"),
                     xlab="Time(days)",pval =T,
                     ggtheme =theme_light())
Fig.3D

##log rank test
surv_diff <- survdiff(Surv(OS.time,OS) ~ group, data = datall)
p = 1 - pchisq(surv_diff$chisq, length(surv_diff$n) - 1)


