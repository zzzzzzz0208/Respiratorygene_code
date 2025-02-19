
###===Figure 6===###

rm (list=ls())
gc()

#====loading data====####

library(data.table)
library(dplyr)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)

####====Fig. 6A====####

dir1 = "D:/TCGA-glycolysis_gene_expression"
file_list1 = list.files(path =dir1, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
dir2 = "D:/TCGA-oxphos_gene_expression"
file_list2 = list.files(path =dir2, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
dir3 = "D:/TCGA-gene_expression"
file_list3 = list.files(path =dir3, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)

for(f in 1:length(file_list1)){
  ##glycolysis
  gly_res<-read.table(file_list1[f],header = T,row.names=1,sep='')
  
  rowmean=apply(gly_res,1,mean) 
  rowsd=apply(gly_res,1,sd)      
  gly_res1=sweep(gly_res,1,rowmean) 
  gly_res2=sweep(gly_res1,1,rowsd,'/') 
  
  gly_colmean<-as.data.frame(colMeans(gly_res2))
  colnames(gly_colmean)[1]<-"score"
  
  exp<-read.table(file_list3[f],header = T,row.names=1,sep = ",")
  del <- which(apply(exp==0, 1, sum) > 0.3*ncol(exp))
  exp <- exp[-del,] 
  
  #计算每个样本得分与每种癌症类型中基因表达之间的相关性
  xiangguan1<-matrix(c("exp","pearson","p-value"),1,3)
  for(j in 1:nrow(exp)){
    cor<-cor.test(as.numeric(gly_colmean[,1]),as.numeric(exp[j,]),method = "pearson")
    if(cor$p.value<0.05){
      xg=c(rownames(exp)[j],cor$estimate,cor$p.value)
      xiangguan1<-rbind(xiangguan1,xg)
    }
  }
  colnames(xiangguan1)<-c("gene","cor","p-value")
  xiangguan1<-xiangguan1[-1,]
  
  setwd("D:/药物/GLYCOLYSIS/相关性结果(p0.05)")
  write.csv(xiangguan1,name,row.names=F)
  
  ##oxphos
  oxp_res<-read.table(file_list2[f],header = T,row.names=1,sep='')
  
  rowmean=apply(oxp_res,1,mean) 
  rowsd=apply(oxp_res,1,sd)      
  oxp_res1=sweep(oxp_res,1,rowmean) 
  oxp_res2=sweep(oxp_res1,1,rowsd,'/') 
  
  oxp_colmean<-as.data.frame(colMeans(oxp_res2))
  colnames(oxp_colmean)[1]<-"score"
  
  #计算每个样本得分与每种癌症类型中基因表达之间的相关性
  xiangguan2<-matrix(c("exp","pearson","p-value"),1,3)
  for(j in 1:nrow(exp)){
    cor<-cor.test(as.numeric(oxp_colmean[,1]),as.numeric(exp[j,]),method = "pearson")
    if(cor$p.value<0.05){
      xg=c(rownames(exp)[j],cor$estimate,cor$p.value)
      xiangguan2<-rbind(xiangguan2,xg)
    }
  }
  colnames(xiangguan2)<-c("gene","cor","p-value")
  xiangguan2<-xiangguan2[-1,]
  
  setwd("D:/药物/OXPHOS/相关性结果(p0.05)")
  write.csv(xiangguan2,name,row.names=F)
}

#提取正、负相关top150基因输入CMap

##
dir = "D:/药物/CMap_result/OXPHOS"
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)

for(f in 1:length(file_list)){
  
  dat<-read.table(file_list[f],header = T,quote="",sep = "\t")
  
  ES<-as.numeric(dat[,2])
  test1<-which(ES>=95)
  up<-as.data.frame(dat[test1,])
  test2<-which(ES<=(-95))
  down<-as.data.frame(dat[test2,])
  
  up<-up[,-c(1,3)]
  setwd("D:/药物/OXPHOS/ES大于95")
  cname<-strsplit(file_list[f],split="\\/")[[1]][2]
  write.table(up,cname,row.names=F,sep="\t",quote=F)
  
  down<-down[,-c(1,3)]
  setwd("D:/药物/OXPHOS/ES小于-95")
  write.table(down,cname,row.names=F,sep="\t",quote=F)
}

##统计药物在多少癌症中相关
###pos###
dir = "D:/药物/OXPHOS/ES大于95" 
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)

merge_ES <- read.table(file_list[1],header = T,quote="",sep="\t")
for(f in 2:length(file_list)){
  
  ES <- read.table(file_list[f],header = T,quote="",sep="\t")
  merge_ES<-rbind(merge_ES,ES)
}
freq<-as.data.frame(table(merge_ES[,2]))
colnames(freq)[1]<-"Drug"
freq<-as.data.frame(freq[order(freq[,2],decreasing = T),]) 
write.table(freq,"D:/药物/OXPHOS/ES大于95/每种药物在多少种癌症中相关.txt",row.names=F,sep="\t",quote=F)

drug_pos<-as.data.frame(freq[freq$Freq>=10,1]) #选择至少在10cancer中相关的药物画图
colnames(drug_pos)<-"ID"

##整理画图数据
geneset_pos<-as.matrix(drug_pos)

dir = "D:/药物/OXPHOS/ES大于95"
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)
ES <- read.table(file_list[1],header = T,quote="",sep="\t")
merge_pos <- ES[ES$ID %in% geneset_pos[,1],]

cname1<-strsplit(file_list[1],split="\\/")[[1]][2]
cname2<-strsplit(cname1,split="\\.")[[1]][1]
merge_pos$cancer<-cname2

for(f in 2:length(file_list)){
  ES <- read.table(file_list[f],header = T,quote="",sep="\t")
  pos <- ES[ES$ID %in% geneset_pos[,1],]
  cname1<-strsplit(file_list[f],split="\\/")[[1]][2]
  cname2<-strsplit(cname1,split="\\.")[[1]][1]
  
  if(nrow(pos)!=0){
    pos$cancer<-cname2 
  }
  merge_pos<-rbind(merge_pos,pos)
}
setwd("D:/药物/OXPHOS/ES大于95")
write.csv(merge_pos,"画图数据(10cancer).csv",row.names=F)

###neg###
dir = "D:/药物/OXPHOS/ES小于-95" 
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)

merge_ES <- read.table(file_list[1],header = T,quote="",sep="\t")
for(f in 2:length(file_list)){
  
  ES <- read.table(file_list[f],header = T,quote="",sep="\t")
  merge_ES<-rbind(merge_ES,ES)
}
freq<-as.data.frame(table(merge_ES[,2]))
colnames(freq)[1]<-"Drug"
freq<-as.data.frame(freq[order(freq[,2],decreasing = T),]) 
write.table(freq,"D:/药物/OXPHOS/ES小于-95/每种药物在多少种癌症中相关.txt",row.names=F,sep="\t",quote=F)

drug_neg<-as.data.frame(freq[freq$Freq>=10,1]) #选择至少在10cancer中相关的药物画图
colnames(drug_neg)<-"ID"

##整理画图数据
geneset_neg<-as.matrix(drug_neg)

dir = "D:/药物/OXPHOS/ES小于-95"
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)
ES <- read.table(file_list[1],header = T,quote="",sep="\t")
merge_neg <- ES[ES$ID %in% geneset_neg[,1],]

cname1<-strsplit(file_list[1],split="\\/")[[1]][2]
cname2<-strsplit(cname1,split="\\.")[[1]][1]
merge_neg$cancer<-cname2

for(f in 2:length(file_list)){
  ES <- read.table(file_list[f],header = T,quote="",sep="\t")
  neg <- ES[ES$ID %in% geneset_neg[,1],]
  cname1<-strsplit(file_list[f],split="\\/")[[1]][2]
  cname2<-strsplit(cname1,split="\\.")[[1]][1]
  
  if(nrow(neg)!=0){
    neg$cancer<-cname2 
  }
  merge_neg<-rbind(merge_neg,neg)
}
setwd("D:/药物/OXPHOS/ES小于-95")
write.csv(merge_neg,"画图数据(10cancer).csv",row.names=F)

##画图
cishu_pos<-as.data.frame(table(merge_pos$Name))
cishu_pos<-as.data.frame(cishu_pos[order(cishu_pos[,2]),]) 
colnames(cishu_pos)[1]<-"Drug"

cishu_neg<-as.data.frame(table(merge_neg$Name))
cishu_neg<-as.data.frame(cishu_neg[order(cishu_neg[,2]),]) 
colnames(cishu_neg)[1]<-"Drug"

library(ggplot2)
group_pos<-as.matrix(cishu_pos[,1]) 
merge_pos$Name <- factor(merge_pos$Name,levels=group_pos) 
group_neg<-as.matrix(cishu_neg[,1]) 
merge_neg$Name <- factor(merge_neg$Name,levels=group_neg) 

Fig.6A <- ggplot(merge_pos,aes(cancer,Name,color=abs(Score)))+
  geom_point(size=2.5)+   #shape=21,colour = "white",stroke = 1.5
  xlab(NULL)+ ylab(NULL)+
  theme(panel.grid.major=element_line(colour='white'))+ 
  theme(axis.text.x=element_text(color = "#424242", size=12, angle=45, hjust=1))+ 
  theme(axis.text.y=element_text(color = "#424242", size=12))+
  theme(axis.ticks.x  = element_blank(),axis.ticks.y  = element_blank())+ 
  scale_colour_gradientn(colors =c("#FDEDEC","#F5B7B1","#E57373"),guide="none")+ #ES>95
  scale_y_discrete(position = "left")

Fig.6A <- ggplot(merge_neg,aes(cancer,Name,color=abs(Score)))+
  geom_point(size=2.5)+   #shape=21,colour = "white",stroke = 1.5
  xlab(NULL)+ ylab(NULL)+
  theme(panel.grid.major=element_line(colour='white'))+ 
  theme(axis.text.x=element_text(color = "#424242", size=12, angle=45, hjust=1))+ 
  theme(axis.text.y=element_text(color = "#424242", size=12))+
  theme(axis.ticks.x  = element_blank(),axis.ticks.y  = element_blank())+ 
  scale_colour_gradientn(colors=c("#F5EEF8","#D1C4E9","#C39BD3"),guide="none")+ #ES<-95 
  scale_y_discrete(position = "left") 
Fig.6A


