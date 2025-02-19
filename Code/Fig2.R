
###===Figure 2===###

rm (list=ls())
gc()

#====loading data====####

library(GSVA)
library(psych)
library(data.table)
library(ggplot2)
library(dplyr)
library(forcats)

####====Fig. 2A====####

geneset<-read.table("D:/Code/Data/Fig2/hallmark_pathway.txt",header = T,sep = "\t")
sig<-read.table("D:/Code/Data/RespiratorySig.txt",header = F,sep="\t")  

dir = "D:/TCGA-gene_expression"  
file_list = list.files(path = dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)      

for(f in 1:length(file_list)){
  
  res <- read.table(file_list[f],header = T,row.names=1,sep=",")
  del <- which(apply(res==0, 1, sum) > 0.3*ncol(res)) 
  res2 <- res[-del,] 
  
  t = as.matrix(res2)  
  gs = as.list(geneset1)  
  gs = lapply(gs, function(x) x[!is.na(x)])
  score = gsva(t, gs, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
  
  setwd("D:/hallmark_pathway-ssgsea_score")
  cname1<-strsplit(file_list1[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split=" _")[[1]][1]
  name1<-paste(cname2,".csv",sep='')
  write.csv(score,name1)
  
  exp<-res[rownames(res) %in% sig[,1],]
  texp<-t(exp)     
  tscore<-t(score)  
  xiangguan<-matrix(c("respiratory_sig","ssgsea_score","pearson","p-value"),1,4)
  for(i in 1:ncol(texp)){    
    for(j in 1:50){        
      cor<-cor.test(texp[,i],tscore[,j],method = "pearson")
      xg=c(colnames(texp)[i],colnames(tscore)[j],cor$estimate,cor$p.value)
      xiangguan<-rbind(xiangguan,xg)
      }
  }
  colnames(xiangguan)<-c("gene","pathway","cor","p-value")
  xiangguan<-xiangguan[-1,]
  setwd("D:/53gene-hallmark_pathway相关性")
  name2<-paste(cname2,".txt",sep='')
  write.table(xiangguan,name2,row.name=F,col.name=T,sep="\t",quote=F)
  
  cor<-as.numeric(xiangguan[,3])
  pvalue<-as.numeric(xiangguan[,4])
  
  test1<-which(cor >= 0.6 & pvalue < 0.01)
  pos<-xiangguan[test1,]
  test2<-which(cor <= (-0.6) & pvalue < 0.01)
  neg<-xiangguan[test2,]
  
  setwd("D:/53gene-hallmark_pathway相关性/pos")
  name3<-paste(cname2,"_pos.txt")
  write.table(pos,name3,sep="\t",row.names = F,quote=F)
  setwd("D:/53gene-hallmark_pathway相关性/neg")
  name4<-paste(cname2,"_neg.txt")
  write.table(neg,name4,sep="\t",row.names = F,quote=F)
}

####====Fig. 2B====####

#统计与每个基因正相关的通路数量
dir = "D:/53gene-hallmark_pathway相关性/pos" 
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)       

merge_data = read.table(file_list[1],header=T,sep="\t") 
for(i in 2:length(file_list)){
  
  data<-read.table(file_list[i],header=T,sep="\t")
  merge_data = rbind(merge_data,data)
}
data<-unique(setDT(merge_data), by = c("gene","pathway"))
pos_freq<-as.data.frame(table(data[,1]))
colnames(pos_freq)<-c("gene","number")
setwd("D:/53gene-hallmark_pathway相关性/画图数据")
write.table(pos_freq,"与每个gene正相关的通路数量.txt",sep="\t",row.names = F,quote=F)

#统计与每个基因负相关的通路数量
dir = "D:/53gene-hallmark_pathway相关性/neg" 
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)       

merge_data = read.table(file_list[1],header=T,sep="\t") 
for(i in 2:length(file_list)){
  
  data<-read.table(file_list[i],header=T,sep="\t")
  merge_data = rbind(merge_data,data)
}
data<-unique(setDT(merge_data), by = c("gene","pathway"))
neg_freq<-as.data.frame(table(data[,1]))
colnames(neg_freq)<-c("gene","number")
setwd("D:/53gene-hallmark_pathway相关性/画图数据")
write.table(neg_freq,"与每个gene负相关的通路数量.txt",sep="\t",row.names = F,quote=F)

#整理画条形图数据
geneset<-read.table("C:/Code/Data/RespiratorySig.txt",header = F,sep="\t")
#pos
rownames(pos_freq)=pos_freq[,1] 
pos_freq=pos_freq[,-1,drop=FALSE] 
pos<-matrix(nr = 53, nc = 1)
rownames(pos)=geneset[,1]  
for (i in 1:nrow(pos)) {
  a<-which(rownames(pos_freq)==rownames(pos)[i])

  if(length(a)!=0) {
    pos[i,1]=pos_freq[which(rownames(pos_freq)==rownames(pos)[i]),1]
  } 
}     
pos[is.na(pos)] <- 0
colnames(pos)<-c("number")
pos<-as.data.frame(pos)
pos$type<-c("Pos") 
pos <- tibble::rownames_to_column(pos, "gene")

#neg
rownames(neg_freq)=neg_freq[,1]
neg_freq=neg_freq[,-1,drop=FALSE]          
neg<-matrix(nr = 53, nc = 1)
rownames(neg)=geneset[,1]  
for (i in 1:nrow(neg)) {
  a<-which(rownames(neg_freq)==rownames(neg)[i])
  if(length(a)!=0) {
    neg[i,1]=neg_freq[which(rownames(neg_freq)==rownames(neg)[i]),1]
  } 
}     
neg[is.na(neg)] <- 0
colnames(neg)<-c("number")
neg<-as.data.frame(neg)
neg$type<-c("Neg")
neg <- tibble::rownames_to_column(neg, "gene")

number<-rbind(pos,neg)

setwd("D:\53gene-hallmark_pathway相关性")
write.table(number, "53gene相关通路数量.txt",sep="\t",row.names = F,quote=F) 

##画图
number[which(number$type=='Neg'),c('number')] <- number[which(number$type=='Neg'),c('number')]*-1
geneset<-read.table("D:/Code/Data/RespiratorySig.txt",header=F,sep="\t")
geneset<-as.matrix(geneset) 
number$gene <- factor(number$gene,levels=geneset)

Fig.2B <- ggplot(data=number, aes(x=gene, y=number,fill=type),fill=group) +  
  geom_col(position = position_dodge(0), width = 1.4) +   
  labs(x = NULL, y = 'Number of pathways') +
  scale_fill_manual(values = c("#2471a3","#c0392b"))+ 
  scale_y_continuous(breaks = seq(-40, 40, 10), 
                     labels = as.character(abs(seq(-40, 40, 10))), 
                     limits = c(-40, 40)) +
  theme_bw()+
  theme(legend.position = "none")+  
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5,vjust = 0.5, size=8,color="#212121"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme(panel.grid=element_blank())+   
  theme(panel.border = element_blank())+  
  theme(axis.ticks.x = element_blank()) 
Fig.2B

####====Fig. 2C====####

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

##画图
Fig.2C <- ggplot(xiangguan3, aes(gly,oxp)) + 
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
Fig.2C


