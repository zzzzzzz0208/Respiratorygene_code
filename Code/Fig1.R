
###===Figure 1===###

rm (list=ls())
gc()

#====loading data====####

library(pheatmap)
library(ggplot2)
library(dplyr)
library(tidyr)
library(limma)
library(VennDiagram)

####====Fig. 1A====####

dat<-read.table("D:/Code/Data/Fig1/venn-geneset.txt",header=F,sep="\t")
colnames(dat)<-c("glycolysis","TCA","oxphos")

Fig.1A <-venn.diagram(list(Glycolysis=dat[,1], TCA=dat[,2],OXPHOS=dat[,3]), 
                        filename = NULL,
                        fill=c("#BBDEFB","#edbb99","#c39bd3"), 
                        alpha=c(0.5,0.5,0.5), cex=2, cat.fontface=4, cat.cex=1.4)
pdf(file="venn.pdf")
grid.draw(Fig.1A)
dev.off()

####====Fig. 1B====####

dir = "D:/TCGA-53gene_expression"
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)

for(f in 1:length(file_list)){
  
  exp<-read.table(file_list[f],header=T,row.names=1,sep=",")
  
  transfer <- function(x){
    x <- as.matrix(x)
    a <- sort(x)
    b <- min(a[which(a!=0)])
    x[which(x==0)] <- b*0.01
    return(x)
  }
  exp <- transfer(exp)
  
  can <- exp[,grep("\\.0[1-9][A-Z]",colnames(exp))]
  nor <- exp[,grep("\\.11[A-Z]",colnames(exp))]
  exper <- cbind(can,nor)
  exper<-as.matrix(exper)
  
  #t test
  out_re<-c()
  for(i in 1:dim(exper)[1]){
    t_test<-t.test(exper[i,1:ncol(can)],exper[i,(ncol(can)+1):ncol(exper)])
    pvalue <- t_test$p.value
    
    normal_ave<-mean(exper[i,(ncol(can)+1):ncol(exper)])
    cancer_ave<-mean(exper[i,1:ncol(can)])
    
    log2foldchange<-cancer_ave - normal_ave
    
    result1<-cbind(rownames(exper)[i],cancer_ave,normal_ave,
                   log2foldchange,pvalue)
    result1<-as.matrix(result1)
    out_re<-rbind(out_re,result1)
  }
  
  out_re<-as.matrix(out_re)
  fdr <- p.adjust(out_re[,5],"BH")
  result <- cbind(out_re[,1:5],fdr)
  colnames(result)<-c("gene","cancar_ave","normal_ave",
                      "log2foldchange","pvalue","fdr")
  
  log2FC<-as.numeric(result[,4])
  fdr<-as.numeric(result[,6])
  test<-which(abs(log2FC)>=2&fdr<0.05)
  diff_gene<-as.data.frame(result[test,])
  
  cname1<-strsplit(file_list[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split='_')[[1]][1]  
  name<-paste(cname2,"_diff.csv")
  setwd('D:/差异表达结果')
  write.csv(result,cname1,row.names=F)
  write.csv(diff_gene,name,row.names=F)
}

data<-read.table("D:/Code/Data/Fig1/53gene-diffexpression_logFC.csv",header=T,row.names=1,sep=",")
data[data > 5] <- 5     
data[data < (-5)] <- (-5)

Fig.1B<-pheatmap(t(data),cellwidth = 8, cellheight = 15,
            cluster_rows = F,cluster_cols= F,fontsize_row = 11,fontsize_col=9, angle_col = "90",
            #legend = FALSE, 
            treeheight_col = 10,border_color="white", 
            color =colorRampPalette(c("#0033CC", "#FFFFFF", "#FF0000"))(100))
Fig.1B

####====Fig. 1C====####

dir = "D:/TCGA-VCAN_expression"
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)      

dat<-read.table(file_list[1],header=T,row.names=1,sep=",")
cname1<-strsplit(file_list[1],split="\\/")[[1]][8]
cname2<-strsplit(cname1,split='\\.')[[1]][1]
name<-trimws(cname2, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
dat$cancer<-name

for(f in 2:length(file_list)){
  
  biaodapu<-read.table(file_list[f],header = T,row.names=1,sep=",")
  cname1<-strsplit(file_list[f],split="\\/")[[1]][8]
  cname2<-strsplit(cname1,split='\\.')[[1]][1]
  name<-trimws(cname2, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  biaodapu$cancer<-name
  
  dat<-rbind(dat,biaodapu)
}

Fig.1C <- ggplot(data = dat,aes(x = cancer, y = value, fill = group))+
  geom_boxplot(aes(colour=group),width=0.6)+  
  theme(axis.text.x = element_blank())+  
  labs(x=NULL, y = NULL)+  
  scale_fill_manual(values = c("#FF6F00","#FFCCBC"))+
  scale_color_manual(values=c("#000000", "#000000"))+ 
  theme_classic()+ 
  theme(legend.position = 'none')+  
  theme(axis.text.x=element_text(size=12))+
  theme(axis.text.y=element_text(size=12))+ 
  theme(axis.title = element_text(size = 14))+
  theme(plot.title = element_text(hjust = 0.5))  
Fig.1C

####====Fig. 1D====####

dir = "D:/TCGA-53gene_expression"
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
###CNV gain###
z<-read.table(file_list[1],header=T,row.names=1,sep=",")
cname1<-strsplit(file_list[1],split="\\/")[[1]][7]
cname2<-strsplit(cname1,split='\\.')[[1]][1]
name<-trimws(cname2, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")

p<-matrix(0,nrow(z),1) 
rownames(p)<-rownames(z)
colnames(p)<-name
for(i in 1:nrow(z)){
  a<-length(which((z[i,]==1)))
  if(a==0){
    p[i]<-"0"
  }else{
    p[i]<-a/ncol(z) 
  }
}
merge_gain<-p

for(f in 2:length(file_list)){
  
  z<-read.table(file_list[f],header=T,row.names=1,sep=",")
  cname1<-strsplit(file_list[f],split="\\/")[[1]][7]
  cname2<-strsplit(cname1,split='\\.')[[1]][1]
  name<-trimws(cname2, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  p<-matrix(0,nrow(z),1) 
  rownames(p)<-rownames(z)
  colnames(p)<-name
  for(i in 1:nrow(z)){
    a<-length(which((z[i,]==1)))
    if(a==0){
      p[i]<-"0"
    }else{
      p[i]<-a/ncol(z) 
    }
  }
  merge_gain = cbind(merge_gain,p)
}
setwd('D:/CNV_result')
write.csv(merge_gain,"cnv_gain.csv") 

###CNV loss###
z<-read.table(file_list[1],header=T,row.names=1,sep=",")
cname1<-strsplit(file_list[1],split="\\/")[[1]][7]
cname2<-strsplit(cname1,split='\\.')[[1]][1]
name<-trimws(cname2, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")

p<-matrix(0,nrow(z),1) 
rownames(p)<-rownames(z)
colnames(p)<-name
for(i in 1:nrow(z)){
  a<-length(which((z[i,]==(-1))))
  if(a==0){
    p[i]<-"0"
  }else{
    p[i]<-a/ncol(z) 
  }
}
merge_loss<-p

for(f in 2:length(file_list)){
  
  z<-read.table(file_list[f],header=T,row.names=1,sep=",")
  cname1<-strsplit(file_list[f],split="\\/")[[1]][7]
  cname2<-strsplit(cname1,split='\\.')[[1]][1]
  name<-trimws(cname2, which = c("both", "left", "right"), whitespace = "[ \t\r\n]")
  
  p<-matrix(0,nrow(z),1) 
  rownames(p)<-rownames(z)
  colnames(p)<-name
  for(i in 1:nrow(z)){
    a<-length(which((z[i,]==(-1))))
    if(a==0){
      p[i]<-"0"
    }else{
      p[i]<-a/ncol(z) 
    }
  }
  merge_loss = cbind(merge_loss,p)
}
setwd('D:/CNV_result')
write.csv(merge_loss,"cnv_loss.csv") 

###画图###
cone1 = data.frame(x = c(0,0,1),                   
                   y = c(0,1,1))
cone1 %>% ggplot(aes(x=x, y=y)) +  geom_polygon(fill="#213c18")
cone2 = data.frame(x = c(0,1,1),                  
                   y = c(0,0,1)) 
cone2 %>%   ggplot(aes(x=x, y=y)) +  geom_polygon(fill = "#668c6f")
cone1 %>%   ggplot(aes(x=x, y=y)) +  geom_polygon(fill="#213c18")+  geom_polygon(data=cone2, fill = "#668c6f")


triangle<-function(a,b,type="up"){
  triangle_down <- function(a,b){
    data.frame(x=c(0,1,1)+a,
               y=c(0,0,1)+b,
               group=paste0(a,"_",b),stringsAsFactors = F)
  }
  triangle_up <- function(a,b){
    data.frame(x=c(0,0,1)+a,
               y=c(0,1,1)+b,
               group=paste0(a,"_",b),stringsAsFactors = F)
  }
  if(type=="up"){
    data <- do.call(rbind,lapply(1:b, function(i){
      do.call(rbind,lapply(1:a, triangle_up,i))
    }))
  }
  if(type=="down"){
    data <- do.call(rbind,lapply(1:b, function(i){
      do.call(rbind,lapply(1:a, triangle_down,i))
    }))
  }
  return(data)
}

updata <- triangle(33,53,"up")
downdata <- triangle(33,53,"down")
updata %>%   ggplot(aes(x=x, y=y)) +  geom_polygon(aes(group=group,fill=group))+  
  theme(legend.position = "none")
downdata %>%   ggplot(aes(x=x, y=y)) +  geom_polygon(aes(group=group,fill=group))+  
  theme(legend.position = "none")

merge_loss <- as.data.frame(merge_loss)
merge_loss <- tibble::rownames_to_column(merge_loss, "X") 
data1<- merge_loss %>% pivot_longer(-X, 
                               names_to = "type", values_to = "exp") %>%   rename(gene=X)

#定义因子顺序
group<-read.table("D:/Code/Data/RespiratorySig.txt",header=F,sep="\t")
group <- as.vector(unlist(group)) 
group <- as.data.frame(rev(group))
group <- as.matrix(group) 
data1$gene <- factor(data1$gene,levels=group)

Amp <- mutate(data1,a=as.numeric(as.factor(type)),
              b=as.numeric(as.factor(gene)))
require(gcookbook)
require(tidyverse)
Amp <- mutate(Amp,group = paste0(a,"_",b))%>%as.data.frame()
Amp1 <-inner_join(Amp,updata,by="group")

#
merge_gain <- as.data.frame(merge_gain)
merge_gain <- tibble::rownames_to_column(merge_gain, "X") 
data2<- merge_gain %>% pivot_longer(-X, 
                               names_to = "type", values_to = "exp") %>%   rename(gene=X)
data2$gene <- factor(data2$gene,levels=group)
DELL <- mutate(data2,a=as.numeric(as.factor(type)),
               b=as.numeric(as.factor(gene)))
require(gcookbook)
require(tidyverse)
DELL <- mutate(DELL,group = paste0(a,"_",b))%>%as.data.frame()
DELL1 <-inner_join(DELL,downdata,by="group")

##画图
Amp1 %>%   ggplot(aes(x=x, y=y)) +
  geom_polygon(aes(group=group,fill=exp))+
  theme(legend.position = "none")

DELL1 %>%   ggplot(aes(x=x, y=y)) +
  geom_polygon(aes(group=group,fill=exp))+
  theme(legend.position = "none")

library(ggnewscale)
xlabels <- Amp1 %>% filter(y==1) %>% distinct(type,.keep_all = T) %>%
  arrange(x) %>%  pull(type)

ylabels <- DELL1 %>%
  filter(x==1) %>%   
  distinct(gene,.keep_all = T) %>%   
  arrange(y) %>%   
  pull(gene) %>%  as.character()

red<-"#B71C1C" 
blue<-"#0D47A1" 
nake<-"#FFFADD" 

Amp1$exp<-as.numeric(Amp1$exp)
DELL1$exp<-as.numeric(DELL1$exp)
Amp1$exp[Amp1$exp > 0.5]<- 0.5
DELL1$exp[DELL1$exp >0.5] <- 0.5

Fig.1D <- ggplot(Amp1,aes(x=x, y=y)) +
  ## 画出上三角
  geom_polygon(data = Amp1,
               aes(x=x, y=y, group=group,fill=exp),
               color="#CACFD2") +  
  ## 上三角用表达值来配色
  scale_fill_gradientn(colours = c(nake,"#1E88E5",blue)) +
  ## 神技能，清空
  new_scale_fill()+ 
  ## 画出下三角
  geom_polygon(data = DELL1,
               aes(x=x, y=y, group=group,fill=exp)) +
  ## 下三角用表达值来配色
  scale_fill_gradientn(colours = c(nake,"#E53935",red))+
  ## 调整x轿
  scale_x_continuous(limits = c(1, NA),
                     expand = c(0,0),breaks = c(1:33)+0.5,
                     labels=xlabels) +
  ## 调整y轿
  scale_y_continuous(limits = c(1, NA),
                     expand = c(0,0),breaks = c(1:52)+0.5,
                     labels=ylabels) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x = element_text(size=12,angle = 45, hjust = 1),
        axis.text.y = element_text(size=12))+
  theme(legend.position="none")+ #bottom
  theme(plot.margin = margin(0.5,0.01,0.5,0.01, "cm"))+
  theme(axis.ticks.x  = element_blank(),axis.ticks.y  = element_blank()) 

Fig.1D

####====Fig. 1E====####

dir = "D:/TCGA-53gene_methylation(GDC_hg38_450k)" 
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)     
for(i in 1:length(file_list)){
  
  exp<-read.table(file_list[i],header = T,row.names=1,sep=",")
  
  cluster<- data.frame(colnames(exp),stringsAsFactors = FALSE)
  cluster[grep("\\.0[1-9][A-Z]",cluster[,1]),2] <- "cancer"
  cluster[grep("\\.11[A-Z]",cluster[,1]),2] <- "normal"
  sample_num<-as.data.frame(table(cluster[,2]))
  
  if(nrow(sample_num)==2 & sample_num[2,2]>5){
    
    group<-cluster
    can <- exp[,grep("\\.0[1-9][A-Z]",colnames(exp))] 
    nor <- exp[,grep("\\.11[A-Z]",colnames(exp))] 
    exper <- cbind(can,nor)
    exper<-as.matrix(exper)
    #limma差异表达
    colnames(group)=c("FileName","Target")  
    lev<-unique(group$Target)           
    f <- factor(group$Target, levels=lev) 
    design <- model.matrix(~0+f)         
    colnames(design) <- lev
    rownames(design) <-group[,1]
    
    eset<-exp
    cont.wt <- makeContrasts("cancer-normal",levels=design) 
    fit <- lmFit(eset, design)
    fit2 <- contrasts.fit(fit, cont.wt) 
    fit2 <- eBayes(fit2) 
    tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)
    tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
    colnames(tT)=c("FDR","P.Value","logFC")
    
    setwd("D:/差异甲基化(limma)")
    cname1<-strsplit(file_list[i],split="\\/")[[1]][8]
    cname2<-strsplit(cname1,split="\\.")[[1]][1]
    name1<-paste(cname2,".csv")
    write.csv(tT,name1)
    
    ##提取差异基因列表
    fdr = 0.05
    logFC = 0.1
    tT[fdr > tT[,"FDR"]  &  tT[,"logFC"] >= logFC, ncol(tT)+1] = "Up"
    tT[fdr > tT[,"FDR"]  & -logFC >= tT[,"logFC"], ncol(tT)] = "Down"
    tT[tT[,"FDR"] >= fdr | logFC > abs(tT[,"logFC"]) , ncol(tT)] = "Normal"
    colnames(tT)[ncol(tT)] = "Regulate"
    deg=tT[which(tT$Regulate!="Normal"),]
    
    name2<-paste(cname2,"_diff.csv")
    write.csv(deg,name2)
  }
}

df <- read.table("D:/Code/Data/Fig1/53gene_logFC.csv",header=T,sep=",")
group <- read.table("D:/Code/Data/RespiratorySig.txt",header=F,sep="\t")
group <- as.matrix(group) 
df$gene <- factor(df$gene,levels=group) 

#气泡图
pos=ifelse(df[,2]<0, "#486e9f", "#ec872d")
Fig.1E<- ggplot(df,aes(gene,cancer))+
  geom_point(aes(size=abs(df[,2])),alpha=0.8,color=pos)+ 
  scale_size(range = c(2.5, 5))+  
  xlab(NULL)+ ylab(NULL)+   
  theme_bw()+  
  theme(panel.grid.major=element_line(colour="#F5F5F5"))+   
  theme(axis.text.x=element_text(color = "#424242", size=12, angle=90,hjust=1))+ 
  theme(axis.text.y=element_text(color = "#424242", size=12))+
  theme(axis.ticks.x  = element_blank(),axis.ticks.y  = element_blank()) 
Fig.1E

#堆叠柱形图
df$type<-NA
for(i in 1:nrow(df)){
  if(any(is.na(df[i,2]))){
    df[i,4]<-NA
  }else{
    if(df[i,2]>0){
      df[i,4]<-"up"
    }else if(df[i,2]<0){
      df[i,4]<-"down"
    }
  }  
}
up<-df[which(df[,4]=="up"),] 
down<-df[which(df[,4]=="down"),] 

#统计每种癌症类型中上/下调呼吸基因数量
up_freq<-as.data.frame(table(up[,3]))
up_freq$type<-"up"
down_freq<-as.data.frame(table(down[,3]))
down_freq$type<-"down"
dat<-rbind(up_freq,down_freq)
colnames(dat)[1]<-"cancer"

library(reshape2)
library(ggpubr)
Fig.1E<-ggplot(dat, aes( x = cancer, weight = Freq, fill = type))+
  scale_fill_manual(values = c("#486e9f","#ec872d"))+ 
  geom_bar( position = "stack",width=0.8)+  
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme_classic()+
  #theme(legend.position = 'none') 
  labs(x=NULL, y = "Number")+
  theme(axis.text.x=element_text(size=12,angle=45,hjust=1))+
  theme(axis.text.y=element_text(size=12))
Fig.1E

####====Fig. 1F====####

dir = "D:/TCGA-gene_mutation" 
file_list = list.files(path =dir, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)

geneset<-read.table("D:/Code/Data/RespiratorySig.txt",header=F,sep="\t")
geneset<-as.matrix(geneset)

mut<-read.table(file_list[1],header=T,sep="\t",fill=T)
exp<-mut[mut$sample %in% geneset[,1],]
rownames(exp)<-exp$sample
exp<-exp[,-1]

cname1<-strsplit(file_list[1],split="\\/")[[1]][2]
cname2<-strsplit(cname1,split='_')[[1]][1]

#计算突变频率
p<-matrix(0,nrow(exp),1) 
rownames(p)<-rownames(exp)
colnames(p)<-cname2
for(i in 1:nrow(exp)){
  a<-length(which((exp[i,]==1)))
  if(a==0){
    p[i]<-"0"
  }else{
    p[i]<-a/ncol(exp) 
  }
}
merge_p<-p

for(f in 2:length(file_list)){
  
  mut<-read.table(file_list[f],header=T,sep="\t",fill=T)
  exp<-mut[mut$sample %in% geneset[,1],]
  rownames(exp)<-exp$sample
  exp<-exp[,-1]
  
  cname1<-strsplit(file_list[f],split="\\/")[[1]][2]
  cname2<-strsplit(cname1,split='_')[[1]][1]
  
  p<-matrix(0,nrow(exp),1) 
  rownames(p)<-rownames(exp)
  colnames(p)<-cname2
  for(i in 1:nrow(exp)){
    a<-length(which((exp[i,]==1)))
    if(a==0){
      p[i]<-"0"
    }else{
      p[i]<-a/ncol(exp) 
    }
  }
  merge_p = cbind(merge_p,p)
}
setwd("D:/突变频率")
write.csv(merge_p,"mutation.csv") 

##画图
library(reshape2)
library(ggpubr)
dat<-read.table('D:/Code/Data/Fig1/53gene_mutation_freq.CSV',header=T,sep=',')
group<-read.table("D:/Code/Data/RespiratorySig.txt",header=F,sep="\t")
group<-as.matrix(group)
dat$gene <- factor(dat$gene,levels=group)

Fig.1F <- ggplot(dat, aes( x = gene, weight = value, fill = group))+
  scale_fill_manual(values = c("#5C6BC0","#7986CB","#9FA8DA","#C5CAE9",
                               "#EC407A","#F06292","#F48FB1","#F8BBD0",
                               "#FF7043","#FF8A65","#FFAB91","#FFCCBC",
                               "#FFA726","#FFB74D","#FFCC80","#FFE0B2",
                               "#FFCA28","#FFD54F","#FFE082","#FFECB3",
                               "#EF5350","#E57373","#EF9A9A","#FFCDD2",
                               "#26A69A","#4DB6AC","#80CBC4","#B2DFDB", 
                               "#26C6DA","#4DD0E1","#80DEEA","#B2EBF2", 
                               "#BBDEFB"))+ 
  geom_bar( position = "stack",width=0.8)+    
  theme_bw()+ 
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  theme_classic()+
  #theme(legend.position = 'none')  
  labs(y = "Percentage of mutation", x=NULL)+   
  theme(axis.text.x=element_text(size=12,angle=90))+
  theme(axis.text.y=element_text(size=12))  
Fig.1F



