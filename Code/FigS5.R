
###===Figure S5===###

rm (list=ls())
gc()

#====loading data====####

library(data.table)
library(clusterProfiler)
library(forcats)
library(tidyverse)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)

####====Fig. S5A====####

dir1 = 'D:/53gene-hallmark_pathway相关性/pos'
file_list1 = list.files(path =dir1, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)
dir2 = 'D:/53gene-hallmark_pathway相关性/neg'
file_list2 = list.files(path =dir2, pattern = "*.txt$",recursive = TRUE,full.names = TRUE)

##正相关
f=13 #KIRP
merge_pos = read.table(file_list1[f],header=T,check.names=F,sep="\t") 

pos_freq<-as.data.frame(table(merge_pos$pathway))
colnames(pos_freq)<-c("pathway","number")
pos_freq<-pos_freq[order(pos_freq$number,decreasing = T),]

##提取gene
for(i in 1:nrow(pos_freq)){
  path <- merge_pos[which(merge_pos[,2]==pos_freq[i,1]),]
  
  if(class(path)=="character"){
    path<-as.matrix(path)
    path<-t(path)
    colnames(path)[1]<-pos_freq[i,1]
    path<-path[,1,drop=F]
  }else{
    path<-as.data.frame(path[,1])
    colnames(path)<-pos_freq[i,1]
  }
  setwd("D:/KIRP-与每个通路相关的gene/pos")
  name<-paste(pos_freq[i,1],".txt")
  write.table(path,name,sep="\t",row.names = F,col.names=T,quote=F)
} 

##负相关
f=13 #KIRP
merge_neg = read.table(file_list2[f],header=T,check.names=F,sep="\t") 

neg_freq<-as.data.frame(table(merge_neg$pathway))
colnames(neg_freq)<-c("pathway","number")
neg_freq<-neg_freq[order(neg_freq$number,decreasing = T),]

##提取gene
for(i in 1:nrow(neg_freq)){
  path <- merge_pos[which(merge_pos[,2]==neg_freq[i,1]),]
  
  if(class(path)=="character"){
    path<-as.matrix(path)
    path<-t(path)
    colnames(path)[1]<-neg_freq[i,1]
    path<-path[,1,drop=F]
  }else{
    path<-as.data.frame(path[,1])
    colnames(path)<-neg_freq[i,1]
  }
  setwd("D:/KIRP-与每个通路相关的gene/neg")
  name<-paste(neg_freq[i,1],".txt")
  write.table(path,name,sep="\t",row.names = F,col.names=T,quote=F)
}

##整理画图数据
dat<-read.table('D:/Code/Data/FigS4/KIRP-与每个通路相关的差异表达基因表达谱.csv',header=T,sep=',')
datnew <- dat[,-2]
group1 <- data.frame(colnames(datnew)[2:ncol(datnew)],stringsAsFactors = F)
group1[grep("\\.0[1-9][A-Z]",group1[,1]),2] <- "cancer"
group1[grep("\\.11[A-Z]",group1[,1]),2] <- "normal"
group1 <- as.data.frame(group1[order(group1[,2]),])

datnew <- datnew[,pmatch(group1[,1],colnames(datnew))] 
datnew$gene<-dat[,1]
datnew <- datnew[,c(ncol(datnew),1:(ncol(datnew)-1))]
datnew <- data.frame(datnew,stringsAsFactors = F)
datnew2<-datnew[,-1]
datnew2<-as.matrix(datnew2)

rownames(group1)<-group1[,1]
group1<-group1[,-1,drop=F] 
colnames(group1)<-'Sample'

##标准化
rowmean=apply(datnew2,1,mean) 
rowsd=apply(datnew2,1,sd)
res1=sweep(datnew2,1,rowmean) 
res2=sweep(res1,1,rowsd,'/')   

res2[res2>2]<-2
res2[res2<(-2)]<-(-2)

group2 = data.frame(Path = factor(rep(c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION","HALLMARK_FATTY_ACID_METABOLISM",
                                        "HALLMARK_BILE_ACID_METABOLISM ","HALLMARK_ADIPOGENESIS",
                                        "HALLMARK_XENOBIOTIC_METABOLISM","HALLMARK_KRAS_SIGNALING_DN",
                                        "HALLMARK_PEROXISOME"), 
                                      c(1,3,3,1,3,1,2))))
ann_colors = list(Path = c(HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION="#2980B9", 
                           HALLMARK_FATTY_ACID_METABOLISM="#F48FB1",
                           HALLMARK_BILE_ACID_METABOLISM="#FBC02D", 
                           HALLMARK_ADIPOGENESIS="#FF8F00",
                           HALLMARK_XENOBIOTIC_METABOLISM="#26C6DA",
                           HALLMARK_KRAS_SIGNALING_DN="#A3E4D7",
                           HALLMARK_PEROXISOME="#D1C4E9"),                 
                  Sample = c(cancer="#EF5350", normal="#00897B"))
#
x <- replicate(14,rnorm(100,0,1),simplify = FALSE)
x <- setNames(x,paste0("A",seq_along(x)))
rownames(group2)<- names(x)
rownames(res2)<-names(x)

library(pheatmap)
Fig.S5A <- pheatmap(res2,
                    cluster_rows = F,   
                    cluster_cols = F,               
                    annotation_col =group1,         
                    annotation_row =group2,
                    annotation_legend=T,              
                    show_rownames = F,       
                    show_colnames = F,             
                    scale = "row",             
                    annotation_colors = ann_colors,
                    color =colorRampPalette(c("#0000FF", "#ffffff","#FF0000"))(100))
Fig.S5A

####====Fig. S5B====####

exp<-read.table("D:/KIRP-gene_expression.csv",header=T,row.names=1,sep=",")
del <- which(apply(exp==0, 1, sum) > 0.3*ncol(exp))
exp <- exp[-del,] 
texp<-as.data.frame(t(exp))
group<-as.data.frame(as.numeric(texp$VCAN > median(texp$VCAN)))

res<-cbind(texp,group[,1])
colnames(res)[ncol(res)]<-'group'
tres<-t(res)
write.csv(tres,"D:/VCAN高低风险组分组/KIRP.csv")

high<-res[which(res[,ncol(res)]=="1"),]  
low<-res[which(res[,ncol(res)]=="0"),]
exper <- rbind(high,low)
exper<-as.matrix(exper[,-(ncol(exper))])

##t test
out_re<-c()
for(i in 1:ncol(exper)){
  
  t_test<-t.test(exper[1:nrow(high),i],exper[(nrow(high)+1):nrow(exper),i])
  pvalue <- t_test$p.value
  
  high_ave<-mean(exper[1:nrow(high),i])
  low_ave<-mean(exper[(nrow(high)+1):nrow(exper),i])
  log2foldchange<-high_ave - low_ave   
  
  result1<-cbind(colnames(exper)[i],high_ave,low_ave,
                 log2foldchange,pvalue)
  result1<-as.matrix(result1)
  out_re<-rbind(out_re,result1)
  
}
out_re<-as.matrix(out_re)
fdr <- p.adjust(out_re[,5],"BH")
result <- cbind(out_re[,1:5],fdr)
colnames(result)<-c("gene","high_ave","low_ave",
                    "log2foldchange","pvalue","fdr")
setwd('D:/VCAN高低风险组分组差异表达')
write.csv(result,"KIRP.csv")

log2<-abs(as.numeric(result[,4]))
fdr<-as.numeric(result[,6])
test<-which(log2>=1&fdr<=0.05)
diff_result<-result[test,]
write.csv(diff_result,"KIRP_difgene.csv")

##富集分析
gene_input<-diff_result[,c('gene','log2foldchange')]
genename <- as.character(gene_input[,1])
gene_map <- bitr(genename, 
                 fromType = "SYMBOL", 
                 toType = c("ENTREZID"), 
                 OrgDb = org.Hs.eg.db) 
colnames(gene_map)[1]<-'gene'

aaa<-merge(gene_map,gene_input,by="gene",all=F)
aaa<-aaa[order(aaa$log2foldchange,decreasing = T),] 
geneList = aaa[,3]
names(geneList) = as.character(aaa[,2])
#
kegmt<-read.gmt("D:/Code/Data/Fig4/h.all.v7.4.entrez.gmt")
hallmark<-GSEA(geneList,TERM2GENE = kegmt,pvalueCutoff=0.05) 

go <- gseGO(geneList, 'org.Hs.eg.db', keyType = "ENTREZID", 
            ont="all", nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
kegg <- gseKEGG(geneList, nPerm = 1000, 
                minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
hallmark_res<-as.data.frame(hallmark) 
go_res<-as.data.frame(go)
kegg_res<-as.data.frame(kegg)

setwd('D:/VCAN富集分析')
write.csv(hallmark_res,"KIRP-HALLMARK_result.csv")
write.csv(go_res,"KIRP-GO_result.csv")
write.csv(kegg_res,"KIRP-KEGG_result.csv")

##画图
GO_upanddown_part<-read.table("D:/Code/Data/FigS4/KIRP-GO_result.csv",header=T,row.names=1,sep=",")
HALLMARK_upanddown_part<-read.table("D:/Code/Data/FigS4/KIRP-HALLMARK_result.csv",header=T,row.names=1,sep=",")
KEGG_upanddown_part<-read.table("D:/Code/Data/FigS4/KIRP-KEGG_result.csv",header=T,row.names=1,sep=",")

CRGs_FEA_plot<-rbind.data.frame(GO_upanddown_part,HALLMARK_upanddown_part,KEGG_upanddown_part,stringsAsFactors = FALSE)
CRGs_FEA_plot<-CRGs_FEA_plot[dim(CRGs_FEA_plot)[1]:1,]
CRGs_FEA_plot$Description<-fct_inorder(as.character(CRGs_FEA_plot$Description))
(scaleFactor <-  max(-log10(CRGs_FEA_plot$pvalue))/max(CRGs_FEA_plot$setSize))

library(ggplot2)
Fig.S5B <- ggplot(CRGs_FEA_plot, aes(x=CRGs_FEA_plot$Description)) + 
  geom_bar(aes(y =  -log10(pvalue)/scaleFactor, fill = ONTOLOGY), stat = "identity", size = 0.8, width = 0.7) +
  scale_fill_manual(values = c("#E21F26","orangered","gold2","#4DD0E1","#2278B3")) + 
  geom_point(aes(y = setSize), shape = 21, fill = "white") +
  geom_line(aes(y = setSize, group = 1)) +
  
  scale_x_discrete(expand = expansion(mult = c(0.015,0.015)), 
                   limits = rev(levels(as.character(CRGs_FEA_plot$Description)))) + 
  scale_y_continuous(name = "Counts",  
                     sec.axis=sec_axis(trans = ~. * scaleFactor, 
                                       name = "-log10(pvalue)",
                                       breaks = c(0,1,2,3,4,5)),
                     expand = c(0, 0)) + 
  labs(x="") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(axis.text.y=element_text(angle=45,size=12))+ 
  coord_flip()
Fig.S5B

####====Fig. S5C====####

library(enrichplot)
Fig.S5C <- gseaplot2(go, "GO:0009893", color = "green", 
                     title='Positive Regulation of Metabolic Process',
                     base_size = 15,
                     pvalue_table = T,subplots=1:2,rel_heights=c(1, .2, .6))
Fig.S5C <- gseaplot2(kegg, "hsa01100", color = "green", 
                     title='Metabolic pathways',
                     base_size = 15,
                     pvalue_table = T,subplots=1:2,rel_heights=c(1, .2, .6))
Fig.S5C


