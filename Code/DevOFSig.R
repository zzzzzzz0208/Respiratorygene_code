
####===Development of respiratory gene signatures===####
library(GSVA)

rm(list=ls())
gc()

##################################################################################
####================================1.导入基因集==============================####
##################################################################################
res_gene<-read.table("D:/Code/Data/收集三个阶段呼吸基因/respiratory gene.txt",header=F,sep="\t")
res_gene<-as.data.frame(res_gene[!duplicated(res_gene[,c(1)]),]) 
colnames(res_gene)<-c("gene")

pathwayGene<-read.table("D:/Code/Data/3functional_pathway.txt",header = T,sep = "\t")

geneset1<-read.table("D:/Code/Data/收集三个阶段呼吸基因/glycolysis.csv",header = F,sep="\t")
geneset2<-read.table("D:/Code/Data/收集三个阶段呼吸基因/TCA cycle.csv",header = F,sep="\t")
geneset3<-read.table("D:/Code/Data/收集三个阶段呼吸基因/oxphos.csv",header = F,sep="\t")

##################################################################################
####===================2.计算肿瘤和正常样本中的差异表达基因===================####
##################################################################################
dir = "D:/TCGA-gene_expression" #17cancer
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)

for(f in 1:length(file_list)){
  
  exp<-read.table(file_list[f],header=T,row.names=1,sep=",")
  
  del <- which(apply(exp==0, 1, sum) > 0.3*ncol(exp))
  exp <- exp[-del,]   
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
  setwd('D:/TCGA差异表达结果')
  write.csv(diff_gene,name,row.names=F)
}

##################################################################################
####=========================3.计算通路活性评分并分组=========================####
##################################################################################
dir = "D:/TCGA-gene_expression" 
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)      
for(f in 1:length(file_list)){
  
  exp<-read.table(file_list[f],header=T,row.names=1,check.names=F,sep=",")
  geneset<-as.matrix(res_gene)
  
  res_exp<-exp[rownames(exp) %in% geneset[,1],]
  
  cname1<-strsplit(file_list[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split='_')[[1]][1]
  name1<-paste(cname2,".csv")
  setwd("D:/TCGA-respiratory_gene_expression")
  write.csv(res_exp,name1)
  
  del<-which(apply(res_exp == 0, 1, sum) > 0.3*ncol(res_exp))
  res_exp<-res_exp[-del,]
  
  ###ssGSEA###
  t = as.matrix(res_exp)
  gs = as.list(pathwayGene)   
  gs = lapply(gs, function(x) x[!is.na(x)])
  score = gsva(t, gs, method = "ssgsea", ssgsea.norm = TRUE, verbose = TRUE)
  setwd("D:/ssgsea_result")
  name2<-paste(cname2,"_ssgsea.csv")
  write.csv(score, name2)

  dat<-as.data.frame(t(score))
  #Glycolysis
  group1<-ifelse(dat$GLYCOLYSIS>median(dat$GLYCOLYSIS),'high','low')
  class1<-cbind(dat[,1],group1)
  rownames(class1)<-rownames(dat)
  colnames(class1)<-c("GLYCOLYSIS","group")
  setwd("D:/根据ssgsea_result对样本分组/GLYCOLYSIS")
  name3<-paste(cname2,"_GLYCOLYSIS.csv")
  write.csv(class1,name3)
  
  #TCA cycle
  group2<-ifelse(dat$TCA_CYCLE>median(dat$TCA_CYCLE),'high','low')
  class2<-cbind(dat[,2],group2)
  rownames(class2)<-rownames(dat)
  colnames(class2)<-c("TCA_CYCLE","group")
  setwd("D:/根据ssgsea_result对样本分组/TCACYCLE")
  name4<-paste(cname2,"_TCA_CYCLE.csv")
  write.csv(class2,name4)
  
  #OXPHOS
  group3<-ifelse(dat$OXPHOS>median(dat$OXPHOS),'high','low')
  class3<-cbind(dat[,3],group3)
  rownames(class3)<-rownames(dat)
  colnames(class3)<-c("OXPHOS","group")
  setwd("D:/根据ssgsea_result对样本分组/OXPHOS")
  name5<-paste(cname2,"_OXPHOS.csv")
  write.csv(class3,name5)
}

##################################################################################
####=======================4.比较高低分组间差异表达基因=======================####
##################################################################################
dir1 = "D:/TCGA-respiratory_gene_expression"
dir2 = "D:/根据ssgsea_result对样本分组/GLYCOLYSIS"
dir3 = "D:/根据ssgsea_result对样本分组/TCACYCLE"
dir4 = "D:/根据ssgsea_result对样本分组/OXPHOS"
file_list1 = list.files(path =dir1, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)     
file_list2 = list.files(path =dir2, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)  
file_list3 = list.files(path =dir3, pattern = "*.csv$",recursive = TRUE,full.names = TRUE) 
file_list4 = list.files(path =dir4, pattern = "*.csv$",recursive = TRUE,full.names = TRUE) 

for(f in 1:17){
  
  exp<-read.table(file_list1[f],header=T,row.names=1,sep=",")
  del<-which(apply(exp==0, 1, sum) > 0.3*ncol(exp))
  exp<-exp[-del,] 
  
  ###GLYCOLYSIS###
  geneset1<-as.matrix(geneset1)
  res1<-exp[rownames(exp)%in% geneset1[,1],]
  
  class1<-read.table(file_list2[f],header=T,row.names=1,sep=",")
  tres1<-as.data.frame(t(res1))
  dat1<-cbind(class1[,2],tres1)
  colnames(dat1)[1]<-'group'
  
  high<-dat1[which(dat1$group=='high'),]
  low<-dat1[which(dat1$group=='low'),]
  exper<-rbind(high,low)
  exper1<-as.matrix(exper[,-1])
  
  #t test
  out_re<-c()
  for(i in 1:ncol(exper1)){
    
    t_test<-t.test(exper1[1:nrow(high),i],exper1[(nrow(high)+1):nrow(exper1),i])
    pvalue <- t_test$p.value
    
    high_ave<-mean(exper1[1:nrow(high),i])
    low_ave<-mean(exper1[(nrow(high)+1):nrow(exper1),i])
    
    log2foldchange <- high_ave - low_ave
    
    result1<-cbind(colnames(exper1)[i],high_ave,low_ave,
                   log2foldchange,pvalue)
    result1<-as.matrix(result1)
    out_re<-rbind(out_re,result1)
  }
  out_re<-as.matrix(out_re)
  fdr <- p.adjust(out_re[,5],"BH")
  result_gly <- cbind(out_re[,1:5],fdr)
  colnames(result_gly)<-c("gene","high_ave","low_ave",
                          "log2foldchange","pvalue","fdr")
  
  result_gly[,4]<-as.matrix(result_gly[,4])
  result_gly[,6]<-as.matrix(result_gly[,6])
  
  log2FC<-as.numeric(result_gly[,4])
  fdr<-as.numeric(result_gly[,6])
  test1<-which(log2FC>=2&fdr<0.05)
  up<-result_gly[test1,]
  test2<-which(log2FC<=(-2)&fdr<0.05)
  down<-result_gly[test2,]
  
  if(length(test1)==1){
    up<-as.data.frame(t(up))
  }
  if(length(test2)==1){
    down<-as.data.frame(t(down))
  }
  
  diff_gene<-rbind(up,down)
  
  setwd("D:/高低组差异表达/GLYCOLYSIS")
  cname1<-strsplit(file_list1[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split='\\.')[[1]][1]
  name1<-paste(cname2,"_up.csv",sep='')
  name2<-paste(cname2,"_down.csv",sep='')
  name3<-paste(cname2,"_diff.csv",sep='')
  write.csv(result_gly,cname1,row.names=F)
  write.csv(up,name1,row.names=F)
  write.csv(down,name2,row.names=F)
  write.csv(diff_gene,name3,row.names=F)
  
  ###TCA CYCLE###
  geneset2<-as.matrix(geneset2)
  res2<-exp[rownames(exp)%in% geneset2[,1],]
  
  class2<-read.table(file_list3[f],header=T,row.names=1,sep=",")
  tres2<-as.data.frame(t(res2))
  dat2<-cbind(class2[,2],tres2)
  colnames(dat2)[1]<-'group'
  
  high<-dat2[which(dat2$group=='high'),]
  low<-dat2[which(dat2$group=='low'),]
  exper<-rbind(high,low)
  exper2<-as.matrix(exper[,-1])

  out_re<-c()
  for(i in 1:ncol(exper2)){
    
    t_test<-t.test(exper2[1:nrow(high),i],exper2[(nrow(high)+1):nrow(exper2),i])
    pvalue <- t_test$p.value
    
    high_ave<-mean(exper2[1:nrow(high),i])
    low_ave<-mean(exper2[(nrow(high)+1):nrow(exper2),i])
    
    log2foldchange<-high_ave - low_ave
    
    result1<-cbind(colnames(exper2)[i],high_ave,low_ave,
                   log2foldchange,pvalue)
    result1<-as.matrix(result1)
    out_re<-rbind(out_re,result1)
  }
  out_re<-as.matrix(out_re)
  fdr <- p.adjust(out_re[,5],"BH")
  result_tca <- cbind(out_re[,1:5],fdr)
  colnames(result_tca)<-c("gene","high_ave","low_ave",
                          "log2foldchange","pvalue","fdr")
  
  result_tca[,4]<-as.matrix(result_tca[,4])
  result_tca[,6]<-as.matrix(result_tca[,6])
  
  log2FC<-as.numeric(result_tca[,4])
  fdr<-as.numeric(result_tca[,6])
  test3<-which(log2FC>=2&fdr<0.05)
  up<-result_tca[test3,]
  test4<-which(log2FC<=(-2)&fdr<0.05)
  down<-result_tca[test4,]
  
  if(length(test3)==1){
    up<-as.data.frame(t(up))
  }
  if(length(test4)==1){
    down<-as.data.frame(t(down))
  }
  
  diff_gene<-rbind(up,down)
  
  setwd("D:/高低组差异表达/TCA CYCLE")
  write.csv(result_tca,cname1,row.names=F)
  write.csv(up,name1,row.names=F)
  write.csv(down,name2,row.names=F)
  write.csv(diff_gene,name3,row.names=F)
  
  ###OXPHOS###
  geneset3<-read.table("D:/geneset/oxphos.csv",header = F,sep="\t")
  geneset3<-as.matrix(geneset3)
  res3<-exp[rownames(exp)%in% geneset3[,1],]
  
  class3<-read.table(file_list4[f],header=T,row.names=1,sep=",")
  tres3<-as.data.frame(t(res3))
  dat3<-cbind(class3[,2],tres3)
  colnames(dat3)[1]<-'group'
  
  high<-dat3[which(dat3$group=='high'),]   
  low<-dat3[which(dat3$group=='low'),]   
  exper<-rbind(high,low)
  exper3<-as.matrix(exper[,-1])
  out_re<-c()
  
  for(i in 1:ncol(exper3)){
    
    t_test<-t.test(exper3[1:nrow(high),i],exper3[(nrow(high)+1):nrow(exper3),i])
    pvalue <- t_test$p.value
    
    high_ave<-mean(exper3[1:nrow(high),i])
    low_ave<-mean(exper3[(nrow(high)+1):nrow(exper3),i])
    
    log2foldchange<-high_ave - low_ave
    
    result1<-cbind(colnames(exper3)[i],high_ave,low_ave,
                   log2foldchange,pvalue)
    result1<-as.matrix(result1)
    out_re<-rbind(out_re,result1)
  }
  out_re<-as.matrix(out_re)
  fdr <- p.adjust(out_re[,5],"BH")
  result_oxp <- cbind(out_re[,1:5],fdr)
  colnames(result_oxp)<-c("gene","high_ave","low_ave",
                      "log2foldchange","pvalue","fdr")
  
  result_oxp[,4]<-as.matrix(result_oxp[,4])
  result_oxp[,6]<-as.matrix(result_oxp[,6])
  
  log2FC<-as.numeric(result_oxp[,4])
  fdr<-as.numeric(result_oxp[,6])
  test5<-which(log2FC>=2&fdr<0.05)
  up<-result_oxp[test5,]
  test6<-which(log2FC<=(-2)&fdr<0.05)
  down<-result_oxp[test6,]
  
  if(length(test5)==1){
    up<-as.data.frame(t(up))
  }
  if(length(test6)==1){
    down<-as.data.frame(t(down))
  }
  
  diff_gene<-rbind(up,down)
  
  setwd("D:/高低组差异表达/OXPHOS")
  write.csv(result_oxp,cname1,row.names=F)
  write.csv(up,name1,row.names=F)
  write.csv(down,name2,row.names=F)
  write.csv(diff_gene,name3,row.names=F)
}

##################################################################################
####=========================5.两部分差异表达基因取交集=======================####
##################################################################################
dir = "D:/TCGA差异表达结果"
dir1 = "D:/高低组差异表达/GLYCOLYSIS"
dir2 = "D:/高低组差异表达/TCA CYCLE"
dir3 = "D:/高低组差异表达/OXPHOS"
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
file_list1 = list.files(path =dir1, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
file_list2 = list.files(path =dir2, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)
file_list3 = list.files(path =dir3, pattern = "*.csv$",recursive = TRUE,full.names = TRUE)

for(f in 1:length(file_list)){
  
  diff_gene<-read.table(file_list[f],header=T,sep=",")
  glycolysis<-read.table(file_list1[f],header=T,sep=",")
  
  jiaoji1<-as.data.frame(intersect(glycolysis[,1],diff_gene[,1]))
  
  cname1<-strsplit(file_list[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split='_')[[1]][1]  
  name1<-paste(cname2,'_glycolysis.csv',sep='')
  setwd("D:/差异表达基因交集")
  write.csv(jiaoji,name1,row.names=F)

  tca_cycle<-read.table(file_list2[f],header=T,sep=",")
  jiaoji2<-as.data.frame(intersect(tca_cycle[,1],diff_gene[,1]))
  name2<-paste(cname2,'_tca_cycle.csv',sep='')
  write.csv(jiaoji2,name2,row.names=F)
  
  oxphos<-read.table(file_list3[f],header=T,sep=",")
  jiaoji3<-as.data.frame(intersect(oxphos[,1],diff_gene[,1]))
  name3<-paste(cname2,'_oxphos.csv',sep='')
  write.csv(jiaoji3,name3,row.names=F)
}


