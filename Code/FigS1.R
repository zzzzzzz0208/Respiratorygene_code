
###===Figure S1===###

rm (list=ls())
gc()

#====loading data====####

library(ggplot2)

####====Fig. S1====####

dir = "D:/GEO-VCAN_expression" 
file_list = list.files(path =dir, pattern = "*.csv$",recursive = TRUE,full.names = TRUE) 

for(f in file_list){
  
  df<-read.table(file_list[f],header = T,row.names=1,sep=",")
  cname1<-strsplit(file_list[f],split="/")[[1]][2]
  cname2<-strsplit(cname1,split='\\.')[[1]][1]
  
  Fig.S1 <- ggplot(data = df,aes(x = group, y = value, fill = group))+
    geom_boxplot(aes(colour=group),width=0.5)+ 
    labs(title=c(cname2))+
    theme_bw()+  
    theme(axis.text.x = element_blank())+ 
    labs(x=NULL, y = NULL)+ 
    theme(axis.text.y=element_text(size=10))+ 
    theme(axis.title = element_text(size = 12))+
    scale_fill_manual(values = c("#FFCDD2","#BBDEFB"))+ 
    scale_color_manual(values=c("#000000", "#000000"))+ 
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+   ##去除网格线
    theme(legend.position = 'none')+   
    theme(panel.border = element_rect(color="#000000"))+  
    theme(plot.margin=unit(rep(1,4),'lines'))+  
    theme(plot.title = element_text(hjust = 0.5))
  Fig.S1
}

