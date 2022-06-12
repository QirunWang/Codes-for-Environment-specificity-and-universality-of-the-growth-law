######GSEA analysis of gene functions
######For figure 4c and figure S3

rm(list=ls())
library(DESeq2)
library(readxl)
library(readr)
temp1<-read_xlsx("Data/riboseq/riboseq_res_minus.xlsx",sheet = 1)
rowname0<-temp1$ORF
temp1<-temp1[,c(-1:-3)]
colnames(temp1)<-c(paste0(rep("wt_glu",8),"_",c(1:8)),
                   paste0(rep("da_glu",8),"_",c(1:8)),
                   paste0(rep("wt_gly",8),"_",c(1:8)),
                   paste0(rep("da_gly",8),"_",c(1:8)))
rownames(temp1)<-rowname0
mycounts<-temp1

temp2<-data.frame(id=colnames(temp1),
                  treatment=factor(c(rep("wt_glu",8),
                                     rep("da_glu",8),
                                     rep("wt_gly",8),
                                     rep("da_gly",8))))
coldata<-temp2
dds <- DESeqDataSetFromMatrix(mycounts, coldata, design= ~ treatment)
dds <- DESeq(dds)
#plotPCA(rlogTransformation(dds), intgroup=c('treatment'))

res_wt<-results(dds, contrast=c("treatment","wt_gly","wt_glu"))
res_da<-results(dds, contrast=c("treatment","da_gly","da_glu"))



###############gsea analysis
library(clusterProfiler)
library(ggplot2)
library(org.Sc.sgd.db)
library(ggthemes)



ki=read.csv("Data/elongation_speed/ki.csv")
colnames(ki)<-c("geneid","ki(Da/s)")
ribo_pro<-read.csv("ribopro_eu.csv")

##for wt

temp1<-as.data.frame(res_wt)
temp1$geneid<-rownames(temp1)
temp2<-setdiff(temp1$geneid,ribo_pro$geneid)
temp3<-data.frame(geneid=temp2)
temp4<-merge(temp3,temp1,by="geneid")
temp5<-temp4[which(temp4$padj<5e-2),]
temp6<-merge(temp5,ki,by="geneid")
non_ribo_genes<-temp6

genes<-as.vector(non_ribo_genes$`ki(Da/s)`)
names(genes)<-non_ribo_genes$geneid
genes<-sort(genes,decreasing = T)
res_gsea_ki<-gseGO(genes,OrgDb="org.Sc.sgd.db", keyType = "ENSEMBL",
                       nPerm = 100000,pAdjustMethod = "BH",
                       pvalueCutoff = 0.25)
temp1<-res_gsea_ki@result
temp2<-temp1[which(temp1$pvalue<0.05),]


genes<-as.vector(non_ribo_genes$log2FoldChange)
names(genes)<-non_ribo_genes$geneid
genes<-sort(genes,decreasing = T)
res_gsea_log2fc<-gseGO(genes,OrgDb="org.Sc.sgd.db", keyType = "ENSEMBL",
                       nPerm = 100000,pAdjustMethod = "BH",
                       pvalueCutoff = 0.25)

temp3<-res_gsea_log2fc@result
temp4<-temp3[which(temp3$pvalue<0.05),]
temp5<-merge(temp2,temp4,by="ID")
gsea_res<-temp5
temp2<-vector()
for (i in 1:nrow(gsea_res)) 
{
  temp2<-c(temp2,paste0(gsea_res$ID[i]," ",gsea_res$Description.x[i])) 
}
gsea_res$Name<-temp2

temp6<-gsea_res[which(gsea_res$NES.x>0&gsea_res$NES.y>0),]
temp6$class<-rep("Response to stress",nrow(temp6))
temp7<-gsea_res[which(gsea_res$NES.x<0&gsea_res$NES.y<0),]
temp7$class<-rep("Response to stress",nrow(temp7))

####bar plot
temp1<-c(log10(gsea_res$p.adjust.x),-log10(gsea_res$p.adjust.y))
temp3<-data.frame(ID=rep(gsea_res$Name,2),log10padj=temp1)
temp3$datasource<-c(rep("ki",8),rep("log2fc",8))
temp3<-temp3[order(temp3$log10padj),]

ggplot()+
  geom_bar(data = temp3,aes(x = reorder(ID,log10padj),y = log10padj,fill = datasource),
           stat = 'identity',position = 'stack')+
  geom_bar(data = temp3,aes(x =reorder(ID,log10padj),y = log10padj,fill = datasource),
          stat = 'identity',position = 'stack')+
  scale_y_continuous(limits = c(-2.5,2.5),breaks = c(-2.5,0,2.5),labels = c(2.5,0,2.5))+
    theme_bw()+
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        axis.text = element_text(size = 25,hjust = 0,color = 'black'),
        axis.title = element_text(size = 25,hjust = 0,color = 'black'),
        legend.text = element_text(size = 20,hjust = 15,color = 'black'),
        legend.title = element_blank())+
    xlab('')+
    ylab('Log')+
    coord_flip() 
ggsave("gsea_wt_padj.eps",width = 20,height = 5)

####scatter plot
ggplot(data=gsea_res, aes(x=NES.x, y=NES.y,color=Name))+
  geom_point(size=5)+
  geom_vline(xintercept=0,lty=2,col="black",lwd=0.8) +
  geom_hline(yintercept =0,lty=2,col="black",lwd=0.8) +
  labs(x="NES of ki",
       y="NES of log2(fold change)")+
  scale_y_continuous(limits = c(-3,3))+
  scale_x_continuous(limits = c(-3,3))+
  theme_bw()+
  theme(axis.text = element_text(size = 25,hjust = 0,color = 'black'),
        axis.title = element_text(size = 25,hjust = 0,color = 'black'),
        legend.text = element_text(size = 20,hjust = 0,color = 'black'),
        legend.title = element_blank())
ggsave("gsea_wt_nes.eps",width = 12
       ,height = 5)


##for da

temp1<-as.data.frame(res_da)
temp1$geneid<-rownames(temp1)
temp2<-setdiff(temp1$geneid,ribo_pro$geneid)
temp3<-data.frame(geneid=temp2)
temp4<-merge(temp3,temp1,by="geneid")
temp5<-temp4[which(temp4$padj<5e-2),]
temp6<-merge(temp5,ki,by="geneid")
non_ribo_genes<-temp6

genes<-as.vector(non_ribo_genes$`ki(Da/s)`)
names(genes)<-non_ribo_genes$geneid
genes<-sort(genes,decreasing = T)
res_gsea_ki<-gseGO(genes,OrgDb="org.Sc.sgd.db", keyType = "ENSEMBL",
                   nPerm = 100000,pAdjustMethod = "BH",
                   pvalueCutoff = 0.25)
temp1<-res_gsea_ki@result
temp2<-temp1[which(temp1$pvalue<0.05),]


genes<-as.vector(non_ribo_genes$log2FoldChange)
names(genes)<-non_ribo_genes$geneid
genes<-sort(genes,decreasing = T)
res_gsea_log2fc<-gseGO(genes,OrgDb="org.Sc.sgd.db", keyType = "ENSEMBL",
                       nPerm = 100000,pAdjustMethod = "BH",
                       pvalueCutoff = 0.25)

temp3<-res_gsea_log2fc@result
temp4<-temp3[which(temp3$pvalue<0.05),]
temp5<-merge(temp2,temp4,by="ID")
gsea_res<-temp5
temp2<-vector()
for (i in 1:nrow(gsea_res)) 
{
  temp2<-c(temp2,paste0(gsea_res$ID[i]," ",gsea_res$Description.x[i])) 
}
gsea_res$Name<-temp2

temp6<-gsea_res[which(gsea_res$NES.x>0&gsea_res$NES.y>0),]
temp6$class<-rep("Response to stress",nrow(temp6))
temp7<-gsea_res[which(gsea_res$NES.x<0&gsea_res$NES.y<0),]
temp7$class<-rep("Response to stress",nrow(temp7))

####bar plot
temp1<-c(log10(gsea_res$p.adjust.x),-log10(gsea_res$p.adjust.y))
temp3<-data.frame(ID=rep(gsea_res$Name,2),log10padj=temp1)
temp3$datasource<-c(rep("ki",6),rep("log2fc",6))


ggplot()+
  geom_bar(data = temp3,aes(x = reorder(ID,log10padj),y = log10padj,fill = datasource),
           stat = 'identity',position = 'stack')+
  geom_bar(data = temp3,aes(x =reorder(ID,log10padj),y = log10padj,fill = datasource),
           stat = 'identity',position = 'stack')+
  scale_y_continuous(limits = c(-3.5,3.5),breaks = c(-3.5,0,3.5),labels = c(3.5,0,3.5))+
  theme_bw()+
  theme(legend.position = 'bottom',
        legend.direction = 'horizontal',
        axis.text = element_text(size = 25,hjust = 0,color = 'black'),
        axis.title = element_text(size = 25,hjust = 0,color = 'black'),
        legend.text = element_text(size = 20,hjust = 15,color = 'black'),
        legend.title = element_blank())+
  xlab('')+
  ylab('Log')+
  coord_flip() 
ggsave("gsea_da_padj.eps",width = 20,height = 5)

####scatter plot
ggplot(data=gsea_res, aes(x=NES.x, y=NES.y,color=Name))+
  geom_point(size=5)+
  geom_vline(xintercept=0,lty=2,col="black",lwd=0.8) +
  geom_hline(yintercept =0,lty=2,col="black",lwd=0.8) +
  labs(x="NES of ki",
       y="NES of log2(fold change)")+
  scale_y_continuous(limits = c(-3,3))+
  scale_x_continuous(limits = c(-3,3))+
  theme_bw()+
  theme(axis.text = element_text(size = 25,hjust = 0,color = 'black'),
        axis.title = element_text(size = 25,hjust = 0,color = 'black'),
        legend.text = element_text(size = 20,hjust = 0,color = 'black'),
        legend.title = element_blank())
ggsave("gsea_da_nes.eps",width = 12
       ,height = 5)

