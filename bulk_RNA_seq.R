library(openxlsx)
library(limma)
####screen for DEGs####
adjPfilter=0.01
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

qx=as.numeric(quantile(data, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  data[data<0]=0
  data=log2(data+1)}
data=normalizeBetweenArrays(data)
data<-data[rowMeans(data)>0,]

sampleName1<-c('NC1','NC2','NC3')
sampleName2<-c('H1','H2','H3')
conData=data[,sampleName1]
treatData=data[,sampleName2]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

allDiff=topTable(fit2,adjust='fdr',number=200000)
logFCfilter=mean(abs(allDiff$logFC))+2*sd(abs(allDiff$logFC))
diffSig=allDiff[with(allDiff, (abs(logFC)>logFCfilter & adj.P.Val < adjPfilter )), ]

library(dplyr)
library(ggplot2)
library(ggrepel)
adj.P.Val.Filter=0.01      

rt = read.table(inputFile, header=T, sep="\t", check.names=F)
rownames(rt)=rt$id
logFCfilter=mean(abs(rt$logFC))+2*sd(abs(rt$logFC))
Sig=ifelse((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter), ifelse(rt$logFC>logFCfilter,"Up","Down"), "Not")

rt = mutate(rt, Sig=Sig)
p = ggplot(rt, aes(logFC, -log10(adj.P.Val)))+
  geom_point(aes(col=Sig))+
  scale_color_manual(values=c("green", "black","red"))+
  labs(title = " ")+
  geom_vline(xintercept=c(-logFCfilter,logFCfilter), col="blue", cex=0.5, linetype=4)+
  geom_hline(yintercept= -log10(adj.P.Val.Filter), col="blue", cex=0.5, linetype=4)+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p1=p+geom_label_repel(data=filter(rt, ((rt$adj.P.Val<adj.P.Val.Filter) & (abs(rt$logFC)>logFCfilter))),
                      box.padding=0.1, point.padding=0.1, min.segment.length=0.05,
                      size=1.8, aes(label=id)) + theme_bw()

####GO enrichment####
library("clusterProfiler")
library(org.Hs.eg.db)
library("enrichplot")
library("ggplot2")
library(GOplot)
library(circlize)
library(RColorBrewer)
library(dplyr)
library(ComplexHeatmap)
rt$entrezID <- mget(rt$Symbol, org.Hs.egSYMBOL2EG, ifnotfound=NA)
rt=rt[is.na(rt[,"entrezID"])==F,]
gene=rt$entrezID

go <- enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =0.05, qvalueCutoff = 0.05,ont="all",readable =T)
GO=as.data.frame(go)
barplot(go, drop = TRUE, showCategory =10,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
enrichplot::dotplot(kk, showCategory=10, orderBy="GeneRatio", split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')

####KEGG enrichment####
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
rt$entrezID <- mget(rt$Symbol, org.Hs.egSYMBOL2EG, ifnotfound=NA) 
rt=rt[is.na(rt[,"entrezID"])==F,]                           
gene=rt$entrezID

#gene=intersect(geneList$`Agrimol B`,geneList$DGEs)
gene=unique(intersect(geneList$`DGEs_Agrimol B`,geneList$Autophagy))
gene=mget(gene, org.Hs.egSYMBOL2EG, ifnotfound=NA) 
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =0.05)   #富集分析
KEGG=as.data.frame(kk)
barplot(kk, drop = TRUE, showCategory = 20)
enrichplot::dotplot(kk, showCategory = 20)
####GSEA analyze####
library(clusterProfiler)
library(enrichplot)

rt<-rt[order(rt$logFC,decreasing = T),]
rt$entrezIDs <- mget(rt$id, org.Hs.egSYMBOL2EG, ifnotfound=NA) 
rt=rt[is.na(rt[,"entrezIDs"])==F,]

geneFC=rt$logFC
gene=rt$entrezIDs
names(geneFC)=gene
kk <- gseKEGG(geneFC, organism = "hsa",pvalueCutoff=1)
sortkk<-kk[order(kk$enrichmentScore,decreasing=T)]

gseaplot2(kk,geneSetID = c('hsa04140'),title = 'Autophagy - animal')
gseaplot2(kk,geneSetID = c('hsa04010'),title = 'MAPK signaling pahtway')

#####GSVA analyze#####
library(GSEABase)
library(GSVA)
library(msigdbr)
library(openxlsx)

rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

qx=as.numeric(quantile(data, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  data[data<0]=0
  data=log2(data+1)
}
data=normalizeBetweenArrays(data)
data<-data[rowMeans(data)>0,]
data=data[which(apply(data,1,function(x){return(sum(x>0))})>ncol(data)*1/2),]

genesets<-c()
genesets <- msigdbr(species = "Homo sapiens",
                    category = "C5") %>%
  dplyr::select(.,"gs_name","gene_symbol")%>% 
  as.data.frame()

genesets <- split(genesets$gene_symbol, genesets$gs_name)

es.max<-gsva(data,genesets,mx.diff=FALSE,verbose=TRUE,parallel.sz=5,method='gsva'
             ,kcdf='Poisson')
groups<-c('Agrimol_B','Contr')
group_list<-factor(c(1,1,1,2,2,2),levels = c(1,2),labels = c(groups))
library(limma)
design<-model.matrix(~0+factor(group_list))
colnames(design)<-levels(factor(group_list))
DEGgroup<-list()
label=c()
contrast.matrix<-makeContrasts('Agrimol_B-Contr',levels = design)
fit<-lmFit(es.max,design)
fit2<-contrasts.fit(fit,contrast.matrix)
fit2<-eBayes(fit2)
res<-decideTests(fit2, p.value=0.05)
summary(res)
tmp<-topTable(fit2,adjust='fdr',number = 200)
DEG<-na.omit(tmp)
DEG$KEGG=rownames(DEG)
a=data.frame(es.max)
a$KEGG=rownames(a)
library(dittoSeq)
Type=c(rep("Controll",3),rep("Agrimol_B",3))
es.max=cbind(es.max[,4:6],es.max[,1:3])
names(Type)=colnames(es.max)
Type=as.data.frame(Type)
pheatmap::pheatmap( es.max, 
                    annotation=Type, 
                    color = colorRampPalette(c("green2", "white", "red2"))(50),
                    annotation_colors = list(Type=c(Controll='#5CB85C',Agrimol_B='#AB201D')),
                    cluster_cols =F,
                    cluster_rows = T,
                    border_color='grey',
                    show_colnames = T,
                    treeheight_row = 0,
                    scale="row",
                    fontsize = 8,
                    fontsize_row=8,
                    fontsize_col=8
                    #,annotation_row = annotation_row
                    #,main = 'Integrin family genes'
)
#####TARGET dataset####

library(limma)
library(ggpubr)

rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp), colnames(exp))
data=matrix(as.numeric(as.matrix(exp)), nrow=nrow(exp), dimnames=dimnames)
rt=avereps(data)

qx=as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC=( (qx[5]>100) || ( (qx[6]-qx[1])>50 && qx[2]>0) )
if(LogC){
  rt[rt<0]=0
  rt=log2(rt+1)}

group=sapply(strsplit(colnames(data),"\\-"),"[",4)
group=sapply(strsplit(group,""),"[",1)
group=gsub("2", "1", group)
conNum=length(group[group=='S'])       
treatNum=length(group[group==0])    

Type=c(rep("con",conNum), rep("treat",treatNum))
geneRT=c('TYMS','DHFR','HSPA8')
data=data[geneRT,,drop=F]

Type=c(rep("Con",conNum), rep("Treat",treatNum))
my_comparisons=list()
my_comparisons[[1]]=levels(factor(Type))

newGeneLists=c()
outTab=data.frame()
for(i in row.names(data)){
  rt1=data.frame(expression=data[i,], Type=Type)
  
  boxplot=ggboxplot(rt1, x="Type", y="expression", color="Type",
                    xlab="",
                    ylab=paste(i, "expression"),
                    legend.title="",
                    palette = c("blue", "red"),
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  
   print(boxplot)
}

library(tidyverse)
library(ggpubr)
library(ggsci)
library(introdataviz)
library(data.table)
library(ggunchained)
df <- melt(data)
names(df)=c('genes','sample','Exp')
groups=gsub('GTEX','normal',sapply(strsplit(as.character(df$sample),"\\-"),"[",1)) 
groups=gsub('TARGET','OS',groups)
df$type2=groups
p=ggplot(df,aes(x = genes,y = Exp,fill = type2)) +
  # split violin
  geom_split_violin(alpha = .5, trim = F,color = NA,width = 1) +
  # mean point
  stat_summary(fun = "mean", geom = "point",position = position_dodge(0.2)) +
  # errorbar
  stat_summary(fun.data = "mean_sd", geom = "errorbar", width = .15,
               size = 0.3,
               position = position_dodge(0.2)) +
  theme_bw(base_size = 16) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 0,color = 'black',hjust = 0.5),
        axis.text.y = element_text(angle = 0,color = 'black',hjust = 0.5),
        legend.position = 'right',legend.direction = 'vertical',legend.title = element_blank()) +
  scale_fill_brewer(palette = 'Set1',direction = -1) +
  #scale_fill_jco(palette = c('blue','red'),alpha = .5) +
  ylim(-.5,18) +labs(x='',y='expression')+
  stat_compare_means(aes(group=type2),
                     symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                      symbols = c("***", "**", "*", "NS")),label = "p.signif",
                     label.y = 17,size = 5,method = 'anova')