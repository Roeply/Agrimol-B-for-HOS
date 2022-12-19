library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(stringr)
OS_dataSt<-function(Xdir,proJ=''){
  counts <- Read10X(data.dir = Xdir)
  OS_data=CreateSeuratObject(counts = counts,project = proJ, names.delim = "_")
  return(OS_data)
}
OS_data.list<-c()
for (i in c('sample1','sample2','sample3','sample4','sample5','sample6')) {
  OS_data1<-OS_dataSt(Xdir = paste0('data/SC/sample/',i),proJ = paste0('OS_',i))
  OS_data.list[[paste0('OS_',i)]]<-OS_data1
}

OS_data<-merge(x = OS_data.list[[1]], 
            y = OS_data.list[c(2:length(OS_data.list))], 
            add.cell.ids = names(OS_data.list)) 

OS_data[["percent.mt"]] <- PercentageFeatureSet(object = OS_data, pattern = "^MT-")
str_subset(string = rownames(OS_data),pattern = "MT-")

VlnPlot(object = OS_data,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident",log = T,pt.size = 0.1)
RidgePlot(object = OS_data,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),# log = T,
          ncol = 1,group.by = "orig.ident")
p1 <- FeatureScatter(object = OS_data,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "orig.ident")
p2 <- FeatureScatter(object = OS_data,feature1 = "nFeature_RNA",feature2 = "percent.mt",group.by = "orig.ident")

OS_data=subset(x = OS_data, subset = nFeature_RNA > 300&nFeature_RNA < 4500 & percent.mt < 10 & percent.mt > -Inf)    #对数据进行过滤

rm(list=ls()[which(ls() != 'OS_data')])
gc()

OS_data.list<-SplitObject(OS_data,split.by = 'orig.ident')
for (i in 1:length(OS_data.list)) {
  a<-OS_data.list[[i]]
  a<-NormalizeData(object = a, normalization.method = "LogNormalize", scale.factor = 10000)
  a <- FindVariableFeatures(object = a, selection.method = "vst", nfeatures = 4000)
  OS_data.list[[i]]<-a
}
AnchorSet<-FindIntegrationAnchors(object.list = OS_data.list,normalization.method = 'LogNormalize',anchor.features = 4000)
AnchorSet
gc()

OS_data<-IntegrateData(anchorset = AnchorSet)

OS_data<-ScaleData(OS_data,features = rownames(OS_data))

VlnPlot(object = OS_data,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),group.by = "orig.ident",log = T,pt.size = 0.1)
RidgePlot(object = OS_data,features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),# log = T,
          ncol = 1,group.by = "orig.ident")
p1 <- FeatureScatter(object = OS_data,feature1 = "nCount_RNA",feature2 = "nFeature_RNA",group.by = "orig.ident")
p2 <- FeatureScatter(object = OS_data,feature1 = "nFeature_RNA",feature2 = "percent.mt",group.by = "orig.ident")

OS_data <- FindVariableFeatures(object = OS_data, selection.method = "vst", nfeatures = 4000)
top10 <- head(x = VariableFeatures(object = OS_data), 10)
plot1 <- VariableFeaturePlot(object = OS_data)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))

library(projectLSI)
OS_data.lsi <- calcLSI(OS_data[["RNA"]]@data[VariableFeatures(OS_data), ],ndims = 50)
OS_data[["pca"]] <- CreateDimReducObject(
  embeddings = OS_data.lsi$matSVD,
  loadings = OS_data.lsi$fLoad,
  assay = "RNA",
  stdev = OS_data.lsi$sdev,
  key = "PC_")
ElbowPlot(OS_data,ndims = 50)#18
DimHeatmap(object = OS_data,reduction = "pca",dims = c(1:4),nfeatures = 30,disp.min = -2.5,balanced = TRUE,projected = FALSE,
           fast = TRUE,raster = TRUE,slot = "scale.data",combine = TRUE)
DimPlot(object = OS_data,reduction = "pca",group.by = "orig.ident",dims = c(1,2),shuffle = TRUE,
        label = TRUE,label.size = 4,label.color = "black",label.box = TRUE,sizes.highlight = 1)
pcSelect=5

OS_data <- RunTSNE(OS_data,dims = 1:pcSelect)

umap.OS_data=uwot::umap(OS_data.lsi$matSVD[,1:pcSelect],
                     n_neighbors = 30,
                     min_dist = 0.3,search_k = 3000,
                     metric = "cosine",
                     ret_model = T,
                     verbose = T)

umap.OS_data.emb <- umap.OS_data$embedding
rownames(umap.OS_data.emb) <- colnames(OS_data)
colnames(umap.OS_data.emb) <- paste0("UMAP_", seq_len(ncol(umap.OS_data.emb)))

OS_data[["umap"]] <- CreateDimReducObject(
  embeddings = umap.OS_data.emb,
  assay = "RNA",
  key = "UMAP_")

OS_data <- FindNeighbors(object = OS_data, dims = 1:pcSelect)   

library(clustree)
DefaultAssay(OS_data)<-'RNA'
obj <- FindClusters(OS_data, resolution = seq(0,0.5,by=0.05))
clustree(obj)
OS_data <- FindClusters(object = OS_data, resolution = 0.15)  #0.3      
table(OS_data@meta.data[["seurat_clusters"]])
UMAPPlot(object = OS_data,pt.size = .4,label = TRUE,group.by = "seurat_clusters")  
TSNEPlot(object = OS_data, pt.size = .4, label = TRUE,group.by = "seurat_clusters")   

logFCfilter=1               
adjPvalFilter=0.05     
OS_data.markers <- FindAllMarkers(object = OS_data,
                               only.pos = FALSE,
                               min.pct = 0.1,
                               logfc.threshold = logFCfilter)
sig.markers=OS_data.markers[(abs(as.numeric(as.vector(OS_data.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(OS_data.markers$p_val_adj))<adjPvalFilter),]
top10 <- OS_data.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

library("Nebulosa")
plot_density(OS_data, showGenes)


OS_data.markers=FindAllMarkers(object = OS_data,
                            only.pos = FALSE,
                            min.pct = 0.25,
                            logfc.threshold = logFCfilter)
sig.cellMarkers=OS_data.markers[(abs(as.numeric(as.vector(OS_data.markers$avg_log2FC)))>logFCfilter & as.numeric(as.vector(OS_data.markers$p_val_adj))<adjPvalFilter),]
write.table(sig.cellMarkers,file="04.cellMarkers.txt",sep="\t",row.names=F,quote=F)

newLabels=c('Myeloid cells','Osteoblasts','NK/T cells'
            ,'Osteoclasts','Myeloid cells','Fibroblasts','Plasma cells','Endothelial cells','B cells','NK/T cells')
names(newLabels)=levels(OS_data@meta.data[["seurat_clusters"]])
celltype=factor( OS_data@meta.data[["seurat_clusters"]],labels  = newLabels)
OS_data@meta.data[["celltype"]]=celltype
OS_data=RenameIdents(OS_data, newLabels)
UMAPPlot(object = OS_data, pt.size = .4, label = TRUE,group.by='celltype') 
TSNEPlot(object = OS_data, pt.size = .4, label = TRUE,group.by='celltype') 

pie(table(OS_data$celltype),labels = paste(names(table(OS_data$celltype)),"-",round(100*table(OS_data$celltype)/sum(table(OS_data$celltype)),2), "%"))

####inferCNV####

inferCNVData=subset(inferCNVData,id=c('Osteoblastic OS cells','Myeloid cells','NK/T cells','Plasma cells','B cells'))

#inferCNVData=OS_data
library(infercnv)
library(AnnoProbe)
geneInfor=annoGene(rownames(inferCNVData),"SYMBOL",'human')
colnames(geneInfor)
geneInfor=geneInfor[with(geneInfor, order(chr, start)),c(1,4:6)]
geneInfor=geneInfor[!duplicated(geneInfor[,1]),]
length(unique(geneInfor[,1]))
head(geneInfor)
dat=inferCNVData[rownames(inferCNVData) %in% geneInfor[,1],]
dat=dat[match( geneInfor[,1], rownames(dat) ),] 
dim(dat)
counts=GetAssayData(dat,slot='counts',assay='RNA')
groupinfo=dat@meta.data
groupinfo[['cell']]=rownames(groupinfo)
groupinfo=cbind(groupinfo$cell,as.character(groupinfo$cell_type))
head(groupinfo)
write.table(groupinfo,file = groupFiles,sep = '\t',quote = F,col.names = F,row.names = F)
head(geneInfor)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=counts,
                                    annotations_file='output_file/groupFiles.txt',
                                    delim="\t",
                                    gene_order_file= 'output_file/geneFile.txt',
                                    chr_exclude = c("chrX", "chrY", "chrM") ,
                                    ref_group_names=c('Myeloid cells','NK/T cells','Plasma cells','B cells'))  
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir='temp/', num_threads = 6,
                             cluster_by_groups=T, save_rds = T,
                             denoise=T,plot_steps = F,
                             HMM=T,output_format = 'pdf')




library(Seurat)
library(projectLSI)
library(patchwork)
bulk.data=read.table('data/genes.txt',row.names = 1,header = T, sep="\t", check.names=F)[,4:6]
colnames(bulk.data)=c('HOS_1','HOS_2','HOS_3')
psudo.all <- psudoSC(bulk.data, n = 1, depth = 3000)
psudo.so <- CreateSeuratObject(psudo.all, project = "bulk")
psudo.so <- NormalizeData(psudo.so)

bulk.matSVD <- projectLSI(psudo.so[["RNA"]]@data, OS_data.lsi)
umap.OS_data1=uwot::umap(OS_data.lsi$matSVD[,1:5],
                      n_neighbors = 30,
                      min_dist = 0.3,search_k = 3000,
                      metric = "cosine",
                      ret_model = T,
                      verbose = T)
umap.OS_data1$embedding=umap.OS_data$embedding

umap.bulk.proj <- uwot::umap_transform(bulk.matSVD[,1:5], umap.OS_data, verbose = T)

DimPlot(OS_data, label = T)

rownames(umap.bulk.proj) <- colnames(psudo.so)
colnames(umap.bulk.proj) <- paste0("UMAP_", seq_len(ncol(umap.bulk.proj)))
psudo.so[["umap"]] <- CreateDimReducObject(
  embeddings = umap.bulk.proj,
  assay = "RNA",
  key = "UMAP_")

DimPlot(psudo.so, label = T)

OS_data.mix <- merge(OS_data, psudo.so)

OS_data.mix[["umap"]] <- CreateDimReducObject(
  embeddings = rbind(OS_data[["umap"]]@cell.embeddings,
                     psudo.so[["umap"]]@cell.embeddings),
  assay = "RNA", key = "UMAP_")

DimPlot(OS_data.mix, label = T)

library(ggunchull)
library(ggsci)
library(ggrepel)

bulk.p=as.data.frame(psudo.so[["umap"]]@cell.embeddings)
rownames(bulk.p)=c('HOS_1','HOS_2','HOS_3')

UMAPPlot(OS_data, label = T,label.size=0,pt.size=.4,group.by='celltype')+#scale_color_npg()+
  geom_point(shape = 23, colour = "black",fill = "red", size = 2,aes(x=UMAP_1,y=UMAP_2),data = bulk.p)+
  geom_text_repel(aes(label=rownames(bulk.p),x=UMAP_1,y=UMAP_2),data = bulk.p,size=3,fontface = 'bold')
UMAPPlot(proC, label = T,label.size=0,pt.size=.4,group.by='celltype')+#scale_color_npg()+
  geom_point(shape = 23, colour = "black",fill = "red", size = 2,aes(x=UMAP_1,y=UMAP_2),data = bulk.p)+
  geom_text_repel(aes(label=rownames(bulk.p),x=UMAP_1,y=UMAP_2),data = bulk.p,size=3,fontface = 'bold')
