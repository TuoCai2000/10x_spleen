library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(cowplot)
library(scater)
library(Hmisc)
tmp=read.table("/Users/wangchenchen/Documents/spleen/pcgenes.txt",header=F,sep="\t")





######preprocess
mi1d1 = Read10X(data.dir = "/share/home/wangchenchen/data/10x/sp/mi1d1_report/outs/filtered_feature_bc_matrix")
class(mi1d1)
mi1d1 = CreateSeuratObject(counts = mi1d1,min.cells = 8,min.features = 200,project = "mi")
mi1d1@assays$RNA@counts = as.matrix(mi1d1@assays$RNA@counts)
mi1d1@assays$RNA@counts = mi1d1@assays$RNA@counts[intersect(row.names(mi1d1@assays$RNA@counts),tmp$V2),]
mi1d2 = Read10X(data.dir = "/share/home/wangchenchen/data/10x/sp/mi1d2_report/outs/filtered_feature_bc_matrix")
class(mi1d2)
mi1d2 = CreateSeuratObject(counts = mi1d2,min.cells = 11,min.features = 200,project = "mi")
mi1d2@assays$RNA@counts = as.matrix(mi1d2@assays$RNA@counts)
mi1d2@assays$RNA@counts = mi1d2@assays$RNA@counts[intersect(row.names(mi1d2@assays$RNA@counts),tmp$V2),]
mi3d1 = Read10X(data.dir = "/share/home/wangchenchen/data/10x/sp/mi3d1_report/outs/filtered_feature_bc_matrix")
class(mi3d1)
mi3d1 = CreateSeuratObject(counts = mi3d1,min.cells = 8,min.features = 200,project = "mi")
mi3d1@assays$RNA@counts = as.matrix(mi3d1@assays$RNA@counts)
mi3d1@assays$RNA@counts = mi3d1@assays$RNA@counts[intersect(row.names(mi3d1@assays$RNA@counts),tmp$V2),]
mi3d2 = Read10X(data.dir = "/share/home/wangchenchen/data/10x/sp/mi3d2_report/outs/filtered_feature_bc_matrix")
class(mi3d2)
mi3d2 = CreateSeuratObject(counts = mi3d2,min.cells = 8,min.features = 200,project = "mi")
mi3d2@assays$RNA@counts = as.matrix(mi3d2@assays$RNA@counts)
mi3d2@assays$RNA@counts = mi3d2@assays$RNA@counts[intersect(row.names(mi3d2@assays$RNA@counts),tmp$V2),]
sham1 = Read10X(data.dir = "/share/home/wangchenchen/data/10x/sp/sham1_report/outs/filtered_feature_bc_matrix")
class(sham1)
sham1 = CreateSeuratObject(counts = sham1,min.cells = 8,min.features = 200,project = "sham")
sham1@assays$RNA@counts = as.matrix(sham1@assays$RNA@counts)
sham1@assays$RNA@counts = sham1@assays$RNA@counts[intersect(row.names(sham1@assays$RNA@counts),tmp$V2),]
sham2 = Read10X(data.dir = "/share/home/wangchenchen/data/10x/sp/sham1_report/outs/filtered_feature_bc_matrix")
class(sham2)
sham2 = CreateSeuratObject(counts = sham2,min.cells = 8,min.features = 200,project = "sham")
sham2@assays$RNA@counts = as.matrix(sham2@assays$RNA@counts)
sham2@assays$RNA@counts = sham2@assays$RNA@counts[intersect(row.names(sham2@assays$RNA@counts),tmp$V2),]
sp <- merge(mi1,y=c(mi2,sham1,sham2),add.cell.ids = c("mi1","mi2","sham1","sham2"),project = "spleen")
sp
table(sp$orig.ident)
rm(sham1)
rm(sham2)
rm(mi1)
rm(mi2)





#####quality control
#VlnPlot(sp, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
sp[["percent.mt"]] <- PercentageFeatureSet(sp, pattern = "^mt-")
size.drop <- scater::isOutlier(sp$nCount_RNA, nmads=5, type="both", log=F,min_diff = NA)
gene.drop <- scater::isOutlier(sp$nCount_RNA, nmads=5, type="both", log=F,min_diff = NA)
mito.drop <- scater::isOutlier(sp$percent.mt, nmads=5, type="both", log=F,min_diff = NA)
sp <- sp[,!(size.drop | gene.drop | mito.drop)]
#####integrate data
sp =as.data.frame(sp@assays$RNA@counts)
sp = CreateSeuratObject(counts = sp,project = "sp")
sp@meta.data$stim[sp@meta.data$orig.ident %in% c("mi1","mi2")] = "mi"
sp@meta.data$stim[sp@meta.data$orig.ident %in% c("sham1","sham2")] = "sham"
######
sp.list <- SplitObject(sp,split.by = "stim")
sp.list <- lapply(X = sp.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
#sp.list$mi@meta.data$sample = substr(row.names(sp.list$mi@meta.data),1,3)
#sp.list$sham@meta.data$sample = substr(row.names(sp.list$sham@meta.data),1,5)
sp.anchors <- FindIntegrationAnchors(object.list = sp.list , anchor.features = 3000,dims = 1:30,  k.anchor = 10)
sp.combined <- IntegrateData(anchorset = sp.anchors, dims = 1:30)
DefaultAssay(sp.combined) <- "integrated"
rm(sp.anchors)
rm(sp.list)

#####dimensional reduction
# Run the standard workflow for visualization and clustering
sp.combined <- ScaleData(sp.combined, verbose = FALSE)
sp.combined <- RunPCA(sp.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
sp.combined <- RunUMAP(sp.combined, reduction = "pca", dims = 1:21)
sp.combined <- FindNeighbors(sp.combined, reduction = "pca", dims = 1:21)
sp.combined <- FindClusters(sp.combined, resolution = 0.245)
# Visualization
p1 <- DimPlot(sp.combined, reduction = "umap",group.by = "orig.ident")
p2 <- DimPlot(sp.combined, reduction = "umap", label = TRUE)
plot_grid(p1)
plot_grid(p2)
DimPlot(sp.combined, reduction = "umap", split.by = "orig.ident")
DimPlot(sp.combined, reduction = "umap", split.by = "stim")





######find markers and celltype annotation
sp.markers <- FindAllMarkers(sp.combined, only.pos = TRUE, min.pct = 0.4, logfc.threshold = 0.5)
top10 <- sp.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
top
DimPlot(sp.combined, label = TRUE)
#FeaturePlot(sp.combined, features = c("Gata2","Itga2b","Gp1bb","Klf1","Gata1","Epor","Ms4a2","Ncr1","C1qb","Flt3"),cols = c("grey","yellow", "red"))
#FeaturePlot(sp.combined, features = c("Irf8","Csf1r","Ly86","Gfi1","Cebpe"),cols = c("grey","yellow", "red"))
#FeaturePlot(sp.combined, features = c("Gata3","Id2","Bcl11b","Ccl5"),cols = c("grey","yellow", "red"))
FeaturePlot(sp.combined, features = c("Procr","Ly6a","Car1","Gata1","Epor","Klf1"),cols = c("grey","yellow", "red"))
FeaturePlot(sp.combined, features = c("Mpo","Ctsg","Elane","H2-Ab1","Cd79a","H2-DMb2","H2-Eb1","H2-Aa"),cols = c("grey","yellow", "red"))
FeaturePlot(sp.combined, features = c("Ncr1","Klrk1","Klrd1","Ncr1","Nkg7","Ms4a2","Gzmb","Cpa3","Mcpt8"),cols = c("grey","yellow", "red"))
FeaturePlot(sp.combined, features = c( "Dntt","Notch1","Ly6d",'Irf8','Ly86','Csf1r',"C1qb","C1qa","Cd68"),cols = c("grey","yellow", "red"))
FeaturePlot(sp.combined, features = c("Pf4","Chil3","Plac8","Saa3", "Arg1","Cd74"),cols = c("grey","yellow", "red"))
#FeaturePlot(sp.combined, features = c("H2-Eb1","H2-Aa","H2-Ab1"),cols = c("grey","yellow", "red"))
#FeaturePlot(sp.combined, features = c("Bcl2a1a","Bcl2a1d","Bcl2a1b"),cols = c("grey","yellow", "red"))
eryp = as.matrix(sp.combined@assays$RNA@counts)[,sp.combined@meta.data$seurat_clusters %in% c("2","3","1","0","9")]
hscmpp = as.matrix(sp.combined@assays$RNA@counts)[,sp.combined@meta.data$seurat_clusters %in% c("5")]
mkp = as.matrix(sp.combined@assays$RNA@counts)[,sp.combined@meta.data$seurat_clusters %in% c("6")]
gmp = as.matrix(sp.combined@assays$RNA@counts)[,sp.combined@meta.data$seurat_clusters %in% c("7")]
clp = as.matrix(sp.combined@assays$RNA@counts)[,sp.combined@meta.data$seurat_clusters %in% c("8")]
write.csv(clp,"/Users/wangchenchen/Documents/spleen/data/clp.csv")
erypd = sp.combined@meta.data[sp.combined@meta.data$seurat_clusters %in% c("2","3","1","0","9"),]
write.csv(erypd,"/Users/wangchenchen/Documents/spleen/data/erypd.csv")
hscdown4= as.matrix(sp.combined@assays$RNA@counts)[,sp.combined@meta.data$seurat_clusters %in% c("5","6","7","8","3","1","0","9")]
write.csv(hscdown4,"/Users/wangchenchen/Documents/spleen/data/hscdown4.csv")
hscdown4pd = sp.combined@meta.data[sp.combined@meta.data$seurat_clusters %in% c("5","6","7","8","3","1","0","9"),]
write.csv(hscdown4pd,"/Users/wangchenchen/Documents/spleen/data/hscdown4pd.csv")
write.csv(sp.combined@assays$RNA@counts,"/Users/wangchenchen/Documents/spleen/data/sp.csv")
write.csv(sp.combined@meta.data,"/Users/wangchenchen/Documents/spleen/data/spmeta.csv")
###cell ratio && degs
sp.combined <- RenameIdents(sp.combined,`0` = "EryP2", `1` = "EryP3", `2` = "MEP", 
    `3` = "EryP4", `4` = "MC/BA", `5` = "HSC/MPP", `6` = "MkP", `7` = "GMP", `8` = "CLP", `9` = "EryP1", 
   `10` = "NK_pro", `11` = "Mono.D_pro",`12`="B_pro",`13`="Macrophage")
##sp.combined <- RenameIdents(sp.combined,`GMP` = "NeuP")
DimPlot(sp.combined, label = TRUE)
Idents(sp.combined) = factor(Idents(sp.combined),levels = c("HSC/MPP", "MkP", "MEP", "EryP1", "EryP2", "EryP3", "EryP4",
                                                           "MC/BA","GMP","Macrophage","CLP","NK_pro","Mono.D_pro","B_pro"))
DotPlot(sp.combined, features = c("Procr","Ly6a","Pf4","Car1","Gata1","Epor","Klf1","Ms4a2","Gzmb","Cpa3","Mcpt8", 
                                  "Mpo","Ctsg","Elane","C1qb","C1qa","Cd68","Dntt","Notch1","Ly6d",
                                  "Ncr1","Klrk1","Klrd1","H2-Ab1","Cd74","H2-DMb2"
                                 ))+ RotatedAxis()
##data quality
ggplot(sp.combined@meta.data,aes(x=orig.ident,y=nFeature_RNA ,color =orig.ident ))+geom_jitter()+
  theme_bw() + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))
sp.combined$celltype =  Idents(sp.combined)
sp.combined$celltype = as.character(sp.combined@meta.data$celltype)
sp.combined@meta.data$celltype[sp.combined@meta.data$celltype %in% c("EryP1","EryP2","EryP3","EryP4")] <- "EryP"

sp.combined$newcelltype <- paste(sp.combined$celltype, sp.combined$stim, sep = "_")
#Idents(sp.combined) = sp.combined$celltype 





#####celltype composition(Ro/e analysis)
roe = read.csv("/Users/wangchenchen/Documents/spleen/fig/last1/sproe.csv")
names(roe)[1] = "cluster"
roe$cluster= factor(roe$cluster,levels = c("HSC/MPP", "MkP", "MEP", "EryP",
                                   "MC/BA","GMP","Macrophage","CLP","NK_pro","Mono.D_pro","B_pro"))
roe = roe[,c(1,4)]
roe =melt(roe)
ggplot(roe,aes(x=cluster,y=value,fill=cluster))+
  geom_bar(stat="identity",width = 0.6)+facet_wrap(~variable)+
  theme_bw() + 
  theme(panel.border = element_blank(),panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ RotatedAxis()

Idents(sp.combined) = sp.combined$newcelltype 
monoddiff <- FindMarkers(sp.combined, ident.1 = "Mono.D_pro_mi", ident.2 = "Mono.D_pro_sham", verbose = TRUE)
write.csv(nkdiff,"/Users/wangchenchen/Documents/spleen/fig/last1/diff/nkdiff.csv")





#####cell cycle analysis
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
s.genes = capitalize(tolower(s.genes))
g2m.genes = capitalize(tolower(g2m.genes))
sp.combined <- CellCycleScoring(sp.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
mimeta = sp.combined@meta.data[sp.combined@meta.data$stim %in% "mi",]
shammeta = sp.combined@meta.data[sp.combined@meta.data$stim %in% "sham",]
ggplot(shammeta,aes(x=stim,fill=Phase))+
  geom_bar(position = "fill",stat = "count")+coord_polar("y")+labs(x="")+scale_fill_brewer(palette="Dark2")+
  theme_minimal()+
  theme(axis.title=element_blank(),axis.ticks=element_blank(),axis.text = element_blank(),legend.title = element_blank())+
  theme(panel.grid=element_blank()) +   
  theme(panel.border=element_blank()) +facet_wrap(~celltype)+ggtitle("sham cell cycle")




######hsc subcluster
library(dplyr)
library(Seurat)
library(patchwork)
library(ggpubr)
library(ggbeeswarm)
library(RColorBrewer)
library(pheatmap)
library(reshape2)
library(Hmisc)
hscmpp = as.matrix(sp.combined@assays$RNA@counts)[,sp.combined$stim %in% c("mi1d","mi3d")&sp.combined$ident %in% c("HSC/MPP")]
hsc <- CreateSeuratObject(counts = hscmpp, project = "hsc")
hsc
hsc[["percent.mt"]] <- PercentageFeatureSet(hsc, pattern = "^mt-")
hsc@meta.data$days = substr(row.names(hsc@meta.data),4,7)
hsc$days[hsc$days %in% "m_sh"]="sham"
hsc$sample = substr(row.names(hsc@meta.data),4,8)
hsc.list <- SplitObject(hsc,split.by = "sample")
hsc.list <- lapply(X = hsc.list, FUN = function(x) {
  x <- subset(x, subset = nFeature_RNA > 500&nFeature_RNA<7000&percent.mt<20) 
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 1000)
})
features <- SelectIntegrationFeatures(object.list = hsc.list)
hsc.anchors <- FindIntegrationAnchors(object.list = hsc.list , anchor.features = features,dims = 1:30,  k.anchor = 10)
hsc.combined <- IntegrateData(anchorset = hsc.anchors, dims = 1:30)
DefaultAssay(hsc.combined) <- "integrated"
rm(hsc.anchors)
rm(hsc.list)
# Run the standard workflow for visualization and clustering
hsc.combined <- ScaleData(hsc.combined, verbose = FALSE)
hsc.combined <- RunPCA(hsc.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
hsc.combined <- RunUMAP(hsc.combined, reduction = "pca", dims = 1:20)
hsc.combined <- FindNeighbors(hsc.combined, reduction = "pca", dims = 1:20)
hsc.combined <- FindClusters(hsc.combined, resolution = 0.15)
DimPlot(hsc.combined, reduction = "umap", label = TRUE)
DimPlot(hsc.combined, reduction = "umap", split.by = "orig.ident")
DimPlot(hsc.combined, reduction = "umap", split.by = "days")
DefaultAssay(hsc.combined) <- "RNA"
hsc.combined <- ScaleData(hsc.combined, verbose = TRUE)
hsc.markers <- FindAllMarkers(hsc.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10 <- hsc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(hsc.combined, features = top10$gene) + NoLegend()
ggplot(hsc.combined@meta.data,aes(orig.ident,fill=seurat_clusters)) + geom_bar(position = "fill",width = 0.6)+
  scale_fill_manual(values=b3)+
  theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),
                   axis.text.x = element_text(size = 14),axis.text.y = element_text(size = 14),
                   panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+
  labs(x="",y="")
hsc.combined@meta.data = hsc.combined@meta.data[,c(1:5,7,9)]
#hsc@meta.data$stim[hsc@meta.data$orig.ident %in% c("mi1","mi2")] = "mi"
#hsc@meta.data$stim[hsc@meta.data$orig.ident %in% c("sham1","sham2")] = "sham"
FeaturePlot(hsc, features = c('Hoxb5','Fgd5','Procr',
      'Slamf1', 'Itga2b', 'Kit', 'Ly6a', 'Bmi1', 'Gata2', 'Hlf', 'Meis1', 
      'Mpl', 'Mcl1', 'Gfi1', 'Gfi1b'),cols = c("grey","yellow", "red"))





######geneset score
#type 1 interferon
hscmpp=subset(sp.combined,subset = ident=="HSC/MPP")
geneset=read.table("type_1_interferon_genelist.txt",sep='\t',header=F)
geneset=geneset$V1
geneset=as.data.frame(geneset)
geneset=as.list(geneset)
hscmpp <- AddModuleScore(hscmpp,features = geneset,name= 'score') 
table(hscmpp@meta.data$stim)
color_panel=c("sham"="#c12727",'mi1d'='#c5a716','mi3d'='#3a5ba8')
library(ggsignif)
library(ggpubr)
library(ggsignif)
library(ggforce)
ggplot(hscmpp@meta.data, aes(x=stim,y=score1,color=stim)) +
  geom_violin()+
  geom_sina(alpha=0.5,size=1)+
stat_summary(fun="mean", geom="point", shape=5, size=8, alpha=1)+
  geom_boxplot(width=0.6,position=position_dodge(0.8),size=1.2,fill="NA")+
  labs(y='',
       x='',
       title='type 1 interferon signaling pathway')+
  theme_classic()+
  scale_color_manual(values=color_panel)+
  #ylim(-0.1,0.3)+
  #scale_x_continuous(breaks=c(1,25,50),labels=custom_x_labels)+
  theme(plot.title=element_text(size=20,colour='black',hjust=0.5), #标题
        axis.title.y=element_text(size=15,color='black',vjust=3,hjust=0.5,angle=90), #y轴标题
        axis.title.x=element_text(size=15,color='black',hjust=0.5,vjust=-1), #x轴标题
        legend.title=element_text(color='black',size=15), #图例标题
        legend.text=element_text(color='black',size=15), #图例文字
        axis.text.x=element_text(size=15,color='black',vjust=0.5,hjust=0.5,angle=0), #x轴坐标轴标题
        axis.text.y=element_text(size=15,color='black',vjust=0.5,hjust=1,angle=0),
        plot.margin=unit(c(0,0,0,0),"in"))+ #y轴坐标轴标题
  coord_fixed(ratio=5)+
  geom_signif(mapping=aes(x=stim,y=score1),
                comparisons=list(c("sham","mi1d"),
                                c("sham","mi3d")),
                map_signif_level=TRUE, #是否显示显著性
                #tip_length=c(0,0,0,0,0,0), #修改显著性线两端的长短
                y_position=c(0.5,0.6), #设置显著性线的位置高度
                size=0.6, #线的粗细
                textsize=6, #修改显著性标记的大小
                test="wilcox.test") #检验的类型
ggsave("type1_interferon_signaling_pathway_box.pdf")
#isg score
geneset=read.table("isg_genelist.txt",sep='\t',header=F)
geneset=geneset$V1
geneset=as.data.frame(geneset)
geneset=as.list(geneset)
hsc.combined <- AddModuleScore(hsc.combined,features = geneset,name= 'isg') 
color_panel=c("0"="#e78411",'1'='#4ca64f','2'='#9793c7','3'='#d86f96')
library(ggsignif)
library(ggpubr)
library(ggsignif)
library(ggforce)
ggplot(hsc.combined@meta.data, aes(x=seurat_clusters,y=isg1,color=seurat_clusters)) +
  geom_violin()+
  geom_sina(alpha=0.5,size=1)+
  stat_summary(fun="mean", geom="point", shape=5, size=8, alpha=1)+
  geom_boxplot(width=0.6,position=position_dodge(0.8),size=1.2,fill="NA")+
  labs(y='',
       x='',
       title='ISG score')+
  theme_classic()+
  scale_color_manual(values=color_panel)+
  #ylim(-0.1,0.3)+
  #scale_x_continuous(breaks=c(1,25,50),labels=custom_x_labels)+
  theme(plot.title=element_text(size=20,colour='black',hjust=0.5), #标题
        axis.title.y=element_text(size=15,color='black',vjust=3,hjust=0.5,angle=90), #y轴标题
        axis.title.x=element_text(size=15,color='black',hjust=0.5,vjust=-1), #x轴标题
        legend.title=element_text(color='black',size=15), #图例标题
        legend.text=element_text(color='black',size=15), #图例文字
        axis.text.x=element_text(size=15,color='black',vjust=0.5,hjust=0.5,angle=0), #x轴坐标轴标题
        axis.text.y=element_text(size=15,color='black',vjust=0.5,hjust=1,angle=0),
        plot.margin=unit(c(0,0,0,0),"in"))+ #y轴坐标轴标题
  coord_fixed(ratio=7)+
  geom_signif(mapping=aes(x=seurat_clusters,y=isg1),
                comparisons=list(c("0","1"),
                                c("0","2"),
                                c("0","3"),
                                c("1","2")),
                map_signif_level=TRUE, #是否显示显著性
                #tip_length=c(0,0,0,0,0,0), #修改显著性线两端的长短
                y_position=c(0.5,0.55,0.6,0.65), #设置显著性线的位置高度
                size=0.6, #线的粗细
                textsize=6, #修改显著性标记的大小
                test="wilcox.test") #检验的类型
ggsave("isg_score_box.pdf")
#isg score projected on umap
hsc.combined@meta.data=hsc.combined@meta.data[order(hsc.combined@meta.data$isg1,decreasing = F),]
ggplot(hsc.combined@meta.data,aes(x = UMAP_1,y = UMAP_2,color = isg1))+
    geom_point(size=2)+
    theme_classic()+
    scale_color_continuous(name="ISG score",
                            limits=c(0,0.46714506910194),
                            breaks=c(0,0.2,0.4),
                            low="#e5e5e5",high='red',na.value="#e5e5e5")+
    theme(plot.title=element_text(size=20,colour='black',hjust=0.5), #标题
        axis.title.y=element_text(size=20,color='black',vjust=3,hjust=0.5,angle=90), #y轴标题
        axis.title.x=element_text(size=20,color='black',hjust=0.5,vjust=-2), #x轴标题
        legend.title=element_text(color='black',size=20), #图例标题
        legend.text=element_text(color='black',size=15), #图例文字
        axis.text.x=element_text(size=15,color='black',vjust=0.5,hjust=0.5,angle=0), #x轴坐标轴标题
        axis.text.y=element_text(size=15,color='black',vjust=0.5,hjust=1,angle=0), #y轴坐标轴标题
        plot.margin = unit(c(0,0,0,0), "in"),
        legend.key=element_blank(),#去除图例后面的灰色背景
        legend.spacing.y=unit(1,'cm'))+
    theme(legend.key.size=unit(30,'pt'))+
    coord_fixed(ratio = 1.6)
    #scale_color_gradient("red",low = "white",high = "red",limits=c(0,0.45))
    #scale_color_gradient2(low = "#e5e5e5",high = "#c21d1f",mid="#e39885",midpoint = 0.35)
ggsave("ISG_score_on_hsc_umap.pdf")
#stemness score
geneset=read.table("stemneds_genelist.txt",sep='\t',header=F)
geneset=geneset$V1
geneset=as.data.frame(geneset)
geneset=as.list(geneset)
hsc.combined <- AddModuleScore(hsc.combined,features = geneset,name= 'STEMNESS_SCORE') 
color_panel=c("0"="#e78411",'1'='#4ca64f','2'='#9793c7','3'='#d86f96')
library(ggsignif)
library(ggpubr)
library(ggsignif)
library(ggforce)
ggplot(hsc.combined@meta.data, aes(x=seurat_clusters,y=STEMNESS_SCORE1,color=seurat_clusters)) +
  geom_violin()+
  geom_sina(alpha=0.5,size=1)+
  geom_boxplot(width=0.6,position=position_dodge(0.8),size=1.2,fill="NA")+
  stat_summary(fun="mean", geom="point", shape=5, size=8, alpha=1)+
  labs(y='',
       x='',
       title='Stem score')+
  theme_classic()+
  scale_color_manual(values=color_panel)+
  #ylim(-0.1,0.3)+
  #scale_x_continuous(breaks=c(1,25,50),labels=custom_x_labels)+
  theme(plot.title=element_text(size=20,colour='black',hjust=0.5), #标题
        axis.title.y=element_text(size=15,color='black',vjust=3,hjust=0.5,angle=90), #y轴标题
        axis.title.x=element_text(size=15,color='black',hjust=0.5,vjust=-1), #x轴标题
        legend.title=element_text(color='black',size=15), #图例标题
        legend.text=element_text(color='black',size=15), #图例文字
        axis.text.x=element_text(size=15,color='black',vjust=0.5,hjust=0.5,angle=0), #x轴坐标轴标题
        axis.text.y=element_text(size=15,color='black',vjust=0.5,hjust=1,angle=0),
        plot.margin=unit(c(0,0,0,0),"in"))+ #y轴坐标轴标题
  coord_fixed(ratio=12)+
  geom_signif(mapping=aes(x=seurat_clusters,y=STEMNESS_SCORE1),
                comparisons=list(c("0","1"),
                                c("0","2"),
                                c("0","3")),
                map_signif_level=TRUE, #是否显示显著性
                #tip_length=c(0,0,0,0,0,0), #修改显著性线两端的长短
                y_position=c(0.45,0.5,0.55,0.6), #设置显著性线的位置高度
                size=0.6, #线的粗细
                textsize=6, #修改显著性标记的大小
                test="wilcox.test") #检验的类型
ggsave("stem_box.pdf")