

library(SingleR)
library(celldex)
library(Seurat)
library(SeuratDisk)
library(dplyr)
library(scater)
library(scRNAseq)
library(parallel)
library(destiny)
library(monocle)

#========================================Multiple Core
detectCores()
detectCores(logical = F)
mc<-makeCluster(getOption("mc.cores",16))
memory.limit(360000)

#========================================Preparation
GSE.number<-"GSE158055"
sc.dir<-"H:/scRNA-seq"
now.dir<-paste(sc.dir,sprintf("%s",GSE.number),sep = "/")
setwd(now.dir)
sc.type<-"10X"

part1<-Read10X("GSE158055_covid19_part1",
               gene.column = 1)
part2<-Read10X("GSE158055_covid19_part2",
               gene.column = 1)

part1.sc<-CreateSeuratObject(counts = part1,project = "part1")
part2.sc<-CreateSeuratObject(counts = part2,project = "part2")

rm(part1)
rm(part2)

sc.annot<-read.csv("GSE158055_cell_annotation.csv")

part1.sc$sampleID<-sc.annot$sampleID
part1.sc$celltype<-sc.annot$celltype
part1.sc$majortype<-sc.annot$majorType

sc.meta<-read.delim("GSE158055_sample_metadata.txt")

sc.meta.pbmc<-sc.meta[grep("PBMC",sc.meta$characteristics..Sample.type),]
sc.meta.pbmc<-sc.meta.pbmc[grep("progression|control",sc.meta.pbmc$characteristics..Sample.time),]

sc.pbmc.grep<-which(part1.sc$sampleID%in%sc.meta.pbmc$Sample.name)

part1.sc<-part1.sc[,sc.pbmc.grep]

part2.sc<-part2.sc[,sc.pbmc.grep]

sc.data<-part1.sc
sc.data@assays$RNA@counts<-sc.data@assays$RNA@counts+part2.sc@assays$RNA@counts
sc.data@assays$RNA@data<-sc.data@assays$RNA@data+part2.sc@assays$RNA@data

sc.data$sampleID<-factor(sc.data$sampleID)
sc.data$celltype<-factor(sc.data$celltype)
sc.data$majortype<-factor(sc.data$majortype)
sc.data$batch<-sc.meta.pbmc$`characteristics...Datasets`[match(sc.data$sampleID,sc.meta.pbmc$Sample.name)]
sc.data$batch<-factor(sc.data$batch)

sc.data.list<-SplitObject(sc.data,split.by = "batch")

# SeuratDisk::Convert(source = "scp_scanpy.gzip.h5ad",
#                     dest="h5Seurat",
#                     overwirte=F)

# x<-SeuratDisk::LoadH5Seurat("scp_scanpy.gzip.h5seurat")
# 
# sc.data.list<-list()
# sc.data.list[[1]]<-CreateSeuratObject(counts = x@assays$RNA@counts,
#                                       min.cells = 3,
#                                       min.features = 200)
# sc.data.list[[1]]$patient<-x$patient
# sc.data.list[[1]]$sort<-x$sort
# sc.data.list[[1]]$cell_type<-x$cell_type
# sc.data.list[[1]]$pheno<-x$pheno
# a<-paste(x$patient,x$pheno,sep='-')
# names(a)<-names(x$patient)
# a<-factor(a)
# sc.data.list[[1]]$orig.ident<-a
# Idents(sc.data.list[[1]])<-sc.data.list[[1]]$orig.ident
# sc.data.list<-list()
# sc.data.list[[1]]<-x

#=======================================Filter
sc.data.list<-lapply(sc.data.list,function(x){a<-x;a[["percent.mt"]]<-PercentageFeatureSet(a, pattern = "^MT-");return(a)})
for(i in 1:length(sc.data.list))
{
  sc.data.list[[i]]$percent.mt[is.nan(sc.data.list[[i]]$percent.mt)]<-0
}
sc.data.list<-lapply(sc.data.list,function(x){subset(x,
                                                     subset = nFeature_RNA>=median(nFeature_RNA)/4&
                                                       nFeature_RNA<=3*median(nFeature_RNA)
                                                       # &percent.mt<=2*median(percent.mt)
                                                     )})


#=======================================Normalization

for(i in 1:length(sc.data.list))
{
  sc.data.list[[i]] <- NormalizeData(sc.data.list[[i]], normalization.method = "LogNormalize", scale.factor = 10000)
}


#=======================================Variable Feature


for(i in 1:length(sc.data.list))
{
  sc.data.list[[i]] <- FindVariableFeatures(sc.data.list[[i]], selection.method = "vst", nfeatures = 2000)
}

#========================================Anchors
sc.anchors<-FindIntegrationAnchors(object.list = sc.data.list)
sc.data<-IntegrateData(anchorset = sc.anchors, dims = 1:30)

DefaultAssay(sc.data)<-"integrated"

#========================================Scale
sc.data<-sc.data[,sc.data$sampleID%in%names(table(sc.data$sampleID)[table(sc.data$sampleID)>500])]

sc.data<-ScaleData(sc.data)

#========================================PCA

sc.data<-RunPCA(sc.data,features = VariableFeatures(sc.data))

VizDimLoadings(sc.data, dims = 1:4, reduction = "pca")

DimPlot(sc.data, reduction = "pca",split.by = NULL)

DimHeatmap(sc.data, dims = 1:15, cells = 500, balanced = TRUE)

#========================================
sc.data <- JackStraw(sc.data, num.replicate = 100)
sc.data <- ScoreJackStraw(sc.data, dims = 1:20)

JackStrawPlot(sc.data, dims = 1:20)

ElbowPlot(sc.data,ndims = 50)

#=======================================Cluster===================

sc.data<-FindNeighbors(sc.data, dims = 1:30)
sc.data<-FindClusters(sc.data, resolution = 1)

#========================================Phylogenetic analysis===================
sc.data<-BuildClusterTree(sc.data,slot = "scale.data")
Tool(object = sc.data, slot = 'BuildClusterTree')
PlotClusterTree(sc.data)


#=======================================Inflection sample==============
sc.data<-CalculateBarcodeInflections(sc.data)
SubsetByBarcodeInflections(sc.data)

#=======================================Dim reduction=================

sc.data<-RunUMAP(sc.data,dims = 1:30)
DimPlot(sc.data,reduction = "umap")

# sc.data.dim.tsne<-RunTSNE(sc.data.cluster)
# DimPlot(sc.data.dim.tsne,reduction = "tsne")

#=======================================Cluster biomaker=================

cluster.markers <- FindAllMarkers(sc.data,only.pos = F,
                                   min.pct = 0.1,logfc.threshold = 0.1
                                   )


cluster.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)

VlnPlot(sc.data, features = c("IFNG"))

FeaturePlot(sc.data,features = c("FCGR3A","CD3D","CD3E","CD3G"))

top10<-cluster.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)
DoHeatmap(sc.data, features = top10$gene) + NoLegend()

#=======================================SingleR annotation===============================

ref1<-celldex::HumanPrimaryCellAtlasData()
 # ref2<-celldex::BlueprintEncodeData()
ref3<-celldex::DatabaseImmuneCellExpressionData()

myref.data<-list()
myref.group<-list()

myref.data[[1]]<-read.delim("F:/myWGCNA/SCRNA Data/GSE128243/GSE128243_logmedian.txt",header = T)
myref.group[[1]]<-colnames(myref.data[[1]])
myref.group[[1]]<-gsub(".*NKT_HS_(.*)[0-9]+$","NKT_\\1ulated",myref.group[[1]])

myref.data[[2]]<-read.delim("F:/myWGCNA/SCRNA Data/GSE124731/GSE124731_Data with annotation.txt",header = T)
myref.data[[2]]<-myref.data[[2]][,grep("CD|NK|MAIT|Vd",colnames(myref.data[[2]]))]
myref.group[[2]]<-read.delim("F:/myWGCNA/SCRNA Data/GSE124731/GSE124731_low_input_rnaseq_meta_data.txt.gz",header = T)
myref.group[[2]]<-myref.group[[2]]$cell_type
myref.group[[2]]<-gsub("CD([48])","CD\\1+_T_cell",myref.group[[2]])
myref.group[[2]]<-gsub("MAIT","T_cell:MAI",myref.group[[2]])
myref.group[[2]]<-gsub("^NK$","NK_cell",myref.group[[2]])
myref.group[[2]]<-gsub("^iNKT$","NKT",myref.group[[2]])
myref.group[[2]]<-gsub("Vd[12]","T_cell:gamma-delta",myref.group[[2]])
myref.data[[2]][grep("gamma-delta",myref.group[[2]])]<-NULL
myref.group[[2]]<- myref.group[[2]][-grep("gamma-delta",myref.group[[2]])]

myref.data[[3]]<-read.delim("F:/myWGCNA/SCRNA Data/GSE128626/GSE128626_data_matrix_sorted_NKT_cells.txt.gz",header = T)
rownames(myref.data[[3]])<-myref.data[[3]][,1]
rownames(myref.data[[3]])<-gsub("\'","",rownames(myref.data[[3]]))
myref.data[[3]]<-myref.data[[3]][,-1]
myref.group[[3]]<-colnames(myref.data[[3]])
myref.group[[3]]<-gsub("NKT_naive.*","NKT_Unstimulated",myref.group[[3]])
myref.group[[3]]<-gsub("NKT_exposed.*","NKT_Stimulated",myref.group[[3]])

myref.data[[4]]<-read.delim("F:/myWGCNA/SCRNA Data/GSE128626/GSE128626_data_matrix_sorted_monocytes.txt.gz",header = T)
rownames(myref.data[[4]])<-myref.data[[4]][,1]
rownames(myref.data[[4]])<-gsub("\'","",rownames(myref.data[[4]]))
myref.data[[4]]<-myref.data[[4]][,-1]
myref.group[[4]]<-colnames(myref.data[[4]])
myref.group[[4]]<-gsub("Monocytes_naive.*","Monocyte:Unstimulated",myref.group[[4]])
myref.group[[4]]<-gsub("Monocytes_exposed.*","Monocyte:Stimulated",myref.group[[4]])

myref.data[[5]]<-read.delim("F:/myWGCNA/Particular cells expression Data/GSE28726/GSE28726_Data with annotation.txt",header = T)
myref.data[[5]]<-myref.data[[5]][,grep("^GSM[0-9]+",colnames(myref.data[[5]]))]
myref.data[[5]]<-log2(myref.data[[5]]+1)
myref.group[[5]]<-read.delim("F:/myWGCNA/Particular cells expression Data/GSE28726/GSE28726_group_mod.txt",header = T)
myref.group[[5]]<-t(myref.group[[5]])
myref.group[[5]]<-myref.group[[5]][,1]
myref.group[[5]]<-gsub(".*CD4 T cell.*","CD4+_T_cell",myref.group[[5]])
myref.group[[5]]<-gsub(".*NKT cell.*resting","NKT_Unstimulated",myref.group[[5]])
myref.group[[5]]<-gsub(".*NKT cell.*stimulated","NKT_Stimulated",myref.group[[5]])
myref.group[[5]]<-gsub(".*CD1d-aGC\\+ Va24- T cell.*resting","dNKT_Unstimulated",myref.group[[5]])
myref.group[[5]]<-gsub(".*CD1d-aGC\\+ Va24- T cell.*stimulated","dNKT_Stimulated",myref.group[[5]])

# myref.data<-read.delim("E:/myWGCNA/SCRNA Data/GSE128243/GSE128243_ReadCounts.txt",
#                      header = T)
# myref.data<-LogNormalize(myref.data)
# myref.group<-read.delim("E:/myWGCNA/SCRNA Data/GSE128243/GSE128243_group_mod.txt",header = T)
# myref.group<-gsub("Human NKT cell ([un]*stimulated) sample [0-9]$","NKT_\\1",myref.group)
# myref.data<-sc.data.nkt@assays$RNA@counts
# myref.data<-as.matrix(myref.data)
# myref.data<-SummarizedExperiment(assays=list(counts=myref.data,logcounts=log2(myref.data+1)))
# myref.group<-sc.data.nkt$orig.ident
myref.group[[1]]<-gsub(".*NKT.*","NKT",myref.group[[1]])
myref.group[[2]]<-gsub(".*NKT.*","NKT",myref.group[[2]])
myref.group[[2]]<-gsub("T_cell:MAI","T_cell_MAIT",myref.group[[2]])
myref.group[[2]]<-gsub(".*CD4.*","T_cell_CD4",myref.group[[2]])
myref.group[[2]]<-gsub(".*CD8.*","T_cell_CD8",myref.group[[2]])
myref.group[[3]]<-gsub(".*NKT.*","NKT",myref.group[[3]])
myref.group[[4]]<-gsub(".*Monocyte.*","Monocyte",myref.group[[4]])
myref.group[[5]]<-gsub(".*NKT.*","NKT",myref.group[[5]])
myref.group[[5]]<-gsub(".*CD4.*","T_cell_CD4",myref.group[[5]])

ref1<-ref1[,ref1$label.main%in%c("B_cell","DC","Erythroblast","Macrophage","Monocyte","Neutrophils","NK_cell","Platelets")]

ref3<-ref3[,ref3$label.main%in%c("B cells","Monocytes","NK cells")]
ref3$label.main<-gsub(".*Monocyte.*","Monocyte",ref3$label.main)
ref3$label.main<-gsub(".*B cells","B_cell",ref3$label.main)
ref3$label.main<-gsub(".*NK cells","NK_cell",ref3$label.main)


# ref.list<-list(ref1,ref2,ref3)
ref.list<-list(ref1,ref3)
# ref.list<-list(ref3)
ref.list<-c(ref.list,myref.data)
# labels.list<-list(ref1$label.main,ref2$label.main,ref3$label.main)
labels.list<-list(ref1$label.main,ref3$label.main)
# labels.list<-list(ref3$label.main)
labels.list<-c(labels.list,myref.group)

Rcpp::sourceCpp(code='
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
IntegerMatrix asMatrix(NumericVector rp,
                       NumericVector cp,
                       NumericVector z,
                       int nrows,
                       int ncols){

  int k = z.size() ;

  IntegerMatrix  mat(nrows, ncols);

  for (int i = 0; i < k; i++){
      mat(rp[i],cp[i]) = z[i];
  }

  return mat;
}')

as_matrix <- function(mat){
  
  row_pos <- mat@i
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])
  
  tmp <- asMatrix(rp = row_pos, cp = col_pos, z = mat@x,
                  nrows =  mat@Dim[1], ncols = mat@Dim[2])
  
  row.names(tmp) <- mat@Dimnames[[1]]
  colnames(tmp) <- mat@Dimnames[[2]]
  return(tmp)
}

sc.data.singleR<-GetAssayData(sc.data,slot = "data")
sc.data.singleR<-as_matrix(sc.data.singleR)

sc.data.annot<-SingleR(test = sc.data.singleR,
                       ref = ref.list,
                       labels = labels.list,
                       de.method = "classic")

rm(sc.data.singleR)

sc.data$annot<-sc.data.annot$pruned.labels
sc.data$cluster_annot<-sc.data.annot$pruned.labels
cluster.group<-unique(sc.data$seurat_clusters)
cluster.group<-cluster.group[order(cluster.group)]
for(i in 1:length(cluster.group))
{
  print(i)
  a<-sc.data.annot$pruned.labels[sc.data$seurat_clusters==cluster.group[i]]
  annot.sum<-table(a)
  annot.sum.ratio<-annot.sum/sum(annot.sum)
  annot.sum.ratio<-annot.sum.ratio[order(annot.sum.ratio,decreasing = T)]
  annot.sum.ratio
  if(max(annot.sum.ratio)>0.7)
  {
    a<-names(annot.sum.ratio)[1]
  } else{
    annot.name.sub<-names(annot.sum.ratio)[annot.sum.ratio>max(annot.sum.ratio)*0.4]
    if(length(annot.name.sub)>1)
    {
      a[!a%in%annot.name.sub]<-names(annot.sum.ratio)[1]
    } else{
      a<-names(annot.sum.ratio)[1]
    }
  }
  sc.data$cluster_annot[sc.data$seurat_clusters==cluster.group[i]]<-a
  i<-i+1
  
}

DimPlot(sc.data, reduction = "umap", label = T,group.by = "cluster_annot")

#========================================marker annotation=========================
setClass(Class = "Cells",
         slots = c(marker="character",
                   positive="logical",
                   subunit="numeric",
                   subgroup="list"),
         sealed = F
         )

setValidity(Class = "Cells",
            method = function(object){
              length(object@positive)==0 ||
                length(object@subunit)==0 ||
                length(object@positive)==length(object@marker)&
                length(object@positive)==length(object@subunit)&
                length(object@marker)==length(object@subunit)
            })

PBMC<-new("Cells")
PBMC@subgroup<-list(HSC=new("Cells",marker=c("CD34","CD38","PTPRC","ITGA2","THY1"),
                            positive=c(T,F,F,T,T)),
                    MPP=new("Cells",marker=c("CD34","CD38","PTPRC","THY1"),
                            positive=c(T,F,F,F)),
                    CLP=new("Cells",marker=c("CD34","CD38","MME","PTPRC"),
                            positive=c(T,T,T,T)),
                    CMP=new("Cells",marker=c("CD34","CD38","CD7","MME","PTPRC","THY1","FLT3"),
                            positive=c(T,T,F,F,F,F,T)),
                    MEP=new("Cells",marker=c("CD34","CD38","CD7","MME","PTPRC","FLT3","IL3RA"),
                            positive=c(T,T,F,F,F,F,F)),
                    GMP=new("Cells",marker=c("CD34","CD38","MME","PTPRC","IL3RA","FLT3"),
                            positive=c(T,T,F,T,T,T)),
                    NK_cell=new("Cells",marker=c("CD3D","CD3E","CD3G","NCAM1"),
                                positive=c(F,F,F,T),
                                subunit=c(1,1,1,0)),
                    T_cell=new("Cells",marker=c("CD3D","CD3E","CD3G"),
                               positive=c(T,T,T),
                               subunit=c(1,1,1)),
                    B_cell=new("Cells",marker=c("CD3D","CD3E","CD3G","CD19","MS4A1"),
                               positive=c(F,F,F,T,T),
                               subunit=c(1,1,1,0,0)),
                    Plasma_cell=new("Cells",marker=c("CD19","SDC1","IL6R","CD52","MZB1"),
                                    positive=c(F,T,T,F,T)),
                    Monocyte=new("Cells",marker=c("CD14"),
                                 positive=c(T)),
                    Macrophage=new("Cells",marker=c("ITGAM","CD68","CD163"),
                                   positive=c(T,T,T)),
                    pDC=new("Cells",marker=c("HLA-DRA","HLA-DRB1","CD209","CLEC4C","IL3RA","LILRA4"),
                            positive=c(T,T,T,T,T,T),
                            subunit=c(1,1,0,0,0,0)),
                    mDC=new("Cells",marker=c("ITGAX","HLA-DRA","HLA-DRB1","CD209","CD1C"),
                            positive=c(F,T,T,T,T),
                            subunit=c(0,1,1,0,0)),
                    Neutrophil=new("Cells",marker=c("ITGAM","CD16","ITGB2","FCGR2A","CD44","CD55","FUT4","ITGA4"),
                                   positive=c(T,T,T,T,T,T,T,F)),
                    Eosinophil=new("Cells",marker=c("PTPRC","IL5RA","CCR3","ADGRE1","ITGAM"),
                                   positive=c(T,T,T,T,T)),
                    Basophil=new("Cells",marker=c("CD19","IL3RA","KIT","ENPP3","FCER1A"),
                                 positive=c(F,T,F,T,T)),
                    Mast_cell=new("Cells",marker=c("FCGR2A","CD33","KIT","ENPP3","FCER1A"),
                                  positive=c(T,T,T,T,T)),
                    Erythroblast=new("Cells",marker=c("GYPA"),
                                     positive=c(T)),
                    Platelets=new("Cells",marker=c("ITGA2B","GP9","GP1BA","ITGB3","PPBP"),
                                  positive=c(T,T,T,T,T))
                    )

subgroup.combine<-function(genes,labels=genes,combine=1:length(genes),positive=rep(T,length(genes)))
{
  if(!all(c(length(genes),length(labels),length(combine),length(positive))==length(genes)))
  {
    stop("Lengths are not equal")
  }
  res.pos<-list()
  res.neg<-list()
  res<-list()
  for(i in unique(combine))
  {
    res.pos<-new("Cells",
                 marker=genes[combine==i],
                 positive=positive[combine==i],
                 subgroup=res)
    res.neg<-new("Cells",
                 marker=genes[combine==i],
                 positive=!positive[combine==i],
                 subgroup=res)
    res<-list(res.pos,res.neg)
    names(res)<-c(paste(labels[combine==i],ifelse(positive[combine==i],"+","-"),sep="",collapse = ""),
                  paste(labels[combine==i],ifelse(positive[combine==i],"-","+"),sep="",collapse = ""))
  }
  res
}

PBMC@subgroup$Platelets@subgroup<-list(active=new("Cells",marker=c("SELP"),
                                                  positive=c(T))
                                       )

PBMC@subgroup$Monocyte@subgroup<-subgroup.combine(genes=c("FCGR3A"),
                                                  labels = c("CD16"))

PBMC@subgroup$NK_cell@subgroup<-subgroup.combine(genes = c("KLRB1"),
                                                 labels = c("CD161"))

PBMC@subgroup$T_cell@subgroup<-list(Th=new("Cells",marker=c("CD4","KLRB1","TRGC1","TRGC2","TRDC"),
                                           positive=c(T,F,F,F,F),
                                           subunit=c(0,0,1,1,1)),
                                    Ts=new("Cells",marker=c("CD8A","KLRB1","TRGC1","TRGC2","TRDC"),
                                           positive=c(T,F,F,F,F),
                                           subunit=c(0,0,1,1,1)),
                                    NKT=new("Cells",marker=c("NCAM1","KLRB1","KLRG1","KLRD1"),
                                            positive=c(T,T,T,T),
                                            subunit=c(0,1,1,1)),
                                    gdT=new("Cells",marker=c("TRGV9","TRDV2"),
                                            positive=c(T,T,T),
                                            subunit=c(1,1,0)),
                                    MAIT=new("Cells",marker=c("SLC4A10","TRAV1-2"),
                                             positive=c(T,T))
                                    )

PBMC@subgroup$T_cell@subgroup$Th@subgroup<-list(active=new("Cells",marker=c("HLA-DRA","HLA-DRB1"),
                                                           positive=c(T,T),
                                                           subunit=c(1,1)),
                                                naive_memory=new("Cells",marker=c("PTPRC"),
                                                                 positive=c(T)),
                                                Treg=new("Cells",marker=c("IL2RA","IL7R"),
                                                         positive=c(T,F))
                                                )

PBMC@subgroup$T_cell@subgroup$Ts@subgroup<-list(active=new("Cells",marker=c("HLA-DRA","HLA-DRB1"),
                                                           positive=c(T,T),
                                                           subunit=c(1,1)),
                                                naive_memory=new("Cells",marker=c("PTPRC"),
                                                                 positive=c(T))
                                                )

PBMC@subgroup$T_cell@subgroup$gdT@subgroup<-subgroup.combine(genes=c("NCAM1","KLRB1","CD8A"),
                                                             labels = c("CD56","CD161","CD8"))
PBMC@subgroup$T_cell@subgroup$NKT@subgroup<-subgroup.combine(genes = c("HAVCR2"),
                                                             labels = c("TIM3"))
cell.marker.annot<-function(marker,log2FC)
{
  if(length(log2FC)!=length(marker))
  {
    stop("The length of marker and log2FC are not equal!")
  }
  cell.type<-"PBMC"
  marker.set<-PBMC
  match.list<-list()
  while (T) {
    match.list[[cell.type]]<-list()
    if(length(marker.set@subgroup)>0)
    {
      temp.cell<-c()
      for(cell.group in names(marker.set@subgroup))
      {
        temp.cell[cell.group]<-0
        group.marker<-rep(NA,length(marker.set@subgroup[[cell.group]]@marker))
        names(group.marker)<-marker.set@subgroup[[cell.group]]@marker
        group.marker<-names(group.marker)%in%marker
        names(group.marker)<-marker.set@subgroup[[cell.group]]@marker
        group.marker[group.marker]<-log2FC[match(names(group.marker)[which(group.marker)],marker)]
        group.log2FC<-group.marker
        group.positive<-ifelse(group.log2FC==0,NA,group.marker>0)
        match.list[[cell.type]][[cell.group]]<-data.frame(marker=marker.set@subgroup[[cell.group]]@marker,
                                                          positive=marker.set@subgroup[[cell.group]]@positive,
                                                          group.marker=ifelse(group.log2FC==0,NA,names(group.log2FC)),
                                                          group.positive=group.positive,
                                                          group.log2FC=group.log2FC)
        if(length(marker.set@subgroup[[cell.group]]@subunit)>0)
        {
          temp.subunit<-marker.set@subgroup[[cell.group]]@subunit
          names(temp.subunit)<-marker.set@subgroup[[cell.group]]@marker
          if(length(temp.subunit[temp.subunit==0])>0)
          {
            if(any(is.na(group.positive[temp.subunit[temp.subunit==0]]))) next
          }
          temp.subunit<-temp.subunit[temp.subunit>0]
          if(length(temp.subunit)>0)
          {
            for(j in 1:length(unique(temp.subunit)))
            {
              group.log2FC[names(temp.subunit[temp.subunit==unique(temp.subunit)[j]])]<-mean(group.log2FC[names(temp.subunit[temp.subunit==unique(temp.subunit)[j]])])
            }
            group.positive<-ifelse(group.log2FC==0,NA,group.log2FC>0)
          }
        }
        
        if(all(group.log2FC*(as.numeric(marker.set@subgroup[[cell.group]]@positive)-0.5)>=0))
        {
          if(all(group.log2FC<0))
          {
            temp.cell[cell.group]<-mean(abs(group.log2FC))
          } else{
            temp.cell[cell.group]<-mean(group.log2FC[group.log2FC>=0])
          }
        }
      }
      
      now.cell<-names(temp.cell)[which(temp.cell==max(temp.cell))][1]
      if(max(temp.cell)>0){
        cell.type<-sprintf("%s_%s",cell.type,now.cell)
      } else{
        break
      }
      marker.set<-marker.set@subgroup[[now.cell]]
    } else{
      break
    }
  }
  cell.type
}

cell.annot<-c()
for(i in unique(cluster.markers$cluster))
{
  marker<-cluster.markers$gene[cluster.markers$cluster==i]
  log2FC<-cluster.markers$avg_log2FC[cluster.markers$cluster==i]
  cell.annot[i]<-cell.marker.annot(marker,log2FC)
}

cell.annot<-gsub("PBMC_","",cell.annot)


sc.data$cluster_marker_annot<-cell.annot[match(as.numeric(sc.data$seurat_clusters),1:length(cell.annot))]

DimPlot(sc.data, reduction = "umap", label = T,group.by = "cluster_marker_annot")


Idents(sc.data)<-paste(as.character(sc.data$seurat_clusters),sc.data$cluster_marker_annot,sep="-")
Idents(sc.data)[grep("NA",Idents(sc.data))]<-"NA"
sc.data<-BuildClusterTree(sc.data,slot = "scale.data")
sc.data@tools$BuildClusterTree$tip.label<-paste(sc.data@tools$BuildClusterTree$tip.label,
                                                1:length(sc.data@tools$BuildClusterTree$tip.label),
                                                sep="-")
Tool(object = sc.data, slot = 'BuildClusterTree')
PlotClusterTree(sc.data,
                type="f",
                node.pos=2,
                no.margin=T)

#========================================NKT annotation========================
sc.data$cluster_final_annot<-sc.data$cluseter_marker_annot
for(i in unique(sc.data$seurat_clusters))
{
  a<-table(sc_data$cluster_annot[sc.data$seurat_clusters==i])
  if(names(a[order(a,decreasing = T)][1])=="NKT") sc_data$cluseter_final_annot[sc.data$seurat_clusters==i]<-"NKT"
}

#========================================pseudo bulk RNA========================
pseudo.bulk.rna<-function(x,...)
{
  UseMethod("pseudo.bulk.rna")
}

pseudo.bulk.rna.list<-function(x,meta=NULL,...){
  if(length(x)==0) 
  {
    stop("list length is 0\n")
  }
  if(class(x[[1]][[1]])!="Seurat")
  {
    stop("list seems not Seurat object list\n")
  }
  common.gene<-rownames(x[[1]]@assays$RNA@counts)
  if(length(x)>1)
  {
    for(i in 2:length(x))
    {
      common.gene<-common.gene[common.gene%in%rownames(x[[i]]@assays$RNA@counts)]
    }
  }
  
  pseudo.bulk<-matrix(0,
                      nrow = length(common.gene),
                      ncol = length(x),
                      dimnames = list(row=common.gene,
                                      col=names(x)
                                      )
                      )
  for(i in 1:ncol(pseudo.bulk))
  {
    pseudo.bulk[,i]<-rowSums(x[[i]]@assays$RNA@counts[rownames(x[[i]]@assays$RNA@counts)%in%rownames(pseudo.bulk),])
  }
  if(!is.null(meta))
  {
    colnames(pseudo.bulk)<-meta[match(colnames(pseudo.bulk),meta[,1]),2]
  } else{
    colnames(pseudo.bulk)<-gsub("^(GSM[0-9]+).*","\\1",colnames(pseudo.bulk))
  }
  pseudo.bulk<-pseudo.bulk[,!is.na(colnames(pseudo.bulk))]
  
  return(pseudo.bulk)
}

pseudo.bulk.rna.Seurat<-function(x,split.by=Idents(x),meta=NULL,...){
  x.list<-SplitObject(x,split.by = split.by)
  pseudo.bulk<-matrix(0,
                      nrow = nrow(x@assays$RNA@counts),
                      ncol = length(x.list),
                      dimnames = list(row=rownames(x@assays$RNA@counts),
                                      col=names(x.list)
                                      )
                      )
  for(i in 1:ncol(pseudo.bulk))
  {
    pseudo.bulk[,i]<-rowSums(x.list[[i]]@assays$RNA@counts)
  }
  if(!is.null(meta))
  {
    colnames(pseudo.bulk)<-meta[match(colnames(pseudo.bulk),meta[,1]),2]
  } else{
    colnames(pseudo.bulk)<-gsub("^(GSM[0-9]+).*","\\1",colnames(pseudo.bulk))
  }
  pseudo.bulk<-pseudo.bulk[,!is.na(colnames(pseudo.bulk))]
  
  return(pseudo.bulk)
}

pseudo.bulk.rna.default<-function(x,...)
{
  cat('You should try a Seurat object or Seurat object list.\n')
}
#=======================================NKT cell==============================

modify_vlnplot<- function(obj,
                          feature,
                          pt.size = 0,
                          plot.margin = unit(c(0, 0, 0, 0), "cm"),
                          ...) {
  p<- VlnPlot(obj, features = feature, pt.size = pt.size, ... )  +
    xlab("") + ylab(feature) + ggtitle("") +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.ticks.y = element_line(),
          axis.title.y = element_text(size = rel(1), angle = 0, vjust = 0.5),
          plot.margin = plot.margin )
  return(p)
}

StackedVlnPlot<- function(obj, features,
                          pt.size = 0,
                          plot.margin = unit(c(0, 0, 0, 0), "cm"),
                          ...) {
  
  plot_list<- purrr::map(features, function(x) modify_vlnplot(obj = obj,feature = x, ...))
  plot_list[[length(plot_list)]]<- plot_list[[length(plot_list)]] +
    theme(axis.text.x=element_text(angle = 45,hjust = 1,vjust = 1), axis.ticks.x = element_line())
  
  p<- patchwork::wrap_plots(plotlist = plot_list, ncol = 1)
  return(p)
}


sc.data.nkt<-sc.data[,grep("T_cell_NKT",sc.data$cluster_final_annot)]

sc.data.nkt<-FindVariableFeatures(sc.data.nkt)

sc.data.nkt<-RunPCA(sc.data.nkt,features = VariableFeatures(sc.data.nkt))

sc.data.nkt <- JackStraw(sc.data.nkt, num.replicate = 100)
sc.data.nkt <- ScoreJackStraw(sc.data.nkt, dims = 1:20)

JackStrawPlot(sc.data.nkt, dims = 1:20)

ElbowPlot(sc.data.nkt,ndims = 50)

pca.dim<-c(2:9,12,13,18)

sc.data.nkt<-FindNeighbors(sc.data.nkt, dims = pca.dim)
sc.data.nkt<-FindClusters(sc.data.nkt, resolution = 0.5)

sc.data.nkt<-BuildClusterTree(sc.data.nkt,slot = "scale.data")
Tool(object = sc.data.nkt, slot = 'BuildClusterTree')
PlotClusterTree(sc.data.nkt)

sc.data.nkt<-RunUMAP(sc.data.nkt,dims = pca.dim)
DimPlot(sc.data.nkt,reduction = "umap")

DimPlot(sc.data.nkt,reduction = "umap",label = T,group.by = "seurat_clusters")

cluster.markers.nkt<-list()

levels.nkt<-levels(factor(sc.data.nkt$seurat_clusters))
cluster.markers.nkt<-list()

for(i in 1:length(cluster.markers.nkt))
{
  cluster.markers.nkt[[levels.nkt[i]]]<-FindMarkers(sc.data.nkt,
                                                    ident.1 = levels.nkt[i],
                                                    only.pos = F,
                                                    min.pct = 0.1,
                                                    logfc.threshold = 0.1)
}

nkt.annot<-c()
nkt.annot[c("0")]<-"NKT_CD8"
nkt.annot[c("7","11")]<-"NKT_CD8_TIM3"
nkt.annot["13"]<-"NKT_CD8_CD62L"
nkt.annot["8"]<-"NKT_CD8_CD62L"
nkt.annot["6"]<-"NKT_CD4_CD40LG" #C0
nkt.annot["16"]<-"NKT_CD8" #C1
nkt.annot[c("15","12","4","14")]<-"NKT_DN_ITGAX" #C1

nkt.annot[c("3")]<-"NKT_CD8"
nkt.annot["14"]<-"NKT_CD8"
nkt.annot[c("1","2","5","9","10")]<-"NKT_CD8"

nkt.annot<-nkt.annot[order(as.numeric(names(nkt.annot)))]

sc.data.nkt$cluster_annot_nkt<-nkt.annot[match(sc.data.nkt$seurat_clusters,names(nkt.annot))]

DimPlot(sc.data.nkt,reduction = "umap",label = T,group.by = "cluster_annot_nkt",repel = T)

sc.data.nkt$cluster_annot_nkt[sc.data.nkt@assays$integrated@data["CD8A",]<1&sc.data.nkt$seurat_clusters=="7"]<-"NKT_CD4_TIM3"


sc.data$cluster_nkt_tim3<-"NA"
sc.data$cluster_nkt_tim3[match(colnames(sc.data.nkt)[sc.data.nkt$cluster_annot_nkt%in%c("NKT_CD8_TIM3_CD62L","NKT_CD4_TIM3_CD62L")],
                               colnames(sc.data))]<-"NKT_TIM3pos"
sc.data$cluster_nkt_tim3[match(colnames(sc.data.nkt)[!sc.data.nkt$cluster_annot_nkt%in%c("NKT_CD8_TIM3_CD62L","NKT_CD4_TIM3_CD62L")],
                               colnames(sc.data))]<-"NKT_TIM3neg"

sc.data.nkt$cluster_nkt_tim3<-"NKT_TIM3neg"
sc.data.nkt$cluster_nkt_tim3[sc.data.nkt$cluster_annot_nkt=="NKT_CD8_TIM3_CD62L"]<-"NKT_TIM3pos"
sc.data.nkt$cluster_nkt_tim3[sc.data.nkt$cluster_annot_nkt=="NKT_CD4_TIM3_CD62L"]<-"NKT_TIM3pos"


#=======================================pseudo time

sc.data.nkt2<-sc.data.nkt[,sc.data.nkt$sampleID%in%sc.meta.cell$Sample.name[sc.meta.cell$.Datasets%in%c("Batch02","Batch07","Batch06","Batch07")]]

cds <- newCellDataSet(as(as_matrix(sc.data.nkt2@assays$RNA@counts),"sparseMatrix"),
                      phenoData = AnnotatedDataFrame(data=sc.data.nkt2@meta.data),
                      featureData = AnnotatedDataFrame(data=data.frame(gene_short_name=rownames(sc.data.nkt2@assays$RNA@counts),
                                                                       row.names = rownames(sc.data.nkt2@assays$RNA@counts))),
                      lowerDetectionLimit = 0.1,
                      expressionFamily = negbinomial.size())
rm(sc.data.nkt2)
cds<-estimateSizeFactors(cds)
cds<-estimateDispersions(cds)


cds<-setOrderingFilter(cds,VariableFeatures(sc.data.nkt))
plot_ordering_genes(cds)

cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)













