DotPlot(major_annotation_ibd_subtype,features = c('CD3E','CD4','CD8A'))
FeaturePlot(major_annotation_ibd_subtype,features = c('CD3E','CD4','CD8A','FOXP3','CXCR5','BATF','IL17A'))
DimPlot(major_annotation_ibd_subtype,group.by = 'cell_subtype',label = T)
DimPlot(major_annotation_ibd_subtype,group.by = 'RNA_snn_res.0.5',label = T)
DotPlot(major_annotation_ibd_subtype,features = c('CD3E','CD4','CD8A','MKI67','TRDC','CXCR5','LYZ','CD19','PTPRC','EPCAM','BATF','MS4A1'),group.by = 'RNA_snn_res.0.5')
table(major_annotation_ibd_subtype$cell_subtype,major_annotation_ibd_subtype$RNA_snn_res.1)
sce.all.filt <- major_annotation_ibd_subtype
celltype=data.frame(ClusterID=0:23,
                    celltype= 0:23) 
celltype[celltype$ClusterID %in% c(7,13,22,11),2]='endothelial_cells'
celltype[celltype$ClusterID %in% c(5,21),2]='Epithelial'
celltype[celltype$ClusterID %in% c(15,0,8),2]='CD4T'
celltype[celltype$ClusterID %in% c(4,16),2]='CD8T'
celltype[celltype$ClusterID %in% c(6,12,23),2]='Stromal'
celltype[celltype$ClusterID %in% c(14),2]='glia'
celltype[celltype$ClusterID %in% c(10),2]='pericyte'
celltype[celltype$ClusterID %in% c(8),2]='Myloid'
celltype[celltype$ClusterID %in% c(19),2]='MAST'
celltype[celltype$ClusterID %in% c(2,1,17,20,18),2]='plasma'
celltype[celltype$ClusterID %in% c(3),2]='B'
celltype[celltype$ClusterID %in% c(9),2]='Cyclin_immune'

sce.all.filt@meta.data$marjor_type = "NA"
for(i in 1:nrow(celltype)){
  sce.all.filt@meta.data[which(sce.all.filt@meta.data$RNA_snn_res.0.5==celltype$ClusterID[i]),'marjor_type'] <- celltype$celltype[i]}
DimPlot(sce.all.filt,group.by = "marjor_type",reduction = 'umap',label = T)
cellRatioPlot(object = sce.all.filt,
              sample.name = "disease",
              celltype.name = "marjor_type")
IBD_immune <- subset(sce.all.filt,marjor_type=='CD4T')

sce.all.filt <- NormalizeData(IBD_immune, 
                              normalization.method = "LogNormalize",
                              scale.factor = 1e4) 
sce.all.filt <- FindVariableFeatures(sce.all.filt)
sce.all.filt <- ScaleData(sce.all.filt)
sce.all.filt <- RunPCA(sce.all.filt, features = VariableFeatures(object = sce.all.filt))
library(harmony)
sce.all.filt <- RunHarmony(sce.all.filt, "orig.ident")

sce.all.filt <- RunUMAP(sce.all.filt,  dims = 1:15, 
                        reduction = "harmony")

sce.all.filt <- FindNeighbors(sce.all.filt, reduction = "harmony",
                              dims = 1:15) 
for (res in c(0.01, 0.05, 0.1, 0.2, 0.3, 0.5,0.8,1)) {
  sce.all.filt=FindClusters(sce.all.filt, #graph.name = "CCA_snn", 
                            resolution = res, algorithm = 1)
}
library(clustree)
p2_tree=clustree(sce.all.filt@meta.data, prefix = "RNA_snn_res.")
p2_tree

p = DotPlot(sce.all.filt, features = unique(c('CD3E','CD4','FOXP3','BATF','IL17A','CXCR5')),
            group.by = "RNA_snn_res.0.8",
            assay='RNA'  ) + coord_flip()
p
DimPlot(sce.all.filt,reduction = "umap",label=T,group.by = 'RNA_snn_res.1' ) 
