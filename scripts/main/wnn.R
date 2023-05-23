

working_dir <- '/data/ShernData/CS032989_Taylor_Kondo/'; setwd(working_dir)

source('scripts/main/load_packages.r')

plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']])

DefaultAssay(object) <- 'integratedrna'
VariableFeatures(object) <- rownames(object)

object <- object %>% 
            ScaleData(verbose=FALSE) %>% RunPCA(verbose=FALSE) %>% 
            FindNeighbors(reduction='pca', dims=1:50) %>% FindClusters(algorithm=3, resolution=0.5) %>% StashIdent(save.name='rna_clusters') %>% 
            
            RunUMAP(reduction='pca', dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_') %>%
            RunUMAP(reduction='integrated_lsi', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_')



object <- object %>% 
            FindMultiModalNeighbors(reduction.list=list('pca', 'integrated_lsi'), dims.list=list(1:50, 2:50)) %>% 
            FindClusters(graph.name='wsnn', algorithm=3, resolution=0.5) %>% StashIdent(save.name='wnn_clusters') %>% 
            RunUMAP(nn.name='weighted.nn', reduction.name='wnn.umap', reduction.key='wnnUMAP_')

object[['seurat_clusters']] <- NULL


saveRDS(object, snakemake@output[['seurat_object']])



