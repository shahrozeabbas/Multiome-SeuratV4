working_dir <- '/data/CARD_singlecell/multiome-test/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


reduction.list <- lapply(snakemake@input[['reduction']], readRDS)
object.list <- future.apply::future_lapply(snakemake@input[['seurat_object']], readRDS)

clustering_algorithm <- 3
clustering_resolution <- 0.3

object <- merge(x=object.list[[1]], y=object.list[-1])

for (counter in seq_along(reduction.list)) {

    assay <- DefaultAssay(reduction.list[[counter]])
    reduction.name <- paste0('harmony_', tolower(assay))
    object[[reduction.name]] <- reduction.list[[counter]]

}

message('calculating rna clusters...')

object <- object %>% 

            FindNeighbors(reduction='harmony_rna', dims=1:50) %>% 
            FindClusters(algorithm=clustering_algorithm, resolution=clustering_resolution) %>% 
            RunUMAP(reduction='harmony_rna', dims=1:50, reduction.name='umap.rna', reduction.key='rnaUMAP_')
            
            
object[['rna_clusters']] <- Idents(object)
            
            
message('calculating atac clusters...')

object <- object %>% 

            FindNeighbors(reduction='harmony_atac', dims=2:50) %>% 
            FindClusters(algorithm=clustering_algorithm, resolution=clustering_resolution) %>% 
            RunUMAP(reduction='harmony_atac', dims=2:50, reduction.name='umap.atac', reduction.key='atacUMAP_') 
            
            
object[['atac_clusters']] <- Idents(object)


message('calculating wnn clusters...')

object <- object %>% 

            FindMultiModalNeighbors(
                    dims.list=list(1:50, 2:50),
                    reduction.list=list('harmony_rna', 'harmony_atac') 
                ) %>% 

            FindClusters(graph.name='wsnn', algorithm=clustering_algorithm, resolution=clustering_resolution) %>% 
            RunUMAP(nn.name='weighted.nn', reduction.name='wnn.umap', reduction.key='wnnUMAP_')

            
object[['wnn_clusters']] <- Idents(object)

object[['seurat_clusters']] <- NULL


m <- copy(object@meta.data)
setDT(m, keep.rownames='cells')

fwrite(x=m, file=snakemake@output[['metadata']])


saveRDS(object, snakemake@output[['seurat_object']])



