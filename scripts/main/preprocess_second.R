
working_dir <- '/data/ShernData/CS032989_Taylor_Kondo/'; setwd(working_dir)

source('scripts/main/load_packages.r')

projects <- fread(snakemake@input[['samples']])[[1]]
filtered.path <- sapply(projects, function(project) paste0(data_path, project, '/outs/filtered_feature_bc_matrix.h5'), simplify=FALSE) 
filtered.input.list <- sapply(filtered.path, Read10X_h5, simplify=FALSE)


object.list <- readRDS(snakemake@input[['seurat_object']])

hub <- AnnotationHub()
reference <- hub[['AH75011']]

frag.files <- sapply(projects, function(project) paste0(data_path, project, '/outs/atac_fragments.tsv.gz'), simplify=FALSE) 

object.list <- sapply(projects, function(dataset) {
    
    object.list[[dataset]] %>% 
        AddChromiumAssay(
            input.data=filtered.input.list[[dataset]], 
            enDB=reference, frag.file=frag.files[[dataset]]) %>% 
            
            NucleosomeSignal(assay='ATAC') %>% TSSEnrichment(assay='ATAC')

}, simplify=FALSE)


object.list <- object.list %>% PlotQC(project.name='CS032989_Taylor_Kondo', outpath='plots/')


saveRDS(object.list, snakemake@output[['seurat_object']])
