
working_dir <- '/data/CARD_singlecell/multiome-test/'; setwd(working_dir)

source('scripts/main/load_packages.r')


object <- readRDS(snakemake@input[['seurat_object']])

hub <- AnnotationHub::AnnotationHub()
reference <- hub[['AH75011']]

project <- Project(object)
batch <- levels(object$batch)

frag.file <- paste0(data_path, batch, '/Multiome/', project, '/outs/atac_fragments.tsv.gz')
counts <- Read10X_h5(paste0(data_path, batch, '/Multiome/', project, '/outs/filtered_feature_bc_matrix.h5'))

object <- object %>% 

    AddChromiumAssay(input.data=counts, enDB=reference, frag.file=frag.file) %>% 
        
    NucleosomeSignal(assay='ATAC') %>% TSSEnrichment(assay='ATAC')


# object.list <- object.list %>% PlotQC(project.name=snakemake@params[['project_name']], outpath=snakemake@params[['plots']])


saveRDS(object, snakemake@output[['seurat_object']])
