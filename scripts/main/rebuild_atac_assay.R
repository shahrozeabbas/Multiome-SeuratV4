working_dir <- '/data/CARD_singlecell/multiome-test/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']])
peaks.use <- readRDS(snakemake@input[['merged_peaks']])

annotations <- get_annotations()
    
project <- Project(object)
batch <- levels(object$batch)

fragments <- Fragments(object[['ATAC']])[[1]]
cells.use <- GetFragmentData(fragments, slot='cells')               
fragpath <- paste0(data_path, batch, '/Multiome/', project, '/outs/atac_fragments.tsv.gz')

atac_counts <- FeatureMatrix( 
    cells=cells.use,
    features=peaks.use,
    fragments=fragments,
    process_n=5000
)

object <- object %>% subset(cells=cells.use)

object[['ATAC']] <- CreateChromatinAssay(
    counts=atac_counts,
    fragments=fragpath,
    annotation=annotations
)


saveRDS(object, snakemake@output[['seurat_object']])
