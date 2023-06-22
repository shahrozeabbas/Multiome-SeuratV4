working_dir <- '/data/CARD_singlecell/multiome-test/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']])


DefaultAssay(object) <- 'ATAC'

object <- object %>% RunTFIDF() %>% FindTopFeatures(min.cutoff='q25')


saveRDS(object, snakemake@output[['seurat_object']])

