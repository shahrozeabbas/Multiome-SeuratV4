working_dir <- '/data/CARD_singlecell/multiome-test/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']])


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

    
DefaultAssay(object) <- 'RNA'

all.genes <- rownames(object)

object <- object %>% 

    NormalizeData() %>% 
    CleanVarGenes(nHVG=2000) %>% 
    ScaleData(features=all.genes, verbose=FALSE) %>% 
    CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes)


saveRDS(object, snakemake@output[['seurat_object']])