
working_dir <- '/data/CARD_singlecell/multiome-test/'; setwd(working_dir)

source('scripts/main/load_packages.r')


object <- readRDS(snakemake@input[['seurat_object']])
        
object <- object %>% subset(
    
    subset=
            nCount_RNA > 300 & nCount_ATAC > 300 &
            nFeature_RNA > 300 & nFeature_ATAC > 300 &
            nucleosome_signal < 2 & TSS.enrichment > 2 &
            percent.mt < 10,
    cells=
            fread(snakemake@input[['metadata']])[
                sample %chin% Project(object) & 
                predicted_gmm_doublets %chin% 'singlet', cells
            ]  
    )

saveRDS(object, snakemake@output[['seurat_object']])

macs_path <- '/data/abbass2/mambaforge/envs/macs/bin/macs3'

peaks <- object %>% 

            CallPeaks(assay='ATAC', macs2.path=macs_path) %>% 
            GenomeInfoDb::keepStandardChromosomes(pruning.mode='coarse') %>% 
            IRanges::subsetByOverlaps(ranges=blacklist_hg38_unified, invert=TRUE)


saveRDS(peaks, snakemake@output[['peaks']])

