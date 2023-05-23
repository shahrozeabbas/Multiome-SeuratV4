
working_dir <- '/data/ShernData/CS032989_Taylor_Kondo/'; setwd(working_dir)

source('scripts/main/load_packages.r')

projects <- fread(snakemake@input[['samples']])[[1]]

object.list <- readRDS(snakemake@input[['seurat_object']])

object.list <- sapply(object.list, function(object) {
        object %>% subset(
            subset=
                    nCount_RNA %between% c(500, 5e4) &
                    nCount_ATAC %between% c(500, 2e5) &
                    nFeature_RNA > 300 & nFeature_ATAC > 300 &
                    nucleosome_signal < 2 & TSS.enrichment > 2 &
                    percent.mt < 50
            )
    }, simplify=FALSE)
    

macs_path <- '/data/abbass2/Apps/conda/envs/macs/bin/macs3'

peaks.list <- sapply(projects, function(project) {
    object.list[[project]] %>% 
        CallPeaks(assay='ATAC', macs2.path=macs_path) %>% 
        keepStandardChromosomes(pruning.mode='coarse') %>% 
        subsetByOverlaps(ranges=blacklist_hg38_unified, invert=TRUE)
})

peaks.use <- suppressWarnings(reduce(x=unlist(as(peaks.list, 'GRangesList'))))


hub <- AnnotationHub()
reference <- hub[['AH75011']]

annotations <- GetGRangesFromEnsDb(ensdb=reference)

genome(annotations) <- 'hg38'
annotations <- renameSeqlevels(annotations, mapSeqlevels(seqlevels(annotations), 'UCSC'))

object.list <- sapply(projects, function(project) {
    object <- object.list[[project]]
    fragments <- Fragments(object[['ATAC']])[[1]]
    cells.use <- GetFragmentData(fragments, slot='cells')               
    
    object <- object %>% subset(cells=cells.use)
    fragpath <- paste0(data_path, project, '/outs/atac_fragments.tsv.gz')
    
    atac_counts <- FeatureMatrix( 
        cells=cells.use,
        features=peaks.use,
        fragments=fragments
    )

    object[['ATAC']] <- CreateChromatinAssay(
        counts=atac_counts,
        fragments=fragpath,
        annotation=annotations
    )

    object
}, simplify=FALSE)


saveRDS(object.list, snakemake@output[['seurat_object']])
