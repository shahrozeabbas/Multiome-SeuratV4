working_dir <- '/data/CARD_singlecell/snakemake_multiome/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object.list <- readRDS(snakemake@input[['seurat_object']])


s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

object.list <- lapply(object.list, function(object) {
    
    DefaultAssay(object) <- 'RNA'; all.genes <- rownames(object)

    object %>% 

        NormalizeData() %>% 
        CleanVarGenes(nHVG=2000) %>% 
        ScaleData(features=all.genes, verbose=FALSE) %>% 
        CellCycleScoring(s.features=s.genes, g2m.features=g2m.genes)

})

top.genes <- object.list %>% SelectIntegrationFeatures(nfeatures=5000)

object.list <- lapply(object.list, function(object) {

    DefaultAssay(object) <- 'ATAC'

    object %>% RunTFIDF() %>% FindTopFeatures(min.cutoff=10)

})

top.peaks <- Reduce(union, lapply(object.list, VariableFeatures))


object <- merge(x=object.list[[1]], y=object.list[-1])

batch <- 'sample'; noise <- c('percent.mt', 'percent.rb', 'nFeature_RNA', 'nCount_RNA', 'doublet_scores', 'G2M.Score', 'S.Score')


DefaultAssay(object) <- 'RNA'

object <- object %>% 

            ScaleData(vars.to.regress=noise, features=top.genes, verbose=FALSE) %>% RunPCA(verbose=FALSE) %>% 
            harmony::RunHarmony(reduction='pca', group.by.vars=batch, assay.use='RNA')

object@reductions$harmony_rna <- object@reductions$harmony


DefaultAssay(object) <- 'ATAC'

object <- object %>% 

            RunSVD(features=top.peaks) %>% 
            harmony::RunHarmony(reduction='lsi', group.by.vars=batch, assay.use='ATAC', project.dim=FALSE, dims.use=2:50)

object@reductions$harmony_atac <- object@reductions$harmony


saveRDS(object, snakemake@output[['seurat_object']])
