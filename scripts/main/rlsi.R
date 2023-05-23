working_dir <- '/data/ShernData/CS032989_Taylor_Kondo/'; setwd(working_dir)

source('scripts/main/load_packages.r')

plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


object <- readRDS(snakemake@input[['seurat_object']])

DefaultAssay(object) <- 'ATAC'
object.list <- object %>% RunTFIDF() %>% SplitObject(split.by='laneID')

object.list <- sapply(object.list, FindTopFeatures, min.cutoff='q5', simplify=FALSE)
features <- Reduce(union, lapply(object.list, VariableFeatures))

object <- object %>% FindTopFeatures(min.cutoff='q0') %>% RunSVD(features=features)
object.list <- sapply(object.list, RunSVD, features=features, simplify=FALSE)

anchors <- object.list %>% FindIntegrationAnchors(anchor.features=features, reduction='rlsi', k.anchor=7, dims=2:50)
object <- IntegrateEmbeddings(anchorset=anchors, reductions=object[['lsi']], new.reduction.name='integrated_lsi', dims.to.integrate=1:50, k.weight=50)

saveRDS(object, snakemake@output[['seurat_object']])






