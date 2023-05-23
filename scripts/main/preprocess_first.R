
working_dir <- '/data/ShernData/CS032989_Taylor_Kondo/'; setwd(working_dir)

source('scripts/main/load_packages.r')

projects <- fread(snakemake@input[['samples']])[[1]]

raw.path <- sapply(projects, function(project) paste0(data_path, project, '/outs/raw_feature_bc_matrix.h5'), simplify=FALSE) 
filtered.path <- sapply(projects, function(project) paste0(data_path, project, '/outs/filtered_feature_bc_matrix.h5'), simplify=FALSE) 

raw.input.list <- sapply(raw.path, Read10X_h5, simplify=FALSE)
filtered.input.list <- sapply(filtered.path, Read10X_h5, simplify=FALSE)


soup_rate <- 0.20
reticulate::source_python('scripts/utility/scrublet_py.py')

object.list <- sapply(projects, function(dataset) {

    adj.matrix <- suppressWarnings(SoupCorrect(raw.input.list[[dataset]], filtered.input.list[[dataset]], contamination_rate=soup_rate))
    object <- CreateSeuratObject(adj.matrix, min.cells=0, min.features=0, project=dataset)

    doublet_rate <- (ncol(object) / 1000) * 0.008
    object[['percent.mt']] <- PercentageFeatureSet(object, pattern='^MT-')
    object[['percent.rb']] <- PercentageFeatureSet(object, pattern='^RP[SL]')

    m <- copy(object@meta.data)
    setDT(m, keep.rownames='cells')

    m[, laneID := rep(dataset, nrow(m))]

    m[,
        `:=` (
            laneID=tstrsplit(laneID, '_')[[1]],
            condition=tstrsplit(laneID, '_')[[2]],
            timepoint=tstrsplit(laneID, '_')[[3]]
            )
    ]

    lane <- m[, laneID]
    condition <- m[, condition]
    timepoint <- m[, timepoint]

    names(lane) <- names(condition) <- names(timepoint) <- m[, cells]

    object %>% 
        
        AddMetaData(metadata=factor(lane), col.name='laneID') %>% 
        AddMetaData(metadata=factor(condition), col.name='condition') %>% 
        AddMetaData(metadata=factor(timepoint), col.name='timepoint') %>%
        
        scrublet(n_prin_comps=30, expected_doublet_rate=doublet_rate) 

}, simplify=FALSE)


saveRDS(object.list, snakemake@output[['seurat_object']])






