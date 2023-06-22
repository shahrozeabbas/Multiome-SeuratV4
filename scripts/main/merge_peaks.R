working_dir <- '/data/CARD_singlecell/multiome-test/'; setwd(working_dir)

source('scripts/main/load_packages.r')

future::plan('multicore', workers=snakemake@threads)
options(future.globals.maxSize=50 * 1000 * 1024^2)


peaks.list <- future.apply::future_lapply(snakemake@input[['peaks']], readRDS)

peaks.use <- suppressWarnings(reduce(x=unlist(as(peaks.list, 'GRangesList'))))

saveRDS(peaks.use, snakemake@output[['merged_peaks']])