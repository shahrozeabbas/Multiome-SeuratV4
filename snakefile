import pandas

num_workers = 4

assays = ['RNA', 'ATAC']
project = 'card_sc_brain_atlas'

input_table = 'input/random_samples.csv'
samples = pandas.read_csv(input_table, header=None).loc[:, 0].tolist()[1:]


rule all:
    input:
        'plots/qc_plot1.pdf', 
        'plots/qc_plot2.pdf',
        'output/seurat_wnn_integrated_cluster_markers.csv',
        'plots/seurat_wnn_seurat_clusters_sample_batch_umap.pdf'


rule preprocess_first:
    input:
        samples=input_table
    output:
        seurat_object='objects/seurat_object_{sample}_rna_preprocessed_01.rds'
    params:
        soup_rate=0.10,
        sample='{sample}'
    conda:
        'envs/multiome.yml'
    script: 
        'scripts/main/preprocess_first.R'

rule preprocess_second:
    input: 
        seurat_object='objects/seurat_object_{sample}_rna_preprocessed_01.rds'
    output:
        seurat_object='objects/seurat_object_{sample}_atac_preprocessed_02.rds'
    conda:
        'envs/multiome.yml'
    script: 
        'scripts/main/preprocess_second.R'

rule call_doublets:
    input:
        seurat_object=expand('objects/seurat_object_{sample}_atac_preprocessed_02.rds', sample=samples)
    output:
        metadata='output/unfiltered_metadata.csv'
    params:
        project_name=project
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/gmm_doublet_calling.R'

rule plot_qc:
    input:
        metadata='output/unfiltered_metadata.csv'
    output:
        plot_1='plots/qc_plot1.pdf', plot_2='plots/qc_plot2.pdf'
    params:
        project_name=project
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/plot_qc_metrics.R'

rule preprocess_third:
    input:
        metadata='output/unfiltered_metadata.csv',        
        seurat_object='objects/seurat_object_{sample}_atac_preprocessed_02.rds'
    output:
        peaks='output/{sample}_peaks.rds',
        seurat_object='objects/seurat_object_{sample}_atac_preprocessed_filtered_03.rds'
    conda:
        'envs/multiome.yml'
    script: 
        'scripts/main/preprocess_third.R'

rule merge_peaks:
    input:
        peaks=expand('output/{sample}_peaks.rds', sample=samples)
    output:
        merged_peaks='output/merged_atac_peaks.rds'
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/merge_peaks.R'

rule rebuild_atac:
    input:
        merged_peaks='output/merged_atac_peaks.rds',
        seurat_object='objects/seurat_object_{sample}_atac_preprocessed_filtered_03.rds'
    output:
        seurat_object='objects/seurat_object_{sample}_preprocessed_filtered_rebuilt_04.rds'
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/rebuild_atac_assay.R'

rule process:
    input:
        seurat_object='objects/seurat_object_{sample}_preprocessed_filtered_rebuilt_04.rds'
    output:
        seurat_object='objects/seurat_object_{sample}_preprocessed_filtered_rebuilt_normalized_{assay}_05.rds'
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/{wildcards.assay}_process.R'

rule harmony:
    input:
        expand('objects/seurat_object_{sample}_preprocessed_filtered_rebuilt_normalized_{assay}_05.rds', sample=samples, assay=assays)
    output:
        'reductions/{assay}_harmony.rds'
    threads:
        num_workers * 2
    shell:
        '''
        source /data/abbass2/mambaforge/bin/activate harmony
        Rscript scripts/main/{wildcards.assay}_harmony.R {input} {output} {threads}
        '''
        
rule wnn:
    input:
        reduction=expand('reductions/{assay}_harmony.rds', assay=assays),
        seurat_object=expand('objects/seurat_object_{sample}_preprocessed_filtered_rebuilt_04.rds', sample=samples)
    output:
        metadata='output/final_metadata.csv',
        seurat_object='objects/seurat_object_preprocessed_filtered_wnn_integrated_06.rds'
    conda:
        'envs/multiome.yml'
    threads:
        num_workers
    script:
        'scripts/main/wnn.R' 

rule markers:
    input:
        seurat_object='objects/seurat_object_preprocessed_filtered_wnn_integrated_06.rds'
    output:
        markers='output/seurat_wnn_integrated_cluster_markers.csv',
        umap='plots/seurat_wnn_seurat_clusters_sample_batch_umap.pdf'
    conda:
        'envs/multiome.yml'
    threads:
        num_workers * 4
    script:
        'scripts/main/markers.R'
