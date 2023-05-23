
rule all:
    input:
        'output/seurat_wnn_integrated_cluster_markers.csv',
        'plots/seurat_wnn_seurat_clusters_lanes_condition_timepoint_umap.pdf'


rule preprocess_first:
    input: 
        samples='input/samples.csv'
    output:
        seurat_object='objects/seurat_object_list_rna_preprocessed_01.rds'
    conda:
        'envs/sc-analysis.yml'
    script: 
        'scripts/main/preprocess_first.R'

rule preprocess_second:
    input: 
        samples='input/samples.csv',
        seurat_object='objects/seurat_object_list_rna_preprocessed_01.rds'
    output:
        seurat_object='objects/seurat_object_list_atac_preprocessed_02.rds'
    conda:
        'envs/sc-analysis.yml'
    script: 
        'scripts/main/preprocess_second.R'


rule preprocess_third:
    input: 
        samples='input/samples.csv',
        seurat_object='objects/seurat_object_list_atac_preprocessed_02.rds'
    output:
        seurat_object='objects/seurat_object_list_rna_atac_preprocessed_filtered_03.rds'
    conda:
        'envs/sc-analysis.yml'
    script: 
        'scripts/main/preprocess_third.R'
        

rule rpca:
    input:
        seurat_object='objects/seurat_object_list_rna_atac_preprocessed_filtered_03.rds'
    output:
        seurat_object='objects/seurat_object_preprocessed_filtered_rpca_integrated_04.rds'
    conda:
        'envs/sc-analysis.yml'
    threads:
        8
    script:
        'scripts/main/rpca.R' 

rule rlsi:
    input:
        seurat_object='objects/seurat_object_preprocessed_filtered_rpca_integrated_04.rds'
    output:
        seurat_object='objects/seurat_object_preprocessed_filtered_rpca_rlsi_integrated_05.rds'
    conda:
        'envs/sc-analysis.yml'     
    threads:
        8   
    script:
        'scripts/main/rlsi.R' 

rule wnn:
    input:
        seurat_object='objects/seurat_object_preprocessed_filtered_rpca_rlsi_integrated_05.rds'
    output:
        seurat_object='objects/seurat_object_preprocessed_filtered_wnn_integrated_06.rds'
    conda:
        'envs/sc-analysis.yml'
    threads:
        8
    script:
        'scripts/main/wnn.R' 

rule markers:
    input:
        seurat_object='objects/seurat_object_preprocessed_filtered_wnn_integrated_06.rds'
    output:
        markers='output/seurat_wnn_integrated_cluster_markers.csv',
        umap='plots/seurat_wnn_seurat_clusters_lanes_condition_timepoint_umap.pdf'
    conda:
        'envs/sc-analysis.yml'
    threads:
        16
    script:
        'scripts/main/markers.R'


# rule subcluster:
#     input:
#         seurat_object='objects/seurat_object_preprocessed_filtered_wnn_integrated_06.rds'
#     output:
#         markers='output/seurat_wnn_integrated_subcluster_markers.csv',
#         seurat_object='objects/seurat_object_preprocessed_filtered_wnn_integrated_subclustered_07.rds'
#     conda:
#         'envs/sc-analysis.yml'
#     threads:
#         16
#     script:
#         'scripts/main/subcluster.R'
