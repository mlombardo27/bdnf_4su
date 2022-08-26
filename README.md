#bdnf-4su

This project analyzes the output of SLAM-sequencing, which labels newly transcribed RNA as well as total RNA. This analysis uses various RNA alignment data outputs. 

## RNA input data:

`grandslam_new_rna` -output from GRAND-SLAM, aligned with HISAT-3N converted reads (T>C conversions)

`feature_counts_full_bdnf` -output from FeatureCounts; aligned with HISAT-3N; provides raw counts of sequencing data

`new_ratio_bayesian_p_de` -uses bayesian model to determine fraction of newly transcribed RNA; aligned with HISAT-3N; used to determine newly transcribed RNA categories (i.e. `bdnf_higher_new_rna` 

`bdnf_lower_new_rna` `bdnf_equals_control` `unclear`)

`salmon_data` -aligned using Salmon; used for IsoformSwitchAnalyzeR

`splicing_data` - aligned using MAJIQ; used for splicing analysis

## Analyses completed and their functions:

### `featureCounts_data_explore_remove_well4` 

Here we analyze the featureCounts data for total rna. Outputs screeplots, biplots, eigencorplots, volcano plot (both with and without plate 4).

This is where we determined the poor RNA quality for plate 4 that we removed from the dataset for the rest of the analysis. 

#### Functions included:

`create_feature_count_table` 

- inputs:
    * file.path to folder of featureCounts .txt files
- outputs a `counts_object` that widens the featureCounts data and puts it in format for pca analysis

* make sure that all the files end in this pattern: `_featureCounts_results.txt`

`make_volcano_plot` 

- inputs:
    * `results` pca table
    * `gene_names` of genes of interest you want labelled 

- outputs a volcano plot showing the log2foldchange and padj of total rna

### `standard_deseq2_total_rna_and_GO`

Here, we run the differential expression analysis and GO analysis on the total RNA without plate 4 included for all timepoints. Outputs volcano plots, p-value charts for differential expression analysis. Outputs gostplots, cnetplots, dotplots for GO analysis. 

#### Functions included:

`run_standard_deseq` 

- inputs:
    * `folder_of_featurecounts` 
    * `base_grep` 
    * `contrast_grep` 
    * `grep_pattern` 
    * `baseName`
    * `contrastName` 

- outputs a deseq object with `results_table` and `deseq_obj`

`label_significant` 

- inputs:
    * deseq `results_table`
    * `log2FoldCut`
    * `log10padj` 
    
- outputs a volcano plot with the significant genes (determined by `log2FoldCut` and `log10padj`) labelled

`filter_significant`

- inputs:
    * deseq `results_table`
    * `padj_cutoff`
    * `log2Fold_cutoff`
    * `direction` 

-outputs a list of ENSG ids that fit your significance inputs

`cp_GO_analysis` 

- inputs:
    * `significant_genes`
    * `bg_ids`
    * `typeGO` 
    
- outputs a GO object

### `deseq_on_grandslam_data`

This is where we run differential expression analysis on GRANDSLAM output (`grandslam_new_rna`), which is of the converted (T>C conversion), or new rna reads. Outputs p-value charts, log2foldchange graph, new rna categorization graph, biplots, screeplots, eigencorplots, volcano plots, venn diagrams of (log2foldchange_total < 0 and log2foldchange_new >0) genes, and density plots.

Here we determined that we couldn't use the converted reads to analyze due to low depth of the reads, which is evident in the poor distribution of the p-value charts. After this, we used a bayesian model to determine the fraction of newly transcribed RNA for the rest of the analysis instead.

#### inputs: 

`grandslam_new_rna` 
`new_ratio_bayesian_p_de`
`full_de_min2` -combined timepoints deseq objects from grandslam_new_rna

#### Functions included: 

`min2_deseq` 

- inputs:
* `filtered_grandslam_data`
* `time` 

- outputs a deseq object with `results_table` and `deseq_obj`

`create_table_1` 

- inputs:
* `filtered_grandslam_data` 
* `time` 

- outputs a formatted counts table

`create_meta_table` 

- inputs: formatted counts table (from create_table_1) and outputs the `metadata` table

`label_significant` 

- inputs:
    * deseq `results_table`
    * `log2FoldCut`
    * `log10padj` 

- outputs a volcano plot with the significant genes (determined by `log2FoldCut` and `log10padj`) labelled

#### Major outputs:

* log2foldchange of new and total rna

* new and total rna categorization

* biplot of all timepoints together for converted reads (new rna)

### `m6a_fraction_total_transcriptome_and_axonal_genes`

This is where we pulled known m6a sites and determined where these sites mostly exist.

We found that 30% of the total transcriptome has an m6a site. Around 77% of axonal genes have an m6a site. Axonal genes only make up 5% of the total transcriptome. 

### `4su_effect_1vs6_1vs2`

We looked at what genes were being differentially expressed due to 4su labelling alone by comparing the 1vs6 hours of the control condition.

We saw that most of the genes being differentially expressed with 4sU were different from those that are differentially expressed with BDNF treatment.

#### Inputs:

`featureCounts`

`new_ratio_bayesian_p`

#### Functions included:

`run_standard_deseq` 
- inputs:
    * `folder_of_featurecounts` 
    * `base_grep`
    * `contrast_grep`
    * `grep_pattern`
    * `baseName`
    * `contrastName` 
    
- outputs a deseq object with `results_table` and `deseq_obj`

    * for `1vs6_dds`: 
    
    `base_grep` = 1 hr
    
    `contrast_grep` = 6hr
    
    `grep_pattern` = CONTROL_.*.*_1|CONTROL_.*.*_6
    
    `baseName` = onehour
    
    `contrastName` = 6hr 

`label_significant`

#### Major outputs:

* volcano plots of 1vs6 hours and 1vs2 hours

* venn diagrams of the overlap of 4su and bdnf for upregulated and downregulated genes

* GO of 4su vs bdnf related genes

* log2foldchange of the total rna for 4su vs bdnf related genes

### `isoform_switch_analysis`

Here we analyzed the salmon alignment data to see which genes were having changes in isoform switching.

We saw that the majority of changes caused by BDNF treatment led to )"unproductive" transcripts, or those with ORF loss, intron retention, coding>noncoding. 

#### Vignette used:
[IsoformSwitchAnalyzer](https://bioconductor.org/packages/devel/bioc/vignettes/IsoformSwitchAnalyzeR/inst/doc/IsoformSwitchAnalyzeR.html)

* the isoformSwitchAnalysisPart2() was nonfunctional on my computer, so I rewrote the code to do step by step

#### Inputs:

`salmon_data` -filepath to salmon aligned data folder

external analyses: 

* this package requires use of webservers to do external analyses that it includes in its output

`CPC2` -identifies coding potential

`PFAM` -identifies protein families

`signalIP` -identifies signal peptides

`IUPred` -identifies disordered protein regions

#### Important functions of the Vignette:

`importIsoformExpression` -inputs the filepath to the salmon data folder and outputs an object with `counts` and `abundance`

`importRdata` -inputs:

    `isoformCountMatrix` :counts from salmon data
    `isoformRepExpression` :abundance from salmon data
    `designMatrix` :dataframe with `sampleID` and `condition` information
    `comparison_to_make` :dataframe identifying which is condition 1 vs condition 2, default would be alphabetical 
    `isoformExonAnnoation` :path to the .gtf file from gencode (I used v40)
    `isoformNtFasta` :path to .fa file 

`isoformSwitchAnalysisPart1` -extracts isoform switches and sequences, does not yet do any analysis

`preFilter` -this package requires a cutoff of 500 genes, so you will have to filter to remove the lowly expressed features

`isoformSwitchTestDEXSeq` -runs differential isoforjm usage tests; controls for FDR and applies effect size cutoffs to correct for confounding effects

`isoformSwitchAnalysisPart2` -nonfunctional, so was replaced by step-by-step running of the source code, ultimately outputs `top_isoform_switches` `consequence_summary` `consequence_enrichment` and `splicing_enrichment` 
    
    * note: 
    for the `analyzeSwitchConsequences()` portion of the part 2 analysis, you select the `consequencesToAnalyze`

    options for consequences to analyze are: 
    "tss","tts","last_exon","isoform_length","intron_structure","intron_retention","isoform_class_code","coding_potential","ORF_genomic","ORF_length","5_utr_length","3_utr_length","isoform_seq_similarity","ORF_seq_similarity","5_utr_seq_similarity","3_utr_seq_similarity","NMD_status","domains_identified","genomic_domain_position","domain_length","signal_peptide_identified","IDR_identified","IDR_length","IDR_type"

#### Major outputs:

most significant switchplots for each timepont: provides visualization of the isoforms available, the overall expression of the gene, the expression of each isoform with bdnf vs control, and the usage of each isoform for bdnf vs control

* note:
    `isoform_expression` is the number of reads each isoform transcript has
    `isoform_usage` is the fraction of each isoform transcript over the total gene expression

categorization of the new and total rna for most significant switch types, such as coding>noncoding, intron retention, ORF loss, more upstream last exon

venn diagram of the overlap between most signficant switch types

### `clustering_bdnf-data`

Here I used clustering techniques on the `new_ratio_bayesian_p_de` data to try and find genes with commonalities in their log2foldchange.

The best functioning cluster methods for the bdnf-4su data was with time series, or partitional, clustering and kmedoid clustering

#### Inputs: 

`new_ratio_bayesian_p` -reorganized this data to filter out the never significant genes, and widen the data so that each gene is its own row with the columns being the time points, and the data being the log2foldchange

#### Clustering methods used:

`kmeans` -most common, uses means for clustering

`part_clust` -partitional, time series clustering

`dbscan` -density-based spatial clustering

`kmed` -uses medians for clustering

`mclust` -machine-based clustering, uses Gaussian mixture modelling
