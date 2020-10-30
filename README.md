# microarray_pipeline
Rscripts for downloading, processing, and findind DE genes using microarray data

## 0) Install packages

```r
install.packages("stringr")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("biomaRt",
                       "limma",
                       "oligo",
                       "GEOquery",
                       "GEOmetadb",
                       "affycoretools",
                       "biomaRt",
                       "qvalue"))
```

## 1) getGEO_build-latest-geometadb.R

This is the first script used to get microarray or RNA-seq data. This will make a local copy of the database for GEO. This does not require a qsub and only takes a few minutes to run.

- args[1] path to the directory you want the database in

```
Rscript getGEO_build-latest-geometadb.R path/to/output_dir
```
## 2) download_raw_microarray_data.R

This is the second script used to analyze microarray data. This script creates a downloaded-files directory into which it downloads GSMs incluced in given GSE(s). The script also outputs a sample_table.csv file with meta data about the GSMs.
 
- args[1] is the file path to the database
- args[2] is the GSE(s) for the experiments containing the GSMs you want to download and compare in small quotes. i.e. for one GSE: "'GSE69078'". For more than one GSE put in double quotes and comma separate: "'GSE69078','GSE69079'".
- args[3] path to the dirctory to save the sample table and downloaded CEL files to. (needs to be high I/O compatiable, don't have trailing forward slash)

```
Rscript download_raw_microarray_data.R path/to/GEOmetadb.sqlite "'GSE69078','GSE69079'" path/to/save/files/in
```
## 3) analyze_microarray.R

This is the third script used to analyze microarray data. Input is the sample_table.csv with the samples you want to process and compare. The script processes data from different platforms separetly and if there is data from more than one platform, it combines the data at the gene level. If data from mouse and human are to be combined, the mouse genes are converted to human genes, and the data from one-to-one orthologs are combined. Script will need to be modified for additional species. 

After the data is combined, DE genes are idenified using limma. Batch effects caused by data originating from different platforms (GPLs) or experiments (GSEs) are added as covariates in the linear model. If more than two conditions are supplied in args[5] or the Condition column of the sample table, then all possible pairwise comparisons are made, and p-values from pair-wise t-tests as well as a one-way ANOVA are supplied.

- args[1] is the path to the directory containing the downloaded-files dir (not to the downloaded-files dir itself)
- args[2] name for the new output folder and output file prefixes
- args[3] is the path to a sample_table.csv file containing rows for the samples to be included in DE analysis
- args[4] is the path to BrainArrayMapping_withGLP.csv (supplied here in the data directory)
- args[5] is  a comma separated list of the groups each sample belongs to [optional]. ie for 6 samples "treated, treated, treated, untreated, untreated, untreated". If args[5] is not included this information must be added to the smaple table indepenently as a column named "Condition".

The script outputs:
- RMA_df.csv containing the combined processed expression values for each GSM
- batch_corrected_mean_log_expression.csv containing the average corrected log expression values compared to id DE genes
- de_anaylsis.csv with the fold change and pvalue for every comparison for every gene 

```
Rscript analyze_microarray.R path/to/dir/containing_the_downloaded-files_folder name_for_this_analysis path/to/sample_table.csv path/to/BrainArrayMapping_withGLP.csv "treated, treated, treated, untreated, untreated, untreated"(optional)
```
