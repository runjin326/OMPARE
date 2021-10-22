
# Detect Outliers using OUTRIDER [@doi:10.1016/j.ajhg.2018.10.025]
BiocManager::install("OUTRIDER")
BiocManager::install("RMariaDB")
suppressPackageStartupMessages({
  library(optparse) 
  library(readr)
  library(OUTRIDER)
  library(RMariaDB)
  library(AnnotationDbi)
  library(tidyverse)
  library(BiocParallel)
})

# Define Directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
data_dir <- file.path(root_dir, "data")
ref_dir <- file.path(root_dir, "references")

analysis_dir <- file.path(root_dir, "code", "outrider_detection")
plots_dir <- file.path(analysis_dir, "plots")
if(!dir.exists(plots_dir)){
  dir.create(plots_dir)
}

results_dir <- file.path(analysis_dir, "results", "outrider_results")
if(!dir.exists(results_dir)){
  dir.create(results_dir, recursive = TRUE)
}

# output directory for FPKM Filter QC plots
fpkm_plot_dir <- file.path(plots_dir, "fpkm_filter_qc")
if(!dir.exists(fpkm_plot_dir)){
  dir.create(fpkm_plot_dir)
}

# output directory for count correlation pre correction
count_cor_pre_dir <- file.path(plots_dir, "count_cor_pre")
if(!dir.exists(count_cor_pre_dir)){
  dir.create(count_cor_pre_dir)
}

# output directory for count correlation post correction
count_cor_post_dir <- file.path(plots_dir, "count_cor_post")
if(!dir.exists(count_cor_post_dir)){
  dir.create(count_cor_post_dir)
}

# output directory for hyper parameter selection figure
hyper_parameter_dir <- file.path(plots_dir, "hyper_param")
if(!dir.exists(hyper_parameter_dir)){
  dir.create(hyper_parameter_dir)
}
  
# Define key parameters
confounding_factors <- c("harmonized_diagnosis", "RNA_library", "tumor_descriptor", "CNS_region")
plotting_group <- c("harmonized_diagnosis", "RNA_library", "tumor_descriptor", "CNS_region")

patients_of_interest<-c("PNOC008-33", "PNOC008-30", "PNOC008-27", "PNOC008-19")

# Read in files necessary for analyses
histology <- readr::read_tsv(file.path(data_dir, "histologies.tsv"), guess_max = 100000)
gene_count <- readRDS(file.path(data_dir, "gene-counts-rsem-expected_count-collapsed.rds"))
clinical_pnoc008 <- readRDS(file.path(data_dir, "pnoc008_clinical.rds"))
gtf <- file.path(ref_dir, "gencode.v27.primary_assembly.annotation.gtf.gz")

########## load txdb for generating steps 
annotation_file <- file.path(ref_dir, "ucsc.knownGenes.db")
mapping_file <- file.path(ref_dir, "mapping_tx_genesymbol.tsv")

# if (!file.exists(mapping_file)) {
#   # establish connection
  con <- dbConnect(MariaDB(), host='genome-mysql.soe.ucsc.edu',
                   dbname="hg19", user='genome')
  map <- dbGetQuery(con, 'select kgId AS TXNAME, geneSymbol from kgXref')
  # save it so that we do not need to reconnect every time
#   readr::write_tsv(map, mapping_file)
  dbDisconnect(con)
# } else {
#   map <- readr::read_tsv(mapping_file)
# }
  
if (!file.exists(annotation_file)) {
  # Define the annotations for the hg38 genome
  txdbUrl <- paste0("https://cmm.in.tum.de/public/",
                    "paper/mitoMultiOmics/ucsc.knownGenes.db")
  # save the file to references so that we don't need to re-run 
  download.file(txdbUrl, annotation_file)
  txdb <- loadDb(annotation_file)
} else {
  txdb <- loadDb(annotation_file)
}

################ Generate histology files with samples of interest 
# Subset files to contain all PBTA samples
histology_pbta <- histology %>% 
  dplyr::filter(cohort == "PBTA",
         sample_type == "Tumor", 
         experimental_strategy == "RNA-Seq")  %>%
  dplyr::group_by(RNA_library, harmonized_diagnosis, primary_site, tumor_descriptor) %>% 
  dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
  # define sampleID for ods assembly step
  dplyr::rename(sampleID=Kids_First_Biospecimen_ID )

# filter to HGG samples based on pathology diagnosis or pathology diagnosis free text
histology_pbta_hgg <- histology %>% 
  dplyr::filter(cohort == "PBTA",
                sample_type == "Tumor", 
                experimental_strategy == "RNA-Seq")  %>%
  dplyr::filter(pathology_diagnosis %in% c("High-grade glioma/astrocytoma (WHO grade III/IV)", 
                                    "Brainstem glioma- Diffuse intrinsic pontine glioma") |
         pathology_free_text_diagnosis =="anaplastic gliomatosis cerebri (who grade 4)") %>%
  dplyr::group_by(RNA_library, harmonized_diagnosis, primary_site, tumor_descriptor) %>% 
  dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
  # define sampleID for ods assembly step
  dplyr::rename(sampleID=Kids_First_Biospecimen_ID )

# filter to GTEx brain region as normal control
histology_gtex <- histology %>% 
  dplyr::filter(cohort == "GTEx",
                 gtex_group == "Brain", 
                 experimental_strategy == "RNA-Seq") 

# use pnoc008 clinical information to get PNOC008 histology df 
clinical_pnoc008_df <- histology %>% 
  dplyr::filter(Kids_First_Biospecimen_ID %in% clinical_pnoc008$Kids_First_Biospecimen_ID)
# merge with GTEx to generate a combined data frame for downstream analysis
histology_gtex_pnoc008 <- bind_rows(histology_gtex, clinical_pnoc008_df)  %>%
  dplyr::group_by(RNA_library, harmonized_diagnosis, primary_site, tumor_descriptor) %>% 
  dplyr::distinct(Kids_First_Participant_ID, .keep_all = TRUE) %>% 
  # define sampleID for ods assembly step
  dplyr::rename(sampleID=Kids_First_Biospecimen_ID)

################ filter the gene count data to contain only protein coding genes
# read gtf and filter to protein coding 
gencode_gtf <- rtracklayer::import(con = gtf)
gencode_gtf <- as.data.frame(gencode_gtf)
gencode_gtf <- gencode_gtf %>%
  dplyr::select(gene_id, gene_name, gene_type) %>%
  filter(gene_type == "protein_coding") %>%
  unique()

# filter expression count file to contain only protein coding gene
gene_count_coding <- gene_count[rownames(gene_count) %in% gencode_gtf$gene_name,]

################ Subset gene counts data to samples of interest  
# select sample of interest 
gene_count_pbta <- gene_count_coding %>% 
  dplyr::select(histology_pbta$sampleID)
#change to integer for next step
gene_count_pbta[] <- lapply(gene_count_pbta, as.integer) 

# select sample of interest 
gene_count_pbta_hgg <- gene_count_coding %>% 
  dplyr::select(histology_pbta_hgg$sampleID)
#change to integer for next step
gene_count_pbta_hgg[] <- lapply(gene_count_pbta_hgg, as.integer) 

# select sample of interest 
gene_count_gtex_pnoc008 <- gene_count_coding %>% 
  dplyr::select(histology_gtex_pnoc008$sampleID)
#change to integer for next  step
gene_count_gtex_pnoc008[] <- lapply(gene_count_gtex_pnoc008, as.integer) 

############ Assemble outrider data sets
# generate pbta
ods_pbta <- OutriderDataSet(countData=gene_count_pbta, colData=histology_pbta)
# generate pbta HGG
ods_pbta_hgg <- OutriderDataSet(countData=gene_count_pbta_hgg, colData=histology_pbta_hgg)
# generate GTEx PNOC008
ods_gtex_pnoc008 <- OutriderDataSet(countData=gene_count_gtex_pnoc008, colData=histology_gtex_pnoc008)

# generate a list so that we can avoid duplicating codes
# ods_list <- list(ods_pbta, ods_pbta_hgg, ods_gtex_pnoc008)
# name_list <- list("ods_pbta", "ods_pbta_hgg", "ods_gtex_pnoc008")

# # for now just use PBTA HGG run as an example
ods_list <- list(ods_pbta_hgg)
name_list <- list("ods_pbta_hgg")

############# calculate FPKM values, label not expressed genes and save plots as QC 
for(i in 1:length(name_list)){
  ods <- ods_list[[i]]
  ods <- filterExpression(ods, txdb, mapping=map,
                               filterGenes=FALSE, savefpkm=TRUE)
  pdf(file.path(fpkm_plot_dir, paste0(name_list[[i]],'_filter.pdf')))
  print(plotFPKM(ods))
  dev.off()
  
  # Filter based on filter labels
  ods <- ods[mcols(ods)$passedFilter,]
  ods_list[[i]]<-ods
}

############# Visualize correlation for confounding factors
for(i in 1:length(name_list)){
  ods<-ods_list[[i]]
  # plot count correlation before the confounder correction
  file_name <- paste0(name_list[[i]], '_count_cor_pre.pdf')
  
  pdf(file.path(count_cor_pre_dir, file_name), width=10, height=6)
  ods <- plotCountCorHeatmap(ods, colGroups=plotting_group,
                             normalized=FALSE, nRowCluster=4)
  dev.off()
  ods_list[[i]]<-ods
}

############# Estimate Size Factors, find encoding dim (q) and do confounding correction
# Define multicore for parallel processing
ncores <- 32
register(MulticoreParam(ncores, ncores*2, progressbar = TRUE))

ods_list <- lapply(ods_list, function(x){
  x<- estimateSizeFactors(x)
  print("finish estimating size")
  # # get encoding dim
  # x <- findEncodingDim(x, BPPARAM=bpparam())
  # # save hyper parameter plots
  # file_name_hyper <- paste0(name_list[[i]],'_hyper_parameter.pdf')
  # pdf(file.path(hyper_parameter_dir, file_name_hyper), width=6, height=6)
  # plotEncDimSearch(x)
  # dev.off()
  # model with confounding factors
  x <- controlForConfounders(x, 
                             BPPARAM=bpparam(), 
                             iterations=1)
})

# plot correlation heatmap after confouding factors corrected
for(i in 1:length(name_list)){
  ods <- ods_list[[i]]
  # plot count correlation AFTER the confounder correction
  file_name <- paste0(name_list[[i]], '_count_cor_post.pdf')
  pdf(file.path(count_cor_post_dir, file_name), width=10, height=6)
  ods <- plotCountCorHeatmap(ods, colGroups=plotting_group,
                             normalized=TRUE, nRowCluster=4)
  dev.off()
}

############# Set exclusion mask, fit binomial model and calculate p value
ods_list <- lapply(ods_list, function(x){
  # set exclusion mask - since they are all independent, all FALSE
  sampleExclusionMask(x) <- FALSE
  # fit negative binomial model
  x <- fit(x)
  # calculate p values
  x <- computePvalues(x, alternative="two.sided", method="BY")
  x <- computeZscores(x)
})


############# Extract results from the table
for(i in 1:length(name_list)){
  ods <- ods_list[[i]]
  # plot count correlation AFTER the confounder correction
  file_name <- paste0(name_list[[i]], '_outrider_res.tsv')
  res <- results(ods)
  readr::write_tsv(res, file.path(results_dir, file_name))
  
  # filter to contain only samples of interest
  sample_interest <- clinical_pnoc008 %>% 
    dplyr::filter(subjectID %in% patients_of_interest) %>% 
    pull(Kids_First_Biospecimen_ID) %>% 
    unique()
  
  # find the subjectID - biospecimen ID match
  match <- clinical_pnoc008 %>% 
    dplyr::select(subjectID, Kids_First_Biospecimen_ID) %>% 
    dplyr::rename(sampleID = Kids_First_Biospecimen_ID)
  
  # annotate subjectID to results
  res_poi <- res %>% 
    dplyr::filter(sampleID %in% sample_interest) %>% 
    dplyr::left_join(match)
  
  # write out results
  file_name <- paste0(name_list[[i]], '_outrider_res_poi.tsv')
  readr::write_tsv(res_poi, file.path(results_dir, file_name))
}



  