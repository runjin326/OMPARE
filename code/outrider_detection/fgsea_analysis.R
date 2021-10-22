# Run fGSEA results on OUTRIDER results to determine pathways sig. diff.

suppressPackageStartupMessages({
  library(optparse) 
  library(readr)
  library(GSEABase)
  library(tidyverse)
})

# define directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
ref_dir <- file.path(root_dir, "references")
module_dir <- file.path(root_dir, "code", "outrider_detection")
input_dir <- file.path(module_dir, "results", "outrider_results")

results_dir <- file.path(module_dir, "results", "fgsea")
if(!dir.exists(results_dir)){
  dir.create(results_dir)
}

# load files
gene_set <- getGmt(file.path(ref_dir, 'c2.cp.reactome.v7.4.symbols.gmt'), collectionType = BroadCollection(), geneIdType = SymbolIdentifier())
gene_set <- geneIds(gene_set)

ods_pbta_outrider <- readr::read_tsv(input_dir, "ods_pbta_outrider_res_poi.tsv")
subject_id_list <- ods_pbta_outrider %>% 
  pull(subjectID) %>% unique()

set.seed(2021)
all_fgsea_results <- data.frame()
for (i in 1:length(subject_id_list)){
    patient_of_interest <- subject_id_list[i]
    ods_pbta_outrider_each <- ods_pbta_outrider %>% 
      dplyr::filter(subjectID == patient_of_interest) %>% 
      arrange(zScore) %>% 
      dplyr::select(geneID, zScore) %>% 
      distinct() %>%
      deframe()
    
    # run fgsea for results >50 entries
    if(length(ods_pbta_outrider_each)>50){
      fgseaRes <- fgsea(pathways=gene_set, 
                        stats=ods_pbta_outrider_each, 
                        minSize = 5, 
                        maxSize = 500, 
                        eps = 0.0
                        )
    }
    fgseaRes <- fgseaRes %>% 
      mutate(subjectID = patient_of_interest)
    
    all_fgsea_results <- bind_rows(all_fgsea_results, fgseaRes)
}

all_fgsea_results %>% arrange(padj, ascending = TRUE) %>%
  write_tsv(file.path(results_dir, "pbta_all_fgsea_results_poi.tsv"))
