# directories
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
module_dir <- file.path(root_dir, "code", "survival_analysis")
output_dir <- file.path(patient_dir, "output", "survival_analysis")
dir.create(output_dir, showWarnings = F, recursive = T)

# source functions
source(file.path(module_dir, "utils", "kaplan_meier.R"))

# survival analysis
fname <- file.path(output_dir, "kaplan_meier_adult.pdf")
kaplan_meier_adult <- kaplan_meier(all_cor = tcga_gbm_pnoc008_nn_table, 
                                   surv_data = tcga_gbm_survival)

# save plot
pdf(fname)
print(kaplan_meier_adult, newpage = FALSE)
dev.off()
