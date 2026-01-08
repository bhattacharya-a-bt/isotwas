## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, echo=F, message=F-------------------------------------------------
library(isotwas)

## ----generate_genetic_data----------------------------------------------------
train_bed = system.file("extdata", "train.bed", package = "isotwas")
snps = bigsnpr::snp_attach(bigsnpr::snp_readBed2(train_bed,
                                                 backingfile = tempfile()))
dim(snps$genotypes[])

isoform_mat = readRDS(system.file("extdata", "isoform_exp.RDS", package = "isotwas"))
dim(isoform_mat)

## ----bootstrap----------------------------------------------------------------
boot_list = lapply(1:10,function(i){jitter(isoform_mat)})
isoform_mat_rep = rlist::list.rbind(boot_list)
rownames(isoform_mat_rep) = rep(paste0('Sample',1:nrow(snps$fam)),10)
colnames(isoform_mat_rep) = paste0('Isoform',1:ncol(isoform_mat_rep))
dim(isoform_mat_rep)

## ----train, message = F, warning=FALSE, eval=F--------------------------------
# snp_mat = as.matrix(snps$genotypes[])
# colnames(snp_mat) = snps$map$marker.ID
# rownames(snp_mat) = snps$fam$sample.ID
# 
# # Train with all methods (recommended for best model selection)
# isotwas_model = compute_isotwas(
#   X = snp_mat,
#   Y = isoform_mat,
#   Y.rep = isoform_mat_rep,
#   R = 10,
#   id = rownames(isoform_mat_rep),
#   omega_est = 'replicates',
#   omega_nlambda = 5,
#   nfolds = 5,
#   verbose = TRUE,
#   seed = 1789,
#   run_all = TRUE,        # Run all 9 methods
#   return_all = TRUE      # Return results from all methods
# )

## ----train_specific, message = F, warning=FALSE, eval=F-----------------------
# # Run only a subset of methods
# isotwas_model = compute_isotwas(
#   X = snp_mat,
#   Y = isoform_mat,
#   Y.rep = isoform_mat_rep,
#   R = 10,
#   id = rownames(isoform_mat_rep),
#   method = c('mrce_lasso', 'multi_enet', 'sgl', 'univariate'),
#   nfolds = 5,
#   verbose = TRUE,
#   seed = 1789,
#   run_all = FALSE
# )

## ----graph_reg, eval=F--------------------------------------------------------
# # Create similarity matrix from GTF (based on shared exons)
# sim_result = similarity_from_gtf(
#   gtf_path = "gencode.v45.annotation.gtf.gz",
#   gene = "BRCA1",
#   method = "jaccard"  # Options: "jaccard", "overlap_coef", "binary"
# )
# 
# # Use in model training
# isotwas_model = compute_isotwas(
#   X = snp_mat,
#   Y = isoform_mat,
#   Y.rep = isoform_mat_rep,
#   R = 10,
#   id = rownames(isoform_mat_rep),
#   similarity_matrix = sim_result$similarity_matrix,
#   method = c('multi_enet', 'graph_reg', 'univariate'),
#   run_all = FALSE,
#   verbose = TRUE,
#   seed = 1789
# )

## ----examine, eval=F----------------------------------------------------------
# # Print summary
# print(isotwas_model)
# 
# # Access the comparison table (R² for each method)
# isotwas_model$comparison
# 
# # Access the best model
# isotwas_model$best_models
# 
# # Access individual transcript results
# isotwas_model$best_models$Isoform1$r2
# isotwas_model$best_models$Isoform1$weights
# 
# # Get the full weight matrix
# weight_mat = get_weight_matrix(isotwas_model$best_models)

## ----sample-------------------------------------------------------------------
model_isotwas = readRDS(system.file("extdata", "model_example.RDS", package = "isotwas"))
class(model_isotwas)
length(model_isotwas)
names(model_isotwas)
model_isotwas$Model[[1]]

## ----convert------------------------------------------------------------------
model_tsv = convert_model(model_isotwas,
                          snp_annot = snps$map,
                          snp_var = 'marker.ID')
model_tsv

