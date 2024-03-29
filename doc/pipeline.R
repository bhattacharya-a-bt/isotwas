## ---- include = FALSE---------------------------------------------------------
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
#  snp_mat = as.matrix(snps$genotypes[])
#  colnames(snp_mat) = snps$map$marker.ID
#  rownames(snp_mat) = snps$fam$sample.ID
#  isotwas_model = compute_isotwas(X = snp_mat,
#                                  Y = isoform_mat,
#                                  Y.rep = isoform_mat_rep,
#                                  R = 10,
#                                  id = rownames(isoform_mat_rep),
#                                  omega_est = 'replicates',
#                                  omega_nlambda = 5,
#                                  method = c('mrce_lasso',
#                                             'multi_enet',
#                                             'univariate',
#                                             'joinet',
#                                             'spls'),
#                                  predict_nlambda = 10,
#                                  family = 'gaussian',
#                                  scale = FALSE,
#                                  alpha = 0.5,
#                                  nfolds = 5,
#                                  verbose = TRUE,
#                                  tx_names = paste0('Isoform',1:3),
#                                  seed = 1789,
#                                  run_all = FALSE,
#                                  return_all = TRUE)

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

