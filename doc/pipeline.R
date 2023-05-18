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

## ----gwas---------------------------------------------------------------------
gwas = read.table(system.file("extdata", "gwas.txt", package = "isotwas"),header = T)
head(gwas)

ld_bed = system.file("extdata", "ld.bed", package = "isotwas")
ld_ref = bigsnpr::snp_attach(bigsnpr::snp_readBed2(ld_bed,
                                                   backingfile = tempfile()))
ld_cor = bigsnpr::snp_cor(ld_ref$genotypes)
colnames(ld_cor) = rownames(ld_cor) = ld_ref$map$marker.ID

## ----runZ---------------------------------------------------------------------
model_tsv$Chromosome = model_tsv$chromosome
model_tsv$Position = model_tsv$physical.pos
model_tsv$A1 = model_tsv$allele1
model_tsv$A2 = model_tsv$allele2
res.df = data.frame(Gene = c(),
                    Transcript = c(),
                    Z = c(),
                    P = c(),
                    permute.P = c(),
                    topSNP = c(),
                    topSNP.P = c())
for (isoform in unique(model_tsv$Transcript)){
  
  model_this = subset(model_tsv,
                      Transcript == isoform)
  gwas_cur = subset(gwas,SNP %in% model_this$SNP)
  gwas_cur = gwas_cur[match(model_this$SNP,
                            gwas_cur$SNP),]
  res.df = rbind(res.df,
                 burdenTest(mod = model_this,
                            ld = ld_cor,
                            gene = 'GeneExample',
                            sumStats = gwas_cur,
                            chr = 'Chromosome',
                            pos = 'Position',
                            a1 = 'A1',
                            a2 = 'A2',
                            Z = 'Z',
                            beta = 'Beta',
                            se = 'SE',
                            R2cutoff = 0.01,
                            alpha = 0.25,
                            nperms = 1e2,
                            usePos = F))
  
}
res.df

## ----focus--------------------------------------------------------------------
Omega = matrix(nrow = length(unique(model_tsv$SNP)),
               ncol = nrow(res.df))
colnames(Omega) = res.df$Transcript
rownames(Omega) = unique(model_tsv$SNP)
V = ld_cor[rownames(Omega),rownames(Omega)]
for (i in 1:ncol(Omega)){
  for (j in 1:nrow(Omega)){
    
    mod_this = subset(model_tsv,Transcript == colnames(Omega)[i])
    Omega[j,i] = ifelse(rownames(Omega)[j] %in% mod_this$SNP,
                        mod_this$Weight[mod_this$SNP == rownames(Omega)[j]],
                        0)
    
    
  }
}

res.df = twas_finemap(res.df = res.df,
                      z = 'Z',
                      Omega = Omega,
                      V = V, 
                      max_enum = 3,
                      cutoff = .9)
res.df

