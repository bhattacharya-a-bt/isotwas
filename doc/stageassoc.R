## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE,warning = FALSE)

## ----lookmod------------------------------------------------------------------
set.seed(1789)
model_test = readRDS(system.file("extdata",
                                 "ENSG00000070831_isoTWAS.RDS",
                                 package = "isotwas"))
head(model_test)

## ----lookLD-------------------------------------------------------------------
LD_test = readRDS(system.file("extdata",
                                 "ENSG00000070831_LDMatrix.RDS",
                                 package = "isotwas"))
head(as.matrix(LD_test)[,1:10])

## ----lookGWAS-----------------------------------------------------------------
gwas = read.table(system.file("extdata",
                              "tutorial_SCZ_2022.hm3_filter.sumstats.gz",
                              package = "isotwas"),
                  header=T)
head(gwas)

## ----burdenTest---------------------------------------------------------------
gene_names = c('ENSG00000070831',
               'ENSG00000162552',
               'ENSG00000230068')
out_df = data.frame(Gene = c(),
                    Feature = c(),
                    Z = c(),
                    P = c(),
                    permute.P = c(),
                    topSNP = c(),
                    topSNP.P = c())

for (gene in gene_names){
  
  model = readRDS(system.file("extdata",
                              paste0(gene,"_isoTWAS.RDS"),
                              package = "isotwas"))
  ld = readRDS(system.file("extdata",
                           paste0(gene,"_LDMatrix.RDS"),
                           package = "isotwas"))
  
  for (tx in unique(model$Feature)){
    
    sumstats.cur = subset(gwas,SNP %in% subset(model, Feature == tx)$SNP)
    tx_df = isotwas::burdenTest(mod = subset(model, Feature == tx),
                       ld = ld,
                       gene = gene,
                       sumStats = sumstats.cur,
                       chr = 'Chromosome',
                       pos = 'Position',
                       a1 = 'A1',
                       a2 = 'A2',
                       a1_mod = 'ALT',
                       a2_mod = 'REF',
                       snpName = 'SNP',
                       Z = 'Z',
                       beta = NULL,
                       se = NULL,
                       featureName = 'Feature',
                       R2cutoff = .01,
                       alpha = 1e-3,
                       nperms = 1000,
                       usePos = F)
    out_df = rbind(out_df,tx_df)
    
  }
}

head(out_df)

## ----stage1-------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))
gene = out_df %>%
  group_by(Gene) %>%
  summarise(Screen.P = isotwas::p_screen(P))
gene = as.data.frame(gene)
head(gene)

alpha1=.05
G = nrow(gene)
gene$Screen.P.Adjusted = p.adjust(gene$Screen.P,method = 'fdr')
R = length(unique(gene$Gene[gene$Screen.P.Adjusted < alpha1]))
alpha2 = (R*alpha1)/G
print(gene)

## ----stage2-------------------------------------------------------------------
isoform_new = as.data.frame(matrix(nrow = 0,
                                   ncol = ncol(out_df)+2))
colnames(isoform_new) = c(colnames(out_df),'Screen.P','Confirmation.P')
gene = gene[order(gene$Screen.P),]
ttt = merge(out_df,
            gene[,c('Gene','Screen.P',
                    'Screen.P.Adjusted')],
            by = 'Gene')
  
isoform_new = ttt %>%
  group_by(Gene) %>%
  summarise(Feature = Feature,
            Confirmation.P = isotwas::p_confirm(P,alpha = alpha2))
isoform_new = merge(isoform_new,ttt,by=c('Gene','Feature'))
isoform_new$Confirmation.P = ifelse(isoform_new$Screen.P.Adjusted < 0.05,
                                    isoform_new$Confirmation.P,
                                    1)
isoform_new = isoform_new[,c('Gene','Feature','Z','P','permute.P',
                             'topSNP','topSNP.P',
                             'Screen.P','Screen.P.Adjusted','Confirmation.P')]
print(isoform_new)
print(subset(isoform_new,Screen.P.Adjusted < alpha1 &
               Confirmation.P < alpha2 &
               permute.P < 0.05))

## ----biomart, eval=F----------------------------------------------------------
#  suppressPackageStartupMessages(require(biomaRt))
#  gene_names = unique(isoform_new$Gene)
#  ensembl <- useEnsembl(biomart = "ensembl",
#                     dataset = "hsapiens_gene_ensembl",
#                     mirror = "useast")
#  bm = getBM(attributes = c('ensembl_gene_id',
#                            'chromosome_name',
#                            'start_position',
#                            'end_position'),
#        filters = 'ensembl_gene_id',
#        values = gene_names,
#        mart = ensembl)
#  colnames(bm) = c('Gene','Chromosome','Start','End')

## ----subsetToOverlap----------------------------------------------------------
bm = data.frame(Gene = c('ENSG00000070831','ENSG00000162552','ENSG00000230068'),
                Chromosome = 1,
                Start = c(22052627,22117313,22059197),
                End = c(22101360,22143969,22064199))
isoform_new = merge(bm,isoform_new,by='Gene')

isoform_new = isoform_new[order(isoform_new$Chromosome,
                                isoform_new$Start,
                                decreasing = F),]
isoform_sig = subset(isoform_new,
                     Screen.P.Adjusted < alpha1 &
                       Confirmation.P < alpha2 &
                       permute.P < 0.05)

keep.isoform = c()
if (nrow(isoform_sig) > 1){
      for (i in 1:(nrow(isoform_sig)-1)){
        if (isoform_sig$End[i] > isoform_sig$Start[i+1] - 1e6){
          keep.isoform = unique(c(keep.isoform,
                              c(isoform_sig$Feature[c(i,i+1)])))
        }
      }
  }
isoform_sig = subset(isoform_sig,Feature %in% keep.isoform)

## ----aggregateModels----------------------------------------------------------
all.snps = c()
omega = c()
pos = c()
gene = c()
snp.chr = c()

### COLLECT WEIGHTS FOR SNPS IN THE MODELS
for (i in 1:nrow(isoform_sig)){
  gene_in = isoform_sig$Gene[i]
  model_in = readRDS(system.file("extdata",
                              paste0(gene_in,"_isoTWAS.RDS"),
                              package = "isotwas"))
  model_in = subset(model_in,
                    Feature == isoform_sig$Feature[i])
  Model = data.frame(SNP = model_in$SNP,
                     Chromosome = model_in$Chromosome,
                     Position = model_in$Position,
                     Effect = model_in$Weight,
                     A1 = model_in$ALT,
                     A2 = model_in$REF)
  Model = subset(Model,Effect!=0)
  Model = Model[!duplicated(Model$SNP),]
  all.snps = c(all.snps,
               as.character(Model$SNP))
  omega = c(omega,
            as.numeric(Model$Effect))
  gene = c(gene,
           rep(isoform_sig$Feature[i],nrow(Model)))
  snp.chr = c(snp.chr,
              as.numeric(Model$Chromosome))
  pos = c(pos,as.numeric(Model$Position))
}

tot.df = data.frame(SNP = all.snps,
                    Gene = gene,
                    Effect = omega,
                    Chromosome = snp.chr)

model.df = as.data.frame(matrix(nrow = length(unique(all.snps)),
                                ncol = nrow(isoform_sig)+1))
colnames(model.df) = c('SNP',isoform_sig$Feature)
model.df$SNP = as.character(unique(all.snps))

for (q in 1:nrow(isoform_sig)){
  cur.tot.df = subset(tot.df,Gene == isoform_sig$Feature[q])
  cur.tot.df$SNP = as.character(cur.tot.df$SNP)
  for (i in 1:nrow(model.df)){
    w = which(cur.tot.df$SNP == model.df$SNP[i])
    model.df[i,q+1] = ifelse(length(w) != 0,
                             cur.tot.df$Effect[w],
                             0)
  }
}

model.df$Chromosome = 2
for (i in 1:nrow(model.df)){
  rrr = subset(tot.df,SNP == model.df$SNP[i])
  model.df$Chromosome[i] = rrr$Chromosome[1]
}

head(model.df)

## ----grabLD-------------------------------------------------------------------
V = readRDS(system.file("extdata",
                        "test_LD_finemapping.RDS",
                        package = "isotwas"))

## ----finemapping--------------------------------------------------------------
V = V[model.df$SNP,model.df$SNP]
Omega = Matrix::Matrix(as.matrix(model.df[,-c(1,ncol(model.df))]))
zscores = isoform_sig$Z
m = length(zscores)

### COMPUTE LD BETWEEN TX ON THE GENETIC LEVEL
wcor = isotwas::estimate_cor(as.matrix(Omega),
                             as.matrix(V),
                             intercept=T)[[1]]
diag(wcor) = 1
wcor[is.na(wcor)] = 0

### COMPUTE LD INTERCEPT BETWEEN ISOFRM ON THE GENETIC LEVEL
swld = isotwas::estimate_cor(as.matrix(Omega),
                             as.matrix(V),
                             intercept=T)[[2]]
        
null_res = m * log(1 - 1e-3)
marginal = m * log(1 - 1e-3)
comb_list = list()
for (n in 1:min(2,length(zscores))){
  comb_list = c(comb_list,
                combn(1:length(zscores),n,simplify=F))
  }

pips = rep(0,length(zscores))

### COMPUTE BAYES FACTORS/LIKELIHOOD AT EACH CAUSAL CONFIG
for (j in 1:length(comb_list)){
  subset = comb_list[[j]]
  local = isotwas::bayes_factor(zscores,
                                idx_set = subset,
                                wcor = wcor)
  
  marginal = log(exp(local) + exp(marginal))
  for (idx in subset){
    if (pips[idx] == 0){
      pips[idx] = local
      } else {
        pips[idx] = log(exp(pips[idx]) + exp(local))
      }
  }
  }

pips = exp(pips - marginal)
null_res = exp(null_res - marginal)
isoform_sig$pip = pips
isoform_sig = isoform_sig[order(isoform_sig$pip,decreasing = T),]
npost = isoform_sig$pip/sum(isoform_sig$pip)
csum = cumsum(npost)
isoform_sig$in_cred_set = F

for (i in 1:nrow(isoform_sig)){
  isoform_sig$in_cred_set[i] = T
  if (i > 1){
    if (csum[i] > .9 & csum[i-1] < .9){
      isoform_sig$in_cred_set[i] = T
      }
    if (csum[i] < .9){
      isoform_sig$in_cred_set[i] = T
      }
    if (csum[i] > .9 & csum[i-1] > .9){
      isoform_sig$in_cred_set[i] = F
    }
  }
  }

print(isoform_sig)

