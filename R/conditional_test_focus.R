#' Conditional gene-level test
#'
#' The function tests the gene association with trait conditional
#' on isoforms
#'
#' @param w_gene vector, weights for gene
#' @param w_iso matrix, weights for significant isoforms
#' @param ld matrix, LD matrix
#' @param z vector, GWAS Z-scores
#' @param gene_name character, gene name
#' @param tx_name vector, character vector of tx names in order of columns of Y
#'
#' @return data frame of conditional Z-scores for all combinations of tx
#'
#' @importFrom MASS ginv
#'
#' @export
conditional_test_focus <- function(w_gene,
                                   w_iso,
                                   ld,
                                   z,
                                   gene_name,
                                   tx_name){

    if (length(w_gene) != nrow(w_iso)){
        stop('Gene weights are not the same length as iso weights')
    }

    if (length(w_gene) != length(z)){
        stop('Gene weights are not the same length as GWAS z-scores')
    }


    if (length(w_gene) != nrow(ld)){
        stop('Gene weights are not the same length as LD')
    }

    make_twas = function(w,z,ld){
        return(list(effect = w %*% z,
                    se = sqrt(w %*% ld %*% w),
                    z = (w %*% z)/(sqrt(w %*% ld %*% w))))
    }

    ### marginal test statistics
    z_twas_gene = make_twas(w_gene,z,ld)$z
    z_twas_gene_effect = make_twas(w_gene,z,ld)$effect
    z_twas_gene_se = make_twas(w_gene,z,ld)$se
    z_twas_iso = apply(w_iso,2,function(w) make_twas(w,z,ld)$z)

    mean.mat = t(cbind(w_gene,w_iso)) %*% ld %*% rep(1,nrow(ld))
    var.mat = t(cbind(w_gene,w_iso)) %*% ld %*% cbind(w_gene,w_iso)
    w_all = cbind(w_gene,w_iso)

    df = data.frame(Test = c(paste0('Overall: ',
                                    c(gene_name,tx_name))),
                    Z = c(z_twas_gene,z_twas_iso),
                    P = 2*pnorm(-abs(c(z_twas_gene,z_twas_iso))))

    if (ncol(w_all) == 2){
        perms = list(as.matrix(2))
    } else {
        perms = list()
        for (l in 1:(ncol(w_all)-1)){
            perms = rlist::list.append(perms,
                                       combn(2:ncol(w_all),
                                             l))
        }}

    for (i in 1:length(perms)){
        ppp = perms[[i]]
        for (c in 1:ncol(ppp)){

           cond.dist = condMVNorm::condMVN(mean = as.numeric(mean.mat),
                                           sigma = var.mat,
                                           dependent.ind = ppp[,c],
                                           given.ind = 1,
                                           X.given = z_twas_iso[ppp[,c]-1])
            df.piece = data.frame(Test = paste0(gene_name,', conditional on ',
                                                paste(tx_name[ppp[,c]-1],
                                                      collapse=', ')),
                                  Z = cond_z_score,
                                  P = 2*pnorm(-abs(cond_z_score)))
            df = rbind(df,df.piece)


        }
    }

    return(df)

}
