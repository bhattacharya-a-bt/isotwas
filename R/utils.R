"%^%" <- function(x, n) {
    with(eigen(x), vectors %*% (values^n * t(vectors)))
    }

.onUnload <- function (libpath) { library.dynam.unload("isoTWAS", libpath)}

calc.mse <- function(obs, pred){
    if(is.vector(obs)) obs <- as.matrix(obs)
    if(is.vector(pred)) pred <- as.matrix(pred)

    n <- nrow(obs)
    rss <- colSums((obs - pred)^2, na.rm = TRUE)
    rss/n
}

calc.r2 <- function(i,obs,pred){

    r2 = summary(stats::lm(obs[,i]~pred[,i]))$adj.r.sq
    ct = stats::cor.test(obs[,i],pred[,i])
    if (is.na(r2)){r2 = 0}
    p = ct$p.value
    return(list(R2 = r2,
                P = p))

}


get_best <- function(list_mods,G=G){
    r2 = matrix(unlist(sapply(list_mods, function(y) sapply(y,function(x) x$R2))),
                nrow = G)
    l_nrows = sapply(list_mods, function(y) sapply(y,function(x) nrow(x$Model)))
    r2[l_nrows == 0] = 0
    r2[is.na(r2)] = 0
    bundle = cbind(apply(r2,1,which.max),1:G)
    return(lapply(1:G,
                  function(i) list_mods[[bundle[i,1]]][[bundle[i,2]]]))
}

cluster_weight = function(x){
  if (length(x) == 0){return(list(sum = 0,
                                  index = NA))}
  if (length(x) == 1){return(list(sum = 0,
                                  index = 1))}
  p = stats::p.adjust(stats::pnorm(
    as.numeric(scale(x)),
    lower.tail=F),
    'fdr')
  if (sum(p < 0.05) > 0){
    return(list(sum = sum(p < 0.05),
                index = which(p < 0.05)))
  }
  if (sum(p < .25) > 0){
    return(list(sum = sum(p < 0.25),
                index = which(p < 0.25)))
  }
  if (sum(p < .25) == 0){
    return(list(sum = length(x),
                index = 1:length(x)))
  }
}

