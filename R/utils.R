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

        r2 = summary(lm(obs[,i]~pred[,i]))$adj.r.sq
        ct = cor.test(obs[,i],pred[,i])
        if (is.na(r2)){r2 = 0}
        p = ct$p.value
        return(list(R2 = r2,
                    P = p))

    }

    PRESS <- function(linear.model) {
      #' calculate the predictive residuals
      pr <- residuals(linear.model)/(1-lm.influence(linear.model)$hat)
      #' calculate the PRESS
      PRESS <- sum(pr^2)

      return(PRESS)
    }

    pred_r_squared <- function(linear.model) {
      #' Use anova() to get the sum of squares for the linear model
      lm.anova <- anova(linear.model)
      #' Calculate the total sum of squares
      tss <- sum(lm.anova$'Sum Sq')
      # Calculate the predictive R^2
      pred.r.squared <- 1-PRESS(linear.model)/(tss)

      return(pred.r.squared)
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
