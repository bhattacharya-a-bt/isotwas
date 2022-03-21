#' Compute screening P-value for one gene
#'
#' The function runs a modified version of P_ACT to test
#' hypothesis of at least one gene-transcript association
#'
#' @param teststat numeric, vector of TWAS Z-scores or P-values
#' @param mode character, is the teststat Z-scores or P-values?
#'
#' @return screen P-value for one gene
#'
#'
#' @importFrom mvtnorm pmvnorm
#'
#'
#' @export
p_screen <- function(teststat,
                     mode = 'P'){


    if (mode == 'Z'){

        P = 2*pnorm(-abs(teststat))

    } else {
        P = teststat
    }

    return(ACAT::ACAT(P))



}
