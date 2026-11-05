rblasso <-
function(s, m, om, nlam, n, B0=NULL, soft=NULL, objective=0, tol=1e-5, maxit=500, quiet=TRUE)
{
  p=dim(s)[1]
  q=dim(om)[1]

  # Handle warm start - convert matrix to vector if provided
 if(is.null(B0)) {
    B0 = as.double(rep(0, p*q))
  } else if(is.matrix(B0)) {
    # Warm start provided as matrix - convert to column-major vector
    B0 = as.double(as.vector(B0))
  } else {
    B0 = as.double(B0)
  }

  # Compute soft-thresholding vector from current B0 if not provided.
  # blasso() never updates an entry's soft value when that entry changes, so
  # soft[r,c] must leave out B[r,c]'s own term: M - S B Omega plus
  # S[r,r] Omega[c,c] B[r,c]. Without the add-back a warm start does not
  # converge to the lasso solution.
  if(is.null(soft)) {
    B_mat = matrix(B0, nrow=p, ncol=q)
    soft = as.double(m - s %*% B_mat %*% om +
                     outer(diag(as.matrix(s)), diag(as.matrix(om))) * B_mat)
  } else {
    soft = as.double(soft)
  }

  objective=as.double(objective)
  s = as.double(s)
  m = as.double(m)
  om = as.double(om)
  nlam=as.double(nlam)
  tol=as.double(tol)
  totalit=0
  mode(n) = "integer"
  mode(p) = "integer"
  mode(q) = "integer"
  mode(maxit) = "integer"
  mode(totalit)="integer"
  dotCoutput=.C("blasso", B=B0, S=s, M=m, Om=om, soft=soft, pin=p,
            qin=q, nin=n, lam=nlam, tol=tol, maxit=maxit, totalit=totalit, objective=objective)
  B = matrix(dotCoutput$B, nrow=p, ncol=q)
  return(list(B=B, obj=dotCoutput$objective, iterations=dotCoutput$totalit))
}

