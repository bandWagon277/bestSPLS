#' A basic implementation of constrained PLSSVD
#' @name svdSPLS
#' @param X description a n*q matrix
#' @param Y description a n*p matrix
#' @return a list of pls components
#' @examples
#' X=matrix(rnorm(200),10,20)
#' Y=matrix(rnorm(100),5,20)
#' fit = sparsePLS(X,Y,2)
#' @importFrom pls plsr
#' @importFrom irlba irlba
#' @export

sparsePLS_0<-function (X, Y, h = 2, fullrank = TRUE, maxiter = 150, lambda1 = 0.5, 
          lambda2 = 0.5) {
  p <- ncol(X)
  q <- ncol(Y)
  X <- scale(X, center = TRUE, scale = TRUE)
  Y <- scale(Y, center = TRUE, scale = TRUE)
  loadingX <- NULL
  loadingY <- NULL
  Xscore <- NULL
  Yscore <- NULL
  weight <- NULL
  constrain <- function(lambda, V) {
    dis = abs(V) - lambda
    dis[dis < 0] = 0
    return(sign(V) * dis)
  }
  for (i in c(1:h)) {
    M <- t(X) %*% Y
    if (fullrank) {
      M_svd <- svd(M)
    }
    else {
      M_svd <- irlba(M, 1)
    }
    u_old = M_svd$u[, 1]
    v_old = M_svd$v[, 1]
    tol = 1e-04
    iter = 1
    u_new = rep(0, p)
    v_new = rep(0, q)
    while (iter < maxiter) {
      u_temp = constrain(lambda1, M %*% v_old)
      u_new = u_temp/sqrt(sum(u_temp^2))
      v_temp = constrain(lambda2, t(M) %*% u_old)
      v_new = v_temp/sqrt(sum(v_temp^2))
      if (max(abs(u_old - u_new)) < tol) {
        break
      }
      u_old = u_new
      v_old = v_new
      iter = iter + 1
    }
    norm_u = sum(u_new^2)
    norm_v = sum(v_new^2)
    
    weight_x <- u_new/norm_u
    loadingX <- cbind(loadingX, weight_x)
    weight_y <- v_new/norm_v
    loadingY <- cbind(loadingY, weight_y)
    
    e = (X %*% u_new)/norm_u
    w = (Y %*% v_new)/norm_v
    
    Xscore = cbind(Xscore, e)
    Yscore = cbind(Yscore, w)
    
    #deflation
    norm_e = sum(e^2)
    norm_w = sum(w^2)
    inner_ew = sum(e * w)
    c = (t(X) %*% e)/norm_e
    d = (t(Y) %*% e)/inner_ew
    X = X - e %*% t(c)
    Y = Y - e %*% t(d)
  }
  return(list(loadingX = loadingX, loadingY = loadingY, 
              Xscore = Xscore, Yscore = Yscore))
}
