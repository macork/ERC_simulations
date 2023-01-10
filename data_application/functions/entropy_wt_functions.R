# What is happening here? Not really sure
invisible(capture.output({
  obj_EB = Rcpp::cppFunction("
                             double obj_EB(Eigen::Map<Eigen::MatrixXd> C,
                             Eigen::Map<Eigen::VectorXd> Z,
                             Eigen::Map<Eigen::VectorXd> Q,
                             Eigen::Map<Eigen::VectorXd> &W,
                             Eigen::Map<Eigen::VectorXd> &Wsum)
                             {
                             W = (-1.0 * C * Z).array().exp() * Q.array();
                             Wsum(0) = W.array().sum();
                             return log(Wsum(0));
                             }", depends = "RcppEigen")
grad_EB = Rcpp::cppFunction("
                            Eigen::VectorXd grad_EB(Eigen::Map<Eigen::MatrixXd> C,
                            Eigen::Map<Eigen::VectorXd> Z,
                            Eigen::Map<Eigen::VectorXd> W,
                            Eigen::Map<Eigen::VectorXd> Wsum)
                            {
                            return (-1.0 * C.transpose() * W).array() / Wsum(0);
                            }", depends = "RcppEigen")
}))

ebal <- function(treatment, # vector or matrix of 1 or several treatments
                 covar_matrix, # matrix of covariates (categories as indicators)
                 t_moments = 1, # balancing moments for treatment skew
                 base_weights = rep(1, length(treatment)),
                 verbose = 0, 
                 reltol = sqrt(.Machine$double.eps), 
                 maxit = 200) 
{
  # Create centered moments of treatment(s)
  treatment <- as.matrix(treatment)
  n <- nrow(treatment)
  p <- ncol(treatment)
  tp <- scale(do.call(cbind, lapply(1:p, function(t) {
    sapply(1:t_moments, function(m) treatment[,t]^m)
  })), center = TRUE, scale = FALSE)
  tp <- scale(tp, center = FALSE, scale = apply(tp, 2, function(x) max(abs(x))))
  
  Q <- base_weights
  W <- rep(0, n)
  Wsum <- c(0)
  C <- cbind(covar_matrix, tp, 
             do.call(cbind, lapply(1:p, function(t) {
               tp[,(t - 1) * t_moments + 1] * covar_matrix })))
  EB <- function(Z) obj_EB(C, Z, Q, W, Wsum)
  gEB <- function(Z) grad_EB(C, Z, W, Wsum)
  o <- optim(par = rep(0, ncol(C)),
             fn = EB, gr = gEB, method = "BFGS",
             control = list(trace = verbose,
                            reltol = reltol,
                            maxit = maxit))
  return(list(weights = W * n / sum(W), opt = o))
}