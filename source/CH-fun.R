#############
## CH-fun.R
#############
##
## Descriptuion:  Supporting functions for all R scripts.
##

##############
## get.data
##############
## 
##  Description: data preprocessing.
##
get.data <- function(variant_dataset){
  
  variant_dataset <- variant_dataset[which(apply(variant_dataset[, c("VAF_baseline", 
                                                                     "VAF_FU",
                                                                     "VAF_thirdvisit")], 1, function(r){sum(!is.na(r))}) >= 2),]
  
  # data preprocessing
  variant_dataset$t0 <- 0
  variant_dataset$t1 <- variant_dataset$date_diff_baseline_FU
  variant_dataset$t2 <- variant_dataset$date_diff_baseline_FU + variant_dataset$date_diff_FU_thirdvisit
  
  t0 <- 0
  t1 <- variant_dataset$date_diff_baseline_FU
  t1 <- t1[!is.na(t1)]
  t2 <- variant_dataset$date_diff_baseline_FU + variant_dataset$date_diff_FU_thirdvisit
  t2 <- t2[!is.na(t2)]
  
  alltps <- unique(sort(c(t0, t1, t2)))
  
  as <- unique(variant_dataset[,c("indiv_ID", "age", "sex")])
  colnames(as) <- c("ID", "a", "s")
  rownames(as) <- as.character(as[,"ID"])
  as <- as[,-1]
  as[as[,"s"] == "Female", "s"] <- 0
  as[as[,"s"] == "Male", "s"] <- 1
  
  
  as2 <- cbind(as.numeric(as[,"a"]),
               as.numeric(as[,"s"]))
  rownames(as2) <- rownames(as)
  
  as <- as2
  rm(as2)
  as <- cbind(1, as)
  
  allGenes <- unique(variant_dataset$gene)
  # allProtCh <- unique(variant_dataset$protein_change)
  allcDNA_change <- unique(variant_dataset$cDNA_change)
  
  IasJKt <- variant_dataset[rep(1:nrow(variant_dataset), each = 3),c("indiv_ID", "age", "sex", "gene", "cDNA_change")]
  IasJKt <- cbind(IasJKt, 
                  time = c(t(variant_dataset[,c("t0", "t1", "t2")])),
                  VAF = c(t(variant_dataset[, c("VAF_baseline", "VAF_FU", "VAF_thirdvisit")])),
                  VAF_01 = c(t(variant_dataset[, c("VAF_baseline", "VAF_FU", "VAF_thirdvisit")]))/100,
                  # reads = c(t(variant_dataset[, c("variant_reads_baseline", "variant_reads_FU", "variant_reads_thirdvisit")])),
                  reads = c(t(variant_dataset[, c("total_reads_baseline", "total_reads_FU", "total_reads_thirdvisit")])),
                  i = sapply(IasJKt$gene, function(gene){which(allGenes == gene)}),
                  j = sapply(IasJKt$cDNA_change, function(cDNAch){which(allcDNA_change == cDNAch)}),
                  k = sapply(IasJKt$indiv_ID, function(id){which(rownames(as) == id)}))
  
  IasJKt <- IasJKt[which(!is.na(IasJKt$VAF_01)),]
  
  IasJKt <- cbind(IasJKt, 
                  t = sapply(IasJKt$time, function(tp){which(alltps == tp)}))
  
  colnames(as) <- c("1", "a", "s")
  rownames(as) <- NULL
  
  IasJKt <- IasJKt[order(IasJKt$i, 
                         IasJKt$j,
                         IasJKt$k,
                         IasJKt$t), ]
  
  tps <- alltps
  
  ijkt <- IasJKt[, c("i", "j", "k", "t")]
  ijk_sampled <- unique(ijkt[, c("i", "j", "k")])
  
  X0 <- cbind(1, as[ijkt$k, c("a", "s")],
              tps[ijkt$t], as[ijkt$k, c("a", "s")]*tps[ijkt$t])
  
  Z0 <- bdiag(lapply(1:nrow(ijk_sampled), function(r){
    cbind(1, tps[ijkt[which(ijkt$i == ijk_sampled[r,"i"] & 
                              ijkt$j == ijk_sampled[r,"j"] &
                              ijkt$k == ijk_sampled[r,"k"]), "t"]])
  }))
  
  ijktyr <- IasJKt[, c("i", "j", "k", "t", "VAF_01", "reads")]
  
  Y0 <- get.logit(ijktyr$VAF_01)
  R0 <- log(ijktyr$reads)
  
  ijkty <- cbind(ijkt, Y0)
  
  B0 <- build_B_sparse(ijk_sampled = ijk_sampled, p = 2)
  B0_lvls <- as.numeric(sapply(unique(colnames(B0)), function(cn){
    length(unlist(strsplit(cn, split = ":", fixed = TRUE)))}))
  
  colnames(X0) <- c("01", "0a", "0s", 
                    "r1", "ra", "rs")  
  
  res <- list(X0 = X0,
              Z0 = Z0,
              Y0 = Y0,
              R0 = R0,
              B0 = B0,
              B0_lvls = B0_lvls,
              tps = tps,
              IasJKt = IasJKt,
              ijktyr = ijktyr,
              ijkty = ijkty,
              ijkt = ijkt,
              ijk_sampled = ijk_sampled,
              as = as
  )
  
  return(res)
}

##############
## get.logit
##############
## 
##  Description: this function calculates the logit (log-odds) of x.
##
get.logit <- function(x){
  return(log(x/(1-x)))
}

##################
## get.logistic
##################
## 
##  Description: this function calculates the logistic of x.
##
get.logistic <- function(x){
  return(exp(x)/(1 + exp(x)))
}

##################
## ldet
##################
## 
##  Description: this function calculates the log-determinant of A.
##
ldet <- function(A){
  return(as.numeric(determinant(A, logarithm = TRUE)$modulus))
}

##################
## tr
##################
## 
##  Description: this function calculates the trace of A.
##
tr <- function(A){
  return(as.numeric(sum(diag(A))))
}

##################
## get.nl
##################
## 
##  Description:  this function calculates the negeative log-likelihood of 
##                the fixed-effects only model.
##
get.nl <- function(psi,
                   X,
                   Y,
                   R){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  
  Sigma <- c(exp(s0[1] + s0[2]*R))
  Sigmai <- Diagonal(n = length(Sigma), Sigma^-1)
  Sigma <- Diagonal(n = length(Sigma), Sigma)
  
  H <- X %*% beta
  nl <- as.numeric(.5*ldet(Sigma) +.5*t(Y - H)%*%Sigmai%*%(Y - H))
  
  if(is.infinite(nl) | is.nan(nl)){
    return(1e8)
  }
  
  return(nl)
}

##################
## get.gnl
##################
## 
##  Description:  this function calculates the gradient of the negeative log-likelihood of 
##                the fixed-effects only model.
##
get.gnl <- function(psi,
                    X,
                    Y,
                    R){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  
  
  Sigma <- c(exp(s0[1] + s0[2]*R))
  Sigmai <- Diagonal(n = length(Sigma), Sigma^-1)
  Sigma <- Diagonal(n = length(Sigma), Sigma)
  
  H <- X %*% beta

  dnl_beta <- as.numeric(-t(X) %*% Sigmai %*% Y + t(X) %*% Sigmai %*% X %*% beta)
  
  dSigma_s0 <- Sigma
  dSigma_s1 <- exp(s0[1] + s0[2]*R)*R
  dSigma_s1 <- Diagonal(n = length(dSigma_s1), dSigma_s1)
  
  dnl_s0 <- as.numeric(.5*nrow(Sigma) - .5*t(Y - H) %*% Sigmai %*% (Y - H))
  dnl_s1 <- as.numeric(.5*tr(Sigmai %*% dSigma_s1) - .5*t(Y - H) %*% Sigmai %*% dSigma_s1 %*% Sigmai %*% (Y - H))
  
  gnl <- c(dnl_beta,
           dnl_s0,
           dnl_s1)
  
  if(sum(is.infinite(gnl)) > 0 | sum(is.nan(gnl)) > 0){
    return(rep(1e8, length(psi)))
  }
  
  return(gnl)
}

##################
## build_B_sparse
##################
## 
##  Description:  this function generates the block matrix, where each block
##                is a selection matrix selecting the proper level of the
##                random effects in the hierarchical model.
##                See the corresponding paper for details.
##
build_B_sparse <- function(ijk_sampled, p = 2) {
  
  if (!all(c("i","j","k") %in% colnames(ijk_sampled))){
    stop("ijk_sampled must contain columns: i, j, k")
  }
  
  if (!is.numeric(p) || p <= 0) stop("p must be a positive integer")
  
  # sort rows
  ijk_sampled <- ijk_sampled[order(ijk_sampled$i,
                                   ijk_sampled$j,
                                   ijk_sampled$k), ]
  
  n_rows <- nrow(ijk_sampled)
  if (n_rows == 0) stop("ijk_sampled has no rows")
  
  # Hierarchical column structure
  
  # i-level
  i_vals <- sort(unique(ijk_sampled$i))
  names_i <- as.character(i_vals)
  
  # ij-level
  ij_vals <- unique(ijk_sampled[, c("i","j")])
  ij_vals <- ij_vals[order(ij_vals$i, ij_vals$j), ]
  names_ij <- paste(ij_vals$i, ij_vals$j, sep=":")
  
  # ijk-level
  ijk_vals <- ijk_sampled
  names_ijk <- paste(ijk_vals$i,
                     ijk_vals$j,
                     ijk_vals$k,
                     sep=":")
  
  base_col_names <- c(names_i, names_ij, names_ijk)
  
  # Column indices
  col_i <- match(ijk_sampled$i, i_vals)
  
  col_ij <- match(
    paste(ijk_sampled$i, ijk_sampled$j),
    paste(ij_vals$i, ij_vals$j)
  ) + length(i_vals)
  
  col_ijk <- seq_len(n_rows) +
    length(i_vals) +
    nrow(ij_vals)
  
  n_cols <- length(base_col_names)
  
  # Sparse design matrix X
  row_indices <- rep(seq_len(n_rows), 3)
  col_indices <- c(col_i, col_ij, col_ijk)
  
  X <- Matrix::sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = 1,
    dims = c(n_rows, n_cols)
  )
  
  # Expand to p×p block matrix
  Ip <- Matrix::Diagonal(p)
  B  <- Matrix::kronecker(X, Ip)
  
  # Column names (block-expanded)
  block_colnames <- as.vector(sapply(base_col_names, function(nm){
    rep(nm, p)
  }))
  
  colnames(B) <- block_colnames
  
  # Row names (also block-expanded)
  row_base <- paste(ijk_sampled$i,
                    ijk_sampled$j,
                    ijk_sampled$k,
                    sep=":")
  
  block_rownames <- as.vector(sapply(row_base, function(nm){
    rep(nm, p)
  }))
  
  rownames(B) <- block_rownames
  
  return(B)
}

####################
## build_L2_sparse
####################
## 
##  Description:  this function generates the block matrix L2, where each block
##                is a selection matrix selecting the proper level of the
##                random effects in the hierarchical model, to recover the 
##                expected mixed effects Euij = L2 %*% Ev from Ev.
##                See the corresponding paper for details.
##
build_L2_sparse <- function(ijk_sampled, p = 2) {
  
  if (!all(c("i","j","k") %in% colnames(ijk_sampled))) {
    stop("ijk_sampled must contain columns: i, j, k")
  }
  
  if (!is.numeric(p) || p <= 0) stop("p must be a positive integer")
  
  # Sort
  ijk_sampled <- ijk_sampled[order(ijk_sampled$i,
                                   ijk_sampled$j,
                                   ijk_sampled$k), ]
  
  # i-level
  i_vals <- sort(unique(ijk_sampled$i))
  
  # ij-level
  ij_vals <- unique(ijk_sampled[, c("i","j")])
  ij_vals <- ij_vals[order(ij_vals$i, ij_vals$j), ]
  
  # ijk-level
  ijk_vals <- ijk_sampled
  
  n_i   <- length(i_vals)
  n_ij  <- nrow(ij_vals)
  n_ijk <- nrow(ijk_vals)
  
  # Column structure identical to B
  n_cols <- n_i + n_ij + n_ijk
  
  # Row structure = number of unique ij
  n_rows <- n_ij
  
  # Column indices
  col_i <- match(ij_vals$i, i_vals)
  col_ij <- seq_len(n_ij) + n_i
  
  row_indices <- rep(seq_len(n_rows), 2)
  col_indices <- c(col_i, col_ij)
  
  X <- Matrix::sparseMatrix(
    i = row_indices,
    j = col_indices,
    x = 1,
    dims = c(n_rows, n_cols)
  )
  
  # Expand blocks
  Ip  <- Matrix::Diagonal(p)
  L2  <- Matrix::kronecker(X, Ip)
  
  return(L2)
}

####################
## build_L1_sparse
####################
## 
##  Description:  this function generates the block matrix L1, where each block
##                is a selection matrix selecting the proper level of the
##                random effects in the hierarchical model, to recover the 
##                expected mixed effects Eui = L1 %*% Ev from Ev.
##                See the corresponding paper for details.
##
build_L1_sparse <- function(ijk_sampled, p = 2) {
  
  if (!all(c("i","j","k") %in% colnames(ijk_sampled))) {
    stop("ijk_sampled must contain columns: i, j, k")
  }
  
  if (!is.numeric(p) || p <= 0) stop("p must be a positive integer")
  
  # Sort
  ijk_sampled <- ijk_sampled[order(ijk_sampled$i,
                                   ijk_sampled$j,
                                   ijk_sampled$k), ]
  
  # i-level
  i_vals <- sort(unique(ijk_sampled$i))
  
  # ij-level
  ij_vals <- unique(ijk_sampled[, c("i","j")])
  ij_vals <- ij_vals[order(ij_vals$i, ij_vals$j), ]
  
  # ijk-level
  ijk_vals <- ijk_sampled
  
  n_i   <- length(i_vals)
  n_ij  <- nrow(ij_vals)
  n_ijk <- nrow(ijk_vals)
  
  n_cols <- n_i + n_ij + n_ijk
  n_rows <- n_i
  
  X <- Matrix::sparseMatrix(
    i = seq_len(n_rows),
    j = seq_len(n_rows),  # η block is first
    x = 1,
    dims = c(n_rows, n_cols)
  )
  
  Ip <- Matrix::Diagonal(p)
  L1 <- Matrix::kronecker(X, Ip)
  
  return(L1)
}

####################
## get.Sigma_v
####################
## 
##  Description:  this function generates the covariance matrix of the
##                random-effects v.
##
get.Sigma_v <- function(lvls,
                        D_l){
  return(bdiag(D_l[lvls]))
}

####################
## get.EVu
####################
## 
##  Description:  this function generates the conditional expectation of the
##                random-effects v, given the current parameters psi.
##
get.EVu <- function(psi,
                    X,
                    Z,
                    Y,
                    R,
                    B_lvls){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1
  
  
  S_v <- c(Dtau[,B_lvls])
  S_u <- S_v
  S_ui <- S_u^-1
  
  Vu <- solve(t(Z) %*% (Z*Sigma_eps_i) + Diagonal(n = length(S_ui), x = S_ui))
  Eu <- Vu %*% t(Z) %*% ((Y - X %*% beta)*Sigma_eps_i)
  
  return(list(Eu = as.numeric(Eu),
              Vu = Vu))
}

####################
## get.nQ
####################
## 
##  Description:  negative E-step function.
##
get.nQ <- function(psi,
                   X,
                   Z,
                   Y,
                   R,
                   B_lvls,
                   Euc,
                   Euc2,
                   dEuc2_rep,
                   dS_u){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1
  
  
  S_v <- c(Dtau[,B_lvls])
  S_u <- S_v
  S_ui <- S_u^-1
  
  nQ <- -as.numeric(-.5*sum(log(Sigma_eps)) - .5*t(Y - X %*% beta) %*% ((Y - X %*% beta)*Sigma_eps_i)
                    -.5*(-2*t(Euc) %*% t(Z) %*% ((Y - X %*% beta)*Sigma_eps_i) + tr(t(Z) %*% (Z*Sigma_eps_i) %*% Euc2))
                    -.5*sum(log(S_u)) -.5*sum(S_ui * diag(Euc2)))
  
  return(nQ)
}

####################
## get.nQ.try
####################
## 
##  Description:  try-error layer for the negative E-step function.
##
get.nQ.try <- function(psi,
                       X,
                       Z,
                       Y,
                       R,
                       B_lvls,
                       Euc,
                       Euc2,
                       dEuc2_rep,
                       dS_u){
  nQ.res <- try(quiet(get.nQ(psi = psi,
                             X = X,
                             Z = Z,
                             Y = Y,
                             R = R,
                             B_lvls = B_lvls,
                             Euc = Euc,
                             Euc2 = Euc2,
                             dEuc2_rep = dEuc2_rep,
                             dS_u = dS_u)), silent = TRUE)
  if(inherits(nQ.res,'try-error') |
     is.nan(nQ.res) |
     is.infinite(nQ.res) |
     is.na(nQ.res)){
    return(1e8)
  }
  return(nQ.res)
}

####################
## get.dSv
####################
## 
##  Description:  partial derivatives of the covariance matrix of the 
##                random-effects v, w.r.t. the parameters tau.
##
get.dSv <- function(lvls){
  return(c(sapply(1:6, function(i){
    tau <- rep(0, 6)
    tau[i] <- 1
    Dtau <- matrix(data = tau,
                   nrow = 2, ncol = 3, byrow = FALSE)
    
    return(c(Dtau[,lvls]))
  })))
}

####################
## get.gnQ
####################
## 
##  Description:  gradient of the negative E-step function w.r.t. the 
##                free/unknown parameters to optimise.
##
get.gnQ <- function(psi,
                    X,
                    Z,
                    Y,
                    R,
                    B_lvls,
                    Euc,
                    Euc2,
                    dEuc2_rep,
                    dS_u){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1
  
  
  S_v <- c(Dtau[,B_lvls])
  S_u <- S_v
  S_ui <- S_u^-1
  
  dSigma_eps_s0 <- Sigma_eps
  dSigma_eps_s1 <- exp(s0[1] + s0[2]*R)*R
  
  nQ_beta <- -as.numeric(t(Y - X %*% beta)%*%(X*Sigma_eps_i) -t(Euc) %*% t(Z) %*% (X*Sigma_eps_i))
  
  nQ_s0 <- .5*as.numeric(
    -t(Y - X %*% beta)%*%((Y - X %*% beta)*Sigma_eps_i^2*dSigma_eps_s0)
    + 2*t(Euc)%*%t(Z)%*%((Y - X %*% beta)*Sigma_eps_i^2*dSigma_eps_s0)
    +tr((-t(Z)%*%(Z*Sigma_eps_i^2*dSigma_eps_s0))%*%Euc2)
    +sum(Sigma_eps_i*dSigma_eps_s0)
  )
  
  nQ_s1 <- .5*as.numeric(
    -t(Y - X %*% beta)%*%((Y - X %*% beta)*Sigma_eps_i^2*dSigma_eps_s1)
    + 2*t(Euc)%*%t(Z)%*%((Y - X %*% beta)*Sigma_eps_i^2*dSigma_eps_s1)
    +tr((-t(Z)%*%(Z*Sigma_eps_i^2*dSigma_eps_s1))%*%Euc2)
    +sum(Sigma_eps_i*dSigma_eps_s1)
  )
  
  S_ui_rep <- rep(S_ui, times = 6)
  
  nQ_tau <- colSums(Matrix(data = .5*(-S_ui_rep^2*dS_u*dEuc2_rep + S_ui_rep*dS_u), 
                           length(S_ui), 6)) * exp(tail(psi, 6))
  
  return(c(nQ_beta,
           nQ_s0,
           nQ_s1,
           nQ_tau))
}

####################
## get.gnQ.try
####################
## 
##  Description:  try-error layer for the gradient of the negative E-step function.
##
get.gnQ.try <- function(psi,
                        X,
                        Z,
                        Y,
                        R,
                        B_lvls,
                        Euc,
                        Euc2,
                        dEuc2_rep,
                        dS_u){
  gnQ.res <- try(quiet(get.gnQ(psi = psi,
                               X = X,
                               Z = Z,
                               Y = Y,
                               R = R,
                               B_lvls = B_lvls,
                               Euc = Euc,
                               Euc2 = Euc2,
                               dEuc2_rep = dEuc2_rep,
                               dS_u = dS_u)), silent = TRUE)
  if(inherits(gnQ.res,'try-error') | 
     sum(is.infinite(gnQ.res)) > 0 | 
     sum(is.nan(gnQ.res)) > 0 | 
     sum(is.na(gnQ.res)) > 0){
    return(rep(1e8, length(psi)))
  }
  return(gnQ.res)
}




# numDeriv::grad(func = get.nQ.try,
#                x = psi,
#                X = X0,
#                Z = Z,
#                Y = Y0,
#                R = R0,
#                B = B,
#                Euc = Euc,
#                Vuc = Vuc,
#                dS_u = dS_u)
# 
# 
# get.gnQ.try(psi = psi,
#             X = X0,
#             Z = Z,
#             Y = Y0,
#             R = R0,
#             B = B,
#             Euc = Euc,
#             Vuc = Vuc,
#             dS_u = dS_u)

####################
## get.EM
####################
## 
##  Description:  E-M algorithm.
##
get.EM <- function(psi,
                   X,
                   Y,
                   Z,
                   R,
                   B,
                   control = list(verbose = TRUE,
                                  maxIT = 10000,
                                  maxTime = 3, # in hours
                                  xtol = 1e-8,
                                  ftol = 1e-8,
                                  local = FALSE),
                   BFGS_ctrl = list(maxit = 10000,
                                    factr = 0,
                                    lmm = 250,
                                    trace = 1,
                                    REPORT = 1)){
  
  B0_lvls <- as.numeric(sapply(unique(colnames(B)), function(cn){
    length(unlist(strsplit(cn, split = ":", fixed = TRUE)))}))
  
  dS_u0 <- get.dSv(B0_lvls)
  
  maxIt <- control$maxIT
  nIT <- 0
  runTime <- 0
  eps_x <- 1e10
  eps_f <- 1e10
  
  psi_c <- psi
  
  EVuc <- get.EVu(psi = psi_c,
                  X = X,
                  Z = Z %*% B0,
                  Y = Y,
                  R = R,
                  B_lvls = B0_lvls)
  
  Vuc0 <- EVuc$Vu
  Euc0 <- EVuc$Eu
  
  Euc20 <- (Vuc0 + Euc0 %*% t(Euc0))
  dEuc20_rep <- rep(as.numeric(diag(Euc20)), times = 6)
  
  f_c <- get.nQ.try(psi = psi_c,
                    X = X,
                    Z = Z %*% B,
                    Y = Y,
                    R = R,
                    B_lvls = B0_lvls,
                    Euc = Euc0,
                    Euc2 = Euc20,
                    dEuc2_rep = dEuc20_rep,
                    dS_u = dS_u0)  
  
  
  while ((eps_x > control$xtol | eps_f > control$ftol) 
         & nIT < maxIt
         & runTime <= (control$maxTime)*60^2) {
    
    stepTime.start <- proc.time()[3]
    
    if(control$verbose){
      cat(paste0("\nIt: ", nIT,
                 "\terr_x: ", eps_x,
                 "\terr_f: ", eps_f,
                 "\n"
      ))
    }
    
    psi_o <- psi_c
    f_o <- f_c
    
    EVuc <- get.EVu(psi = psi_c,
                    X = X,
                    Z = Z %*% B0,
                    Y = Y,
                    R = R,
                    B_lvls = B0_lvls)
    
    Vuc0 <- EVuc$Vu
    Euc0 <- EVuc$Eu
    
    Euc20 <- (Vuc0 + Euc0 %*% t(Euc0))
    dEuc20_rep <- rep(as.numeric(diag(Euc20)), times = 6)
    
    opt.res_mixed <- optim(par = psi_c,
                           fn = get.nQ.try,
                           gr = get.gnQ.try,
                           method = "L-BFGS-B",
                           X = X,
                           Z = Z %*% B0,
                           Y = Y,
                           R = R,
                           B_lvls = B0_lvls,
                           Euc = Euc0,
                           Euc2 = Euc20,
                           dEuc2_rep = dEuc20_rep,
                           dS_u = dS_u0,
                           control = BFGS_ctrl)
    
    
    psi_c <- opt.res_mixed$par
    f_c <- opt.res_mixed$value
    
    eps_x <- norm(psi_c - psi_o, type = "2")/norm(psi_o, type = "2")
    eps_f <- norm(f_c - f_o, type = "2")/norm(f_o, type = "2")
    
    nIT <- nIT + 1
    
    stepTime.end <- proc.time()[3]
    stepTime <- stepTime.end - stepTime.start
    
    runTime <- runTime + stepTime
    
    # ###
    # plot(psi.true,
    #      psi_c,
    #      xlim = range(psi.true,
    #                   psi_c),
    #      ylim = range(psi.true,
    #                   psi_c))
    # lines(range(psi.true,
    #             psi_c),
    #       range(psi.true,
    #             psi_c))
    # cbind(psi.true,
    #       psi_c)  
    # ###
  }
  
  res <- list()
  res$x <- psi_c
  res$f <- f_c
  res$EVuc <- EVuc
  res$eps_x <- eps_x
  res$eps_f <- eps_f
  
  return(res)
}

####################
## get.dfc
####################
## 
##  Description:  Calculate the degrees of freedom for the conditional AIC/BIC.
##                See paper for details.
##
get.dfc <- function(psi,
                    X,
                    Y,
                    Z,
                    R,
                    B){
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1
  Sigma_eps <- Diagonal(n = length(Sigma_eps), 
                        x = Sigma_eps)
  
  B_lvls <- as.numeric(sapply(unique(colnames(B)), function(cn){
    length(unlist(strsplit(cn, split = ":", fixed = TRUE)))}))
  
  S_v <- c(Dtau[,B_lvls])
  S_u <- S_v
  S_ui <- S_u^-1
  S_u <- Diagonal(n = length(S_u),
                  x = S_u)
  
  V <- Z %*% S_u %*% t(Z) + Sigma_eps
  V_i <- solve(V)
  
  df <- (tr(X %*% solve(t(X) %*% V_i %*% X) %*% t(X) %*% V_i) 
         + tr(Z %*% S_u %*% t(Z) %*% V_i) 
         - tr(Z %*% S_u %*% t(Z) %*% V_i %*% X %*% solve(t(X) %*% V_i %*% X) %*% t(X) %*% V_i))
  return(df)
}

####################
## get.cll
####################
## 
##  Description:  Calculate the conditional log-likelihood.
##                See paper for details.
##
get.cll <- function(psi,
                    X,
                    Z,
                    Y,
                    R,
                    Euc
){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1

  cll <- as.numeric(-.5*(sum(log(Sigma_eps)) 
                         + t(Y - X %*% beta - Z %*% Euc) %*% ((Y - X %*% beta - Z %*% Euc) * Sigma_eps_i) 
                         + length(Y)*log(2*pi)))
  return(cll)
}

####################
## get.cAIC
####################
## 
##  Description:  Calculate the conditional AIC.
##                See paper for details.
##
get.cAIC <- function(psi,
                     X,
                     Y,
                     Z,
                     R,
                     B,
                     Euc){
  
  dfc <- get.dfc(psi = psi,
                 X = X,
                 Y = Y,
                 Z = Z,
                 R = R,
                 B = B)
  
  cll <- get.cll(psi = psi,
                 X = X,
                 Z = Z,
                 Y = Y,
                 R = R,
                 Euc = Euc)
  
  cAIC <- -2*cll + 2*dfc
  
  return(cAIC)
}

####################
## get.cBIC
####################
## 
##  Description:  Calculate the conditional BIC.
##                See paper for details.
##
get.cBIC <- function(psi,
                     X,
                     Y,
                     Z,
                     R,
                     B,
                     Euc){
  
  dfc <- get.dfc(psi = psi,
                 X = X,
                 Y = Y,
                 Z = Z,
                 R = R,
                 B = B)
  
  cll <- get.cll(psi = psi,
                 X = X,
                 Z = Z,
                 Y = Y,
                 R = R,
                 Euc = Euc)
  
  cBIC <- -2*cll + log(length(Y))*dfc
  
  return(cBIC)
}

####################
## get.mll
####################
## 
##  Description:  Calculate the marginal log-likelihood.
##                See paper for details.
##
get.mll <- function(psi,
                    X,
                    Z,
                    Y,
                    R,
                    B){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1
  
  B_lvls <- as.numeric(sapply(unique(colnames(B)), function(cn){
    length(unlist(strsplit(cn, split = ":", fixed = TRUE)))}))
  
  S_v <- c(Dtau[,B_lvls])
  S_v <- Diagonal(n = length(S_v), x = S_v)
  
  V <- Z %*% S_v %*% t(Z) + Diagonal(n = length(Sigma_eps), x = Sigma_eps)
  V_i <- solve(V)
  
  mll <- as.numeric(-.5*(ldet(V) 
                         + t(Y - X %*% beta) %*% V_i %*% (Y - X %*% beta)
                         + length(Y)*log(2*pi)))
  
  return(mll)
}

####################
## get.mll.grad
####################
## 
##  Description:  Calculate the gradient of the marginal log-likelihood w.r.t. psi.
##                See paper for details.
##
get.mll.grad <- function(psi,
                         X,
                         Z,
                         Y,
                         R,
                         B){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1
  
  B_lvls <- as.numeric(sapply(unique(colnames(B)), function(cn){
    length(unlist(strsplit(cn, split = ":", fixed = TRUE)))}))
  
  S_v <- c(Dtau[,B_lvls])
  S_v <- Diagonal(n = length(S_v), x = S_v)
  S_u <- S_v
  S_ui <- S_u^-1
  
  dS_u0 <- get.dSv(B_lvls)
  
  V <- Z %*% S_v %*% t(Z) + Diagonal(n = length(Sigma_eps), x = Sigma_eps)
  V_i <- solve(V)
  
  dSigma_eps_s0 <- Sigma_eps
  dSigma_eps_s1 <- exp(s0[1] + s0[2]*R)*R
  
  dV_s0 <- Diagonal(n = length(dSigma_eps_s0), x = dSigma_eps_s0)
  dV_s1 <- Diagonal(n = length(dSigma_eps_s1), x = dSigma_eps_s1)
  
  dmll.dbeta <- as.numeric(t(X) %*% V_i %*% (Y - X %*% beta))
  
  dmll.s0 <- as.numeric(-.5*(-t(Y - X %*% beta) %*% V_i %*% dV_s0 %*% V_i %*% (Y - X %*% beta) + sum(V_i * t(dV_s0))))
  dmll.s1 <- as.numeric(-.5*(-t(Y - X %*% beta) %*% V_i %*% dV_s1 %*% V_i %*% (Y - X %*% beta) + sum(V_i * t(dV_s1))))
  
  dmll.tau <- sapply(1:6, function(i) {
    tau <- rep(0, 6)
    tau[i] <- 1
    Dtau <- matrix(data = tau, nrow = 2, ncol = 3, byrow = FALSE)
    
    dS_v <- Diagonal(n = length(c(Dtau[, B_lvls])), x = c(Dtau[, B_lvls]))
    dV_dtau <- Z %*% dS_v %*% t(Z)
    
    dmll.tauj <- as.numeric(-.5*(-(t(Y - X %*% beta) %*% V_i) %*% (dV_dtau %*% (V_i %*% (Y - X %*% beta))) + sum(V_i * t(dV_dtau))))
    
    return(dmll.tauj)
  })
  
  dmll.tau <- dmll.tau*exp(tail(psi, 6))
  
  return(c(dmll.dbeta,
           dmll.s0,
           dmll.s1,
           dmll.tau))
}

####################
## get.mll.hess
####################
## 
##  Description:  Calculate the hessian of the marginal log-likelihood w.r.t. psi.
##                See paper for details.
##
get.mll.hess <- function(psi,
                         X,
                         Z,
                         Y,
                         R,
                         B){
  
  beta <- psi[1:ncol(X)]
  s0 <- psi[(ncol(X) + 1):(ncol(X) + 2)]
  Dtau <- matrix(data = exp(tail(psi, 6)),
                 nrow = 2, ncol = 3, byrow = FALSE)
  
  Sigma_eps <- c(exp(s0[1] + s0[2]*R))
  Sigma_eps_i <- Sigma_eps^-1
  
  B_lvls <- as.numeric(sapply(unique(colnames(B)), function(cn){
    length(unlist(strsplit(cn, split = ":", fixed = TRUE)))}))
  
  S_v <- c(Dtau[,B_lvls])
  S_v <- Diagonal(n = length(S_v), x = S_v)
  S_u <- S_v
  S_ui <- S_u^-1
  
  dS_u0 <- get.dSv(B_lvls)
  
  V <- Z %*% S_v %*% t(Z) + Diagonal(n = length(Sigma_eps), x = Sigma_eps)
  V_i <- solve(V)
  
  dSigma_eps_s0 <- Sigma_eps
  dSigma_eps_s1 <- exp(s0[1] + s0[2]*R)*R
  
  dV_s0 <- Diagonal(n = length(dSigma_eps_s0), x = dSigma_eps_s0)
  dV_s1 <- Diagonal(n = length(dSigma_eps_s1), x = dSigma_eps_s1)
  
  d2V_s01 <- dV_s1
  d2V_s00 <- dV_s0
  d2V_s11 <- exp(s0[1] + s0[2]*R)*R^2
  d2V_s11 <- Diagonal(n = length(d2V_s11), x = d2V_s11)
  
  d2mll.beta <- -t(X)%*%V_i%*%X
  
  d2mll.s00 <- as.numeric(-.5*(t(Y - X %*% beta) %*% V_i %*% dV_s0 %*% V_i %*% dV_s0 %*% V_i%*% (Y - X %*% beta)
                               - t(Y - X %*% beta) %*% V_i %*% d2V_s00 %*% V_i %*% (Y - X %*% beta)
                               + t(Y - X %*% beta) %*% V_i %*% dV_s0 %*% V_i %*% dV_s0 %*% V_i%*% (Y - X %*% beta)
                               -sum((V_i %*% dV_s0)*t(V_i %*% dV_s0)) + sum(V_i * t(d2V_s00))))
  
  
  d2mll.s01 <- as.numeric(-.5*(t(Y - X %*% beta) %*% V_i %*% dV_s0 %*% V_i %*% dV_s1 %*% V_i%*% (Y - X %*% beta)
                               - t(Y - X %*% beta) %*% V_i %*% d2V_s01 %*% V_i %*% (Y - X %*% beta)
                               + t(Y - X %*% beta) %*% V_i %*% dV_s1 %*% V_i %*% dV_s0 %*% V_i%*% (Y - X %*% beta)
                               -sum((V_i %*% dV_s0)*t(V_i %*% dV_s1)) + sum(V_i * t(d2V_s01))))
  d2mll.s10 <- d2mll.s01
  
  d2mll.s11 <- as.numeric(-.5*(t(Y - X %*% beta) %*% V_i %*% dV_s1 %*% V_i %*% dV_s1 %*% V_i%*% (Y - X %*% beta)
                               - t(Y - X %*% beta) %*% V_i %*% d2V_s11 %*% V_i %*% (Y - X %*% beta)
                               + t(Y - X %*% beta) %*% V_i %*% dV_s1 %*% V_i %*% dV_s1 %*% V_i%*% (Y - X %*% beta)
                               -sum((V_i %*% dV_s1)*t(V_i %*% dV_s1)) + sum(V_i * t(d2V_s11))))
  
  d2mll.s2 <- matrix(data = c(d2mll.s00, d2mll.s01,
                              d2mll.s10, d2mll.s11), 
                     nrow = 2, ncol = 2, byrow = TRUE)
  
  dV_dtauis_V_i_times_dV_dtaui <- lapply(1:6, function(i){
    tau <- rep(0, 6)
    tau[i] <- 1
    Dtau <- matrix(data = tau, nrow = 2, ncol = 3, byrow = FALSE)
    dS_v <- Diagonal(n = length(c(Dtau[, B_lvls])), x = c(Dtau[, B_lvls]))
    dV_dtaui <- Z %*% dS_v %*% t(Z)
    
    V_i_times_dV_dtaui <- (V_i %*% dV_dtaui)
    
    return(list = c(dV_dtaui = dV_dtaui,
                    V_i_times_dV_dtaui = V_i_times_dV_dtaui))
  })
  
  d2mll.tau2 <- matrix(data = NA, nrow = 6, ncol = 6)
  for (i in 1:6) {
    
    dV_dtaui <- dV_dtauis_V_i_times_dV_dtaui[[i]]$dV_dtaui
    V_i_times_dV_dtaui <- dV_dtauis_V_i_times_dV_dtaui[[i]]$V_i_times_dV_dtaui
    
    for(j in i:6){
      
      dV_dtauj <- dV_dtauis_V_i_times_dV_dtaui[[j]]$dV_dtaui
      V_i_times_dV_dtauj <- dV_dtauis_V_i_times_dV_dtaui[[j]]$V_i_times_dV_dtaui
      
      d2mll.tau2[i,j] <- as.numeric(-.5*(t(Y - X %*% beta) %*% V_i %*% dV_dtaui %*% V_i %*% dV_dtauj %*% V_i %*% (Y - X %*% beta)
                                         +t(Y - X %*% beta) %*% V_i %*% dV_dtauj %*% V_i %*% dV_dtaui %*% V_i %*% (Y - X %*% beta)
                                         -sum(V_i_times_dV_dtaui*t(V_i_times_dV_dtauj))))*exp(tail(psi, 6)[i] + tail(psi, 6)[j])
      d2mll.tau2[j,i] <- d2mll.tau2[i,j]
    }
  }
  
  dmll.tau <- sapply(1:6, function(i) {
    dV_dtau <- dV_dtauis_V_i_times_dV_dtaui[[i]]$dV_dtaui
    dmll.tauj <- as.numeric(-.5*(-(t(Y - X %*% beta) %*% V_i) %*% (dV_dtau %*% (V_i %*% (Y - X %*% beta))) + sum(V_i * t(dV_dtau))))
    return(dmll.tauj)
  })
  
  dmll.tau <- dmll.tau*exp(tail(psi, 6))
  d2mll.tau2 <- d2mll.tau2 + diag(dmll.tau)
  
  d2mll.tau_sigma <- matrix(data = NA, nrow = 6, ncol = 2)
  for (i in 1:6){
    dV_dtaui <- dV_dtauis_V_i_times_dV_dtaui[[i]]$dV_dtaui
    V_i_times_dV_dtaui <- dV_dtauis_V_i_times_dV_dtaui[[i]]$V_i_times_dV_dtaui
    
    for(j in 1:2){
      if(j == 1){
        dV_sj <- dV_s0
      }else{
        dV_sj <- dV_s1
      }
      V_i_times_dV_sj <- (V_i %*% dV_sj)
      
      d2mll.tau_sigma[i,j] <- as.numeric(-.5*(t(Y - X %*% beta) %*% V_i %*% dV_dtaui %*% V_i %*% dV_sj %*% V_i %*% (Y - X %*% beta)
                                              +t(Y - X %*% beta) %*% V_i %*% dV_sj %*% V_i %*% dV_dtaui %*% V_i %*% (Y - X %*% beta)
                                              -sum(V_i_times_dV_dtaui*t(V_i_times_dV_sj))))
    }
  }
  
  d2mll.tau_sigma <- d2mll.tau_sigma*exp(tail(psi, 6))
  
  d2mll.tau_beta <- matrix(data = NA, nrow = ncol(X), ncol = 6)
  for (i in 1:6) {
    dV_dtaui <- dV_dtauis_V_i_times_dV_dtaui[[i]]$dV_dtaui
    d2mll.tau_beta[,i] <- as.numeric(-t(X) %*% V_i %*% dV_dtaui %*% V_i %*% (Y - X %*% beta))*exp(tail(psi, 6)[i])
  }
  
  d2mll.s_beta <- matrix(data = NA, nrow = ncol(X), ncol = 2)
  for (j in 1:2) {
    if(j == 1){
      dV_sj <- dV_s0
    }else{
      dV_sj <- dV_s1
    }
    
    d2mll.s_beta[,j] <- as.numeric(-t(X) %*% V_i %*% dV_sj %*% V_i %*% (Y - X %*% beta))
  }
  
  H <- rbind(cbind(d2mll.beta, d2mll.s_beta, d2mll.tau_beta),
             cbind(t(d2mll.s_beta), d2mll.s2, t(d2mll.tau_sigma)),
             cbind(t(d2mll.tau_beta), d2mll.tau_sigma, d2mll.tau2))
  
  return(H)
}

####################
## get.mBIC
####################
## 
##  Description:  Calculate the marginal BIC.
##                See paper for details.
##
get.mBIC <- function(psi,
                     X,
                     Z,
                     Y,
                     R,
                     B){
  
  mll <- get.mll(psi,
                 X,
                 Z,
                 Y,
                 R,
                 B)
  
  mBIC <- -2*mll + length(psi)*log(length(Y))
  
  return(mll)
}


# get.nQ.try(psi = psi_c,
#            X = X,
#            Z = Z %*% B,
#            Y = Y,
#            R = R,
#            B_lvls = B0_lvls,
#            Euc = Euc0,
#            Euc2 = Euc20,
#            dEuc2_rep = dEuc20_rep,
#            dS_u = dS_u0) 
# 
# 
# get.gnQ.try(psi = psi_c,
#             X = X,
#             Z = Z %*% B,
#             Y = Y,
#             R = R,
#             B_lvls = B0_lvls,
#             Euc = Euc0,
#             Euc2 = Euc20,
#             dEuc2_rep = dEuc20_rep,
#             dS_u = dS_u0) 
# 
# 
# numDeriv::grad(func = get.nQ.try,
#                x = psi,
#                X = X,
#                Z = Z %*% B,
#                Y = Y,
#                R = R,
#                B_lvls = B0_lvls,
#                Euc = Euc0,
#                Euc2 = Euc20,
#                dEuc2_rep = dEuc20_rep,
#                dS_u = dS_u0)



