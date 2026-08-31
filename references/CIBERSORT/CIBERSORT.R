#' Main functions
#'
#' The Main function of CIBERSORT
#' @param sig_matrix  sig_matrix file path to gene expression from isolated cells, or a matrix of expression profile of cells.
#'
#' @param mixture_file mixture_file file path to heterogenous mixed expression file, or a matrix of heterogenous mixed expression
#'
#' @param perm Number of permutations
#' @param QN Perform quantile normalization or not (TRUE/FALSE)
#' @param mixture_scale "linear", "log2" (without a pseudocount), or legacy "auto".
#' @import utils
#' @importFrom preprocessCore normalize.quantiles
#' @importFrom stats sd
#' @export
#' @examples
#' \dontrun{
#'   ## example 1
#'   sig_matrix <- system.file("extdata", "LM22.txt", package = "CIBERSORT")
#'   mixture_file <- system.file("extdata", "exampleForLUAD.txt", package = "CIBERSORT")
#'   results <- cibersort(sig_matrix, mixture_file)
#'   ## example 2
#'   data(LM22)
#'   data(mixed_expr)
#'   results <- cibersort(sig_matrix = LM22, mixture_file = mixed_expr)
#' }
cibersort <- function(sig_matrix, mixture_file, perm = 0, QN = TRUE,
                     mixture_scale = c("auto", "linear", "log2")){

  mixture_scale <- match.arg(mixture_scale)

  #read in data
  if (is.character(sig_matrix)) {
    X <- read.delim(sig_matrix, header=T, sep="\t", row.names=1, check.names = F)
    X <- data.matrix(X)
  } else {
    X <- sig_matrix
  }

  if (is.character(mixture_file)) {
    Y <- read.delim(mixture_file, header=T, sep="\t", row.names=1, check.names = F)
    Y <- data.matrix(Y)
  } else {
    Y <- mixture_file
  }


  #order
  X <- X[order(rownames(X)), , drop = FALSE]
  Y <- Y[order(rownames(Y)), , drop = FALSE]

  P <- perm #number of permutations

  # Preserve the legacy heuristic only for direct callers that have not
  # declared their scale. The template wrapper always supplies linear data.
  if (mixture_scale == "log2" || (mixture_scale == "auto" && max(Y) < 50)) {
    Y <- 2^Y
  }

  #quantile normalization of mixture file
  if(QN == TRUE){
    if (!requireNamespace("preprocessCore", quietly = TRUE)) {
      stop("Package 'preprocessCore' is required for quantile normalization. ",
           "Install it via BiocManager::install('preprocessCore').")
    }
    tmpc <- colnames(Y)
    tmpr <- rownames(Y)
    Y <- preprocessCore::normalize.quantiles(Y)
    colnames(Y) <- tmpc
    rownames(Y) <- tmpr
  }

  #intersect genes
  Xgns <- row.names(X)
  Ygns <- row.names(Y)
  YintX <- Ygns %in% Xgns
  Y <- Y[YintX, , drop = FALSE]
  XintY <- Xgns %in% row.names(Y)
  X <- X[XintY, , drop = FALSE]

  #standardize sig matrix
  X <- (X - mean(X)) / sd(as.vector(X))

  #empirical null distribution of correlation coefficients
  if (P > 0) {
    nulldist <- sort(doPerm(P, X, Y)$dist)
  }

  #print(nulldist)

  header <- c('Mixture',colnames(X),"P-value","Correlation","RMSE")
  #print(header)

  output <- matrix()
  itor <- 1
  mixtures <- dim(Y)[2]
  pval <- NA_real_ # No permutation p-value is available when perm = 0.

  #iterate through mixtures
  while (itor <= mixtures) {

    y <- Y[,itor]

    #standardize mixture
    y <- (y - mean(y)) / sd(y)

    #run SVR core algorithm
    result <- CoreAlg(X, y)

    #get results
    w <- result$w
    mix_r <- result$mix_r
    mix_rmse <- result$mix_rmse

    #calculate p-value
    if (P > 0) {
      pval <- 1 - (which.min(abs(nulldist - mix_r)) / length(nulldist))
    }

    #print output
    out <- c(colnames(Y)[itor],w,pval,mix_r,mix_rmse)
    if(itor == 1) {
      output <- out
    } else {
      output <- rbind(output, out)
    }

    itor <- itor + 1

  }

  #return matrix object containing all results
  obj <- rbind(header,output)
  obj <- obj[, -1, drop = FALSE]
  obj <- obj[-1, , drop = FALSE]
  obj <- matrix(as.numeric(unlist(obj)),nrow=nrow(obj))
  rownames(obj) <- colnames(Y)
  colnames(obj) <- c(colnames(X),"P-value","Correlation","RMSE")
  obj
}

# Core algorithm: nu-SVR with three nu values, returning the best model by RMSE.
CoreAlg <- function(X, y) {
  svn_itor <- 3
  res <- function(i) {
    if (i == 1) nus <- 0.25
    if (i == 2) nus <- 0.5
    if (i == 3) nus <- 0.75
    model <- e1071::svm(X, y, type = "nu-regression", kernel = "linear",
                        nu = nus, scale = FALSE)
    model
  }

  enableParallel()
  if (Sys.info()["sysname"] == "Windows") {
    out <- furrr::future_map(1:svn_itor, res)
  } else {
    if (svn_itor <= future::availableCores() - 2) {
      enableParallel(nThreads = svn_itor)
    } else {
      enableParallel()
    }
    out <- furrr::future_map(1:svn_itor, res)
  }

  nusvm <- rep(0, svn_itor)
  corrv <- rep(0, svn_itor)
  t <- 1
  while (t <= svn_itor) {
    weights <- t(out[[t]]$coefs) %*% out[[t]]$SV
    weights[which(weights < 0)] <- 0
    w <- weights / sum(weights)
    u <- sweep(X, MARGIN = 2, w, "*")
    k <- apply(u, 1, sum)
    nusvm[t] <- sqrt(mean((k - y)^2))
    corrv[t] <- cor(k, y)
    t <- t + 1
  }

  rmses <- nusvm
  mn <- which.min(rmses)
  model <- out[[mn]]
  q <- t(model$coefs) %*% model$SV
  q[which(q < 0)] <- 0
  w <- (q / sum(q))
  mix_rmse <- rmses[mn]
  mix_r <- corrv[mn]

  list(w = w, mix_rmse = mix_rmse, mix_r = mix_r)
}

# Empirical null distribution of correlation coefficients by permutation.
doPerm <- function(perm, X, Y) {
  itor <- 1
  Ylist <- as.list(data.matrix(Y))
  dist <- matrix()
  itorect <- function(Ylist, X) {
    yr <- as.numeric(Ylist[sample(length(Ylist), dim(X)[1])])
    yr <- (yr - mean(yr)) / sd(yr)
    result <- CoreAlg(X, yr)
    mix_r <- result$mix_r
    return(mix_r)
  }
  if (perm == 1) {
    dist <- itorect(Ylist = Ylist, X = X)
  } else {
    dist <- purrr::map(1:perm, ~itorect(Ylist = Ylist, X = X)) |>
      purrr::reduce(rbind)
  }
  list(dist = dist)
}

# Small helper to enable future parallel backend.
enableParallel <- function(nThreads = NULL) {
  if (is.null(nThreads)) {
    future::plan("multisession", workers = max(1, future::availableCores() - 2))
  } else {
    future::plan("multisession", workers = nThreads)
  }
}
