require(MASS)
require(data.table)
require(magrittr)

#' Samples cells for donors using a multivariate normal a mis-clustering confusion matrix A and cell type distributions for cases and controls. 
#' @param ndonors
#' @param mean_case vector of means for cluster observations in cases. Will be normalized to add up to 1 
#' @param cov_case covariance matrix for vector of means for cluster observations in cases.
#' @param mean_ctrl same as above, but for controls. 
#' @param cov_ctrl same as above, but for controls
#' @param A confusion matrix describing probability of seeing a cell in cluster i (column) given that it is of type j (row). Must be non-negative and rowSum to 1. 
#' @param N A noise matrix describing the probability of seeing BAD (quality = 0) cells in cluster i (column) given that it is of type j (row). Must be non-negative and rowSum to 1. 
#' @param ncells_fixed should number of cells per donor be fixed or sampled with a Poisson? 
#' @param ncells_means either number of cells per donor or mean number of cells per donor (Poisson rate parameter)
#' @param status_proportion_fixed should number of cases and controls be fixed or sampled with a binomial?
#' @param status_proportion ratio of cases to controls
#' @param status_str character strings defining control and case status. Always put control first. 
simulateCells <- function(ndonors, mean_case, cov_case, mean_ctrl, cov_control, A, N,  ncells_fixed = FALSE, 
                            ncells_mean = 100, status_proportion_fixed = FALSE, status_proportion = .5, 
                            status_str = c("case", "ctrl")) {
  ## Check inputs
  if (any(abs(rowSums(A) - 1) > 1e-10)) {
    print("WARNING: Type confusion matrix may not be row normalized. Row normalizing it here.")
    A <- sweep(A, 1, rowSums(A), "/")
  }
  if (any(abs(rowSums(N) - 1) > 1e-10)) {
    print("WARNING: Noise matrix may not be row normalized. Row normalizing it here.")
    N <- sweep(N, 1, rowSums(A), "/")
  }
  if (all(dim(A) != dim(N))) {
    print("ERROR: Type confusion matrix and noise matrix are not the same shape")
    stop()
  }
  
  ## get number of cells per donor
  if (ncells_fixed) {
    ncells <- rep(ncells_mean, ndonors)
  } else {
    ncells <- rpois(ndonors, ncells_mean)
  }
  
  ## get donor status
  if (status_proportion_fixed) {
    status <- c(rep(1, round(ndonors * status_proportion)), 
                rep(2, ndonors - round(ndonors * status_proportion)))
  } else {
    status <- 1 + 1 * rbernoulli(ndonors, status_proportion)
  }
  
  ## sample a theta for each donor
  alpha_prior <- matrix(c(mean_case, mean_ctrl), length(mean_case), 2) %>% t
  lambda_prior <- log(alpha_prior)
  cov_prior <- list(cov_case, cov_ctrl)
  # alpha_prior <- matrix(c(alpha_case, alpha_ctrl), length(alpha_case), 2) %>% t
  # alpha_post <- alpha_prior %*% A
  status_counts <- table(status)
  theta <- do.call(rbind, lapply(1:nrow(lambda_prior), function(status_value) {
    exp(mvrnorm(n = status_counts[status_value],
                mu = lambda_prior[status_value,],  
                Sigma = cov_prior[[status_value]]))
  }))
  # row-normalize theta
  theta <- sweep(theta, 1, rowSums(theta), "/")
  status <- rep.int(as.integer(names(status_counts)), status_counts) ## reorder status to match theta
  
  ## Generate cell quality scores
  cell_quality <- do.call(c, lapply(1:ndonors, function(i) {
    rbeta(ncells[i], 8, 2)}))
  
  ## sample cells (cluster_ids) within each donor
  cluster_ids <- do.call(c, lapply(1:ndonors, function(i) {
    do.call(c, lapply(1:ncells[i], function (j) {
      # Use correct cell quality score depending on donor
      if (i == 1) {idx <- j } else { idx <- sum(ncells[1:(i-1)]) + j}
      cell_confusion <- cell_quality[idx] * A + (1 - cell_quality[idx]) * N
      cell_theta <- theta[i,] %*% cell_confusion
      sample(1:ncol(theta), 1, prob = cell_theta, replace = T)
    }))
  }))
  
  data.table(donor = rep.int(1:ndonors, ncells), 
             status = rep.int(status, ncells), 
             quality = cell_quality,
             cluster = cluster_ids)[
               data.table(status_str) %>% tibble::rowid_to_column("status"), 
               on = "status"
               ][
                 , `:=`(status = status_str, status_str = NULL)
                 ] %>% 
    data.table() ## for weird behavior with := assignment
}

#' Generates type-type confusion matrix
#' @param type integer values for cell types
#' @param subtype (optional) integer values for cell subtypes
#' @param cluster (optional) integer values for cell clusters
#' @param delta_type probability of confusing two cell types
#' @param delta_subtype probabilit of confusing two cell subtypes
#' @return A row normalized transition matrix from types (rows) to clusters (columns)
## NOTE: type, subtype, and cluster must be contain only consecutive integers (though not necessarily in order)
## NOTE: cluster is useless right now, but may be interesting later when 1 cluster can have 2 (sub)types.
makeAMat <- function(type, subtype=NULL, cluster=NULL, delta_type = 0.01, delta_subtype = 0.10) {
  type_order <- order(type)
  type <- type[type_order]
  if (is.null(subtype)) subtype <- unlist(sapply(unique(type), function(t) 1:(table(type)[t])))
  else subtype <- subtype[type_order]
  if (is.null(cluster)) cluster <- 1:length(type)
  else cluster <- cluster[type_order]
  
  if (any(data.table(type, subtype, cluster) %>%
          tidyr::unite(type, type, subtype) %>%
          with(table(type, cluster)) %>% 
          colSums() > 1)) {
    print(data.table(type, subtype, cluster) %>%
            tidyr::unite(type, type, subtype) %>%
            with(table(type, cluster)))
    stop("Cluster cannot contains more than one (sub)type")
  }
  if (length(setdiff(1:max(type), type)) > 0) stop("Types must be consecutive integers")
  if (length(setdiff(1:max(subtype), subtype)) > 0) stop("Subtypes must be consecutive integers")
  if (length(setdiff(1:max(cluster), cluster)) > 0) stop("Clusters must be consecutive integers")
  
  A <- tidyr::crossing(data.table(type, subtype), cluster = cluster) %>% unique %>%
    merge(data.table(type, subtype, cluster),
          by = c("cluster"), suffixes = c("", "_main")) %>% 
    dplyr::mutate(value = ifelse(type != type_main, delta_type, 
                                 ifelse(subtype != subtype_main, delta_subtype, 0))) %>% 
    dplyr::select(cluster, type, subtype, value) %>% 
    tidyr::unite(type, type, subtype) %>% 
    tidyr::spread(cluster, value)
  row.names(A) <- A$type
  A$type <- NULL
  
  ## distribute the rest of the prob density equally
  A <- apply(A, 1, function(A_row) {
    p <- (1 - sum(A_row)) / sum(A_row == 0)
    ifelse(A_row == 0, p, A_row)
  }) %>% t
  if (any(A < 0)) {
    stop("Matrix A cannot be negative, adjust delta values")
  }
  return(as.matrix(A))
}

