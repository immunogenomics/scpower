#!/usr/bin/env/Rscript
source("/PHShome/cyf0/Projects/scpower/src/simulation.libs.R")

require("optparse")

option_list = list(
  make_option(c("--ndonors"), type = "numeric", default = NULL,
              help = "number of donors to simulate", metavar = "numeric"),
  make_option(c("--ncells"), type = "numeric", default = NULL,
              help = "number of cells per donor", metavar = "numeric"),
  make_option(c("--foldchange"), type = "numeric", default = NULL,
                help = "foldchange increase in cluster 2 proportion for cases", metavar = "numeric"),
  make_option(c("-o", "--out"), type = "character", default = getwd(), 
              help = "output directory [default = %default]", metavar = "character"),
  make_option(c("--filename"), type = "character", default = "default",
              help = "output filename [default = %default]", metavar = "character"),
  make_option(c("-v", "--verbose"), action = "store_true", default = TRUE,
              help = "Should the program print extra stuff out? [default %default]")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

# Create output filename
if (opt$filename == "default") {
  filename <- paste("simulation_", opt$ncells, "cells_", opt$ndonors,
                    "donors_", round(opt$foldchange, 2), "foldchange.rds", sep = "")
  } else { 
  filename <- opt$filename
}

# Read back options to user if verbose
if (opt$v) {
  cat(paste("Simulating dataset with", opt$ncells, "cells for", 
            opt$ndonors, "donors at", round(opt$foldchange, 2), "foldchange"))
  cat(paste("\nSaving to:", file.path(opt$o, filename), "\n"))
}

# mean_case, cov_case, mean_ctrl, cov_ctrl, A, and N are set directly in the 
# script for now
mean_case <- c(0.15, 0.2, 0.4, 0.15, 0.1)
mean_ctrl <- c(0.15, 0.2, 0.4, 0.15, 0.1)
cov_case <- diag(0.2, nrow = 5, ncol = 5)
cov_ctrl <- diag(0.2, nrow = 5, ncol = 5)
# Create type confusion matrix (all cell types should be classified correctly)
A <- matrix(0, 5, 5)
diag(A) <- 1
rownames(A) <- paste(1:5, 1, sep = "_")
colnames(A) <- as.character(1:5)
# Generate noise matrix
N <- matrix(1/5, 5, 5)

# Adjust mean_case for foldchange
cluster2 <- mean_case[2] * opt$foldchange
otherclusters <- mean_case[c(1, 3:5)] - (cluster2 - mean_case[2])/(length(mean_case) - 1)
new_case <- c(otherclusters[1], cluster2, otherclusters[2:4])

dataset <- simulateCells(ndonors = opt$ndonors, mean_case = new_case, mean_ctrl = mean_ctrl, 
                         cov_case = cov_case, cov_control = cov_ctrl, A = A, N = N, ncells_fixed = TRUE, 
                         ncells_mean = opt$ncells, status_proportion_fixed = TRUE, status_proportion = 0.5)
# Save simulation dataset
saveRDS(object = dataset, file = file.path(opt$out, filename))

