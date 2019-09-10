# Common utility functions for scRNA-seq analysis and plotting

read10x <- function(run, suffix) {
#     barcode.loc <- file.path(run, "barcodes.tsv")
#     gene.loc <- file.path(run, "genes.tsv")
#     matrix.loc <- file.path(run, "matrix.mtx")
    barcode.loc <- list.files(run, pattern = 'barcodes.tsv', full.names = TRUE)
    gene.loc <- list.files(run, pattern = 'features.tsv|genes.tsv', full.names = TRUE)
    matrix.loc <- list.files(run, pattern = 'matrix.mtx', full.names = TRUE)

    data <- readMM(file = matrix.loc) %>% as("dgCMatrix")
    cell.names <- readLines(barcode.loc)
    cell.names <- gsub("-1$", "", cell.names)
    if (!missing(suffix)) {
        cell.names %<>% paste(suffix, sep = "_")
    }
    
    gene.names <- fread(gene.loc, header = FALSE)$V2
    row.names(data) <- gene.names
    colnames(data) <- cell.names

    
    return(as(data, "dgCMatrix"))
#     return(as(sumOverRowNames(data), "dgCMatrix"))
}

BuildSNNSeurat <- function (data.use, k.param = 30, prune.SNN = 1/15, nn.eps = 0) {
    my.knn <- nn2(data = data.use, k = k.param, searchtype = "standard", eps = nn.eps)
    nn.ranked <- my.knn$nn.idx
    snn_res <- ComputeSNN(nn_ranked = nn.ranked, prune = prune.SNN)
    rownames(snn_res) <- row.names(data.use)
    colnames(snn_res) <- row.names(data.use)
    return(snn_res)
}
environment(BuildSNNSeurat) <- asNamespace("Seurat")

fig.size <- function (height, width) {
    options(repr.plot.height = height, repr.plot.width = width)
}

# Generic theme settings for ggplot
GenericTheme <- theme(legend.text = element_text(size = 18),
                      legend.text.align = 0,
                      legend.title = element_text(size = 18, face = 'bold'),
                      legend.key = element_rect(fill = NA),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.background = element_blank(),
                      axis.line.x = element_line(size = 0.5, color = 'black'),
                      axis.line.y = element_line(size = 0.5, color = 'black'),
                      strip.text = element_text(size = 24, face = 'bold'),
                      axis.text = element_text(size = 18, colour = 'black'),
                      axis.title = element_text(size = 24, face = 'bold'),
                      plot.title = element_text(size = 32, face = 'bold'),
                      plot.subtitle = element_text(size = 28))

# Theme setting for creating SNE plots
sneTheme <- theme(axis.text = element_blank(),
                  axis.title = element_blank(),
                  axis.ticks = element_blank(),
                  axis.line = element_blank(),
                  plot.title = element_text(size = 28, hjust = 0.5),
                  panel.border = element_rect(fill = NA, size = 0.5))

# Theme settings for publication, larger font size, etc
PubTheme <- theme(legend.text = element_text(size = 24),
                  legend.text.align = 0,
                  legend.title = element_text(size = 24, face = 'bold'),
                  legend.key=element_rect(fill = NA),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line.x = element_line(size = 0.5, color = 'black'),
                  axis.line.y = element_line(size = 0.5, color = 'black'),
                  strip.text = element_text(size = 32, face = 'bold'),
                  axis.text = element_text(size = 28, colour = 'black'),
                  axis.title = element_text(size = 32, face = 'bold'),
                  plot.title = element_text(size = 44, face = 'bold'),
                  plot.subtitle = element_text(size = 36))

markerColorbar <- function(marker) {
  # Calculate order of magnitude by using the 99.99% value to prevent a few, 
  # real outlier signals from completely changing the colorscale
  max <- quantile(marker, probs = 0.9999)
  order.of.mag <- ceiling(log10(max))
  break.intervals <- c(0, 10^seq(0, order.of.mag))
  colorscale <- scale_colour_gradientn(trans = asinh_trans(), breaks = break.intervals, limits = c(0, 10^order.of.mag), colors = stain.cols,
                                       labels = expression(0, 10^0, 10^1, 10^2, 10^3, 10^4, 10^5, 10^6)[1:(order.of.mag + 2)],
                                       guide = guide_colorbar(barheight = 8, nbin = 1000))
  return(colorscale)
}


# Fancy formatting for scientifc numbers
fancy_scientific <- function(l, digit_precision = 3, hide_significand = FALSE,
                             return_char = FALSE) {
  # turn in to character string in scientific notation 
  l <- format(l, digits = digit_precision, scientific = TRUE)
  # Fix 0 case
  l <- gsub("0e\\+00","0",l)
  # quote the part before the exponent to keep all the digits 
  l <- gsub("^(.*)e", "'\\1'e", l)
  # remove + after exponent, if exists. E.g.: (3x10^+2 -> 3x10^2) 
  l <- gsub("e\\+","e",l)
  # turn the 'e+' into plotmath format 
  l <- gsub("e", "%*%10^", l)
  # convert 1x10^ or 1.000x10^ -> 10^
  if (hide_significand == TRUE) {
    l <- gsub("\\'1[\\.0]*\\'\\%\\*\\%", "", l)
  }
  # return as an expression or as a character vector if specified
  if (return_char == TRUE) {
    return(l)
  } else {
    return(parse(text = l))
  }
}

# Get density of points in 2 dimensions.
# @param x A numeric vector.
# @param y A numeric vector.
# @param n Create a square n by n grid to compute density.
# @return The density within each square.
get_density <- function(x, y, n = 100) {
  dens <- MASS::kde2d(x = x, y = y, n = n)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

insertElems = function(vect, pos, elems) {
  l = length(vect)
  j = 0
  for (i in 1:length(pos)){
    if (pos[i]==1)
      vect = c(elems[j+1], vect)
    else if (pos[i] == length(vect)+1)
      vect = c(vect, elems[j+1])
    else
      vect = c(vect[1:(pos[i]-1+j)], elems[j+1], vect[(pos[i]+j):(l+j)])
    j = j+1
  }
  return(vect)
}

minMax <- function(x)
{
    return((x- min(x)) /(max(x)-min(x)))
}
