
# GFA ---------------------------------------------------------------------


read_segment_sizes_from_gfa <- function(filepath) {
  con <- file(filepath, "r")
  out <- integer()
  
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    if(stringr::str_starts(line, 'S')){
      line <- strsplit(line, '\t')[[1]]
      segment_name <- line[2]
      segment_length <- nchar(line[3])
      out <- c(out, set_names(segment_length, segment_name))
    }
  }
  
  close(con)
  return(out)
}

read_links_from_gfa <- function(filepath) {
  con <- file(filepath, "r")
  out <- character()
  
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    
    if(stringr::str_starts(line, 'L')){
      out <- c(out, line)
    }
  }
  
  close(con)
  return(out)
}



# Files -------------------------------------------------------------------


basename_no_ext <- function(x) {
  tools::file_path_sans_ext(basename(x))
}


# Tidying -----------------------------------------------------------------

get_prcomp_plotdata <- function(pca, rownames_col='rname') {
  out <- as_tibble(pca$x)
  pc_cols <- names(out)
  out[, rownames_col] <- rownames(pca$x)
  
  return(out[, c(rownames_col, pc_cols)])
}

# Predicates --------------------------------------------------------------


is_duplicate_pair <- function(x, y) {
  stopifnot(length(x) == length(y))
  pairs <- lapply(seq_along(x), function(i) sort(c(x[i], y[i])))
  return(duplicated(pairs))
}

all_are_unique <- function(x) {
  length(unique(x)) == length(x)
}

all_are_identical <- function(x) {
  length(unique(x)) == 1
}

all_have_same_length <- function(...) {
  list(...) %>%
    lengths() %>%
    all_are_identical()
}

# Colors, Plotting --------------------------------------------------------

invert_hex <- function(hex_code) {
  # Convert hex code to RGB values
  rgb_vals <- col2rgb(hex_code)
  
  # Invert RGB values
  inv_rgb_vals <- 255 - rgb_vals
  
  # Convert inverted RGB values to hex code
  inv_hex <- rgb(inv_rgb_vals[1, ], inv_rgb_vals[2, ], inv_rgb_vals[3, ], maxColorValue = 255)
  
  return(inv_hex)
}


# Distances ---------------------------------------------------------------


cosine_similarity_ <- function(x, y, min_overlaps=5, scale=TRUE) {
  stopifnot(all_have_same_length(x, y))
  
  overlapping_ix <- !is.na(x) & !is.na(y)
  
  if(sum(overlapping_ix) < min_overlaps) {
    return(NA)
  }
  
  x <- x[overlapping_ix]
  y <- y[overlapping_ix]
  
  if(scale) {
    x <- x / sqrt(sum(x^2))
    y <- y / sqrt(sum(y^2))
  }
  
  
  norm_x <- sqrt(sum(x^2))
  norm_y <- sqrt(sum(y^2))
  cosine_sim <- (x %*% y) / (norm_x * norm_y)
  
  return(cosine_sim)
}

cosine_similarity <-function(mat, ...) {
  similarity_mat <- 
    matrix(nrow=nrow(mat), ncol=nrow(mat))
  
  dimnames(similarity_mat) <- list(rownames(mat), rownames(mat))
  n <- length(rownames(mat))
  
  for(i in 1:n) {
    for(j in i:n) {
      sim <- cosine_similarity_(mat[i, ], mat[j, ])
      similarity_mat[i, j] <- sim
      similarity_mat[j, i] <- sim
    }
  }
  
  return(similarity_mat)
}

scale_vector <- function(x) {x / sqrt(sum(x^2))}


# Misc --------------------------------------------------------------------

pull_distinct <- function(x, col) {
  x %>% 
    distinct({{col}}) %>% 
    pull()
}

cut_range_n <- function(x, n_bins, max_value=max(x)) {
  
  stopifnot(max_value > n_bins)
  binwidth <- max_value %/% n_bins # ~floor(x/n), related to `%%`
  breaks <- binwidth * seq_len(n_bins-1)
  
  cuts <-
    map_int(x, function(xx) sum(xx > breaks)) + 1
  
  return(cuts)
  
}