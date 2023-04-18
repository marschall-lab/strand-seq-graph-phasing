
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
