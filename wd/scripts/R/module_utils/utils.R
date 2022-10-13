
# GFA ---------------------------------------------------------------------


read_rnames_from_gfa <- function(x) {
  grep('^S',
       readLines(x),
       value = TRUE) %>%
    strsplit("\t") %>%
    vapply(function(x) x[[2]], character(1)) # grab the second field in every line
}

read_segment_lengths_from_gfa <- function(x) {
  split_lines <-
    grep('^S',
         readLines(x),
         value = TRUE) %>%
    strsplit("\t")
  
  tig_lengths <-
    split_lines %>%
    vapply(function(x) nchar(x[[3]]), integer(1))
  
  tig_names <-
    split_lines %>%
    vapply(function(x) x[[2]], character(1))
  
  names(tig_lengths) <- tig_names
  
  return(tig_lengths)
}


# Files -------------------------------------------------------------------


basename_no_ext <- function(x) {
  tools::file_path_sans_ext(basename(x))
}

dir_create_if_does_not_exist <- function(x) {
  if (!file.exists(x)) {
    dir.create(x)
    return(invisible(TRUE))
  } else {
    return(invisible(FALSE))
  }
}



# MAtrices, Dataframes, etc -----------------------------------------------


subset_rownames <- function(m, x) {
  
  rn <- rownames(m)
  idx <- rn %in% x
  out <- m[idx, ,drop=FALSE]
  rownames(out) <- rn[idx]
  return(out)
}

convert_dt_col_to_rownames <- function(dt, col) {
  rn <- dt[[col]]
  cols_to_keep <- names(dt) != col
  dt <- subset(dt, select=cols_to_keep)
  rownames(dt) <- rn
  return(dt)
}


dt_from_named_vector <- function(x, value='value', name='name') {
  # Probably a less dumb looking way to do this
  out <-
    data.table(value = x, name = names(x)) %>%
    setnames(old = c('value', 'name'), new = c(value, name))
  
  return(out)
}

rownamed_matrix_to_dt <- function(x, rownames_name='rownames') {
  
  rownames <- rownames(x)
  x <- as.data.table(x, keep.rownames = )
  if(!is.null(rownames)) {
    x[[rownames_name]] <- rownames
  }
  
  return(x)
}
