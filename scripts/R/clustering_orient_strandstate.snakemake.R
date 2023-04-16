# Functions ---------------------------------------------------------------

strip_range <-
  function(x) {
    stringr::str_remove(x, ':.*$')
  }


is_duplicate_pair <- function(x, y) {
  stopifnot(length(x) == length(y))
  pairs <- lapply(seq_along(x), function(i) sort(c(x[i], y[i])))
  return(duplicated(pairs))
}

make_wc_matrix <- function(watson, crick, lib, unitig) {
  stopifnot(length(unique(lengths(list(watson, crick, lib, unitig))))==1)
  
  # factors can cause issue ~ is the integer index or name used? When making
  # dimnames, it appears the name, when filling the matrix, it appears the
  # integer index is used.
  stopifnot(is.character(lib), is.character(unitig)) 
  stopifnot(is.integer(watson), is.integer(crick)) 
  stopifnot(!is_duplicate_pair(lib, unitig))
  
  w_frac <- (watson - crick) / (watson + crick)
  
  dimnames <-
    list(sort(unique(unitig)), sort(unique(lib)))
  
  mat <-
    matrix(
      nrow = length(dimnames[[1]]),
      ncol = length(dimnames[[2]]),
      dimnames = dimnames
    )
  
  for (i in seq_along(w_frac)) {
    mat[unitig[i], lib[i]] <- w_frac[i]
  }
  
  mat[is.na(mat)] <- 0
  
  return(mat)
  
}

extract_exact_matches <- function(x) {
  lines <- readLines(x)
  
  # There is a weird thing where all SQ tags following the first are preceded
  # by a '//'. No idea why
  
  # Filter valid lines. There are some line(s) that are just rogue '//'
  lines <-
    lines %>% 
    stringr::str_replace('^//SQ', 'SQ')
  
  is_header_line <- stringr::str_starts(lines,'SQ') 
  is_match_line <- stringr::str_starts(lines,'EM')
  
  valid_lines <-
    is_header_line | is_match_line
  
  lines <- lines[valid_lines]
  
  # Filter Exact Matches
  groupings <-
    cumsum(stringr::str_starts(lines,'SQ'))
  
  lines <- 
    split(lines, groupings)
  
  lines <-
    lines %>% 
    keep(function(x) length(x) == 2) # one header and one match line
  
  # SQ tags designate a particular read and the read length. 
  # SQ - Read Name - Read Length
  
  # EM tags designate an exact match for a segment of the read
  # EM - Start pos on read - End pos on Read - Number of matches - Matched unitig(s).
  
  # The number of matched unitigs that are printed in the output is controlled
  # by the -w parameter. A '*' in the matched unitigs position indicates that
  # there are more matches than the -w value used in the fastmap call.
  # Otherwise, up to -w matches are printed.
  
  # Dataframes
  
  header <-
    lines %>%
    map_chr(1)
  
  header <-
    readr::read_tsv(I(header),
                    col_names = c('tag', 'qname', 'qlen'),
                    col_types = '_c_')
  
  header <-
    header %>%
    mutate(group = names(lines))
  
  matches <-
    lines %>%
    map_chr(2)
  
  matches <-
    readr::read_tsv(
      I(matches),
      col_names = c('tag', 'lpos', 'rpos', 'n_matches', 'map_info'),
      col_types = '___ic'
    )
  
  matches <-
    matches %>%
    mutate(group = names(lines)) %>% 
    filter(n_matches == 1)
  
  # Tidy map info
  matches <-
    matches %>%
    tidyr::separate(
      col = map_info,
      into = c('unitig', 'strand_info'),
      sep = ':'
    ) %>%
    tidyr::separate(
      strand_info,
      into = c(
        'strand',
        "I don't know this field's meaning fastmap is poorly documented"
      ),
      sep = 1
    ) %>%
    select(unitig, strand, group)
  
  out <-
    inner_join(header,  matches, by='group') %>% 
    select(-group)
  
  return(out)
  
}


# Command Line ------------------------------------------------------------


## Args --------------------------------------------------------------------


args <- commandArgs(trailingOnly = FALSE)
setwd("~/Documents/wd")

args <-
  c(
    "/opt/conda/lib/R/bin/exec/R"                                                                                                                                   ,
    "--no-echo"                                                                                                                                                     ,
    "--no-restore"                                                                                                                                                  ,
    "--vanilla"                                                                                                                                                     ,
    "--file=/Users/henglinm/Documents/GitHub/strand-seq-graph-phasing/scripts/R/clustering_orient_strandstate.snakemake.R"                ,
    "--args"                                                                                                                                                        ,
    "--bam"                                                                                                                                                         ,
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20407.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20450.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20315.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20378.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20414.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20405.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20479.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20352.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20491.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20308.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20425.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20339.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20351.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20340.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20331.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20364.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20439.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20345.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20466.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20419.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20417.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20357.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20328.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20454.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20467.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20303.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20403.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20325.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20490.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20456.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20377.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20334.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20388.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20431.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20335.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20359.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20473.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20318.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20422.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20494.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20368.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20488.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20361.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20392.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20376.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20329.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20332.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20483.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20387.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20350.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20353.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20391.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20428.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20341.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20362.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20327.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20486.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20355.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20306.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20452.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20307.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20484.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20492.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20319.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20487.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20337.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20336.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20393.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20432.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20313.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20469.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20381.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20424.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20358.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20363.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20464.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20434.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20485.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20374.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20301.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20342.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20356.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20416.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20305.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20440.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20421.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20347.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20389.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20379.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20458.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20430.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20343.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20445.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20435.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x01PE20367.mdup.bam",
    "unmerged_bwa_ss_unitigs/unmerged_compressed_ss/NA24385_good_frac-50_rep-1_test/HWVTJAFXY_HG002x1xHG002x2_19s005136-1-1_Hasenfeld_lane1HG002x02PE20443.mdup.bam",
    "--gfa"                                                                                                                                                         ,
    "NA24385_good_frac-50_rep-1_test/assembly.gfa"                                                                                                                  ,
    "--homology"                                                                                                                                                    ,
    "mashmap/NA24385_good_frac-50_rep-1_test_mashmap_unitig-matches.mashssv"                                                                                        ,
    "--output-prefix"                                                                                                                                               ,
    "clustering_orientation_strandstate/NA24385_good_frac-50_rep-1_test/"                                                                                           ,
    "--threads"                                                                                                                                                     ,
    "8"                                                                                                                                                             ,
    "--segment-length-threshold"                                                                                                                                    ,
    "200000"                                                                                                                                                        ,
    "--log"                                                                                                                                                         ,
    "log/cos_NA24385_good_frac-50_rep-1_test.log"         
  )
## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    ## Input
    '--bam',
    '--gfa',
    '--homology',
    ## Output
    '--output-prefix',
    '--log',
    
    ## Params
    '--threads',
    '--segment-length-threshold'
  )

# Have to handle multiple inputs per tag
arg_idx <- sort(which(args %in% expected_args))
arg_idx <- c(arg_idx, length(args) + 1) # edge case of last tag

get_values <- function(arg, singular=TRUE){
  idx <- which(args == arg)
  stopifnot(length(idx) == 1)
  
  next_idx <- arg_idx[which.max(arg_idx > idx)]
  values <- args[(idx + 1):(next_idx - 1)]
  
  # More than one value? return list. One value ~ remove from list structure. It
  # is probably bad practice to return two different types like that, but most
  # arguments have single values and it makes the code less annoying to look at.
  if(singular) {
    stopifnot(length(values)==1)
    values <- values[[1]]
  } else {
    # stopifnot(length(values)>1)
  }
  
  return(values)
}


get_script_dir <- function() {
  needle <- '--file='
  script_path <- gsub(needle, '',  args[grepl(needle, args)])
  return(dirname(script_path))
}


print(args)

# Library -----------------------------------------------------------------

library(dplyr)
library(purrr)
# library(data.table)
library(Rsamtools)
library(doParallel)
library(matrixStats)
library(contiBAIT)

source(file.path(get_script_dir(), "module_utils/utils.R"))

reference <-
  pafr::read_paf('NA24385_T2Tv11_hg002Yv2_chm13_ref-aln.paf', tibble=TRUE, include_tags=FALSE) %>% 
  group_by(qname) %>% 
  slice_max(alen, with_ties = FALSE)
# Parsing -----------------------------------------------------------------

## Input

input.alignment.files <- get_values("--bam", singular=FALSE)
gfa <- get_values("--gfa", singular=TRUE)
mashmap <- get_values("--homology", singular=TRUE)

## Parameters
numCPU <- as.numeric(get_values("--threads"))
segment_length_threshold <- as.numeric(get_values('--segment-length-threshold'))

## Output
outputfolder <- get_values('--output-prefix')

parsed_args <-
  list(
    input.alignment.files = input.alignment.files,
    gfa=gfa,
    numCPU = numCPU,
    segment_length_threshold=segment_length_threshold,
    outputfolder = outputfolder
  )

print("Parsed Args:")
print(parsed_args)


# Import ------------------------------------------------------------------

## Graph Explosion ---------------------------------------------------------

# Other Script lol
source('~/Documents/GitHub/strand-seq-graph-phasing/scripts/R/explode_largest_component.R')


## Load Headers ------------------------------------------------------------

# All bam headers should be the same right? Only need one
unitig_lengths <- scanBamHeader(input.alignment.files[[1]], what='targets')[[1]]$targets
unitig_lengths_df <- tibble::enframe(unitig_lengths, name='unitig', value='length')

long_unitigs_df <- 
  unitig_lengths_df %>% 
  filter(length >= segment_length_threshold) %>% 
  select(unitig)

# contibait names
# Concatenate unitig with its length. Necessary for contiBAIT functions to not break.
contibait_names_df <- 
  unitig_lengths_df %>% 
  mutate(unitig_range = paste0(unitig, ':1-', length)) %>% 
  select(-length)


## Load Mashmap ------------------------------------------------------------
homology <- readr::read_delim(mashmap, delim='\t', col_names = c('unitig_1', 'unitig_2'))
# deduplicate



# filter out identical pairs, just with unitig_1 and unitig_2 switched
homology <-
  homology %>%
  filter(!is_duplicate_pair(unitig_1, unitig_2))

# tidy
homology <-
  homology %>% 
  mutate(bubble = paste0('bubble_', 1:n())) %>% 
  tidyr::pivot_longer(cols = c('unitig_1', 'unitig_2'), names_to = 'bubble_arm', values_to ='unitig') %>% 
  mutate(bubble_arm = stringr::str_replace(bubble_arm, '^unitig_', 'arm_'))

# Remove unititgs that are homologous to more than one other unitig.
homology <-
  homology %>% 
  group_by(unitig) %>% 
  filter(n() == 1) %>% 
  ungroup()

homology <-
  homology %>% 
  group_by(bubble) %>% 
  filter(n() == 2) %>% 
  ungroup()
## Count Alignments ---------------------------------------------------------


### Count BWA Alignments ------------------------------------------------------
print('counting w/c reads...\n')
lib.names <- map_chr(input.alignment.files, function(x) gsub('.mdup.bam$', '', basename(x)))

cl <- makeCluster(numCPU)
registerDoParallel(cl)

counts <- foreach(bam=input.alignment.files, .packages=c('Rsamtools', 'dplyr')) %dopar%{
  cat('counting directional reads in', basename(bam), '\n')
  aln = scanBam(file = bam,
                param = ScanBamParam(
                  what = c('qname', 'rname', 'strand', 'pos', 'qwidth', 'mrnm', 'isize', 'mapq'),
                  flag = scanBamFlag(
                    isSupplementaryAlignment = FALSE,
                    isSecondaryAlignment = FALSE,
                    isDuplicate = FALSE,
                    # for the purpose of determining qname/direction mapping, having both mates simply provides redundant information?
                    isFirstMateRead = TRUE,
                    isProperPair = TRUE
                  )
                ))
  
  aln <- as_tibble(aln[[1]])
  
  # filter out alignments to short unitigs

  
  aln <- 
    aln %>%  
    mutate(rname = as.character(rname)) %>% # default is as.factor import?
    filter(rname %in% long_unitigs_df$unitig)
  
  
  # Keep only reads that successfully aligned
  aln <- 
    aln %>% 
    filter(strand %in% c('+','-'))
  
  # filter any ss_reads that map to too many unitigs
  # dplyr::filter with lots of groups can be very slow -> duplicated is faster
  duplicated_qnames <-
    with(aln, qname[duplicated(qname)])
  
  if(any(duplicated_qnames)) {
    aln <-
      aln %>%
      filter(!(qname %in% duplicated_qnames))
  }
  
  # to simplify, only keep alignments where both mates land on the same rname
  aln <- 
    aln %>% 
    filter(rname == mrnm)
  
  # Counting
  out <-
    aln %>%
    # filter(mapq > 0) %>% 
    group_by(rname) %>%
    summarise(c = sum(strand == '+'), w = sum(strand == '-')) %>%
    ungroup()
  
  return(out)
  
  
}
stopCluster(cl)

# TODO experiment with filtering on the total/average/something number of
# alignments?
counts <-
  counts %>%
  set_names(lib.names) %>% 
  bind_rows(.id = 'lib') %>%
  dplyr::rename(unitig = rname) 

counts <-
  counts %>%
  tidyr::complete(lib, unitig, fill=list(c=0, w=0))

counts <-
  counts %>% 
  mutate(n = c+w)




# semi join on long unitigs?
 
# contiBait Clustering ----------------------------------------------------

## W Fraction Matrix -------------------------------------------------------

# Name alignments with range. Necessary for contiBAIT functions to not break.
wfrac.matrix <-
  counts %>% 
  left_join(contibait_names_df, by='unitig') %>% 
  with(make_wc_matrix(w, c, lib, unitig_range))

# 
# 
# strandFrequencyList <- strandSeqFreqTable(
#   input.alignment.files,
#   qual = 0,
#   pairedEnd = TRUE,
#   BAITtables = TRUE
# )
# 
# plotWCdistribution(strandFrequencyList$strandTable)

## ContiBAIT QC ------------------------------------------------------------

strand.freq <- StrandFreqMatrix(wfrac.matrix)


# getMethod(plotWCdistribution, "StrandFreqMatrix")
# debugonce(plotWCdistribution, signature = "StrandFreqMatrix")
plotWCdistribution(strand.freq)
# 


# debugonce(preprocessStrandTable, signature = 'StrandFreqMatrix')
strand.states <- preprocessStrandTable(strand.freq, filterThreshold = 0.7)

# filterThreshold = 0.7 is a stricter values than the default (filterThreshold =
# 0.8) because, with larger, higher quality initial assemblies, I expected that
# more unitigs will be more well behaved

# TODO experiment with setting filterThreshold even more strictly (maybe ~ .65,
# .6) and splitting large unitigs that fail QC to see if segments will pass WC

# TODO consider setting filter threshold according to the data. EG, given the
# distribution of aligned counts per unitig, what is the expected variance in
# unitig fraction when the strand state is WC? Set the parameter according to
# that variance.

# TODO The counts fo watson and crick following binomial distribution for when
# the state is WC appears to be underdispered. Experiment with overdispersed
# bionmials (beta binomials?)


## Contibait Chromosome Clustering -----------------------------------------

# TODO incorporate  mapping quality as well? Maybe not, filtering based on
# read quality seems to produce weird results with contiBAIT, so maybe
# clustering with read quality is also not a good idea

# weight the unitigs that have more alignments to be more likely to be
# selected earlier by the contiBAIT clustering algorithm.
mean_coverage <- 
  counts %>% 
  group_by(unitig) %>% 
  summarise(coverage = mean(w+c))

weights <-
  mean_coverage %>% 
  left_join(contibait_names_df, by='unitig') %>% 
  with(setNames(coverage, unitig_range))

# arrange weights by order in wfrac matrix
weights <- weights[rownames(strand.states$strandMatrix)]

# getMethod(clusterContigs, "StrandStateMatrix")
# debugonce(clusterContigs, signature = 'StrandStateMatrix')
clust <-
  clusterContigs(strand.states$strandMatrix,
                 recluster = 100,
                 randomWeight = weights,
                 clusterBy = 'hetero',
                 verbose = FALSE)

## Detect Haploid Chromosomes ----------------------------------------------

# debugonce(findSexGroups, signature = c('LinkageGroupList', 'StrandStateMatrix'))
clust <- findSexGroups(clust, strand.states$strandMatrix)

if(length(grep('^sex', names(clust), value = TRUE)) > 1) {
  # TODO what if multiple groups of haploid detected chromosomes?
  warning('More than 1 cluster has been identified as a haploid chromosome cluster')
}


## Enframe -----------------------------------------------------------------

cluster_df <-
  set_names(clust@.Data, clust@names) %>% 
  map(~tibble(unitig = strip_range(.x))) %>% 
  bind_rows(.id = 'cluster')

cluster_df <-
  cluster_df %>% 
  right_join(long_unitigs_df, by='unitig')


# PAR Detection -----------------------------------------------------------

haploid_components <-
  cluster_df %>%
  left_join(components_df, by='unitig') %>% 
  group_by(component) %>% 
  filter(any(grepl('^sex', cluster))) %>% 
  pull(component) %>% 
  unique()

haploid_component_unititgs <-
  components_df %>% 
  filter(component %in% haploid_components) %>% 
  pull(unitig)

par_clusters <- 
  cluster_df %>% 
  filter(unitig %in% haploid_component_unititgs) %>% 
  filter(!grepl('^sex', cluster)) %>% 
  filter(!is.na(cluster)) %>%
  pull(cluster) %>% 
  unique()

if(length(par_clusters) > 1) {
  warning("more than 1 PAR cluster, haven't thought about what will happen in this case")
}

par_unitigs <-
  cluster_df %>% 
  filter(cluster %in% par_clusters) %>% 
  pull(unitig)

# Call WC Libraries -------------------------------------------------------

strand_state_df <-
  strand.states$strandMatrix@.Data %>% 
  as_tibble(rownames = 'unitig_range')

strand_state_df <-
  strand_state_df %>% 
  tidyr::pivot_longer(-unitig_range, names_to = 'lib', values_to = 'state')

strand_state_df <-
  strand_state_df %>% 
  mutate(unitig = strip_range(unitig_range)) %>% 
  left_join(cluster_df, by='unitig')

# 2 ~ heterozygous call
het_fracs <-
  strand_state_df %>% 
  group_by(cluster, lib) %>% 
  summarise(het_frac = mean(state==2)) %>%  
  ungroup()

# Chosen based on looking at a plot
wc_threshold <- 0.9

wc_libraries_df <-
  het_fracs %>% 
  dplyr::filter(het_frac >= wc_threshold) 

ww_libraries_df <-
  het_fracs %>% 
  anti_join(wc_libraries_df, by=c('cluster', 'lib'))


# Label Propagation -------------------------------------------------------

# TODO add a check on the proportion of a component that is clustered? EG. A
# comopnent with only one small node clustered. Can that one node cluster the
# whole component? Or should a component be "mostly clustered" in order too add
# to the cluster in this way?

# TODO, on a related note, if there is a very small cluster attached to a very
# large one, should the largest one take over the smaller ones? Is this better
# handled with some sort of counts threshold, to better control rogue clusters
# in the telomeres and centromeres? ---  I have taken a look and unfortunately a
# countsfilter  is unlikely to work, as some cases will have > 200 alignments
# for many libraries and still cluster separately. 

# TODO What to do about components that are both large and have unassigned unitigs? Should I still attempt to assign them somehow?

### Remove Micro Clusters ---------------------------------------------------

cluster_component_fractions <-
  components_df %>% 
  left_join(cluster_df, by='unitig') %>% 
  left_join(unitig_lengths_df, by = 'unitig')

cluster_component_fractions <-
  cluster_component_fractions %>% 
  group_by(component, cluster) %>% 
  summarise(length = sum(length)) %>% 
  group_by(component) %>% 
  mutate(perc_length = length/sum(length)) %>% 
  ungroup()

# clusters that are contained on only one component
one_component_clusters <-
  cluster_component_fractions %>% 
  filter(!is.na(cluster)) %>% 
  group_by(cluster) %>% 
  filter(n() == 1) %>% 
  pull(cluster) %>% 
  unique()

# Arbitrary threshold
component_fraction_threshold <- 0.05
  
  # TODO filter to only components with multiple clusters? What about a
  # component consisting of only micro clusters? Worry about this stuff later is
  # things start breaking lol. One case maybe to worry about: the PAR taking
  # over both X/Y chromosomes, if no haploid chromosomes detected?

isolated_micro_clusters <-
  cluster_component_fractions %>%
  filter(!is.na(cluster)) %>%
  filter(!(component %in% haploid_components))  %>% # Don't want to accidentally remove the PAR
  filter(cluster %in% one_component_clusters) %>% 
  filter(perc_length <= component_fraction_threshold) %>% 
  pull(cluster)

cluster_df <-
  cluster_df  %>% 
  mutate(cluster = ifelse(cluster %in% isolated_micro_clusters, NA, cluster))

### Propagate One-Cluster-Component Clusters --------------------------------

one_cluster_components_df <-
  cluster_df %>% 
  left_join(components_df, by = 'unitig') %>% 
  filter(!is.na(cluster)) %>% 
  filter(!(component %in% haploid_components))

one_cluster_components_df <-
  one_cluster_components_df %>% 
  distinct(component, cluster) %>% 
  group_by(component) %>% 
  filter(n() == 1) %>% 
  ungroup()

one_cluster_component_unitigs <-
  components_df %>% 
  inner_join(one_cluster_components_df, by='component')

# lookup vector
one_cluster_component_unitigs <-
  with(one_cluster_component_unitigs, set_names(cluster, unitig))

# coalesce
cluster_df <-
  cluster_df %>%
  mutate(cluster = ifelse(
    unitig %in% names(one_cluster_component_unitigs),
    one_cluster_component_unitigs[unitig],
    cluster
  )) 


## Propagate PAR -----------------------------------------------------------


xy_cluster_label <- 'LGXY'

if (length(par_clusters) > 0) {
  cluster_df <-
    cluster_df %>%
    mutate(cluster = ifelse(
      unitig %in% c(haploid_component_unititgs, par_unitigs), xy_cluster_label, cluster))
  
  wc_libraries_df <-
    wc_libraries_df %>% 
    mutate(cluster = ifelse(cluster %in% par_clusters, xy_cluster_label, cluster))
  
  ww_libraries_df <-
    ww_libraries_df %>% 
    mutate(cluster = ifelse(cluster %in% par_clusters, xy_cluster_label, cluster))
}



## Split and recluster if any low quality? -------------------------------_

# TODO, EG ala makeChrTable(splitBy = 1e6), though this would be tiresome to
# implement...


# Orientation Detection w/ Inverted Unitigs -------------------------------

# Add inverted version of every unitig to dataset. Guarantees that there will
# unitigs in both orientations when clustering 

# TODO need to double check /rerun the strand orientation clustering. If The
# haploid and PAR clusters are oriented separately ~ chance that they could
# become unsynchronized. Is there a way to fix this? Do it at the haploid clustering step?

## Invert Unitigs ----------------------------------------------------------

counts <-
  counts %>% 
  mutate(unitig_dir = unitig)

counts <-
  counts %>% 
  bind_rows(
    counts %>% 
      mutate(unitig_dir = paste0(unitig, '_inverted'),
             c = n-c,
             w = n-w)
  )

# cluster_counts_df <- 
#  bind_rows(cluster_counts, inverted_cluster_counts)

## Orientation Detection -------------------------------------------------

# Use WW and CC libraries only for this step? Or does it not really matter? I
# guess the first principal component is picking out all the variation from the
# WW/CC libraries, and filtering to only the WW/CC libraries will therefore only
# have a minimal effect? Do it anyways, it is just a semi_join? This could be a
# vital step for clustering the haploid segments? It appears that for the
# haploid clusters it works well to improve the explained variance of the first
# PC, and makes it comparable to a diploid cluster.

# Is working in a continuous space (eg, principal components of W-fraction)
# instead of a discretized space (eg, k-modes clustering on strand states)
# better for tasks when there are strong prior assumptions?

prcomps <-
  counts %>% 
  left_join(cluster_df, by='unitig') %>% 
  semi_join(ww_libraries_df, by=c('lib', 'cluster')) %>% 
  split(.$cluster) %>% 
  map(function(x) with(x, make_wc_matrix(w,c,lib,unitig_dir))) %>% 
  map(prcomp)


# It appears that I need to do clustering only on the first PC, which
# corresponds to orientation, as otherwise, other PCs can influence the results
# and lead to grouping of a unitig and its inversion together.

strand_orientation_clusters_df <-
  map(prcomps, 'x') %>%
  map(function(x)
    tibble(unitig_dir = rownames(x), strand_cluster = sign(x[, 'PC1']))) %>%
  bind_rows(.id = 'cluster')

strand_orientation_clusters_df <-
  strand_orientation_clusters_df %>% 
  mutate(unitig = stringr::str_remove(unitig_dir, '_inverted$'))

# Warning that checks that every unitig and its invert are in opposite clusters
 bad <-
  strand_orientation_clusters_df %>% 
  group_by(unitig, strand_cluster) %>% 
  filter(n() != 1)

 # TODO, gather warnings in a log somewhere or something?
if(any(bad)) {
  warning('A warning about inversions clustered together or something')
}



# Phase Diploid Chromosomes -----------------------------------------------

# TODO filter haploid chromosomes out of this. Or not WTF

### Count fastmap Alignments ------------------------------------------------

# TODO fastmap input parameter
dir <- 'exact_match/NA24385_good_frac-50_rep-1_test/'
exact_match_files <- list.files(dir, full.names = TRUE)
lib_names <- sapply(exact_match_files, function(x) gsub('_maximal_unique_exact_match.tsv$', '', basename(x)))



cl <- makeCluster(numCPU)
registerDoParallel(cl)

exact_counts <- foreach(x=exact_match_files, .packages=c('purrr', 'dplyr')) %dopar%{
  out <- extract_exact_matches(x)
  
  out <-
    out %>%
    filter(unitig %in% long_unitigs_df$unitig) %>% 
    group_by(unitig) %>%
    summarise(c = sum(strand == '+'), w = sum(strand == '-')) %>%
    ungroup()
  
  return(out)
}
stopCluster(cl)


exact_counts <-
  exact_counts %>%
  set_names(lib_names) %>% 
  bind_rows(.id = 'lib') %>%
  tidyr::complete(lib, unitig, fill=list(c=0, w=0))

exact_counts <-
  exact_counts %>% 
  mutate(n = c+w)

## Orient Counts -----------------------------------------------------------


# concatenate with inverted, then filter based on strand orientation

exact_counts <-
  exact_counts %>% 
  mutate(unitig_dir = unitig)

exact_counts <-
  exact_counts %>% 
  bind_rows(
    exact_counts %>% 
      mutate(unitig_dir = paste0(unitig, '_inverted'),
             c = n-c,
             w = n-w)
  )

exact_counts <-
  exact_counts %>%
  semi_join(
    filter(strand_orientation_clusters_df, strand_cluster == 1),
    by = c('unitig', 'unitig_dir')
  ) 

exact_counts <-
  exact_counts %>%
  select(-unitig_dir)


## Filter to WC libraries for each cluster ---------------------------------

bubble_coverage <-
  exact_counts %>% 
  left_join(cluster_df, by='unitig') %>% 
  semi_join(wc_libraries_df, by=c('lib', 'cluster')) %>% 
  # TODO not innerjoin here? Need to left join and check that all clusters have
  # homology? Maybe make cluster a factor to handle this or something?
  inner_join(homology, by='unitig')

bubble_coverage <-
  bubble_coverage %>% 
  filter(stringr::str_detect(cluster, '^sex', negate=TRUE))

## Count Exact Alignments to Bubbles ---------------------------------------


mutate_coverage <- function(x, coverage_minimum = 5, coverage_ratio_threshold = 0.75 ){
  x %>% 
    mutate(crick_coverage_ratio = ifelse(n >= coverage_minimum, c / n, NA)) %>%
    mutate(
      covered_by = case_when(
        crick_coverage_ratio >= coverage_ratio_threshold ~ 'crick',
        crick_coverage_ratio <= 1-coverage_ratio_threshold ~ 'watson',
        TRUE ~ NA
      )
    )
}

bubble_coverage <-
  bubble_coverage %>%
  mutate_coverage()

# This step will also filter out libraries that are evaluated as WC but that
# cover no bubbles
bubble_coverage <-
  bubble_coverage %>%
  filter(!is.na(covered_by))


## Create StrandphaseR Arrays --------------------------------------------

strandphaser_arrays <-
  bubble_coverage %>% 
  split(.$cluster) %>% 
  map(function(x) {
    
    dimensions <-
      with(x,
           list(
             lib = sort(unique(lib)),
             bubble = sort(unique(bubble)),
             orientation = c('watson', 'crick')
           ))
    
    bubble_coverage_array <-
      array(dim = map_int(dimensions, length), dimnames = dimensions)
    
    x %>% 
      distinct(lib, bubble, covered_by, unitig) %>% 
      pwalk(function(lib, bubble, covered_by, unitig){
        bubble_coverage_array[lib, bubble, covered_by] <<- unitig
      })
    
    return(bubble_coverage_array)
  })


## Sorting StrandphaseR Arrays ---------------------------------------------

calc_concensus_margin <- function(x, ...) {
  if(all(is.na(x))) {
    return(NA)
  }
  # lower is better
  counts <- table(x, ...)
  return(sum(counts) - max(counts)) # if no names(counts) %in% values, then warning and function returns Inf
}

swap_bubbles <- function(phaser_array, ix) {
  tmp <- phaser_array[ix, , 'watson']
  phaser_array[ix, , 'watson'] <- phaser_array[ix, , 'crick']
  phaser_array[ix, , 'crick'] <- tmp
  
  return(phaser_array)
}

strandphaser_sort_array <- function(phaser_array) {
  # x[lib, bubble, watson/crick] <- unitig
  n_libs <-
    dim(phaser_array)[1]
  
  lib_swapped <- 
    logical(length = n_libs) %>% 
    set_names(dimnames(phaser_array)[[1]])
  
  for(i in seq_len(n_libs)){
    
    concensus_score <-
      phaser_array %>% 
      apply(c(2,3), calc_concensus_margin) %>% 
      sum(na.rm = TRUE)
    
    swapped_concensus_score <-
      swap_bubbles(phaser_array, i) %>% 
      apply(c(2,3), calc_concensus_margin) %>% 
      sum(na.rm = TRUE)
    
    if(swapped_concensus_score < concensus_score) {
      print(swapped_concensus_score)
      phaser_array <- swap_bubbles(phaser_array, i)
      lib_swapped[i] <- TRUE
    } else{
      print(concensus_score)
    }

  }
  
  # rownames(phaser_array) <-  ifelse(lib_swapped,
  #                                   paste0('swapped_', rownames(phaser_array)),
  #                                   rownames(phaser_array))
  
  return(lib_swapped)
  
  
}


lib_swaps <- 
  strandphaser_arrays %>% 
  map(strandphaser_sort_array)


## Haplotype Marker Counts -------------------------------------------------


lib_swaps_df <-
  lib_swaps %>%
  map_dfr(function(x)
    tibble::enframe(x, name = 'lib', value = 'swapped'),
    .id = 'cluster')

marker_counts <-
  exact_counts %>% 
  left_join(cluster_df, by='unitig') %>% 
  semi_join(wc_libraries_df, by=c('lib', 'cluster')) %>% 
  left_join(lib_swaps_df, by=c('lib', 'cluster')) 

# NA values indicate WC libraries that mapped to no bubbles, or clusters with no bubbles

# effectively treat NA as swapped=FALSE for now?
marker_counts <-
  marker_counts %>%
  filter(!is.na(swapped)) %>% 
  mutate(c = ifelse(swapped & !is.na(swapped), n - c, c), 
         w = ifelse(swapped & !is.na(swapped), n - w, w))

marker_counts <-
  marker_counts %>% 
  group_by(unitig) %>% 
  summarise(c = sum(c), w=sum(w)) %>% 
  ungroup() 

marker_counts <-
  marker_counts %>% 
  left_join(cluster_df, by='unitig') %>% 
  left_join(homology, by='unitig') %>% 
  arrange(cluster, bubble, bubble_arm)


# 
# # Phase Haploid Chromosomes -----------------------------------------------
# 
# 
# # If no haploid, need an empty vector when looking at excluded unitigs
# haploid_unitigs <- c()
# 
# if(length(haploid_components) > 0) {
#   
#   ## Detect PAR --------------------------------------------------------------
#   
#   # This process is contingent on successfully detecting the PAR untiigs that
#   # are connected to the haploid X and Y portions of the graph.
#   
#   haploid_component_unititgs <-
#     components_df %>% 
#     filter(component %in% haploid_components) %>% 
#     pull(unitig)
#   
#   par_clusters <- 
#     cluster_df %>% 
#     filter(unitig %in% haploid_component_unititgs) %>% 
#     filter(!grepl('^sex', cluster)) %>% 
#     filter(!is.na(cluster)) %>%
#     pull(cluster) %>% 
#     unique()
#   
#   if(length(par_clusters) > 0) {
#     
#     # Filter WC Libraries -----------------------------------------------------
#     
#     par_wc_libs <- 
#       wc_libraries_df %>% 
#       filter(cluster %in% par_clusters) %>% 
#       pull(lib) %>% 
#       unique()
#     
#     
#     par_unitigs <-
#       cluster_df %>% 
#       filter(cluster %in% par_clusters) %>% 
#       pull(unitig)
#     
#     # For the haploid unititgs specifically, it appears that the evidence
#     # ratios (w to C) is roughly the same whether exact or non-exact alignments
#     # are used. For the diploid PAR unitigs, the evidence ratio is much higher
#     # with exact alignments. Therefore, can use exact alignments for the PAR
#     # unititgs, but non-exact alignments for the haploid unititgs to better
#     # ensure that small unitigs will get any counts.
#     
#     # This step is probably excessive, and you could get away with just using
#     # the exact counts, but, I already wrote this code investigating this so
#     # lol.
#     
#     # No you dumbass, the PAR is already phased by the regular procedure as a
#     # diploid cluster. Doing both can be useful for your own personal debugging
#     # to ensure that the PAR results are the same in both procedures, to ensure
#     # synchronization between the PAR and haploid parts.
#     
#     # Exact Counts for PAR ----------------------------------------------------
#     
#     
#     hap_exact_counts <-
#       exact_counts %>% 
#       filter(lib %in% par_wc_libs) %>% 
#       filter(unitig %in% haploid_component_unititgs | unitig %in% par_unitigs)
#     
#     hap_lib_swaps <-
#       lib_swaps_df %>% 
#       filter(cluster %in% par_clusters)
#     
#     # use PAR results for all unitigs
#     hap_exact_counts <-
#       hap_exact_counts %>% 
#       left_join(hap_lib_swaps, by=c('lib')) 
#     
#     hap_exact_counts <-
#       hap_exact_counts %>%
#       filter(!is.na(swapped)) %>% 
#       mutate(c = ifelse(swapped & !is.na(swapped), n - c, c), 
#              w = ifelse(swapped & !is.na(swapped), n - w, w))
#     
#     hap_exact_marker_counts <-
#       hap_exact_counts %>% group_by(unitig) %>% summarise(c=sum(c), w=sum(w))
#     
#     
#     # Counts for Haploid ------------------------------------------------------
#     
# 
#     hap_counts <-
#       counts %>%
#       filter(lib %in% par_wc_libs) %>%
#       filter(unitig %in% haploid_component_unititgs | unitig %in% par_unitigs)
# 
#     # TODO need to double check /rerun the strand orientation clustering. If The
#     # haploid and PAR clsuters are oriented separately ~ chance that they could
#     # become unsynchronized.
#     hap_counts <-
#       hap_counts %>%
#       semi_join(
#         filter(strand_orientation_clusters_df, strand_cluster == 1),
#         by = c('unitig', 'unitig_dir')
#       )
# 
#     hap_counts <-
#       hap_counts %>%
#       left_join(hap_lib_swaps, by=c('lib'))
# 
#     hap_counts <-
#       hap_counts %>%
#       filter(!is.na(swapped)) %>%
#       mutate(c = ifelse(swapped & !is.na(swapped), n - c, c),
#              w = ifelse(swapped & !is.na(swapped), n - w, w))
#     
#     
#     hap_marker_counts <-
#       hap_counts %>% 
#       group_by(unitig) %>% 
#       summarise(c = sum(c), w = sum(w))
#     
# 
#     # Merge Marker Counts -----------------------------------------------------
#     
#     hap_marker_counts <-
#       hap_marker_counts %>% 
#       filter(!(unitig %in% par_unitigs))
#       bind_rows(
#         hap_exact_marker_counts %>% 
#           filter(unitig %in% par_unitigs)
#       )
#     
#     haploid_cluster_df <-
#       cluster_df %>%
#       left_join(contibait_names_df, by='unitig') %>% 
#       filter(stringr::str_starts(cluster, 'sex')) 
#     
#     # 
#     # haploid_strand_states <-
#     #   strand.states$strandMatrix
#     # 
#     # haploid_strand_states@.Data <- haploid_strand_states@.Data[haploid_cluster_df$unitig_range, par_wc_libs, drop=FALSE]
#     # 
#     # 
#     # # Invert Strand States ----------------------------------------------------
#     # inverted_unitig_range <-
#     #   strand_orientation_clusters_df %>%
#     #   left_join(contibait_names_df, by='unitig') %>% 
#     #   filter(unitig %in% haploid_cluster_df$unitig) %>% 
#     #   filter(strand_cluster == 1) %>% 
#     #   filter(grepl('inverted$', unitig_dir)) %>% 
#     #   pull(unitig_range) 
#     # 
#     # 
#     # if(length(inverted_unitig_range) > 0) {
#     #   # Invert the homozygous strand states for the inverted unitigs
#     #   # dumb trick to swap 3s and 1s in the matrix lol
#     #   is_1_or_2 <- haploid_strand_states@.Data[inverted_unitig_range, , drop = FALSE] %in% c(1, 2)
#     #   is_1 <-  haploid_strand_states@.Data[inverted_unitig_range, , drop = FALSE] == 1
#     #   
#     #   haploid_strand_states@.Data[inverted_unitig_range, ] <- is_1_or_2 + is_1 + 1
#     # }
#     # 
#     # 
#     # # Clustering on Homozygous Strand State -----------------------------------
#     # weights <-
#     #   mean_coverage %>%
#     #   left_join(contibait_names_df, by='unitig') %>% 
#     #   filter(unitig_range %in% rownames(haploid_strand_states))
#     # 
#     # weights <-
#     #   with(weights, set_names(coverage, unitig_range))
#     # 
#     # weights <- weights[rownames(haploid_strand_states)]
#     # 
#     # hap_cl <-
#     #   clusterContigs(haploid_strand_states,
#     #                  recluster = 100, # only a small number so why not a bunch? 
#     #                  randomWeight = weights,
#     #                  clusterBy = 'homo',
#     #                  verbose = FALSE)
#     # 
#     # # Cluster Selection -------------------------------------------------------
#     # 
#     # # If there are more than two clusters, keep only the most homozygous
#     # # clusters. If there are two or fewer clusters, nothings happens.
#     # 
#     # # Other ways to retain clusters? The ones with the most coverage? The
#     # # largest?
#     # clusters <-
#     #   set_names(hap_cl@.Data, hap_cl@names)
#     # 
#     # het_scores <-
#     #   map_dbl(clusters, function(x) {
#     #     
#     #     # Interesting bug. If you do mat %in% c(1,3), the result is a vector.
#     #     # However, if you do mat != 2, then the result is still a matrix lol
#     #     # what.
#     #     (haploid_strand_states[x, ,drop=FALSE] != 2) %>% 
#     #       apply(2, mean) %>% 
#     #       mean()
#     #   })
#     # 
#     # het_scores <-
#     #   het_scores %>% 
#     #   sort(decreasing = TRUE)
#     # 
#     # n_clusters_to_keep <-
#     #   ifelse(length(het_scores) == 1, 1, 2)
#     # 
#     # clusters_to_keep <-
#     #   names(het_scores)[seq_len(n_clusters_to_keep)]
#     # 
#     # clusters <- clusters[clusters_to_keep]
#     
#     # Phase -------------------------------------------------------------------
#     
#     # TODO Use exact counts to rescue highWC parts of the sex chromosomes? That
#     # would require changing how orientation is handled, as it is currently only
#     # calculated on not high-wc unititgs.
#   }  
#   
#   
#   # remove haploid clusters from rest of process.
#   # clust[haploid_unitigs] <- NULL
# }
# 
# # Excluded Unitigs --------------------------------------------------------
# 
# # At this point, all untigs that will be filtered out, should already be
# # filtered out. Can now record reason why each unitig was excluded
# short_rnames <-
#   names(rlengths[rlengths < segment_length_threshold])
# 
# high_wc_rnames <-as.character(strand.states$AWCcontigs@seqnames@values)
# 
# low_wc_rnames <- haploid_unitigs
# 
# excluded_unitigs <-
#   list(
#     tibble(unitig_name = short_rnames, exclusion_reason = paste('Length less than threshold:', as.character(segment_length_threshold))),
#     tibble(unitig_name = high_wc_rnames, exclusion_reason = 'Too many WC SSlib'),
#     tibble(unitig_name = low_wc_rnames, exclusion_reason = 'Too few WC SSlib')
#   ) %>%
#   bind_rows() %>% 
#   filter(!is.na(unitig_name))
# 
# stopifnot(with(excluded_unitigs, length(unitig_name) == length(unique(unitig_name))))
# 
# 
# 
# # Exclusion Check ---------------------------------------------------------
# clustered_rnames <- 
#   clusters_dt[invert==0, rname] %>% 
#   unique()
# 
# excluded_rnames <- excluded_unitigs$unitig_name
# stopifnot(setequal(rnames, c(clustered_rnames, excluded_rnames)))
# 
# 
# # Export Preparation ------------------------------------------------------
# # Matching old output 
# 
# old_clust <- unique(clusters_dt$clust)
# 
# clusters_dict <-
#   seq_along(old_clust) %>% 
#   setNames(old_clust) 
# 
# ## Unitig Clustering
# # names: rname, first_clust, second_clust
# unitig_clusters <- 
#   clusters_dt %>%
#   filter(invert==0) %>% 
#   distinct(rname, clust) %>% 
#   left_join(strand_orientation_clusters, by = c("clust", "rname"))
# 
# ## Cluster ID table
# # names: first_clust, second_clust, chrom_clust
# new_cluster_ids <-
#   unitig_clusters %>% 
#   distinct(clust) %>%
#   mutate(strand_clust = list(c(1,2))) %>% 
#   tidyr::unnest(cols = c(strand_clust)) %>% 
#   mutate(first_clust = 1:n()) %>% 
#   mutate(second_clust = ifelse(first_clust %% 2 == 0, first_clust-1, first_clust+1)) %>% 
#   mutate(chrom_clust = as.integer(as.factor(clust)))
# 
# unitig_clusters <-
#   unitig_clusters %>% 
#   left_join(new_cluster_ids, by = c("clust", "strand_clust")) %>% 
#   select(rname, first_clust, second_clust, chrom_clust)
# 
# ## List of WC libraries for each cluster
# ##### names: lib, thetawc, clust.forward, clust.backward
# temp <-
#   new_cluster_ids %>% 
#   distinct(clust, chrom_clust, first_clust,second_clust) %>% 
#   group_by(clust) %>% 
#   dplyr::slice(1) %>% 
#   ungroup()
# 
# wc_libraries <-
#   wc_libraries %>% 
#   left_join(temp, by = "clust")  %>% 
#   select(library, het_frac, first_clust, second_clust)
# 
# 
# 
# # Export ------------------------------------------------------------------
# 
# #Create a master output directory
# outputfolder <- file.path(outputfolder)
# dir_create_if_does_not_exist(outputfolder)
# 
# # TODO to refactor these output, it appear that one would have to learn/ remove
# # the whatshap_split step, as changing the variable names/order in these tables
# # can lead to unclear errors, and the scripts is very confusing. Consider
# # eliminating that step all together lol.
# 
# ## Unitig Clustering
# # names: rname, first_clust, second_clust, chrom_clust
# unitig_clusters %>% 
#   dplyr::select(rname, first_clust, chrom_clust) %>% 
#   dplyr::rename(`#rname` = rname) %>% 
#   fwrite(file=file.path(outputfolder, 'unitig_clusters.tsv'), 
#          sep='\t',
#          quote = F,
#          row.names = F)
# 
# ## Cluster ID table
# # names: first_clust, second_clust, chrom_clust
# new_cluster_ids  %>% 
#   select(first_clust, second_clust, chrom_clust) %>% 
#   filter(first_clust < second_clust) %>% 
#   fwrite( 
#     file=file.path(outputfolder, 'clust_partners.tsv'), 
#     sep='\t',
#     quote = F,
#     row.names = F)
# 
# ## List of WC libraries for each cluster
# ##### names: lib, thetawc, clust.forward, clust.backward
# fwrite(wc_libraries, 
#        file=file.path(outputfolder, 'wc_libraries.tsv'), 
#        sep='\t',
#        quote = F,
#        row.names = F)
# 
# 
# 
# 
# # excluded unitigs:
# excluded_unitigs %>% 
#   fwrite(
#     file=file.path(outputfolder, 'excluded_unitigs.tsv'), 
#     sep='\t',
#     quote = F,
#     row.names = F
#   )




