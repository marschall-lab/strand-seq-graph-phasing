
setwd("~/Documents/wd")



base_args <- c(
  "/opt/conda/lib/R/bin/exec/R"                                                                                                           ,
  "--no-echo"                                                                                                                             ,
  "--no-restore"                                                                                                                          ,
  "--vanilla"                                                                                                                             ,
  "--file=~/Documents/GitHub/strand-seq-graph-phasing/scripts/R/clustering_orient_strandstate.snakemake.R",
  "--args"
)
# NA19129 -----------------------------------------------------------------


args <-
  c(base_args,
    "--mem-counts"                                                                                                                          ,
    "sseq_alignment_counts/NA19129_sseq_mem_counts.csv"                                                                                     ,
    "--fastmap-counts"                                                                                                                      ,
    "sseq_alignment_counts/NA19129_sseq_fastmap_counts.csv"                                                                                 ,
    "--connected-components"                                                                                                                ,
    "gfa/ccs/NA19129_exploded_ccs.tsv"                                                                                                      ,
    "--segment-length-threshold"                                                                                                            ,
    "250000"                                                                                                                                ,
    "--expect-XY-separate"                                                                                                                  ,
    "True"                                                                                                                                  ,
    "--threads"                                                                                                                             ,
    "1"                                                                                                                                     ,
    "--output-marker-counts"                                                                                                                ,
    "haplotype_marker_counts/NA19129_haplotype_marker_counts.csv"                                                                           ,
    "--output-lib"                                                                                                                          ,
    "library_weights/NA19129_library_weights.csv"
  )

# Functions ---------------------------------------------------------------


export_clusters_test <- function() {
  cluster_palette <-
    cluster_df %>%
    distinct(cluster) %>%
    mutate(Color = rainbow(n())) %>%
    mutate(Color = ifelse(is.na(cluster), 'black', Color))

  ref <-
    pafr::read_paf(
      'reference_alignments/T2Tv11_hg002Yv2_chm13/NA19129_T2Tv11_hg002Yv2_chm13_ref-aln.paf',
      tibble = TRUE,
      include_tags = FALSE) %>%
    group_by(qname) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::rename(unitig = qname)

  cluster_df %>%
    left_join(cluster_palette, by = 'cluster') %>%
    left_join(homology_df) %>%
    left_join(ref, by='unitig') %>%
    select(unitig, everything()) %>%
    readr::write_csv('test_colors.csv')

  invisible(TRUE)
}

