
setwd("~/Documents/wd")



base_args <- c(
  "/opt/conda/lib/R/bin/exec/R"                                                                                                           ,
  "--no-echo"                                                                                                                             ,
  "--no-restore"                                                                                                                          ,
  "--vanilla"                                                                                                                             ,
  "--file=~/Documents/GitHub/strand-seq-graph-phasing/scripts/R/clustering_orient_strandstate.snakemake.R",
  "--args"
)
# HG03248 -----------------------------------------------------------------


args <-
  c(base_args,
    "--mem-counts"                                                                                                                          ,
    "sseq_alignment_counts/HG03248_sseq_mem_counts.csv"                                                                                     ,
    "--fastmap-counts"                                                                                                                      ,
    "sseq_alignment_counts/HG03248_sseq_fastmap_counts.csv"                                                                                 ,
    "--connected-components"                                                                                                                ,
    "gfa/ccs/HG03248_exploded_ccs.tsv"                                                                                                      ,
    '--intermediate-output-dir',
    'intermediate_output/HG03248/',
    '--included-libs',
    # 'intermediate_output/HG03248/included_libraries.tsv',
    '--initial-clusters',
    # 'intermediate_output/HG03248/initial_clusters.tsv',
    '--final-clusters',
    # 'intermediate_output/HG03248/final_clusters.tsv',
    '--unitig-orientation',
    # 'intermediate_output/HG03248/unitig_orientation.tsv',
    '--counting-methods',
    # 'intermediate_output/HG03248/counting_methods.tsv',
    "--segment-length-threshold"                                                                                                            ,
    "250000"                                                                                                                                ,
    "--cluster-PAR-with-haploid"                                                                                                                  ,
    "True"                                                                                                                                  ,
    "--threads"                                                                                                                             ,
    "1"                                                                                                                                     ,
    "--output-marker-counts"                                                                                                                ,
    "haplotype_marker_counts/HG03248_haplotype_marker_counts.csv"                                                                           ,
    "--output-lib"                                                                                                                          ,
    "library_weights/HG03248_library_weights.csv"
  )



# WHarvey -----------------------------------------------------------------
# 
# args <-
#   c(base_args,
#     "--mem-counts"                                                                                                                          ,
#     '/Users/henglinm/Downloads/NA19983_count_haplotypes/sseq_alignment_counts/NA19983_sseq_mem_counts.csv'                                                                                     ,
#     "--fastmap-counts"                                                                                                                      ,
#     '/Users/henglinm/Downloads/NA19983_count_haplotypes/sseq_alignment_counts/NA19983_sseq_fastmap_counts.csv'                                                                                  ,
#     "--connected-components"                                                                                                                ,
#     '/Users/henglinm/Downloads/NA19983_count_haplotypes/gfa/ccs/NA19983_exploded_ccs.tsv'                                                                                                      ,
#     "--segment-length-threshold"                                                                                                            ,
#     "250000"                                                                                                                                ,
#     "--expect-XY-separate"                                                                                                                  ,
#     "True"                                                                                                                                  ,
#     "--threads"                                                                                                                             ,
#     "1"                                                                                                                                     ,
#     "--output-marker-counts"                                                                                                                ,
#     "haplotype_marker_counts/NA18936_haplotype_marker_counts.csv"                                                                           ,
#     "--output-lib"                                                                                                                          ,
#     "library_weights/NA18936_library_weights.csv"
#   )

# Functions ---------------------------------------------------------------


export_clusters_test <- function() {
  cluster_palette <-
    cluster_df %>%
    distinct(cluster) %>%
    mutate(Color = rainbow(n())) %>%
    mutate(Color = ifelse(is.na(cluster), 'black', Color))

  ref <-
    pafr::read_paf(
      'reference_alignments/T2Tv11_hg002Yv2_chm13/HG03248_T2Tv11_hg002Yv2_chm13_ref-aln.paf',
      tibble = TRUE,
      include_tags = FALSE) %>%
    group_by(qname) %>%
    slice_head(n = 1) %>%
    ungroup() %>%
    dplyr::rename(unitig = qname)

  cluster_df %>%
    left_join(cluster_palette, by = 'cluster') %>%
    left_join(ref, by='unitig') %>%
    select(unitig, everything()) %>%
    readr::write_csv('test_colors.csv')

  invisible(TRUE)
}

