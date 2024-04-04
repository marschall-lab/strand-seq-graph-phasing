
# functions ---------------------------------------------------------------
set_tname_factor <- function(x) {
  factor(x, stringr::str_sort(unique(x), numeric=TRUE))
}

set_sample_factor <- function(x) {
  factor(x, stringr::str_sort(samples, numeric=TRUE))
}


# Command Line ------------------------------------------------------------
## Args --------------------------------------------------------------------

args <- commandArgs(trailingOnly = FALSE)

## Parsing Helper ----------------------------------------------------------
## Expected Args
expected_args <-
  c(
    # Plot Options
    '--included-diploid-chroms-regex',
    '--included-haploid-chroms-regex',
    '--included-hemiploid-chroms-regex',
    '--hta-chroms-regex',
    # Output Options
    '--plot-width',
    '--plot-height',
    # Parameters
    '--samples',
    '--segment-length-threshold',
    # Input
    '--haplotype-marker-counts',
    '--rukki-paths',
    '--ref-alignments',
    '--output'
  )

# Have to handle multiple inputs per tag
arg_idx <- sort(which(args %in% expected_args))
arg_idx <- c(arg_idx, length(args) + 1) # edge case of last tag

get_values <- function(arg, singular=TRUE){
  idx <- which(args == arg)
  stopifnot(length(idx) == 1)
  
  next_idx <- arg_idx[which.max(arg_idx > idx)]
  values <- args[idx:next_idx]
  # remove needles:
  values <- values[-c(1, length(values))]
  
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
library(pafr)
library(ggplot2)
# library(ggokabeito)


# Parsing -----------------------------------------------------------------

# samples <- c('HG00731', 'HG00732', 'HG00733', 'NA19434', 'NA24385')
samples <- get_values('--samples', singular = FALSE)

# slt <- 250e3
slt <- as.numeric(get_values('--segment-length-threshold', singular = TRUE))

# Input
# haplotype_marker_counts <- glue::glue('haplotype_marker_counts/{samples}_fudged_haplotype_marker_counts.csv')
haplotype_marker_counts <- get_values('--haplotype-marker-counts', singular=FALSE)

# rukki_paths <- glue::glue("rukki/{samples}/{samples}_rukki_paths.tsv")
rukki_paths <- get_values('--rukki-paths', singular=FALSE)

# ref_alns <-  '.none.'
# ref_alns <-  glue::glue('reference_alignments/T2Tv11_hg002Yv2_chm13/{samples}_T2Tv11_hg002Yv2_chm13_ref-aln.paf')
ref_alns <- get_values('--ref-alignments', singular=FALSE)

if(length(ref_alns) == 0) {
  ref_alns <- vector(mode='list', length=length(haplotype_marker_counts))
}

# output <- 'test_output_plots.pdf'
output <- get_values('--output', singular=TRUE)

out_width <- as.numeric(get_values('--plot-width', singular=TRUE))
out_height <- as.numeric(get_values('--plot-height', singular=TRUE))

included_dip_regex <- get_values('--included-diploid-chroms-regex', singular=TRUE)
included_hap_regex <- get_values('--included-haploid-chroms-regex', singular = TRUE)
included_hem_regex  <- get_values('--included-hemiploid-chroms-regex', singular = TRUE)
hta_regex <- get_values('--hta-chroms-regex', singular = TRUE)

# TODO condition handling on empty reference alignments
arg_lengths <-
  map_int(list(haplotype_marker_counts, rukki_paths, ref_alns), length) 

stopifnot(all(arg_lengths >= 1), n_distinct(arg_lengths) == 1) 
 
# Import Marker Counts ----------------------------------------------------

haplotype_marker_df <- 
  haplotype_marker_counts %>% 
  set_names(samples) %>% 
  map(readr::read_csv, show_col_types = FALSE, .progress = TRUE) %>% 
  list_rbind(names_to = 'sample') %>% 
  mutate(sample = set_sample_factor(sample))

haplotype_marker_df <-
  haplotype_marker_df %>% 
  mutate(color = case_when(length < slt ~ 'Threshold',
                           hap_1_counts == 0 & hap_2_counts == 0 ~ 'Uncounted',
                           TRUE ~ 'Attempted')) 

# Import Rukki Paths ------------------------------------------------------

import_rukki_paths <- function(x) {
  paths <- readr::read_tsv(x, col_types = 'ccc')
  
  # tidying paths
  paths <-
    paths %>%
    tidyr::separate_longer_delim(path, delim = ',') %>%
    group_by(name) %>%
    mutate(order = 1:n()) %>%
    ungroup() 
  
  paths <-
    paths %>%
    dplyr::rename(unitig = path) %>%
    mutate(unitig = gsub("[+-]+$", '', unitig))
  
  # run-marking
  paths <-
    paths %>%
    mutate(is_gap = grepl('\\[', unitig)) %>%
    group_by(name) %>%
    mutate(run = cumsum(is_gap)) %>%
    ungroup()
  
  paths <-
    paths %>%
    mutate(gap_size = ifelse(is_gap, stringr::str_extract(unitig, '[0-9]+'), NA)) %>%
    mutate(gap_size = as.integer(gap_size))
  
  return(paths)
}



rukki_paths_df <-
  rukki_paths %>% 
  set_names(samples) %>%
  map(import_rukki_paths, .progress = TRUE) %>% 
  list_rbind(names_to = 'sample') %>% 
  mutate(sample = set_sample_factor(sample))

# 
# rukki_paths_df <-
#   rukki_paths_df %>%
#   left_join(unitig_lengths_df, by = c('sample', 'unitig')) %>%
#   mutate(node_len = coalesce(qlen, gap_size)) %>%
#   select(sample, name, unitig, assignment, order, is_gap, run, node_len)

# 
# rukki_paths <-
#   rukki_paths %>%
#   mutate(is_not_gap = !is_gap)

## Rukki Hom and Issue -----------------------------------------------------


hom_nodes <-
  rukki_paths_df %>%
  filter(!is_gap) %>%
  group_by(sample, unitig) %>%
  mutate(n = n()) %>%
  filter('HAPLOTYPE1' %in% assignment & 'HAPLOTYPE2' %in% assignment) %>%
  ungroup() %>%
  filter(n == 2) %>%
  distinct(sample, unitig) %>%
  mutate(is_hom = TRUE)

issue_nodes <-
  rukki_paths_df %>%
  filter(is.na(assignment)) %>%
  distinct(sample, unitig) %>%
  mutate(is_issue = TRUE)

stopifnot(nrow(semi_join(hom_nodes, issue_nodes)) == 0)

rukki_paths_df <-
  rukki_paths_df %>%
  left_join(hom_nodes) %>%
  left_join(issue_nodes) %>%
  mutate(assignment = case_when(
    is_hom ~ 'HOM', 
    # is_issue ~ 'RUKKI_NA',
    TRUE ~ assignment
  )) %>% 
  mutate(assignment = as.factor(assignment))


# Import Reference Alignments ---------------------------------------------

# pafr is slow unfortunately
import_ref_aln <- function(x) {
  
  if(!is_null(x)) {
    ref_aln <- pafr::read_paf(x, tibble = TRUE, include_tags = FALSE)
  } else {
    ref_aln <-
      tibble(
        sample=character(),
        qname=character(),
        qlen=double(),
        qstart=double(),
        qend=double(),
        strand=character(),
        tname=set_tname_factor(character()),
        tlen=double(),
        tstart=double(),
        tend=double(),
        nmatch=double(),
        alen=double(),
        mapq=double(),
      )
  }

  
  ref_aln <-
    ref_aln %>% 
    dplyr::rename(unitig = qname)
  
  return(ref_aln)
}


ref_aln_df <-
  # ref_alns[0] %>% 
  ref_alns %>%
  set_names(samples) %>% 
  map(import_ref_aln, .progress = TRUE) %>% 
  list_rbind(names_to = 'sample') %>% 
  mutate(tname = set_tname_factor(tname)) %>%
  mutate(sample = set_sample_factor(sample))


## Included Chromosomes ----------------------------------------------------

tnames <-
  ref_aln_df %>% 
  distinct(tname) %>% 
  pull(tname)

# TODO make these arguments passable from the command line?
included_diploid_chrs <-
  stringr::str_subset(tnames, included_dip_regex) %>% 
  stringr::str_sort(numeric = TRUE)

# sex chrs ususlly haploid
included_haploid_chrs <-
  stringr::str_subset(tnames, included_hap_regex) %>% 
  stringr::str_sort(numeric = TRUE)

included_hemiploid_chrs <-
  stringr::str_subset(tnames, included_hem_regex) %>% 
  stringr::str_sort(numeric = TRUE)

included_chrs <-
  c(included_diploid_chrs, included_haploid_chrs, included_hemiploid_chrs)

hard_to_assemble_chrs <-
  stringr::str_subset(tnames, hta_regex) %>% 
  stringr::str_sort(numeric = TRUE)

included_chrs_df <-
  tibble(tname = set_tname_factor(included_chrs)) %>% 
  mutate(ploidy = case_when(
    tname %in% included_diploid_chrs ~ 'Diploid',
    tname %in% included_haploid_chrs ~ 'Haploid',
    TRUE  ~ 'Hemiploid')) %>% 
  # mutate(is_diploid = ploidy == 'diploid') %>% 
  mutate(is_hard_to_assemble = tname %in% hard_to_assemble_chrs)

# Derived DFs -------------------------------------------------------------


all_unitigs_df <-
  haplotype_marker_df %>%
  distinct(sample, unitig)

unitig_lengths_df <-
  haplotype_marker_df %>% 
  distinct(sample, unitig, length) %>% 
  rename(qlen = length)

chrom_lengths_df <-
  ref_aln_df %>% 
  distinct(sample, tname, tlen)

# Take the first alignment listed each time. I think there is an order to them
# as output by minimap so this seems okay? Or should this be done by, eg total
# alignment length? Weighted by quality? [This is also slow.]
unitig_to_tname_df <-
  ref_aln_df %>%
  group_by(sample, unitig) %>%
  slice_head(n=1) %>% 
  ungroup() %>% 
  distinct(sample, unitig, tname)

cluster_df <-
  haplotype_marker_df %>% 
  distinct(sample, unitig, cluster)
# Palette -----------------------------------------------------------------

# Stupid fucking conda doesn't have ggokabeito availible on the cluster.

okabeito_palette <-
  c("#E69F00",
    "#56B4E9",
    "#009E73",
    "#F0E442",
    "#0072B2",
    "#D55E00",
    "#CC79A7",
    "#999999",
    "#000000")

# Qlen Histograms ---------------------------------------------------------


## Plain -------------------------------------------------------------------


plots <- list()

p <-
  unitig_lengths_df %>% 
  ggplot() +
  geom_vline(xintercept=log10(slt)) +
  geom_histogram(aes(log10(qlen))) +
  facet_wrap(~sample, scale='free_y', drop=FALSE) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  scale_x_log10() +
  xlab('Log10 Unitig Length') +
  ylab('Count') +
  theme_bw() +
  theme(axis.title.y = element_blank()) +
  ggtitle('Qlen Histogram')

plots[['qlen_histogram']] <- p

## Rukki-Colored -----------------------------------------------------------


p <-
  rukki_paths_df %>%
  filter(!is_gap) %>%
  full_join(unitig_lengths_df, by=c('sample', 'unitig')) %>% 
  ggplot() +
  geom_vline(aes(xintercept = slt)) +
  geom_histogram(aes(qlen, fill=assignment)) +
  facet_wrap(~sample, scales = 'free_y', drop=FALSE) +
  # scale_fill_okabe_ito(name=NULL, na.value='grey50') +
  scale_fill_manual(values=okabeito_palette, name=NULL, na.value='grey50') +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  xlab('Unitig Length') +
  ylab('Count') +
  ggtitle('Rukki-Colored Unitig Length Histogram') +
  theme_bw()

plots[['rukki_colored_qlen_histogram']] <- p


## Weighted & Rukki-Colored -----------------------------------------------------------


p <-
  rukki_paths_df %>%
  filter(!is_gap) %>%
  full_join(unitig_lengths_df, by=c('sample', 'unitig')) %>% 
  ggplot() +
  geom_vline(aes(xintercept = slt)) +
  geom_histogram(aes(qlen, fill=assignment, weight=qlen/1e6)) +
  facet_wrap(~sample, scales = 'free_y', drop=FALSE) +
  # scale_fill_okabe_ito(name=NULL, na.value='grey50') +
  scale_fill_manual(values=okabeito_palette, name=NULL, na.value='grey50') +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  scale_x_log10() +
  xlab('Unitig Length') +
  ylab('Mbp') +
  ggtitle('Rukki-Colored Unitig Length Weighted Histogram') +
  theme_bw()

plots[['rukki_colored_qlen_histogram_weighted']] <- p

# Assembly Alignment and Assignment Fractions ------------------------------

plot_data_alignment <-
  ref_aln_df %>% 
  count(sample, tname, tlen, wt=alen) %>% 
  inner_join(included_chrs_df, by='tname') %>% 
  mutate(frac =n/tlen) %>% 
  select(sample, tname, ploidy, frac)

# Mostly for Debugging ~ should not differ too greatly from assembly alignment
# fractions.
plot_data_assignment <-
  unitig_to_tname_df %>% 
  full_join(unitig_lengths_df, by=c('sample', 'unitig')) %>% 
  full_join(chrom_lengths_df, by=c('sample', 'tname')) %>% 
  count(sample, tname, tlen, wt=qlen) %>% 
  inner_join(included_chrs_df, by='tname') %>% 
  mutate(frac_assign =n/tlen) %>% 
  select(sample, tname, ploidy, frac_assign)

plot_data <-
  full_join(plot_data_alignment, plot_data_assignment)

p <-
  plot_data %>% 
  ggplot() +
  geom_linerange(aes(x=tname, ymin=0, ymax=frac, color=ploidy)) +
  geom_point(aes(x = tname, y=frac, fill=ploidy, color=ploidy), shape=21) +
  facet_wrap( ~ sample, scales = 'free_y', drop=FALSE)  +
  xlab('Chromosome') +
  # scale_fill_okabe_ito(name='Ploidy') +
  scale_fill_manual(values=okabeito_palette, name='Ploidy') +
  # scale_color_okabe_ito(name='Ploidy') +
  scale_color_manual(values=okabeito_palette, name='Ploidy') +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) 

# TODO refactor so that the horizontal lines are behind the lolipops
## Alignment ---------------------------------------------------------------


plots[['assembly_alignment_fractions']] <- 
  p +
  geom_hline(yintercept=1) +
  geom_hline(yintercept=2) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  ylab('Percent of Chromosome') +
  ggtitle('Assembly Alignment Fractions')


## Assignment --------------------------------------------------------------

plots[['assembly_assignment_fractions']] <- 
  (p %+% mutate(plot_data, frac = frac_assign)) +
  geom_hline(yintercept=1) +
  geom_hline(yintercept=2) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  ylab('Percent of Chromosome') +
  ggtitle('Assembly Assignment Fractions')


## Alignment-Assignment ----------------------------------------------------


plots[['assembly_alignment-assignment_fraction_differences']] <- 
  (p %+% mutate(plot_data, frac = frac-frac_assign)) +
  ylab('Percent of Chromosome (Alignment-Assignment)') +
  ggtitle('Assembly Alignment-Assignment Fraction Differences')

# Largest Cluster Size vs tlen by Chrom ----------------------------------

plot_data <-
  haplotype_marker_df %>% 
  full_join(unitig_to_tname_df, by=c('sample', 'unitig')) %>%  
  count(sample, cluster, tname, color, wt=length, name='qlen') 

plot_data <-
  plot_data %>%
  group_by(sample, tname) %>%
  summarise(largest_cluster = max(qlen[color=='Attempted'], 0), .groups='drop')

plot_data <-
  plot_data %>% 
  inner_join(included_chrs_df, by='tname') %>% 
  left_join(chrom_lengths_df) %>% 
  mutate(frac=largest_cluster/tlen)

p <-
  plot_data %>%
  ggplot() +
  geom_hline(yintercept=1) +
  geom_hline(yintercept=2) +
  geom_linerange(aes(x=tname, ymin=0, ymax=frac, color=ploidy)) +
  geom_point(aes(x = tname, y=frac, fill=ploidy, color=ploidy), shape=21) +
  facet_wrap(~sample, drop = FALSE) +
  # scale_fill_okabe_ito(name='Ploidy') +
  scale_fill_manual(values=okabeito_palette, name='Ploidy') +
  # scale_color_okabe_ito(name='Ploidy') +
  scale_color_manual(values=okabeito_palette, name='Ploidy') +
  ylab('Mbp') +
  ggtitle('Largest Cluster as Percent of Length of Ref Chromosome') +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
  ) 

plots[['largest_cluster_div_tlen']] <- p


# PCA Variances -----------------------------------------------------------

plot_data <-
  haplotype_marker_df %>% 
  group_by(sample, cluster, PC1_pvar, PC2_pvar) %>% 
  summarise(length = sum(length))

p <-
  plot_data %>% 
  ggplot() +
  geom_point(aes(x=PC1_pvar, y=PC2_pvar, size=length/1e6), alpha=0.25, shape=21, fill='darkgrey') +
  facet_wrap(~sample) +
  scale_size_area(name = 'Total Cluster Length (Mbp)') +
  xlab('PC1 %Var') +
  ylab('PC2 %Var') +
  ggtitle('Principal Components Variation Explained') +
  theme_linedraw()

plots[['pca_variances']] <- p


p <-
  plot_data %>% 
  ggplot() +
  geom_text(aes(x=PC1_pvar, y=PC2_pvar, label=cluster), alpha=0.33) +
  facet_wrap(~sample) +
  xlab('PC1 %Var') +
  ylab('PC2 %Var') +
  ggtitle('Principal Components Variation Explained') +
  theme_linedraw()

plots[['pca_variances_label']] <- p


# 0 Marker Long Unitigs --------------------------------------------------------


plot_data <-
  haplotype_marker_df %>% 
  filter(length >= segment_length_threshold) %>% 
  filter(hap_1_counts == 0 & hap_2_counts == 0) %>% 
  mutate(clustered = !is.na(cluster)) %>% 
  arrange(desc(length)) %>% 
  mutate(unitig = factor(unitig, levels = unique(unitig)))

p <-
  ggplot(plot_data) +
  geom_histogram(aes(length, fill = clustered), color = 'black') +
  scale_fill_manual(name = 'Has Cluster', values=okabeito_palette) +
  facet_wrap( ~ sample, axes = 'all_x') +
  scale_y_continuous(expand = expansion(c(0, 0.1))) +
  scale_x_log10(limits = c(slt, NA)) +
  xlab('Log10 Length') +
  ylab('Count') +
  ggtitle(paste0('0 Marker Unitigs >= ', slt/1e6, 'Mbp')) +
  theme_classic() +
  theme(panel.border = element_rect(fill = NA, color='black', linewidth = 1))

plots[['0_marker_unitigs_hist']] <- p


p <-
  ggplot(plot_data) +
  geom_point(aes(y = unitig, x = length, fill= clustered), shape=21) +
  scale_fill_manual(name = 'Has Cluster', values=okabeito_palette) +
  facet_wrap( ~ sample, scales = 'free_y', axes = 'all_x') +
  scale_x_log10(limits = c(slt, NA)) +
  xlab('Log10 Length') +
  ylab('Unitig') +
  ggtitle(paste0('0 Marker Unitigs >= ', slt/1e6, 'Mbp')) +
  theme_linedraw() +
  theme(panel.border = element_rect(fill = NA, color='black', linewidth = 1))

plots[['0_marker_unitigs_point']] <- p


# Cluster Error By Chromosome ---------------------------------------------
merge_hap_tnames <- function(x, hap_chrs) {
  stopifnot(is.factor(x))
  x <- forcats::fct_relabel(x, function(y) if_else(y %in% hap_chrs, 'chrHAP', y))

  return(x)
}

to_merge <- c(included_haploid_chrs, included_hemiploid_chrs)

modified_inlcuded_chrs <-
  included_chrs_df %>% 
  mutate(tname = merge_hap_tnames(tname, to_merge)) %>% 
  mutate(
    ploidy = ifelse(tname == 'chrHAP', 'Haploid', ploidy)
  ) %>% 
  group_by(tname) %>% 
  # TODO when there are no references, this function returns 0 rows, which is deprecated for summarise (needs exactly 1 row)
  summarise(ploidy = unique(ploidy), is_hard_to_assemble = any(is_hard_to_assemble)) %>% 
  ungroup()

modified_tlengths <-
  chrom_lengths_df %>% 
  mutate(tname = merge_hap_tnames(tname, to_merge)) %>% 
  group_by(sample, tname) %>% 
  summarise(tlen = sum(tlen)) %>% 
  ungroup()

modified_u2t <- 
  unitig_to_tname_df %>% 
  mutate(tname = merge_hap_tnames(tname, to_merge))

# Major chromosome for each cluster:
cluster_chrom_counts <-
  haplotype_marker_df %>% 
  left_join(modified_u2t, by = c('sample', 'unitig')) %>% 
  left_join(modified_tlengths,  by = c('sample', 'tname')) %>% 
  count(sample, cluster, tname, wt=length)

major_alignments <-
  cluster_chrom_counts %>% 
  filter(!is.na(cluster)) %>% 
  group_by(sample, cluster) %>% 
  slice_max(n) %>% 
  ungroup() %>% 
  select(-n)

plot_data <-
  haplotype_marker_df %>% 
  left_join(modified_u2t, by = c('sample', 'unitig')) %>% 
  anti_join(major_alignments, by = c('sample', 'cluster', 'tname')) 

plot_data <-
  plot_data %>%
  count(sample, tname, color, wt=length) %>% # to split bar colors, group by cluster
  left_join(modified_tlengths, by = c('sample', 'tname')) %>% 
  mutate(perc = n / tlen) %>%
  filter(!is.na(perc))

plot_data <-
  plot_data %>% 
  inner_join(modified_inlcuded_chrs, by='tname')

label_data <-
  plot_data%>%
  group_by(sample, tname) %>%
  summarise(perc=sum(perc)) %>%
  group_by(sample) %>%
  slice_max(order_by = perc) %>%
  ungroup()

p <-
  ggplot(mapping = aes(x = tname, y = perc)) +
  geom_vline(
    xintercept = hard_to_assemble_chrs,
    alpha = 0.33,
    linetype = 'dashed'
  ) +
  geom_col(aes(fill = color), color='black', linewidth=0.2, data = plot_data) +
  geom_text(aes(label = paste0(round(100 * perc, 1), '%')), size=2, nudge_y = 0.05, data =
              label_data) +
  facet_wrap( ~ sample, scales = 'free_y', drop=FALSE) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  ylab('Percent of Ref Chromosome') +
  xlab('Chromosome') +
  ggtitle('Cluster Error by Percent of Length of Ref Chromosome') +
  # scale_fill_okabe_ito(name=NULL) +
  scale_fill_manual(values=okabeito_palette, name=NULL) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )


plots[['clust_chrom_error']] <- p



# Rukki Assignment Sizes --------------------------------------------------


## Aggregate ---------------------------------------------------------------


#TODO lots of lines visible in plots in pdf ~ count all in each assignment before
#plotting to remove lines?
plot_data <-
  rukki_paths_df %>%
  filter(!is_gap) %>%
  full_join(unitig_lengths_df, by=c('sample', 'unitig')) #%>% 
  # count(sample, assignment, wt=qlen, name='qlen') 

p <-
  plot_data %>% 
  ggplot() +
  geom_col(aes(x=assignment, y=qlen/1e6, fill=assignment)) +
  facet_wrap(~ sample, drop=FALSE) + 
  ylab('Size (Mbp)') +
  ggtitle('Rukki Assignment Sizes') +
  xlab(NULL) +
  # scale_fill_okabe_ito(name=NULL, na.value='grey50') +
  scale_fill_manual(values=okabeito_palette, name=NULL, na.value='grey50') +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )

plots[['rukki_assignment_sizes']] <- p

## NA Only -----------------------------------------------------------------
# Rukki Assignment Sizes NA only
p <-
  plot_data %>% 
  filter(is.na(assignment)) %>% 
  ggplot() +
  geom_col(aes(x=sample, y=qlen/1e6, fill=assignment), color='black', linewidth = 0.05) +
  ylab('Size (Mbp)') +
  ggtitle('Rukki NA Assignment Sizes') +
  xlab(NULL) +
  # scale_fill_okabe_ito(name=NULL, na.value='grey50') +
  scale_fill_manual(values=okabeito_palette, name=NULL, na.value='grey50') +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )

plots[['rukki_assignment_sizes_na']] <- p


## By Cluster --------------------------------------------------------------------
# Rukki Assignment Sizes by connected component?
# TODO
p <-
  plot_data %>% 
  left_join(cluster_df, by = c('sample', 'unitig')) %>% 
  ggplot()  +
  geom_col(aes(x=cluster, y=qlen/1e6, fill=assignment))+
  facet_wrap(~sample, scales = 'free_x') +
  ylab('Size (Mbp)')  +
  xlab('Cluster') +
  ggtitle('Rukki Assignment Sizes by Cluster') +
  scale_fill_manual(values=okabeito_palette, name=NULL, na.value='grey50') +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
  )

plots[['rukki_assignment_sizes_clust']] <- p 


p <-
  plot_data %>% 
  left_join(cluster_df, by = c('sample', 'unitig')) %>% 
  group_by(sample, cluster) %>% 
  mutate(qlen = qlen/sum(qlen, na.rm = TRUE)) %>% 
  ggplot()  +
  geom_col(aes(x=cluster, y=qlen, fill=assignment))+
  facet_wrap(~sample, scales = 'free_x') +
  ylab('Fraction of Cluster')  +
  xlab('Cluster') +
  ggtitle('Rukki Assignment Fractions by Cluster') +
  scale_fill_manual(values=okabeito_palette, name=NULL, na.value='grey50') +
  scale_y_continuous(expand = expansion(mult=c(0,0)), breaks = seq(0, 1, 0.1)) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
  )

plots[['rukki_assignment_fractions_clust']] <- p
 
## By Ref Chromosome ---------------------------------------------------------

# Rukki Assignment Sizes by Ref Chrom
plot_data <-
  rukki_paths_df %>%
  filter(!is_gap) %>%
  full_join(unitig_lengths_df, by=c('sample', 'unitig')) %>% 
  full_join(unitig_to_tname_df)  %>% 
  count(sample, tname, assignment, wt=qlen, name='qlen') 

plot_data <-
  plot_data  %>% 
  left_join(chrom_lengths_df) 

plot_data <-
  plot_data %>% 
  tidyr::complete(sample, tname, assignment, fill=list(qlen=0, tlen=1)) %>% 
  inner_join(included_chrs_df)
p <-
  plot_data %>% 
  ggplot() +
  # geom_linerange(aes(x=tname, ymin=0, ymax=qlen/tlen, color=assignment),  position = position_dodge(width=0.75, preserve = 'total')) +
  # geom_point(aes(x=tname, y=qlen/tlen, fill=assignment, color=assignment), shape=21,  position = position_dodge(width=0.75, preserve = 'total')) +
  geom_col(aes(x=tname, y=qlen/tlen, fill=assignment),color='black', linewidth=0.1) +
  facet_wrap(~ sample, drop=FALSE) + 
  ylab('Fraction of Ref Chromosome') +
  ggtitle('Rukki Assignment Sizes by Fraction of Ref Chromosome') +
  xlab(NULL) +
  # scale_fill_okabe_ito(name=NULL, na.value='grey50', drop=FALSE) +
  scale_fill_manual(values=okabeito_palette, na.value='grey50', drop=FALSE) +
  # scale_color_okabe_ito(name=NULL, na.value='grey50', drop=FALSE) +
  scale_color_manual(values=okabeito_palette, na.value='grey50', drop=FALSE) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
  )

plots[['rukki_assignment_sizes_by_ref_chrom_frac']] <- p

# TODO A plot with the homozygous nodes added to hap 1 and 2?


# Haploid Phasing -----------------------------------------------------

plot_data <-
  haplotype_marker_df %>% 
  left_join(unitig_to_tname_df, by = c('sample', 'unitig')) %>%
  filter(tname %in%  c(included_haploid_chrs, included_hemiploid_chrs))

plot_data <-
  plot_data %>%
  mutate(hap_1_counts = hap_1_counts + 1,
         hap_2_counts = hap_2_counts + 1) %>%
  mutate(ssf = (hap_1_counts - hap_2_counts)/(hap_1_counts+hap_2_counts))

plot_data <-
  plot_data %>% 
  mutate(`Length (Mbp)` = length/1e6)

p <- 
  plot_data %>%
  ggplot() +
  geom_point(aes(
    x = ssf,
    y = tname,
    size = `Length (Mbp)`,
    fill = tname,
    group = sample
  ),
  shape = 21) +
  facet_wrap( ~ sample, drop=FALSE) +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  scale_size_area()  +
  # scale_fill_okabe_ito(name=NULL) +
  scale_fill_manual(values=okabeito_palette, name=NULL) +
  theme_bw() +
  ylab('Chromosome') +
  xlab('Strand State Frequency') +
  ggtitle('Haploid Chromosome SSFs') +
  # theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
  )

plots[['hap_marker_ratios']] <- p


# Contiguity --------------------------------------------------------------


unitig_lengths <-
  unitig_lengths_df %>% 
  rename(length = qlen) %>% 
  mutate(metric = 'Unitig')

cluster_lengths <-
  haplotype_marker_df %>% 
  count(sample, cluster, wt=length, name='length') %>% 
  mutate(metric = 'Cluster')

.path_lengths <-
  rukki_paths_df %>% 
  left_join(unitig_lengths_df, by=c('sample', 'unitig')) %>% 
  mutate(node_len = coalesce(gap_size, qlen)) 

scaftig_lengths <-
  .path_lengths %>%
  count(sample, name, run, wt=node_len * (!is_gap), name='length') %>% 
  mutate(metric = 'Scaftig')

scaffold_lengths <-
  .path_lengths %>%
  count(sample, name, wt=node_len, name='length') %>% 
  mutate(metric = 'Scaffold')

ref_lengths <-
  ref_aln_df %>% 
  distinct(sample, tname, tlen) %>% 
  rename(length = tlen) %>% 
  mutate(metric = 'Ref')

# TODO handle the the columns that vary between them ~ no harm but messy?
all_lengths <-
  bind_rows(
    unitig_lengths,
    cluster_lengths,
    scaftig_lengths,
    scaffold_lengths,
    ref_lengths
  )

## Nx Plots ----------------------------------------------------------------
# TODO
pal <-
  okabeito_palette[c(9,8,1,2,3,7)] %>% 
  set_names(c('Ref', '2 x Ref', 'Scaffold', 'Scaftig', 'Unitig', 'Cluster'))

  
calc_cum_frac <- function(x, decreasing=TRUE, na.rm=TRUE) {
  sort_order<- order(x, decreasing = decreasing)
  unsort_order <- order(sort_order)
  
  # browser()
  cum_length <- cumsum(x[sort_order])
  cum_frac <- cum_length/max(cum_length, na.rm=na.rm)
  out <- cum_frac[unsort_order]
  return(out)
}

plot_data <-
  all_lengths %>% 
  group_by(sample, metric) %>% 
  arrange(desc(length)) %>% 
  mutate(cum_frac = cumsum(length)/sum(length))

left_points <-
  plot_data %>% 
  slice_head(n=1) %>% 
  mutate(cum_frac=0)

plot_data <-
  plot_data %>% 
  bind_rows(left_points) %>% 
  ungroup()


### Rukki Focused -----------------------------------------------------------


plot_data_rukki <-
  plot_data %>% 
  filter(metric != 'Cluster')

p <-
  plot_data_rukki %>% 
  ggplot() +
  geom_step(aes(
    x = cum_frac,
    y = length / 1e6,
    color = metric
  ),
  direction = 'vh') +
  facet_wrap(~sample, scales='free_y', drop=FALSE) +
  scale_color_manual(values=pal, name=NULL) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  xlab('Cumulative Fraction') +
  ylab('Mbp') +
  ggtitle('Rukki Nx Plots') +
  theme_bw() 

plots[['rukki_nx_curves']] <- p



### Cluster Focused ---------------------------------------------------------


plot_data_cluster <-
  plot_data %>% 
  filter(metric %in% c('Cluster', 'Unitig', 'Ref')) %>% 
  mutate(length = ifelse(metric == 'Ref', 2*length, length)) %>% 
  mutate(metric = ifelse(metric == 'Ref', '2 x Ref', metric))

p <-
  plot_data_cluster %>% 
  ggplot() +
  geom_step(aes(
    x = cum_frac,
    y = length / 1e6,
    color = metric
  ),
  direction = 'vh') +
  facet_wrap(~sample, scales='free_y', drop=FALSE) +
  scale_color_manual(values=pal, name=NULL) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  xlab('Cumulative Fraction') +
  ylab('Mbp') +
  ggtitle('Cluster Nx Plots') +
  theme_bw() 

plots[['cluster_nx_curves']] <- p



### Cluster vs Scaffold Focused ---------------------------------------------------------


plot_data_cluster_scaffold <-
  plot_data %>% 
  filter(metric %in% c('Scaffold', 'Ref')) %>% 
  bind_rows(plot_data_cluster)


p <-
  plot_data_cluster_scaffold %>% 
  ggplot() +
  geom_step(aes(
    x = cum_frac,
    y = length / 1e6,
    color = metric
  ),
  direction = 'vh') +
  facet_wrap(~sample, scales='free_y', drop=FALSE) +
  scale_color_manual(values=pal, name=NULL) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  xlab('Cumulative Fraction') +
  ylab('Mbp') +
  ggtitle('Cluster and Scaffold Nx Plots') +
  theme_bw() 

plots[['cluster_scaffold_nx_curves']] <- p

## auN ---------------------------------------------------------------------

calc_aun <- function(x, na.rm=FALSE) {
  sum(x * x/sum(x, na.rm=na.rm), na.rm=na.rm)
}

auns <-
  all_lengths %>% 
  group_by(sample, metric) %>% 
  summarise(aun = calc_aun(length)) %>% 
  ungroup()

ref_auns <-
  auns %>% 
  filter(metric == 'Ref')

auns <-
  auns %>% 
  anti_join(ref_auns)

common <-
  list(
    facet_wrap(~sample, scales = 'free_y', drop=FALSE),
    scale_y_continuous(expand = expansion(mult=c(0,0.1))),
    ylab('auN (Mbp)'),
    xlab(NULL),
    theme_bw(),
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
    )
  )
# 
# p <-
#   auns %>%
#   filter(metric %in% c('Unitig', 'Scaftig', 'Scaffold')) %>% 
#   ggplot() +
#   geom_col(aes(x=metric, y=aun/1e6)) +
#   geom_hline(aes(yintercept=aun/1e6), data=ref_aun, alpha=0.75, linetype='62') +
#   ggtitle('Rukki auNs') +
#   common
# 
# plots[['rukki_auns']] <- p
# 
# 
# p <-
#   auns %>%
#   filter(metric %in% c('Unitig', 'Cluster')) %>% 
#   ggplot() +
#   geom_col(aes(x=metric, y=aun/1e6)) +
#   geom_hline(aes(yintercept=2*aun/1e6), data=ref_aun, alpha=0.75, linetype='62') +
#   ggtitle('Cluster auNs') +
#   common
# 
# plots[['cluster_auns']] <- p

# TODO add lines related to connected component sizes?

p <-
  auns %>% 
  ggplot() +
  geom_hline(aes(yintercept=aun/1e6), data=ref_auns, alpha=0.75, linetype='62') +
  geom_col(aes(x=metric, y=aun/1e6)) +
  geom_hline(aes(yintercept=2*aun/1e6), data=ref_auns, alpha=0.75, linetype='62') +
  ggtitle('auNs') +
  common

plots[['all_auns']] <- p



# SSFs -------------------------------------------------------------------

plot_data <-
  haplotype_marker_df %>% 
  left_join(rukki_paths_df, by=c('sample','unitig'))

plot_data <-
  plot_data %>%
  mutate(hap_1_counts = hap_1_counts + 1,hap_2_counts = hap_2_counts + 1) %>% 
  mutate(ssf = (hap_1_counts - hap_2_counts)/(hap_1_counts+hap_2_counts))

## SSF Histograms --------------------------------------------------------



p <-
  plot_data %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_histogram(aes(ssf, weight=length, fill=assignment)) +
  facet_wrap(~sample, drop=FALSE) +
  # scale_fill_okabe_ito(name='Rukki Assignment', na.value='grey50') +
  scale_fill_manual(values=okabeito_palette, name='Rukki Assignment', na.value='grey50') +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  xlab('Strand State Frequency') +
  ylab(NULL) +
  ggtitle('Strand State Frequency Histogram') +
  theme_bw()

plots[['ssf_histogram']] <- p


## SSF Scatter -----------------------------------------------------------

p <-
  plot_data %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_point(aes(x=ssf, y=assignment, fill=assignment, size=length/1e6), shape=21) +
  facet_wrap(~sample, drop=FALSE) +
  # scale_fill_okabe_ito(name='Rukki Assignment', na.value='grey50') +
  scale_fill_manual(values=okabeito_palette, name='Rukki Assignment', na.value='grey50') +
  scale_size_area(name='Mbp') +
  xlab('Strand State Frequency') +
  ylab(NULL) +
  ggtitle('Strand State Frequency Scatter') +
  theme_bw()

plots[['ssf_scatter']] <- p

## Unexpected SSF Scatter ------------------------------------------------

p <-
  plot_data %>% 
  filter(!(assignment == 'HAPLOTYPE1' & ssf > 0) | is.na(assignment)) %>% 
  filter(!(assignment == 'HAPLOTYPE2' & ssf < 0) | is.na(assignment)) %>% 
  ggplot() +
  geom_vline(xintercept = 0, linetype = 'dotted') +
  geom_point(aes(x=ssf, y=assignment, fill=assignment, size=length/1e6), shape=21) +
  facet_wrap(~sample, drop=FALSE) +
  # scale_fill_okabe_ito(name='Rukki Assignment', na.value='grey50') +
  scale_fill_manual(values=okabeito_palette, name='Rukki Assignment', na.value='grey50') +
  scale_size_area(name='Mbp') +
  xlab('Strand State Frequency') +
  ylab(NULL) +
  ggtitle('Unexpected Strand State Frequency Scatter') +
  theme_bw()

plots[['ssf_scatter_unexpected']] <- p


# Marker Blob Plot --------------------------------------------------------


plot_data <-
  haplotype_marker_df %>% 
  left_join(rukki_paths_df, by=c('sample','unitig'))  # %>% 
  # mutate(assignment = forcats::fct_explicit_na(assignment, na_level='NA'))

symmetry_data <-
  plot_data %>% 
  mutate(tmp = hap_1_counts) %>% 
  mutate(hap_1_counts = hap_2_counts) %>% 
  mutate(hap_2_counts = tmp) 

# plot_data <-
#   bind_rows(plot_data, invisible_data)
## Blob --------------------------------------------------------------------


p <-
  plot_data %>%
  ggplot() +
  geom_abline(intercept = 0, slope=1, alpha=0.75) +
  geom_point(aes(x=hap_1_counts, y=hap_2_counts, fill=assignment, size=length/1e6), shape=21) +
  geom_point(aes(x=hap_1_counts, y=hap_2_counts), color=NA, data=symmetry_data) +
  # facet_wrap(~sample+assignment, drop=FALSE, scales='free') +
  facet_wrap(~sample, drop=FALSE, scales='free') +
  # scale_fill_okabe_ito(name='Rukki Assignment', order=c(1,2,3,8)) +
  scale_fill_manual(values=okabeito_palette[c(1,2,3,8)], name='Rukki Assignment', na.value = 'grey50') +
  scale_size_area(name='Mbp') +
  xlab('Haplotype 1 Marker Counts') +
  ylab('Haplotype 2 Marker Counts') +
  ggtitle('Haplotype Marker Blob Plots') +
  theme_bw()

plots[['hap_marker_blob']] <- p

## Log Blob ---------------------------------------------------------------

# NOTE: The straight lines that can appsear in the log-log blob plot are
# consequences of the +1 pesudocount added for low marker-count unitigs.
p <-
  plot_data %>%
  ggplot() +
  geom_abline(intercept = seq(-10,10,1), slope=1, alpha=0.25 ) +
  geom_point(aes(x=log10(hap_1_counts+1), y=log10(hap_2_counts+1), fill=assignment, size=length/1e6), shape=21) +
  geom_point(aes(x=log10(hap_1_counts+1), y=log10(hap_2_counts+1)), color=NA, data=symmetry_data) +
  facet_wrap(~sample, drop=FALSE, scales='free') +
  # scale_fill_okabe_ito(name='Rukki Assignment', order=c(1,2,3,8)) +
  scale_fill_manual(values=okabeito_palette[c(1,2,3,8)], name='Rukki Assignment', na.value = 'grey50') +
  # scale_x_continuous(trans='exp') +
  # scale_y_continuous(trans='exp') +
  scale_size_area(name='Mbp') +
  xlab('Log10 Haplotype 1 Marker Counts') +
  ylab('Log10 Haplotype 2 Marker Counts') +
  ggtitle('Log-Log Haplotype Marker Blob Plots') +
  theme_bw()

plots[['hap_marker_blob_log']] <- p

# Rukki-Marker Agreement --------------------------------------------------

# Rukki HOM calls will almost certainly disagree with the marker count calls,
# which require exact equality. Therefore, disagreement is only really measured
# for the non-HOM nodes. The HOM nodes as called by RUKKI, the issue nodes as
# called by RUKKI, are labeled, but not disagreement-quantified.

marker_calls <-
  haplotype_marker_df %>%
  mutate(
    marker_assignment =
      case_when(hap_1_counts > hap_2_counts ~ 'HAPLOTYPE1',
                hap_2_counts > hap_1_counts ~ 'HAPLOTYPE2',
                hap_1_counts == hap_2_counts ~'HOM', 
                TRUE ~ NA)
  ) %>% 
  select(sample, unitig, marker_assignment)



# rukki_paths_df <-
#   rukki_paths_df %>% 
#   mutate(
#     assignment = ifelse(unitig %in% hom_nodes$unitig, 'HOM', assignment)
#   ) %>% 
#   mutate(
#     assignment = ifelse(unitig %in% issue_nodes$unitig, 'ISSUE', assignment)
#   )


plot_data <-
  rukki_paths_df %>%
  filter(!is_gap) %>%
  full_join(unitig_lengths_df, by = c('sample', 'unitig')) %>% 
  left_join(marker_calls, by = c('sample', 'unitig'))  %>%
  mutate(
    disagreement = ifelse(assignment != marker_assignment, qlen, 0)
  ) 


plot_data <-
  plot_data %>%
  left_join(hom_nodes) %>%
  left_join(issue_nodes) %>%
  mutate(fill = case_when(
    is_hom ~ 'RUKKI_HOM', 
    is_issue ~ 'RUKKI_NA',
    disagreement > 0 ~ 'DISAGREEMENT',
    disagreement == 0 ~ 'AGREEMENT',
    TRUE ~ 'UNCLASSIFIED'
  ))


## Aggregate ---------------------------------------------------------------

p <-
  plot_data %>%  
  count(sample, fill, wt=qlen) %>% 
  ggplot() +
  geom_col(aes(x=sample, y=n/1e6, fill=fill), color = 'black', linewidth=0.2) +
  ylab('Mbp') +
  xlab(NULL) +
  ggtitle('Rukki - Marker Agreement') +
  # scale_fill_okabe_ito(order = c(3,1,2,4,5,6,7,8), name=NULL) +
  scale_fill_manual(values=okabeito_palette[c(3,1,2,4,5,6,7,8)], name=NULL) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plots[['rukki_marker_agreement']] <- p


## By CC -------------------------------------------------------------------

# TODO

## By Reference Chromosome -------------------------------------------------

plot_data <-
  plot_data %>%
  left_join(unitig_to_tname_df, by = c('sample', 'unitig')) %>% 
  semi_join(included_chrs_df, by='tname') %>% 
  count(sample, tname, fill, wt=qlen) 

p <- 
  plot_data %>%
  ggplot() +
  geom_col(aes(x=tname, y=n/1e6, fill=fill), color = 'black', linewidth=0.2) +
  facet_wrap(~sample, scales = 'free_x', drop=FALSE) +
  ylab('Mbp') +
  xlab(NULL) +
  ggtitle('Rukki - Marker Agreement by Ref Chromosome') +
  # scale_fill_okabe_ito(order = c(3,1,2,4,5,6,7,8), name=NULL) +
  scale_fill_manual(values=okabeito_palette[c(3,1,2,4,5,6,7,8)], name=NULL) +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


plots[['rukki_marker_agreement_by_chrom']] <- p



# Misjoins ----------------------------------------------------------------

tmp_rukki_data <-
  rukki_paths_df %>% 
  left_join(unitig_to_tname_df) %>% 
  inner_join(included_chrs_df) 

paths_with_misjoins <-
  tmp_rukki_data %>% 
  distinct(sample, name, tname) %>% 
  count(sample, name) %>% 
  filter(n > 1) %>% 
  # mutate(n_misjoins = n-1) %>% 
  select(-n)

plot_data <-
  tmp_rukki_data %>% 
  inner_join(paths_with_misjoins) %>% 
  group_by(sample, name) %>% 
  summarise(all_hard_to_assemble = all(is_hard_to_assemble), .groups='drop') %>% 
  count(sample, all_hard_to_assemble)

p <-
  ggplot(plot_data) +
  geom_col(aes(x=sample, y=n, fill=all_hard_to_assemble)) +
  # scale_fill_okabe_ito(name='Aligns Only to Hard-to-Assemble Chromosomes') +
  scale_fill_manual(values=okabeito_palette, name='Aligns Only to\nHard-to-Assemble Chromosomes') +
  scale_y_continuous(expand = expansion(mult=c(0,0.1))) +
  xlab('Sample') +
  ylab('N') +
  ggtitle('Paths Aligning to Multiple Chromosomes') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

plots[['rukki_path_multichrom_alignments']] <- p
# Gaps --------------------------------------------------------------------
# TODO, or already handled with scaftig aun and Nx plots?

# End-to-end alignments ---------------------------------------------------
# TODO

# Cluster Size Barplots ---------------------------------------------------
# TODO

# ECDF --------------------------------------------------------------------
# TODO

# Contiguity - Lx/LGx ------------------------------------------------------
# TODO

# Export ------------------------------------------------------------------

# TODO for reference-free version, remove plots that have no content
pdf(output, width = out_width, height = out_height)
walk(plots, print)
dev.off()
