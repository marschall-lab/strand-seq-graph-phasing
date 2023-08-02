setwd("~/Documents/wd")
# Functions ---------------------------------------------------------------

calculate_n50 <- function(lengths) {
  sorted_lengths <- sort(lengths, decreasing = TRUE)
  total_length <- sum(sorted_lengths)
  half_length <- total_length / 2
  cumulative_length <- cumsum(sorted_lengths)
  n50_index <- min(which(cumulative_length >= half_length))
  n50 <- sorted_lengths[n50_index]
  return(n50)
}
# Library -----------------------------------------------------------------


library(dplyr)
library(purrr)
library(pafr)
library(ggplot2)
library(ggokabeito)
library(progressr)
# Import ------------------------------------------------------------------


## Reference Alignments ----------------------------------------------------

reference_alignments <- list.files('reference_alignments/T2Tv11_hg002Yv2_chm13/', full.names = TRUE)



reference_alignments <-
  list.files('reference_alignments/T2Tv11_hg002Yv2_chm13/', full.names = TRUE) %>% 
  # reference_alignments %>% 
  set_names(hutils::trim_common_affixes(.))


reference_alignments <-
  reference_alignments %>%
  map(pafr::read_paf,
      tibble = TRUE,
      include_tags = FALSE) %>% 
  bind_rows(.id='sample')


# take the first alignment listed each time. I think there is an order to them
# as output by minimap so this seems okay?
reference_alignments <-
  reference_alignments %>% 
  group_by(sample, qname) %>% 
  slice_head(n=1) %>% 
  dplyr::rename(unitig = qname)

reference_alignments <-
  reference_alignments %>% 
  mutate(tname = factor(tname, levels = c(paste0('chr', 1:22), 'chrX', 'chrY', 'chrM')))
## Haplotype Marker Counts -------------------------------------------------


haplotype_marker_counts <- list.files('haplotype_marker_counts/', full.names = TRUE, pattern = 'fudged_haplotype_marker_counts.csv') %>% 
  set_names(hutils::trim_common_affixes(.))


haplotype_marker_counts <-
  haplotype_marker_counts %>% 
  map(readr::read_csv) %>% 
  bind_rows(.id='sample')

haplotype_marker_counts <-
  haplotype_marker_counts %>% 
  filter(!grepl('good_frac', sample) | grepl('frac-100', sample)) 



## Connected Components ----------------------------------------------------

ccs <- list.files('gfa/ccs', full.names = TRUE) %>% 
  set_names(hutils::trim_common_affixes(.))

ccs <-
  ccs %>% 
  map(readr::read_tsv) %>% 
  bind_rows(.id='sample') 

ccs <- 
  ccs %>% mutate(component = ifelse(member_largest_component, 0, component))


## Rukki Paths -------------------------------------------------------------------


rukki_paths <-
  list.files('rukki',recursive = TRUE,  full.names = TRUE, pattern = '*rukki_paths.tsv') %>% 
  set_names(function(x) hutils::trim_common_affixes(basename(x)))

rukki_paths <-
  rukki_paths %>% 
  map_dfr(readr::read_tsv, .id='sample')

rukki_paths <-
  rukki_paths %>% 
  tidyr::separate_longer_delim(path, delim = ',') %>% 
  group_by(sample, name) %>% 
  mutate(order = 1:n()) %>% 
  ungroup() %>% 
  dplyr::rename(unitig = path) %>% 
  mutate(unitig = gsub("[+-]+$", '', unitig))


rukki_paths <-
  rukki_paths %>% 
  mutate(is_gap = grepl('\\[', unitig)) %>% 
  group_by(sample, name) %>% 
  mutate(run = cumsum(is_gap)) %>% 
  ungroup()


rukki_paths <-
  rukki_paths %>%
  mutate(gap_size = ifelse(is_gap, stringr::str_extract(unitig, '[0-9]+'), NA)) %>%
  mutate(gap_size = as.integer(gap_size)) %>%
  left_join(reference_alignments, by = c('sample', 'unitig')) %>%
  mutate(node_len = coalesce(qlen, gap_size)) %>%
  select(sample, name, unitig, assignment, order, is_gap, run, node_len)



rukki_paths <-
  rukki_paths %>% 
  mutate(is_not_gap = !is_gap) 

## Rukki Annotations -------------------------------------------------------

rukki_annotations <-
  list.files('rukki',recursive = TRUE,  full.names = TRUE, pattern = '*out_final_ann.csv') %>% 
  set_names(function(x) hutils::trim_common_affixes(basename(x)))

rukki_annotations <-
  rukki_annotations %>% 
  map_dfr(readr::read_tsv, .id='sample')


rukki_annotations <-
  rukki_annotations %>% 
  select(sample, node, assignment, info) %>% 
  dplyr::rename(unitig = node, annotation = assignment)


# library_weights ---------------------------------------------------------


lib_weights <-
  list.files('library_weights/',recursive = TRUE,  full.names = TRUE) %>% 
  set_names(function(x) hutils::trim_common_affixes(basename(x)))

lib_weights <-
  map_dfr(lib_weights, readr::read_csv, .id='sample')

# Join --------------------------------------------------------------------

left_join(haplotype_marker_counts, reference_alignments, by = c('sample', 'unitig')) %>% 
  # filter(!grepl('^Length less than', exclusion)) %>% 
  split(.$sample) %>% 
  iwalk(function(x, nm) {
    x %>% 
      select(unitig, everything()) %>% 
      readr::write_csv(glue::glue('rmc/{nm}_ref_marker_counts.csv'))
  })


# Size-Ratio Plot ---------------------------------------------------------
plot_data <-
  haplotype_marker_counts %>% 
  filter(!(hap_1_counts + hap_2_counts == 0)) %>% 
  mutate(
    hap_1_counts = hap_1_counts + 1,
    hap_2_counts = hap_2_counts + 1
  )


plot_data <-
  plot_data %>% 
  mutate(rat = hap_1_counts/hap_2_counts) %>% 
  mutate(log_rat = log(hap_1_counts/hap_2_counts)) %>% 
  mutate(primary_rat = exp(abs(log_rat))) %>% 
  mutate(abs_log_rat = abs(log_rat),
         log10_length =log10(length))


mod <- lm(abs_log_rat ~ log10_length, data= plot_data)

plot_data %>% 
  filter(unitig %in% c('utig4-417', 'utig4-418', 'utig4-2736', 'utig4-2737')) %>%
  filter(sample == 'HG02769') %>%
  filter(log10(length) >= 6) %>% 
  ggplot(aes(x = log10(length), abs(log_rat))) +
  geom_point() +
  geom_smooth(method = 'lm') 


# To look at: -------------------------------------------------------------

#1.5 ~ 93%
#1 ~ 98%

to_look_at <- 
  plot_data %>% 
  filter(abs(log_rat) <= 1) %>% 
  select(sample,unitig)

samples_to_look_at <-
  to_look_at %>% 
  pull_distinct(sample)


# Raw Counts Import -------------------------------------------------------


raw_counts <-
  list.files('sseq_alignment_counts/', full.names = TRUE, pattern='fastmap_raw.csv') %>% 
  set_names(function(x) hutils::trim_common_affixes(basename(x)))

raw_counts <- raw_counts['HG02769']

raw_counts_df <-
  map(raw_counts, readr::read_csv) %>% 
  bind_rows(.id = 'sample')

raw_counts_df <- 
  raw_counts_df %>% 
  mutate(tstart = as.double(tstart)) %>% 
  mutate(
    y = ifelse(strand == '+', 1, -1),
    tstart = as.double(tstart)
  ) 


# Unitig Information ------------------------------------------------------


unitig_info <-
  haplotype_marker_counts %>% 
  distinct(sample, unitig, cluster, length, unitig_orientation)

raw_counts_df <- 
  raw_counts_df %>% 
  left_join(unitig_info, by = c('sample', 'unitig')) %>% 
  filter(!is.na(cluster))


# Lib Weights -------------------------------------------------------------


raw_counts_df <-
  raw_counts_df %>% 
  left_join(lib_weights, by=c('sample', 'lib', 'cluster')) %>% 
  filter(!is.na(ww_weight_mem))


# Plots -------------------------------------------------------------------


pd <-
  raw_counts_df %>% 
  group_by(sample, unitig) %>% 
  mutate(tstart = ifelse(unitig_orientation == 1, tstart, max(tstart) -tstart)) %>% 
  ungroup()
  
pd <- 
  pd %>% 
    mutate(color = y * wc_weight_fastmap * unitig_orientation) 

pd <-
  pd %>% 
  semi_join(to_look_at)

unitigs_to_plot <-
  pd %>% 
  distinct(unitig) %>% 
  sample_n(40)

bw <- 100e3
pd %>% 
  filter(length >= 5e5) %>% 
  # semi_join(unitigs_to_plot) %>% 
  ggplot() +
  geom_histogram(aes(tstart, weight=color, group=sign(color), fill = as.factor(sign(color))), binwidth = bw) +
  facet_wrap(~unitig, scale='free') +
  theme_grey(base_size=14)


# Scratch -----------------------------------------------------------------


samples <-
  tibble(sample = list.files('bwa_alignments/fastmap/')) %>% 
  filter(!grepl('hifiasm', sample)) %>% 
  filter(!grepl('NA24385', sample))


n_threads <- 10
if(n_threads > 1) {
  library(furrr)
  plan(multisession, workers=n_threads)
  import_mapper <- furrr::future_map  
} else {
  import_mapper <- purrr::map
}

with_progress({
  p <- progressor(steps = length(samples$sample))
  exact_match_counts <- 
    samples$sample %>% 
    set_names() %>% 
    import_mapper(function(x) {
      p()
      fastmap_alignment_files <- list.files(paste0('bwa_alignments/fastmap/', x), full.names = TRUE)
      
      lib_names <-
        map_chr(fastmap_alignment_files, function(x){
          gsub('_maximal_unique_exact_match.tsv$', '', basename(x))
        })
      
      # There is a common warning stating "incomplete final line found on". I
      # think it is connected to the source of the extra // that appear in the files
      emc_df <- map(fastmap_alignment_files, function(x) {
        cat(paste('counting bwa-fastmap alignments in', basename(x), '\n'))
        out <- extract_exact_matches(x)
        return(out)
      })
      
      
      
      emc_df <-
        emc_df %>%
        set_names(lib_names) %>%
        bind_rows(.id = 'lib') 
      
      return(emc_df)
      
    })
})




if(n_threads > 1) {
  # close workers
  plan(sequential)
}


  

