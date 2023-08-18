
# WD ----------------------------------------------------------------------

setwd("~/Documents/wd")
# Library -----------------------------------------------------------------

library(dplyr)
library(purrr)
library(furrr)
library(ggplot2)
library(ggokabeito)
library(pafr)



## Reference Alignments ----------------------------------------------------

reference_alignments <-
  list(
    HG03683stage4='reference_alignments/T2Tv11_hg002Yv2_chm13/HG03683stage4_T2Tv11_hg002Yv2_chm13_ref-aln.paf',
    HG03683='reference_alignments/T2Tv11_hg002Yv2_chm13/HG03683_T2Tv11_hg002Yv2_chm13_ref-aln.paf'
  )

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

# Rukki -------------------------------------------------------------------


rukki_paths <-
  list(
    HG03683stage4='rukki/HG03683stage4/HG03683stage4_rukki_paths.tsv',
    HG03683='rukki/HG03683/HG03683_rukki_paths.tsv'
    ) %>% 
  map(readr::read_tsv) %>% 
  bind_rows(.id='sample')



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


# Plot_data ---------------------------------------------------------------

plot_data <-
  rukki_paths %>% 
  left_join(reference_alignments, by = c('sample', 'unitig'))

plot_data <-
  plot_data %>% 
  mutate(plot_name = paste0(sample, '_', name))


plot_data %>% 
  filter(!is.na(assignment)) %>% 
  filter(!grepl('unused', name)) %>% 
  filter(qlen >= 250000) %>% 
  arrange(desc(sample), assignment, desc(tstart)) %>%
  mutate(name = factor(name, levels = unique(name))) %>%
  ggplot() +
  geom_segment(aes(x = tstart, xend=tstart + alen, y=name, yend=name, color = paste(sample, assignment)), linewidth=3) +
  facet_wrap(~tname, scales='free') +
  scale_color_okabe_ito(order=c(5,2,6,1,2,3,4)) +
  theme(axis.text.y = element_blank())






# Scratch -----------------------------------------------------------------

ref_aln <-
  pafr::read_paf(
    'correspondance_test/HG03683.ps-sseq.wd/assembly-hpc_ref-aln.paf',
    tibble = TRUE,
    include_tags = FALSE
  )


plot_data <-
  ref_aln %>% 
  rename(unitig = qname) %>% 
  mutate(sample = 'HG03683') %>% 
  filter(qlen >= 250e3) %>% 
  select(sample, everything()) %>%
  arrange(sample, unitig)

plot_data %>% 
  ggplot() +
  geom_histogram(aes(log10(alen/qlen)))

plot_data %>% 
  ggplot() +
  geom_histogram(aes(log10(alen/qlen), weight=alen)) +
  scale_x_continuous(n.breaks = 20)



betweenrs <- 
  plot_data %>% filter(between(alen/qlen, 0.2, 0.8)) %>% 
  distinct(sample, unitig)


plot_data %>% semi_join(betweenrs)
# %>% 
#   group_by(unitig) %>% 
#   filter(n() > 1) %>% 
#   ungroup()

plot_data2 <-
  ref_aln %>% 
  rename(unitig = qname) %>% 
  mutate(sample = 'HG03683') %>% 
  group_by(unitig) %>% 
  # mutate(cum_perc = cumsum(alen/qlen)) %>% 
  # filter(cum_perc <= 0.99) %>% 
  slice_head(n=3) %>%
  filter(qlen >= 250e3) %>% 
  select(sample, everything()) 

plot_data2 <-
  rukki_paths %>% 
  filter(!is_gap) %>% 
  group_by(sample, unitig) %>% 
  filter(n() == 1) %>% 
  ungroup() %>% 
  right_join(plot_data2) %>% 
  filter(!is.na(tname)) 


qlens <-
  plot_data2 %>% 
  group_by(sample, name) %>% 
  distinct(unitig, qlen) %>% 
  summarise(qlen = sum(qlen), .groups='drop') 

plot_data2 %>% 
  group_by(sample, name, tname) %>% 
  summarise(alen = sum(alen)) %>% 
  left_join(qlens) %>% 
  ggplot() +
  geom_histogram(aes(alen/qlen, weight=alen)) +
  scale_x_continuous(n.breaks = 20)

  

plot_data2 %>% 
  group_by(assignment, tname, name) %>% 
  summarise(alen = sum(alen)) %>% 
  ggplot() +
  geom_col(aes(x = name, y = alen, fill=tname), color='black') +
  scale_fill_discrete(guide='none') +
  facet_wrap(~assignment, nrow = 3, scales='free')


plot_data2 %>% 
  # filter(assignment %in% c('HAPLOTYPE1', 'HAPLOTYPE2')) %>% 
  group_by(assignment, tname, name) %>% 
  summarise(alen = sum(alen)) %>% 
  ggplot() +
  geom_col(aes(x = tname, y = alen, fill=name), color='black') +
  scale_fill_discrete(guide='none') +
  facet_wrap(~assignment, nrow = 3, scales='free') +
  theme(axis.text.x = element_text(angle = 270))

# %>% 
#   right_join(rukki_paths)
