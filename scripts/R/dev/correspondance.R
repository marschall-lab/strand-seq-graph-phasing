
# WD ----------------------------------------------------------------------

setwd("~/Documents/wd")
# Library -----------------------------------------------------------------
library(magrittr)
library(dplyr)
library(purrr)
library(furrr)
library(ggplot2)
library(gridExtra)
library(pafr)
library(ggokabeito)



# Pafs --------------------------------------------------------------------


reference_alignments <- 
  list(
    final_to_self= "correspondance_test/HG03683.ps-sseq.wd/assembly-hpc_ref-aln.paf",
       final_to_ref = "reference_alignments/T2Tv11_hg002Yv2_chm13/HG03683_T2Tv11_hg002Yv2_chm13_ref-aln.paf",
       fixed_to_self= "correspondance_test//HG03683.ps-sseqfix.wd/assembly-hpc_ref-aln.paf",
       fixed_to_ref = "reference_alignments/T2Tv11_hg002Yv2_chm13/HG03683stage4_T2Tv11_hg002Yv2_chm13_ref-aln.paf",
       fixed_to_final = "correspondance_test/HG03683/assembly_ref-aln.paf"
    )
  

reference_alignments <-
  reference_alignments %>%
  map(pafr::read_paf,
      tibble = TRUE,
      include_tags = FALSE) 

SAVE <- reference_alignments


reference_alignments
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
  list(HG02059='rukki/HG02059/HG02059_rukki_paths.tsv') %>% 
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



rukki_paths <-
  rukki_paths %>% 
  mutate(is_not_gap = !is_gap) 



# Scratch -----------------------------------------------------------------


# Hom_Check ---------------------------------------------------------------

# HPC to Reference --------------------------------------------------------



hpc_to_ref <- 
  list(
    HG02059= "correspondance_test/ref_aln/HG02059.ps-sseq.wd/assembly-hpc_ref-aln.paf"#,
    # HG00733= "correspondance_test/ref_aln/HG00733.ps-sseq.wd/assembly-hpc_ref-aln.paf"
  )


hpc_to_ref <-
  hpc_to_ref %>%
  map(pafr::read_paf,
      tibble = TRUE,
      include_tags = TRUE) %>% 
  bind_rows(.id= 'sample')

hpc_to_ref <-
  hpc_to_ref %>%
  group_by(sample, qname) %>%
  # slice_head(n=1) %>%
  dplyr::rename(unitig = qname)

# hpc_to_ref <-
#   hpc_to_ref %>%
#   rename(unitig = qname)


SAVE <- hpc_to_ref


hpc_to_ref2 <-
  hpc_to_ref %>%
  group_by(sample, unitig) %>%
  slice_head(n=2) 
hpc_to_ref
# take the first alignment listed each time. I think there is an order to them
# as output by minimap so this seems okay?



  
rukki_paths <-
  list(
    HG02059='rukki/HG02059/HG02059_rukki_paths.tsv',
    HG00733 = 'rukki/HG00733/HG00733_rukki_paths.tsv') %>% 
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


hom_nodes_rukki <-
  rukki_paths %>% 
  filter(!is_gap) %>% 
  group_by(sample, unitig) %>% 
  mutate(n = n()) %>% 
  filter('HAPLOTYPE1' %in% assignment & 'HAPLOTYPE2' %in% assignment) %>% 
  ungroup() %>%
  filter(n == 2) %>% 
  distinct(sample, unitig) %>% 
  mutate(is_hom = TRUE)

# HG02059
het_unitigs <-
  c('utig4-2367', 'utig4-1106', 'utig4-788', 'utig4-789', 'utig4-790', 'utig4-791')

hom_unitigs <-
  c('utig4-982', 'utig4-966', 'utig4-1389','utig4-1105', 'utig4-787', 'utig4-780', 'utig4-2411', 'utig4-800', 'utig4-779', 'uti4-722'
  )

## chr7
hom_unitigs <-
  c('utig4-799', 'utig4-2408', 'utig4-1492','utig4-747', 'utig4-1294', 'utig4-2903', 'utig4-819', 'utig4-816', 'utig4-1473'
  )

## HG00733
# het_unitigs <- 
#   c('utig4-1754', 'utig4-1755','utig4-1900', 'utig4-1901')
# 
# hom_unitigs <-
#   c('utig4-1923', 'utig4-1935', 'utig4-1916','utig4-536', 'utig4-1913', 'utig4-594', 'utig4-1938', 'utig4-1920'
#   )




hpc_to_ref2 %>% 
  filter(sample == 'HG02059') %>% 
  # filter(tp == 'P') %>% 
  filter(unitig %in% c(hom_unitigs)) %>% 
  arrange(unitig) %>% 
  print(n=35)


rukki_paths %>% 
  filter(sample == 'HG02059') %>% 
  filter(unitig %in% c(hom_unitigs)) %>% 
  arrange(unitig)

# 
# 
# 
# hpc_to_ref %>% 
#   filter(sample == 'HG00733') %>% 
#   filter(!grepl('unassigned', tname)) %>% 
#   # filter(alen >= 1e5) %>% 
#   filter(qname %in% hom_unitigs) %>% 
#   arrange(qname) 
# 
#  
# hpc_to_ref %>% 
#   filter(sample == 'HG00733') %>% 
#   # filter(!grepl('unassigned', tname)) %>% 
#   filter(alen >= 1e5) %>% 
#   # filter(tname %in% c('haplotype1-0000004', 'haplotype2-0000167')) %>% 
#   filter(qname %in% hom_unitigs) %>% 
#   group_by(qname) %>% 
#   slice_head(n=2) %>% 
#   arrange(qname) # %>% 
#   # filter(!(qname %in% c('utig4-1105'))) 
# 
# 
# hpc_to_ref %>% 
#   filter(alen >= 1e5) %>% 
#   group_by(tname) %>% 
#   filter(any(qname %in% c(het_unitigs, hom_unitigs))) %>% 
#   ggplot() +
#   geom_segment(aes(x=tstart, xend=tend, y=tname, yend=tname, group=tp, color=qname %in% het_unitigs)) +
#   facet_wrap(~qname, scales = 'free')
# 
# 
# hom_df <-
#   hpc_to_ref2 %>% 
#   semi_join(hom_nodes) %>% 
#   filter(qlen >= 1e6) %>% 
#   filter(!grepl('unassigned', tname)) %>% 
#   arrange(sample, unitig)
# 
# hom_df
# 
# 
# rukki_paths <-
#   rukki_paths %>%
#   mutate(gap_size = ifelse(is_gap, stringr::str_extract(unitig, '[0-9]+'), NA)) %>%
#   mutate(gap_size = as.integer(gap_size)) %>%
#   left_join(hpc_to_ref, by = c('sample', 'unitig')) %>%
#   mutate(node_len = coalesce(qlen, gap_size)) %>%
#   select(sample, name, unitig, assignment, order, is_gap, run, node_len)



rukki_paths <-
  rukki_paths %>% 
  mutate(is_not_gap = !is_gap) 


# Homology Check ----------------------------------------------------------



reference_alignments <- 
  list(
    HG02059= "correspondance_test/self_homology/HG02059.ps-sseq.wd/assembly-homology.paf"
  )


reference_alignments <-
  reference_alignments %>%
  map(pafr::read_paf,
      tibble = TRUE,
      include_tags = TRUE) %>% 
  bind_rows(.id= 'sample')

# reference_alignments <-
#   reference_alignments %>% 
#   group_by(sample, qname) %>% 
#   slice_head(n=1) %>% 
#   dplyr::rename(unitig = qname)

# reference_alignments <-
#   reference_alignments %>%
#   rename(unitig = qname)


SAVE <- reference_alignments


reference_alignments2 <-
  reference_alignments %>%
  filter(tp == 'P') %>% 
  filter(alen >= 250e3) %>% 
  filter(tlen >= 1e6)

tlens <-
  reference_alignments2 %>% 
  distinct(tname, tlen)

anchors <-
  reference_alignments2 %>% 
  group_by(tname) %>% 
  arrange(desc(alen)) %>% 
  slice(c(1,1)) %>% 
  mutate(tstart = unique(c(0, tlen)), tend = unique(c(0, tlen)))

reference_alignments2 %>% 
  # bind_rows(anchors) %>% 
  mutate(color = case_when(
    qname == tname ~ 'identity', 
    grepl('haplotype1', qname) ~ 'hap1',
    grepl('haplotype2', qname) ~ 'hap2',
    grepl('unassigned', qname) ~ 'unassigned',
    TRUE ~ NA)) %>% 
  ggplot() +
  geom_segment(aes(x=tstart, xend=tend, y=qname, yend=qname, color=color, alpha=mapq), linewidth=1.25) +
  geom_vline(xintercept=0, alpha=0) +
  geom_vline(aes(xintercept=tlen), alpha=0, data=tlens) +
  facet_wrap(~tname, scales = 'free') +
  scale_color_okabe_ito(order = c(1,2,3,5,6,7,8,4)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_line(color = 'darkgrey', linewidth = 0.1),
    panel.background = element_blank()
  ) +
  ggtitle('HG02059')



# VCF -----------------------------------------------------------------

vcfs <-
  list(
    het = vcfR::read.vcfR('~/Documents/snp_hom__wd/output/evaluation/qv_estimation/variant_calls/30-split-gtype/HG02059_f_f1-VN061.ERR3988981_short_map-to_assembly.snvs.het.vcf.bgz'),
    hom = vcfR::read.vcfR('~/Documents/snp_hom__wd/output/evaluation/qv_estimation/variant_calls/30-split-gtype/HG02059_f_f1-VN061.ERR3988981_short_map-to_assembly.snvs.hom.vcf.bgz')
  )

vcf_data <- 
  vcfs %>% 
  map(function(vcf) 
    as_tibble(vcf@fix) %>% 
      bind_cols(as_tibble(vcf@gt)) %>% 
      mutate(rnum = 1:n())
  ) %>% 
  bind_rows(.id='snp_type') %>% 
  mutate(POS = as.integer(POS)) %>% 
  dplyr::rename(tname=CHROM) %>%
  mutate(qname = tname)

vcf_data %<>% 
  mutate(genotype = map_chr(HG02059, function(x) stringr::str_split_1(x, ':')[1]))

pd <-
  vcf_data %>%
  filter(!grepl('unassigned', tname)) %>%
  filter(!grepl('unassigned', qname)) %>% 
  group_by(tname) %>% 
  filter(max(POS) >= 10e6) %>% 
  ungroup() %>% 
  mutate(snp_sign = ifelse(snp_type == 'het', 1, 1))

bin_size_kbp <- 250
per_window_kbp <- 100
ymax <- 5

ggplot(pd) +
  geom_histogram(
    aes(
      POS,
      y = after_stat(count) / after_stat(width) * per_window_kbp*1e3,
      weight=snp_sign,
      group = genotype,
      fill=genotype,
    ),
    # color='black',
    position='dodge',
    # bins=20,
    binwidth = bin_size_kbp*1e3,
    linewidth=0.1
    
  ) +
  facet_wrap(~ tname, scales = 'free_x', ncol=5) +
  ylab(glue::glue("SNVs per {per_window_kbp}Kbp")) +
  xlab(glue::glue("Position ({bin_size_kbp}Kbp bins)")) +
  ggtitle('HG02059, (ymax:5)') +
  coord_cartesian(ylim=c(0, ymax)) +
  scale_fill_okabe_ito(order = c(1,2,3))
  # scale_shape_manual(values=c(het=1, hom=21))




vcf_data <-
  vcf_data %>% 
  semi_join(reference_alignments2, by='tname')
#vcf2@fix
# 
# contigs <-
#   c('haplotype2-0000161', 'haplotype2-0000167')



# VCF, Homology -----------------------------------------------------------


contigs <-
  c(
    'haplotype1-0000003',
    'haplotype1-0000006',
    'haplotype1-0000008',
    'haplotype1-0000009',
    'haplotype1-0000013',
    'haplotype1-0000015',
    'haplotype1-0000021',
    'haplotype1-0000023',
    'haplotype1-0000033',
    'haplotype2-0000127',
    'haplotype2-0000130',
    'haplotype2-0000137',
    'haplotype2-0000140',
    'haplotype2-0000143',
    'haplotype2-0000147',
    
    'haplotype1-0000015'

  )

contigs <-
  c(
    'haplotype1-0000032',
    'haplotype2-0000150'
  )

contigs <-
  c(
    'haplotype1-0000017',
    'haplotype2-0000124'
  )
contigs_to_keep <-
  reference_alignments2 %>% 
  filter(tlen >= 10e6) %>% 
  distinct(qname) %>% 
  pull()

plot_data <-
reference_alignments2 %>% 
  filter(tname %in% contigs) %>% 
  bind_rows(filter(anchors, tname %in% contigs)) %>%
  # filter(tname %in% contigs) %>% 
  mutate(color = case_when(
    qname == tname ~ 'identity', 
    grepl('haplotype1', qname) ~ 'hap1',
    grepl('haplotype2', qname) ~ 'hap2',
    grepl('unassigned', qname) ~ 'unassigned',
    TRUE ~ NA)) 

 
  ggplot(plot_data) +
  # geom_rug(aes(x=POS), data=filter(vcf_data, tname %in% contigs)) +
  geom_histogram(aes(POS,  y = after_stat(count) / after_stat(width) * per_window_kbp*1e3, fill=genotype), data = filter(vcf_data, tname %in% plot_data$tname), binwidth=bin_size_kbp*1e3) +
  geom_segment(aes(x=tstart, xend=tend, y=qname, yend=qname, alpha=mapq),  color='black',linewidth=1.75) +
  geom_segment(aes(x=tstart, xend=tend, y=qname, yend=qname, color=color, alpha=mapq), linewidth=1.25) +
  facet_wrap(~tname, scales = 'free') +
  scale_color_okabe_ito(order = c(1,2,3,5,6,7,8,4)) +
  scale_fill_okabe_ito(order = c(1,2,3)) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(breaks = 1:100) +
  theme_classic() +
  theme(
    axis.text.y = element_blank(),
    panel.grid.major.y = element_line(color = 'darkgrey', linewidth = 0.5),
    panel.background = element_blank()
  ) +
    ggtitle('HG02059')


tmp <- pafr::read_paf("reference_alignments/T2Tv11_hg002Yv2_chm13/HG02587_T2Tv11_hg002Yv2_chm13_ref-aln.paf", tibble = TRUE)


