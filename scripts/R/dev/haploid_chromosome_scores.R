

## Haploid Seeding --------------------------------------------------------- #
#This is to handle the unusual case if there is a haploid chromosome present as
#only one unitig. Can a haploid score be calculated to better create one-unitig
#clusters from complete haploid chromosomes? (Haploid chromosomes present in one
#unitig in the graph)

exact_match_counts_df %>% 
  rename(em_c=c, em_w=w, em_n=n) %>%
  left_join(counts_df) %>% 
  mutate(rat = (em_n+1)/(n+1)) %>% mutate(rat = (em_n+1)/(n+1)) %>% 
  group_by(unitig) %>% 
  summarise(rat = mean(rat, na.rm=TRUE), n=sum(n, na.rm=TRUE), em_n=sum(em_n, na.rm=TRUE)) %>% 
  mutate(hap_unitig = unitig %in% paste0('utig4-', c(1122,282,456,371,1123,13))) %>%
  ggplot(aes(x=log10(em_n), y=log10(n))) + 
  geom_point(aes(x=log10(em_n), y=-log2(rat), color=hap_unitig)) 

geom_point(aes(color=hap_unitig)) +
  geom_smooth(method = 'lm')

pd <- exact_match_counts_df %>% 
  rename(em_c=c, em_w=w, em_n=n) %>%
  left_join(counts_df) %>% 
  left_join(unitig_lengths_df) %>%
  mutate(rat = (em_n+1)/(n+1)) %>% 
  mutate(ssf = ifelse(n >= 4, (w-c)/n, NA), ssf_em = ifelse(em_n >= 4, (em_w-em_c)/em_n, NA)) %>% 
  group_by(unitig) %>% 
  summarise(sim = pairwise_complete_cosine_similarity_(ssf, ssf_em), n=sum(n, na.rm=TRUE), em_n=sum(em_n, na.rm=TRUE), rat = mean(rat, na.rm=TRUE) ) %>% 
  ungroup() %>%
  arrange(desc(sim)) %>%
  mutate(hap_unitig = unitig %in% paste0('utig4-', c(1122,282,456,371,1123,13))) %>% 
  left_join(unitig_lengths_df)

pd %>% 
  ggplot(aes(x=log10(em_n), y=log10(n))) + 
  geom_point(aes(x=log10(em_n),y=-log2(rat), color=hap_unitig)) + 
  geom_smooth(method = 'lm') + 
  scale_fill_viridis_c()

pd %>% ggplot() + # geom_point(aes(x = log10(n), y = sim, color=hap_unitig))
  geom_point(aes(x = log10(rat), y = sim, color=hap_unitig))
