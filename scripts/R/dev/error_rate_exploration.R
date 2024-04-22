library(tidyverse)

n_libs <- 96
plot_data <- 
  tibble(n_ww = 1:n_libs, n_wc = n_libs - n_ww) %>% 
  mutate(prob = dbinom(n_ww, n_libs, 0.5)) %>% 
  mutate(eta = list(seq(0, 1, 0.01))) %>% 
  tidyr::unnest(eta)

make_vect <- function(n_ww, n_wc, eta) {
  c(rep(1, n_ww), rep(eta, n_wc))
}

compare <- function(n_ww, n_wc, eta) {
  v_no_noise <- make_vect(n_ww, n_wc, eta=0)
  v_w_noise <- make_vect(n_ww, n_wc, eta=eta)
  # if(eta > 0) browser()
  cosine_similarity_(v_no_noise, v_w_noise)
}

labeller <- function(x) {
  (100 + 100*x)/(100-100*x)
}

col_labeller <- function(x) {
  round(acos(x) * 180/pi, 1)
}
plot_data %>% 
  mutate(cos_sim = pmap_dbl(list(n_ww, n_wc, eta), compare)) %>% 
  ggplot() +
  geom_raster(aes(x = n_ww, y=eta, fill=cos_sim, alpha=prob)) +
  scale_fill_viridis_c(n.breaks = 20) +
  scale_y_continuous(n.breaks = 10, labels = labeller) +
  scale_x_continuous(n.breaks=20)

acos(0.9) * 180/pi
