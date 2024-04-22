
x <- make_vector(haploid=TRUE)
y <- make_vector(haploid=TRUE)
x_par <- ifelse(x != y, 0, x)

uv(x) %*% uv(x_par)
uv(y) %*% uv(x_par)

make_vector <- function(n_libs=96, haploid=FALSE) {
  if(!haploid) {
    sample(c(-1, 0, 1), n_libs, replace=TRUE, prob = c(0.25, 0.5, 0.25))
  } else {
    sample(c(-1, 1), n_libs, replace=TRUE, prob = c(0.5, 0.5))
  }
} 

simulate_unitigs <- function(n_unitigs=20, n_libs=96, noise_sd = 0.1, haploid=FALSE) {
  core_vector <- make_vector(n_libs, haploid)

  map(1:n_unitigs, function(unused) {
    noises <-
      rnorm(length(core_vector), mean = 0, sd=noise_sd)
    
    out <- case_when(
      core_vector == -1 ~ core_vector + abs(noises),
      core_vector == 0 ~ noises,
      core_vector == 1 ~ core_vector - abs(noises)
    )
    
    return(out)
  }) %>% 
    reduce(rbind)
}

simulate_experiment <- function(n_libs=96, n_dip=23, n_hap=2, n_unitigs=10, noise_sd=0.1) {
  
  c(
    map(1:n_hap,function(unused) simulate_unitigs(n_unitigs, n_libs, noise_sd, haploid = TRUE)),
    map(1:n_dip,function(unused) simulate_unitigs(n_unitigs, n_libs, noise_sd))
  )%>% 
    reduce(rbind) %>% 
    pairwise_complete_cosine_similarity() %>% 
    (function(x) x[upper.tri(x, diag = FALSE)]) 
}
  
simulate_perfect_experiment <- function(n_libs=96, n_dip=23, n_hap=2){
  c(
    map(1:n_hap,function(unused) make_vector(n_libs, haploid = TRUE)),
    map(1:n_dip,function(unused) make_vector(n_libs))
  )%>% 
    reduce(rbind) %>% 
    pairwise_complete_cosine_similarity() %>% 
    (function(x) x[upper.tri(x, diag = FALSE)]) 
}

results <-
  map(1:500, function(unused) simulate_perfect_experiment(24, 22, 2))

plot_data <-
  results %>% 
  map(function(x) tibble(sim = x)) %>% 
  bind_rows(.id='replication') 

plot_data %>% 
  ggplot() +
  geom_histogram(aes(abs(sim)), bins=100) +
  geom_vline(xintercept = 0.5)


results <-
  map(1:10, function(unused) simulate_experiment(30, 22, 2, n_unitigs=10, noise_sd=0.2))

results %>% 
  map(function(x) tibble(sim = x)) %>% 
  bind_rows(.id='replication') %>% 
  ggplot() +
  geom_histogram(aes(abs(sim)), bins=100) +
  geom_vline(xintercept = 0.5)
