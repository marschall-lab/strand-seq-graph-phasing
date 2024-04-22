
# Library -----------------------------------------------------------------


library(magrittr)
library(dplyr)
library(ggplot2)
library(patchwork)

# Functions ---------------------------------------------------------------


# score_f <- function(x, y, lambda1=1, lambda2=1) {
#   stopifnot(length(x) == length(y))
#   
#   objective_score <-
#     sum(abs(x)) %>% 
#     as.vector()
#   
#   constraint_1_score <-
#     -lambda1 * (x %*% y)^2 %>% 
#     as.vector()
#     
#   constraint_2_score <-
#     -lambda2 * (x %*% x - 1)^2 %>%
#     as.vector()
# 
#   out <- objective_score + constraint_1_score + constraint_2_score
#   return(out/length(x))
# }
# 
# score_grad <- function(x, y, lambda1=1, lambda2=1) {
#   stopifnot(length(x) == length(y))
#   
#   obj_grad <-
#     sign(x)
#   
#   constraint_1_grad <-
#     -lambda1 * (2 * (x%*%y) * y)
#   
#   constraint_2_grad <-
#     -lambda2 * (2 * (x %*% x) * 2 * x - 4 * x + 1)
#   
#   out <- obj_grad + constraint_1_grad + constraint_2_grad
#   return(out/length(x))
# }
# data <- matrix(c(1, 1 , 1, -1, 1, -1), nrow=2)

score_f <- function(x, data, lambdao=1, lambda1=1) {
  stopifnot(length(x) == nrow(data))

  projections <-
    x %*% data
  
  groups <- split(projections, sign(projections))
  
  group_means <-
    map_dbl(groups, mean)
  
  objective_score <- 
    sum(group_means)^2

  group_sums <-
    map_dbl(groups, sum) 
  
  constraint_1_score <-
     sum(abs(group_sums))
  
  out <- lambdao * objective_score + lambda1 * constraint_1_score
  
  return(list(all = out/length(x), obj = lambdao * objective_score, cn1 = lambda1 * constraint_1_score))
}


score_grad <- function(x, data,lambdao=1, lambda1=1) {
  stopifnot(length(x) == nrow(data))

  projections <-
    x %*% data
  
  signs <-
    as.vector(sign(projections))
  
  neg_group <- signs == -1
  pos_group <- signs == 1
  
  n_neg <- sum(neg_group)
  n_pos <- sum(pos_group)
  
  pos_mean <- 
    (data[, pos_group, drop=FALSE]/n_pos) %>% 
    apply(1, sum)
  
  neg_mean <- 
    (data[, neg_group, drop=FALSE]/n_neg) %>% 
    apply(1, sum)
  
  objective_grad <- as.vector(2 * x %*% (pos_mean + neg_mean)^2)
  
  constraint_1_grad <- as.vector(data %*% abs(signs))
  
  out <- lambdao * objective_grad + lambda1*constraint_1_grad
  return(list(all = out/length(x), obj = lambdao * objective_grad, cn1= lambda1 * constraint_1_grad ))
}


gradient_desc <-
  function(data,
           lambdao = 1,
           lambda1 = 1,
           # lambda2 = 1,
           lambda1_growth_rate = 1 + 1e-4,
           # lambda1_max = 5e2,
           learn_rate = .005,
           conv_threshold = 1e-7,
           max_iter = 5e3) {
    
  # x <- c(0,1)
  x <- rnorm(nrow(data)) # + runif(length(y), min = -0.1, max=0.1)
  x <- x / sqrt(sum(x^2))
  score <- score_f(x, data, lambda1)
  converged <- FALSE
  iterations <- 0
  xs <- list(x)
  scores <- as_tibble(score)
  grads <- tibble()
  while(!converged) {
    # if(lambda1 <= lambda1_max) lambda1 <- lambda1 * lambda1_growth_rate
    
    ## Implement the gradient descent algorithm
    grad <- score_grad(x, data, lambdao, lambda1)
    x <- x - learn_rate * grad$all
    x <- x / sqrt(sum(x^2))

    score_new <- score_f(x, data, lambdao, lambda1)

    xs <- append(xs, list(x))
    scores <- bind_rows(scores, as_tibble(score_new))
    # if(iterations > 100 && abs(score - score_new) <= conv_threshold) {
    #   converged <- TRUE
    #   print('converged!')
    #   return(x)
    # }
    iterations = iterations + 1
    print(iterations/max_iter)
    if(iterations > max_iter) { 
      converged <- TRUE
      print('max_iters reached')
      print(lambda1)
      plot_data <-
        bind_cols(scores, tibble(x=xs))
      return(plot_data)
    }
  }
}


# Test --------------------------------------------------------------------
uv <- function(x) x/sqrt(sum(x^2))

ww_basis <- sample(c(-1, 0, 1), size = 96, replace = TRUE, prob = c(0.25, 0.5, 0.25))
wc_basis <- sample(c(-1, 1), size=96, replace = TRUE) * (ww_basis == 0)

hap1 <- ww_basis + wc_basis
hap2 <- ww_basis - wc_basis
hap_hom <- ww_basis

data <- matrix(c(hap1, hap2, hap_hom), ncol=3, byrow = FALSE)
data <- apply(data, 2, function(x) {x / sqrt(sum(x^2, na.rm=TRUE))}) 
plot_data <- gradient_desc(data, lambdao = 1, lambda1 = 1, max_iter = 20000, learn_rate = 0.01)


p1 <-
  plot_data %>% 
  mutate(n = 1:n()) %>% 
  mutate(wc_cos_sim = map_dbl(x, function(xx) abs(sum(xx * uv(wc_basis))))) %>% 
  select(-x) %>% 
  tidyr::pivot_longer(cols = -n) %>% 
  ggplot() +
  geom_point(aes(x=n, y=value)) +
  facet_wrap(~name, scales = 'free') +
  scale_y_continuous(limits = c(0, NA))

p1

x <- plot_data$x[[nrow(plot_data)]] 

# cube_projection
cube_y <- y / max(abs(y))
cube_z <- 1-cube_y
z <- cube_z / (sqrt(as.vector(cube_z %*% cube_z)))

modified_gram_schmidt <- function(v1, v2) {
  # Normalize the first vector
  v1_norm <- v1 / sqrt(sum(v1^2))
  
  # Project the second vector onto the first
  proj <- sum(v1_norm * v2) * v1_norm
  
  # Subtract the projection from the second vector
  v2_orthogonal <- v2 - proj
  
  # Normalize the orthogonalized second vector
  v2_orthogonal_norm <- v2_orthogonal / sqrt(sum(v2_orthogonal^2))
  
  # Return the normalized vectors
  return(v2_orthogonal_norm)
}

 
# z <- 1/y
# z <- z / sqrt(as.vector(z%*%z))
pd <- tibble(x=x, y=y, z=z, g=modified_gram_schmidt(y, z),  s=modified_gram_schmidt(y, x)) %>% 
  mutate(n=1:n()) %>% 
  tidyr::pivot_longer(c('x', 'y', 'z', 'g', 's'))

p2 <- 
  ggplot(pd) + geom_point(aes(x=n, y=value, color =name), alpha=0.4)

p1 + p2

# TODO: don't try to find an orthogonal basis, simply project each vector onto the WW vector, and then take the negative of the perependicular component.
# What vector to project onto? Onto every other vector and then average? Or project onto the average vector?

project_through <- function(x, mirror_line) {
  # Project a a vector over another vector using gram-schmidt decomposition
  stopifnot(length(mirror_line) == length(x))

  mirror_line <- mirror_line/sqrt(sum(mirror_line^2, na.rm=TRUE))
  x_proj <- sum(x * mirror_line, na.rm=TRUE) * mirror_line
  x_ortho <- x - x_proj
  x_new <- x_proj - x_ortho

  return(x_new)
}



x <- c(1, 0, NA, 1)
y <- c(1, 1, 1, NA)

mirror_line <- c(1,0.2)
x <- c(1, -1)
project_through(x, mirror_line)
plot(c(1,1,1), c(0, 1, -1))

     