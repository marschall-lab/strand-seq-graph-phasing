
# Library -----------------------------------------------------------------

library(magick)
library(purrr)

# Import ------------------------------------------------------------------

pdf_files <-
  list.files('intermediate_output/', recursive = TRUE, pattern = '*plots.pdf', full.names = TRUE) %>% 
  set_names(function(x) basename(dirname(x)))

# Output Layout -----------------------------------------------------------

n_pages <- 
  image_read(pdf_files[[1]]) %>% 
  image_info() %>% 
  nrow()

n_pages <- 4
# 
# width <-ceiling(sqrt(length(pdf_files)))
# arrange_factor <- rep(seq_len(width), each=width)[seq_along(pdf_files)]

# Merge Plots -------------------------------------------------------------

plots <-
  map(seq_len(n_pages), function(i) {
    # Annotate with sample name
    images <- 
      map(pdf_files, image_read_pdf, pages=i, density=100) %>% 
      imap(image_annotate, gravity='north') %>% 
      reduce(c)
    
    return(images)
    
    # p <- 
    #   split(images, arrange_factor) %>% 
    #   map(reduce, function(x, y) image_append(c(x,y))) %>% 
    #   reduce(function(x, y) image_append(c(x, y), stack=TRUE))
    # 
    # return(p)
  })

# plots <-
#   reduce(plots, c)

# Export ------------------------------------------------------------------

iwalk(plots, function(x, i){
  image_write(x, paste0('check', i, '.pdf'), format='pdf')
})
# image_write(plots, 'check.pdf', format='pdf')
