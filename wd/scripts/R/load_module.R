library(here)

get_script_dir <- function() {
  wd <- getwd()
  write(wd, stdout())
  return(wd)
}


all.sources <-
  c(
    "module_saarclust/calcProbs.R",
    "module_saarclust/countDirectionalReads.R",
    "module_saarclust/dumped_functions.R",
    "module_saarclust/EMclust.R",
    "module_saarclust/export.R",
    "module_saarclust/findClusterPartners.R",
    "module_saarclust/hardClust.R",
    "module_saarclust/helperFuctions.R",
    "module_saarclust/import.R",
    "module_saarclust/importReads.R",
    "module_saarclust/SaaRclust_evaluation_plots.R",
    "module_saarclust/SaaRclust.R",
    "module_saarclust/timedMessage.R",
    "module_saarclust/utils.R",
    "module_saarclust/wrapper_parallel.R",
    "module_saarclust/wrapper.R"
  )

# Path Handling
all.sources <- paste0(get_script_dir(), '/wd/scripts/R/', all.sources)

invisible(sapply(all.sources, source))
