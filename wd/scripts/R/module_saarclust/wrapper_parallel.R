
#Load the function below into R if you want to run all steps in one command [ONGOING]

library("foreach")
library("doParallel")

inputfolder <- "/media/daewoooo/WORK/SS2PacBio_alignment_HG00733/Test_cluster_chr21&chr22/Minimap_out_chunks/"

runSaaRclust_parallel <- function(inputfolder=NULL, outputfolder="./SaaRclust_results", num.clusters=44, EM.iter=100, chunk.size=20000, alpha=0.1, numCPU=2) {

  
  data.store <- file.path(outputfolder, 'data')
  file.list <- list.files(path = inputfolder, pattern = "mimimapChunk")
  
  if (numCPU > 1) {
    message("Using ",numCPU," CPUs")
    cl <- makeCluster(numCPU)
    doParallel::registerDoParallel(cl)
    
    ptm <- startTimedMessage("Running clustering ...") 
    temp <- foreach (file = file.list) %dopar% {
    #temp <- foreach (file = file.list, .packages=c('SaaRclust')) %dopar% {  
      tC <- tryCatch({
        tab.in <- importTestData(infile = file.path(inputfolder, file))
      }, error = function(err) {
        stop(file,'\n',err)
      })
    }
    stopCluster(cl)
    stopTimedMessage(ptm)
  }
}
