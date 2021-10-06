library(permute)
library(ChAMP)
library(glue)

correct_idats_path <- function(path_idats){
  if (endsWith(path_idats, "/")){
    substr(path_idats, 1, nchar(path_idats) - 1)
  }
  else{
    return(path_idats)
  }
}

delete_temp_sample_sheet <- function(path_idats){
  temp_file = glue(path_idats, "temp_sample_sheet.csv")
  if (file.exists(temp_file)) {
    file.remove(temp_file)
  }
}


check_if_file_exists <- function(path){
  if (!(file.exists(path))) {
    cat("Dir", path, "does not exists", "\n")
  }
}

count_files <- function(path, file_type, file_format){
  files = length(list.files(path, file_format))
  cat("Founded", file_type, "no.:", files, "\n")
  return(length(files))
}

check_if_necessary <- function(chunk_size, samples_number){
  if(chunk_size > samples_number){
    cat("Chunk size", chunk_size, ">", samples_number, "samples number", "\n")
    stop("Stopped.")
  }
}

create_dir <- function(path, name, extension){
  path = glue(path, name, "_", extension)
  
  if ((file.exists(path))) {
    cat("Dir", path, "already exists", "\n")
    stop("Stopped.")
  }
  
  dir.create(path)
  return(path)
}

overlap_cpgs <- function(list_of_mynorms_mynorms){
  
  cpgs_per_mynorm <- list()
  idx <- 1
  
  for (mynorm in list_of_mynorms_mynorms){
    cpgs_per_mynorm[[idx]] <- row.names(mynorm)
    idx <- idx + 1
  }
  common_cpgs <- Reduce(intersect, cpgs_per_mynorm)
  cat("Found: ", length(common_cpgs), "common across chunks.", "\n")
  return(common_cpgs)
}

concate_mynorms <- function(list_of_mynorms, cpgs_common){
  
  idx <- 1
  for (mynorm in list_of_mynorms){
    list_of_mynorms[[idx]] <- mynorm[cpgs_common, ]
    idx <- idx + 1
  }
  
  myNorm <- do.call("cbind", list_of_mynorms)
  return(myNorm)
}

run_champ <- function(path_idats, QC_path, Norm_path, array_type, force, norm_type, cores){
  
  myLoad <- champ.load(directory = path_idats,
                       method="champ",
                       methValue="B",
                       autoimpute=TRUE,
                       filterDetP=TRUE,
                       ProbeCutoff=0,
                       SampleCutoff=0.1,
                       detPcut=0.01,
                       filterBeads=TRUE,
                       beadCutoff=0.05,
                       filterNoCG=TRUE,
                       filterSNPs=TRUE,
                       population=NULL,
                       filterMultiHit=TRUE,
                       filterXY=TRUE,
                       force=force,
                       arraytype=array_type)
  
  champ.QC(beta = myLoad$beta,
           pheno=myLoad$pd$Sample_Group,
           mdsPlot=TRUE,
           densityPlot=TRUE,
           dendrogram=TRUE,
           PDFplot=TRUE,
           Rplot=TRUE,
           resultsDir=QC_path)
  
  myNorm <- champ.norm(beta=myLoad$beta,
                       rgSet=myLoad$rgSet,
                       mset=myLoad$mset,
                       resultsDir=Norm_path,
                       method=norm_type,
                       plotBMIQ=TRUE,
                       arraytype=array_type,
                       cores = cores)
  
  return(myNorm)
}

load_chunks <- function(path, pattern){
  
  chunks <- list()
  files <- list.files(path, pattern)

  idx <- 1
  for (file in files){
    chunk_path = glue(path, file)
    cat("Loading: ", chunk_path, "\n")
    
    chunk_df <- data.table::fread(chunk_path, data.table = FALSE)
    chunk_df <- data.frame(chunk_df, row.names = 1)
    chunks[[idx]] <- chunk_df 
    idx <- idx + 1
  }
  return(chunks)
}

delete_temp_files <- function(path, pattern){
  temp_files <- list.files(path, pattern)
  
  for (file in temp_files){
    chunk_path = glue(path, file)
    cat("Deleting: ", chunk_path, "\n")
    unlink(chunk_path, recursive = FALSE, force = FALSE)
  }
}

