library(permute)
library(ChAMP)
library(glue)

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

check_if_necessary <- function(batch_size, samples_number){
  if(batch_size > samples_number){
    cat("Batch size", batch_size, ">", samples_number, "samples number", "\n")
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
  
  return(common_cpgs)
  }

concate_mynorms <- function(list_of_mynorms_mynorms, cpgs_common){
  
  idx <- 1
  for (mynorm in list_of_mynorms_mynorms){
    list_of_mynorms_mynorms[[idx]] <- mynorm[cpgs_common, ]
    idx <- idx + 1
  }
  
  myNorm <- do.call("cbind", list_of_mynorms_mynorms)
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

run_distributed_champ <- function(path_ss, path_idats, output, array_type = "EPIC", force = TRUE, norm_type = "BMIQ", cores = 1, chunk_size = 100){
  
  check_if_file_exists(path_idats)
  check_if_file_exists(path_ss)
  n_csv_files <- count_files(path_idats, "SampleSheet", ".csv")
  
  if (n_csv_files > 1){stop("Remove Sample Sheet from idat directory")}
  
  count_files(path_idats, "IDATs", ".idat")
  check_if_file_exists(output)
  
  sample_sheet <- read.csv(path_ss, row.names = 1)

  n_rows <- length(row.names(sample_sheet))
  check_if_necessary(chunk_size, n_rows)
  
  randomized_samples <- shuffle(rownames(sample_sheet))
  mynorms <- list()
  
  for (i in 1:ceiling(n_rows / chunk_size)){

      # Crete QC path
      QC_path <- create_dir(output, "QC", i)
      
      # # Norm path
      Norm_path <- create_dir(output, "Norm", i)

      cat("Run", i, "/", ceiling(n_rows / batch_size), "Each up to: ", batch_size ,"samples.", "\n")
      
      idx_start <- chunk_size * (i - 1) + 1
      idx_end <- (chunk_size * i)
      if (idx_end > n_rows){idx_end <- n_rows}
  
      temp_samples <- randomized_samples[idx_start:idx_end]
      temp_sample_sheet <- sample_sheet[temp_samples, ]
      
      temp_sample_sheet["Sample_Name"] <- rownames(temp_sample_sheet)
      write.csv(temp_sample_sheet, glue(path_idats, "temp_sample_sheet.csv"), sep=",")
  
      mynorm_temp <- run_champ(path_idats = path_idats, QC_path = QC_path, Norm_path = Norm_path, array_type = array_type, force = force, norm_type = norm_type, cores = cores)
      delete_temp_sample_sheet(path_idats)
      mynorms[[i]] <- mynorm_temp
  
  }
  cat("Looking for CpGs common across batches.", "\n")
  common_cpg <- overlap_cpgs(mynorms)
  
  cat("Concating batches into myNorm.", "\n")
  myNorm <- concate_mynorms(mynorms, common_cpg)
  
  cat("Saving myNorm.", "\n")
  path = glue(output, "myNorm", ".csv")
  write.csv(myNorm, path, sep=",")
}
