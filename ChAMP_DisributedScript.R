source("utils/utils.R")
set.seed(101) 

run_distributed_champ <-
  function(path_ss,
           path_idats,
           output,
           method = "champ",
           array_type = "EPIC",
           force = TRUE,
           norm_type = "BMIQ",
           cores = 1,
           chunk_size = 50) {

    check_if_file_exists(path_idats)
    check_if_file_exists(path_ss)
    n_csv_files <- count_files(path_idats, "SampleSheet", ".csv")
    
    if (n_csv_files > 1) {
      stop("Please remove SampleSheet from idat directory!")
    }
    
    count_files(path_idats, "IDATs", ".idat")
    check_if_file_exists(output)
    
    sample_sheet <- read.csv(path_ss, row.names = 1)
    
    n_rows <- length(row.names(sample_sheet))
    check_if_necessary(chunk_size, n_rows)
    
    randomized_samples <- shuffle(rownames(sample_sheet))
    
    for (i in 1:ceiling(n_rows / chunk_size)) {
      
      # Crete QC path
      QC_path <- create_dir(output, "QC", i)
      
      # norm path
      norm_path <- create_dir(output, "Norm", i)
      
      cat(
        "Run",
        i,
        "/",
        ceiling(n_rows / chunk_size),
        "Each up to: ",
        chunk_size ,
        "samples.",
        "\n"
      )
      
      idx_start <- chunk_size * (i - 1) + 1
      idx_end <- (chunk_size * i)
      
      if (idx_end > n_rows) {
        idx_end <- n_rows
      }
      
      temp_samples <- randomized_samples[idx_start:idx_end]
      temp_sample_sheet <- sample_sheet[temp_samples,]
      
      temp_sample_sheet["Sample_Name"] <-
        rownames(temp_sample_sheet)
      
      write.csv(temp_sample_sheet,
                file.path(path_idats, "temp_sample_sheet.csv"),
                sep = ",")
      
      mynorm_temp <-
        run_champ(
          path_idats = path_idats,
          QC_path = QC_path,
          Norm_path = norm_path,
          method = method,
          array_type = array_type,
          force = force,
          norm_type = norm_type,
          cores = cores
        )
      
      delete_temp_sample_sheet(path_idats)
      
      path = file.path(output, glue(i, "_temp_chunk", ".parquet"))
      
      mynorm_temp$CpG <- rownames(mynorm_temp)
      write.parquet(mynorm_temp, path)
      message("Saved chunk no. ", i, "\n")
    }
    
    message("Loading chunks ...", "\n")
    mynorms <- load_chunks(output, "_temp_chunk.parquet")
    
    message("Looking for CpGs common across batches ...", "\n")
    common_cpg <- overlap_cpgs(mynorms)
    
    message("Concating batches into myNorm ...", "\n")
    myNorm <- concate_mynorms(mynorms, common_cpg)
    
    mynorm_path = glue(output, "myNorm.csv")
    cat("Saving myNorm ...", "\n")
    
    write.csv(myNorm, mynorm_path)
    
    cat("Removing temporary files ...", "\n")
    delete_temp_files(output, "temp_chunk.csv")
    
    message("DONE!")
    
  }
