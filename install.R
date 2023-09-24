install.packages(c("glue", "permute", "data.table", "arrow"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChAMP")
