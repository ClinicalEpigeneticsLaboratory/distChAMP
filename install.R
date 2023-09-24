install.packages(c("glue", "permute", "arrow"))

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ChAMP")
