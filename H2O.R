## Load H2O library
library(h2o)

## Connect to H2OCluster
h2o.init(nthreads = -1)

##Define data path
base_path = normalizePath("~/kratschie/LMU/Binf/gobi/Blockteil/")
gene_counts_path = paste0(base_path, "/")
junction_counts_path = paste0(base_path, "/")


gene_counts_hex = h2o.importFile(PATH = gene_counts_path, destination_frame = "gene_counts_hex")
junction_counts_hex = h2o.importFile(PATH = junction_counts_path, destination = "junction_counts_hex")


## Summary
h2o.describe(gene_counts_hex)
h2o.describe(junction_counts_hex)

## plot 
h2o.hist()