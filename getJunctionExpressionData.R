library(stats)
project_accession <- "SRP019994"
pheno_keyword <- "subtypes"

##############################
##########Download############
##############################
if(!dir.exists(project_accession)){
  dir.create(project_accession)
}
setwd(project_accession)
url<-"http://duffel.rail.bio/recount"
###Gene Counts

filename <- "counts_gene.tsv.gz"
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
                quiet = FALSE, mode = "w", cacheOK = TRUE)
}

gene_counts <- read.table(file =filename, sep = '\t', header = TRUE)

### Junction Counts
filename <- "counts_jx.tsv.gz"
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
              quiet = FALSE, mode = "w", cacheOK = TRUE)
}
junction_counts <- read.table(file = filename, sep = '\t', header = TRUE)

### Phenotype
filename <- paste(project_accession, "tsv", sep=".")
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
              quiet = FALSE, mode = "w", cacheOK = TRUE)
}
phenotype <- read.table(file = filename, sep = '\t', header = TRUE)
if(length(which(phenotype$reads_downloaded==0))!=0){
  phenotype <- phenotype[-which(phenotype$reads_downloaded==0),]
}


### Junction Positions
filename <- paste(project_accession, "junction_id_with_transcripts.bed.gz", sep=".")
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
              quiet = FALSE, mode = "w", cacheOK = TRUE)
}
junction_positions <-read.table(file = filename, sep = '\t', header = FALSE)
##############################
##############################
##############################
print("Download successfully completed!")

stopifnot(ncol(gene_counts)==ncol(junction_counts) && ncol(gene_counts) == nrow(phenotype)&&
          nrow(junction_counts)==nrow(junction_positions))

concat_row_names = paste (junction_positions[,1], junction_positions[,2], junction_positions[,3], sep="_")
rownames(junction_counts) <- concat_row_names

phenotype_labels <- data.frame(char=phenotype$characteristics)
rownames(phenotype_labels) <- phenotype$run
phenotype_labels$char <- tolower(phenotype_labels$char)
phenotype_labels$char <- gsub(paste(".*", pheno_keyword, ":", sep=""), "", phenotype_labels$char)
phenotype_labels$char <- gsub("(,|\\)$).*", "", phenotype_labels$char)
phenotype_labels <- transform(phenotype_labels, char = as.integer(factor(char, unique(char))))-1

gene_counts <- log(gene_counts+1)
gc_row_sums <- rowSums(gene_counts)

hist(gc_row_sums, breaks=seq(from=0, to=max(gc_row_sums)+20, by=1))
abline(v = median(gc_row_sums),
       col = "red",
       lwd = 2)
abline(v = mean(gc_row_sums),
       col = "green",
       lwd = 2)
abline(v = quantile(gc_row_sums, 0.10),
       col = "blue",
       lwd = 2)

#Filter unsignificant entries with PCA
gene_counts_filtered <- gene_counts[-which(gc_row_sums==0),]
gene_counts_filtered <- t(gene_counts_filtered)
gene_counts_pca <- prcomp(gene_counts_filtered, center = TRUE, scale = FALSE)

jc_row_sum <- unname(rowSums(junction_counts))
junction_counts_filtered <- junction_counts[-which(jc_row_sum==1),]

#Take random subsample
gene_counts_random <-gene_counts_pca$x[,1:50]
junction_counts_random <- junction_counts[sample(nrow(junction_counts), 2000),]

save(gene_counts_random, junction_counts_filtered, phenotype_labels, file = "Christian.RData")

print("RData generated!")
