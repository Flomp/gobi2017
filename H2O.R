project_accession <- "SRP032775"
pheno_keyword <- "infectious agent"

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


#Reduce dimensions with PCA
gene_counts_filtered <- gene_counts[!gc_row_sums==0,]
gene_counts_filtered <- t(gene_counts_filtered)
gene_counts_pca <- prcomp(gene_counts_filtered, center = TRUE, scale = FALSE)
plot (gene_counts_pca, type="l")

#jc_row_sum <- unname(rowSums(junction_counts))
junction_counts_filtered <- junction_counts[rowSums(junction_counts>5)>(0.5*ncol(junction_counts)),]

#Take random subsample
junction_counts_random <- junction_counts[sample(nrow(junction_counts), 5000),]
junction_counts_random <- junction_counts_random[rowSums(junction_counts_random)!=0,]
gc_pca <-gene_counts_pca$x
save(gc_pca, junction_counts_filtered, phenotype_labels, file = paste(project_accession, ".RData", sep = ""))

print("RData generated!")



library(h2o)

## Connect to H2OCluster
h2o.init(nthreads = -1)

# make h2o objects from data.frames
gene_counts.hex <- as.h2o(gene_counts, destination_frame="gene_counts.hex")
junction_counts.hex <- as.h2o(junction_counts, destination_frame="junction_counts.hex")
phenotype_labels.hex <- as.h2o(phenotype_labels, destination_frame="phenotype_labels.hex")


## Summary
h2o.describe(gene_counts_hex)
h2o.describe(junction_counts_hex)

## plot 
h2o.hist()

## Deep-Learning
###RunPredictions <- function(dat){
# load(dat)

#aframe<-data.frame(outcome=as.character(phenotype_labels[,1]),t(data.matrix(gene_counts_random)))
#train<-sample(nrow(aframe),round(nrow(aframe)*0.6))
#test<-setdiff(1:nrow(aframe),train)
#write.table(aframe[test,],file="~/test.tmp.txt",sep="\t",row.names=F,quote=F)
#write.table(aframe[train,],file="~/train.tmp.txt",sep="\t",row.names=F,quote=F)

# test <- gene_counts_hex
# train <- gene_counts_hex
#test<-h2o.importFile(path = normalizePath("~/test.tmp.txt"))
#train<-h2o.importFile(path = normalizePath("~/train.tmp.txt"))

# hidden <- c(200,200)
# genes.pred <- h2o.deeplearning(,1,training_frame = train,validation_frame = test,epochs = 200, hidden=hidden)

####
#aframe<-data.frame(outcome=as.character(phenotype_labels[,1]),t(data.matrix(junction_counts_random)))

#train<-sample(nrow(aframe),round(nrow(aframe)*0.6))
#test<-setdiff(1:nrow(aframe),train)

#write.table(aframe[test,],file="~/test.tmp.txt",sep="\t",row.names=F,quote=F)
#write.table(aframe[train,],file="~/train.tmp.txt",sep="\t",row.names=F,quote=F)

# test <- junction_counts_hex
# train <- junction_counts_hex
#test<-h2o.importFile(path = normalizePath("~/test.tmp.txt"))
#train<-h2o.importFile(path = normalizePath("~/train.tmp.txt"))

# hidden<-c(200,200)
# junctions.pred<-h2o.deeplearning(,1,training_frame = train,validation_frame = test, epochs = 200, hidden=hidden)
# return(list(genes.pred,junctions.pred))
#}

# dat<-list.files("data",full.names = T)


## GBM
## split data_set

gene_all <- data.frame(gene_outcome = as.character(phenotype_labels[,1]), data.matrix(gc_pca))
gene_all.hex <- as.h2o(gene_all, destination_frame="gene_all.hex")
gene_independent <- gene_all.hex[,-1]
gene_dependent <- gene_all.hex[,1]

junction_all <- data.frame(outcome2 = as.character(phenotype_labels[,1]), t(data.matrix(junction_counts_random)))
junction_independent <- junction_all[,-1]
junction_dependent <- data.frame(junction_all[,1])


gene_gbm1 <- h2o.gbm(y = gene_dependent, x = gene_independent, data = gene_all,
        n.trees = 10, interaction.depth = 3,
        n.minobsinnode = 2, shrinkage = 0.2, distribution= "gaussian")

junction_gbm1 <- h2o.gbm(y = dependent, x = independent, data = .hex,
                     n.trees = 10, interaction.depth = 3,
                     n.minobsinnode = 2, shrinkage = 0.2, distribution= "gaussian")

## GLM

#h2o.glm(y = "CAPSULE", x = c("AGE","RACE","PSA","DCAPS"), data =
#          prostate.hex, family = "binomial", nfolds = 10, alpha = 0.5)
#myX = setdiff(colnames(prostate.hex), c("ID", "DPROS", "DCAPS", "VOL"))
#h2o.glm(y = "VOL", x = myX, data = prostate.hex, family = "gaussian", nfolds = 5, alpha = 0.1)
