library("h2o")

project_accession <- "SRP032775"
pheno_keyword <- "infectious agent"

##############################
##########Download############
##############################
dir.create(project_accession)
setwd(project_accession)
url<-"http://duffel.rail.bio/recount"
filename <- "counts_gene.tsv.gz"
###Gene Counts
download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
              quiet = FALSE, mode = "w", cacheOK = TRUE)
gene_counts <- read.table(file =filename, sep = '\t', header = TRUE)
gene_counts <- gene_counts[-which(rowSums(gene_counts)== 0),]

### Junction Counts
filename <- "counts_jx.tsv.gz"
download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
              quiet = FALSE, mode = "w", cacheOK = TRUE)
junction_counts <- read.table(file = filename, sep = '\t', header = TRUE)

### Phenotype
filename <- paste(project_accession, "tsv", sep=".")
download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
              quiet = FALSE, mode = "w", cacheOK = TRUE)
phenotype <- read.table(file = filename, sep = '\t', header = TRUE)
phenotype <- phenotype[-which(phenotype$reads_downloaded == 0),]


### Junction Positions
filename <- paste(project_accession, "junction_id_with_transcripts.bed.gz", sep=".")
download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
              quiet = FALSE, mode = "w", cacheOK = TRUE)
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
phenotype_labels <- transform(phenotype_labels, char = as.integer(factor(char, unique(char))))

gene_counts <- log10(gene_counts+1)
gc_row_sums <- rowSums(gene_counts)


list_histo <- hist(gc_row_sums, breaks=seq(from=0, to=max(gc_row_sums)+20, by=1))
abline(v = median(gc_row_sums),
       col = "red",
       lwd = 2)
abline(v = mean(gc_row_sums),
       col = "green",
       lwd = 2)
abline(v = quantile(gc_row_sums, 0.10),
       col = "blue",
       lwd = 2)
abline(v = quantile(gc_row_sums, 0.95),
       col = "blue",
       lwd = 2)


#Filter unsignificant entries with PCA
gene_counts_filtered <- gene_counts[-which(gc_row_sums==0),]
gene_counts_filtered <- t(gene_counts_filtered)
gene_counts_pca <- prcomp(gene_counts_filtered, center = TRUE, scale = FALSE)

jc_row_sum <- unname(rowSums(junction_counts))
junction_counts <- junction_counts[-which(jc_row_sum==1),]

#Take randon subsample
gene_counts <-gene_counts[sample(nrow(gene_counts), 2000),]
junction_counts <- junction_counts[sample(nrow(junction_counts), 2000),]

save(gene_counts, junction_counts, phenotype_labels, file = "Rebecca.RData")

print("RData generated!")


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

model <- h2o.deeplearning(x = phenotype_labels.hex,  #A vector containing the names of the predictors in the model
                          y = phenotype_labels.hex,   # column number for label
                          data = gene_counts.hex, # data in H2O format
                          activation = "TanhWithDropout", # or 'Tanh'
                          input_dropout_ratio = 0.2, # % of inputs dropout
                          hidden_dropout_ratios = c(0.5,0.5,0.5), # % for nodes dropout
                          balance_classes = TRUE, 
                          hidden = c(50,50,50), # three layers of 50 nodes
                          epochs = 100) # max. no. of epochs



RunPredictions <- function(dat){
  load(dat)
  
  #aframe<-data.frame(outcome=as.character(phenotype_labels[,1]),t(data.matrix(gene_counts_random)))
  #train<-sample(nrow(aframe),round(nrow(aframe)*0.6))
  #test<-setdiff(1:nrow(aframe),train)
  #write.table(aframe[test,],file="~/test.tmp.txt",sep="\t",row.names=F,quote=F)
  #write.table(aframe[train,],file="~/train.tmp.txt",sep="\t",row.names=F,quote=F)
          
  test <- gene_counts_hex
  train <- gene_counts_hex
  #test<-h2o.importFile(path = normalizePath("~/test.tmp.txt"))
  #train<-h2o.importFile(path = normalizePath("~/train.tmp.txt"))
  
  hidden <- c(200,200)
  genes.pred <- h2o.deeplearning(,1,training_frame = train,validation_frame = test,epochs = 200, hidden=hidden)
  
  ####
  #aframe<-data.frame(outcome=as.character(phenotype_labels[,1]),t(data.matrix(junction_counts_random)))
        
  #train<-sample(nrow(aframe),round(nrow(aframe)*0.6))
  #test<-setdiff(1:nrow(aframe),train)
            
  #write.table(aframe[test,],file="~/test.tmp.txt",sep="\t",row.names=F,quote=F)
  #write.table(aframe[train,],file="~/train.tmp.txt",sep="\t",row.names=F,quote=F)
 
  test <- junction_counts_hex
  train <- junction_counts_hex
  #test<-h2o.importFile(path = normalizePath("~/test.tmp.txt"))
  #train<-h2o.importFile(path = normalizePath("~/train.tmp.txt"))
    
  hidden<-c(200,200)
  junctions.pred<-h2o.deeplearning(,1,training_frame = train,validation_frame = test, epochs = 200, hidden=hidden)
  return(list(genes.pred,junctions.pred))
}

  dat<-list.files("data",full.names = T)