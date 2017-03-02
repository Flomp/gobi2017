library(rgl)
library(stats)
library(raster)
library(ggplot2)
project_accession <- "SRP057500"
pheno_keyword <- "cancer type"

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
phenotype <- phenotype[!phenotype$reads_downloaded==0,]
phenotype <- phenotype[!is.na(phenotype$auc),]



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

####### Für Binär: Transformation auskommentieren
####### Für 7 Kategorien: ifelse auskommmentieren
#phenotype_labels <- transform(phenotype_labels, char = as.integer(factor(char, unique(char))))-1
phenotype_labels$char <- ifelse(grepl('hc',phenotype_labels$char),phenotype_labels$char <- 0, phenotype_labels$char <- 1)

gene_counts_unlogged <- gene_counts
junction_counts_unlogged <- junction_counts


gc_row_sums <- rowSums(gene_counts)



# Reduce dimension gene_counts
gene_counts_filtered <- gene_counts[!gc_row_sums==0,]

Filter2<-function(dat, threshold){
  gene_counts_filtered_healthy <- dat[,phenotype_labels$char == 0]
  gene_counts_filtered_unhealthy <- dat[,phenotype_labels$char == 1]
  
  gc_healthy_means <- apply(gene_counts_filtered_healthy, 1, mean)
  gc_unhealthy_means <-  apply(gene_counts_filtered_unhealthy, 1, mean)
  
  log_fold <- log(((gc_unhealthy_means)+1)/((gc_healthy_means)+1))
  
  gc_means <- apply(gene_counts_filtered,1,mean) 
  gc_sds <- apply(gene_counts_filtered,1,sd) 
  gc_cv <- gc_sds/gc_means
  
  scatter.smooth(log_fold,gc_cv, lpars = list(col = "blue", lwd = 3, lty = 2))
  good <- which(log_fold>threshold | log_fold< -threshold)
  points(log_fold[good], gc_cv[good], col="red")
  
  return (dat[good,])
}

gene_counts_filtered <- Filter2(gene_counts_filtered,2)

Filter<-function(dat, threshhold){
  gc_means <- apply(dat,1,mean) 
  gc_sds <- apply(dat,1,sd) 
  gc_ok <- gc_means>threshhold
  dat <- dat[gc_ok,]

  gc_means_ok <- gc_means[gc_ok] 
  gc_sds_ok<-gc_sds[gc_ok] 
  gc_cv <- sqrt(gc_sds_ok/gc_means_ok)
  
  gc_afit<-loess(gc_cv~gc_means_ok) 
  gc_resids<-gc_afit$residuals
  plot(density(gc_resids))
  good<-which(gc_resids >= quantile(gc_resids,0.90)) 
  dat <- dat[good,]
  
  #plots
  scatter.smooth(gc_means_ok,gc_cv, lpars = list(col = "blue", lwd = 3, lty = 2))
  points(gc_means_ok[good],gc_cv[good],col="red",pch=19) 
  return (dat)
}

gene_counts <- log(gene_counts+1)
junction_counts <- log(junction_counts+1)

gene_counts_filtered <- Filter(gene_counts_filtered, 1)
junction_counts_filtered <- Filter(junction_counts, 0.1)

gene_counts_unlogged_filtered <- Filter(gene_counts_unlogged,1)
junction_counts_unlogged_filtered <- Filter(junction_counts_unlogged,0.1)


###########Write data for deep net###########
filename <- paste(project_accession, "_python", ".csv", sep="")
write.csv(t(gene_counts_unlogged_filtered), file=paste("gc_", filename, sep=""), row.names = FALSE)
write.csv(t(junction_counts_unlogged_filtered), file=paste("jc_", filename, sep=""), row.names=FALSE)
write.csv(phenotype_labels$char, file=paste("labels_", filename, sep=""), row.names=FALSE)


################################
#####Reduce dimensions with PCA

print("Performing PCA...")


gene_counts_filtered <- t(gene_counts_filtered)
gene_counts_pca <- prcomp(gene_counts_filtered, center = TRUE, scale = FALSE)
#plot (gene_counts_pca, type="l")


#junction_counts_filtered <- junction_counts[rowSums(junction_counts>5)>(0.2*ncol(junction_counts)),]
junction_counts_pca <- t(junction_counts_filtered)
junction_counts_pca <- prcomp(junction_counts_pca, center = TRUE, scale = FALSE)

###benötigte Komponenten für Random Forest
g_counts_pca <- gene_counts_pca$x
j_counts_pca <- junction_counts_pca$x

###Plots für 2 PCA Komponenten für Gene und Junctions
png("PCA_Gene.png")
gc_plot <- ggbiplot(gene_counts_pca, choices = 1:2, obs.scale = 1, var.scale = 1, groups = as.factor(phenotype_labels$char), ellipse = TRUE, 
              circle = FALSE, arrow = 0.0,  labels = NULL,labels.size = 0, var.axes = FALSE)
gc_plot <- gc_plot + labs(color=("Patientengruppen"))
gc_plot <- gc_plot + ggtitle("PCA Component Plot (Gene)")
gc_plot <- gc_plot + theme(legend.direction = 'vertical', legend.position = 'right')
print(gc_plot)
dev.off()

png("PCA_Junctions.png")
jc_plot <- ggbiplot(junction_counts_pca, choices = 1:2, obs.scale = 1, var.scale = 1, groups = as.factor(phenotype_labels$char), ellipse = TRUE, 
                    circle = FALSE, arrow = 0.0,  labels = NULL,labels.size = 0, var.axes = FALSE)
jc_plot <- jc_plot + labs(color=("Patientengruppen"))
jc_plot <- jc_plot + ggtitle("PCA Component Plot (Junctions)")
jc_plot <- jc_plot + theme(legend.direction = 'vertical', legend.position = 'right')
print(jc_plot)
dev.off()
