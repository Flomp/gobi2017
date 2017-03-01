library(rgl)
library(stats)
library(raster)
library(openxlsx)
library(ggplot2)
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
phenotype_labels <- transform(phenotype_labels, char = as.integer(factor(char, unique(char))))-1

gene_counts_unlogged <- gene_counts
junction_counts_unlogged <- junction_counts

gene_counts <- log(gene_counts+1)
gc_row_sums <- rowSums(gene_counts)
junction_counts <- log(junction_counts+1)
#jc_row_sum <- unname(rowSums(junction_counts))

# Reduce dimension gene_counts
gene_counts_filtered <- data.matrix(gene_counts[!gc_row_sums==0,])

gc_means <- apply(gene_counts_filtered,1,mean) 
gc_sds <- apply(gene_counts_filtered,1,sd) 


gc_ok <- which(gc_means > 1) 
gc_means_ok <- gc_means[gc_ok] 
gc_sds_ok<-gc_sds[gc_ok] 
gc_cv <- gc_sds_ok/gc_means_ok
plot(gc_means_ok, gc_cv) 
plot(gc_means_ok,gc_sds_ok/gc_means_ok) 
plot(gc_means_ok,sqrt(gc_sds_ok/gc_means_ok)) 
gc_afit<-loess(sqrt(gc_sds_ok/gc_means_ok)~gc_means_ok) 
gc_resids<-gc_afit$residuals
plot(density(gc_resids)) 
good<-which(gc_resids > 0.1) 

#plots
scatter.smooth(gc_means_ok,sqrt(gc_sds_ok/gc_means_ok), lpars =
                 list(col = "blue", lwd = 3, lty = 2))
points(gc_means_ok[good],sqrt(gc_sds_ok[good]/gc_means_ok[good]),col="red",pch=19) 

##unlogged
#gene_counts_filtered_ul <- gene_counts_unlogged[!gc_row_sums==0,]
#gene_counts_filtered_stat_ul <- gene_counts_filtered_ul
#gene_counts_filtered_stat_ul$mean <- apply(gene_counts_filtered_ul, 1, mean)
#gene_counts_filtered_stat_ul$sd <- apply(gene_counts_filtered_ul, 1, sd)
#write.xlsx(gene_counts_filtered_stat_ul, "mydata_ul.xlsx")
####calculate sd/mean in excel, cause of memory
#gene_counts_filtered_stat2 <- read.csv2("~/LMU/Binf/gobi/Blockteil/mydata_gc.csv", header= TRUE)

#gene_counts_filtered_stat2 <- read.csv2("~/LMU/Binf/gobi/Blockteil/mydata.csv", header= TRUE)
#plot(x = gene_counts_filtered_stat2$mean  ,y = gene_counts_filtered_stat2$cv ,type = "p")

##################################
# Reduce dimension junction_counts
junction_counts_filtered_stat <- junction_counts
junction_counts_filtered_stat$mean <- apply(junction_counts, 1, mean)
junction_counts_filtered_stat$sd <- apply(junction_counts, 1, sd)
#gene_counts_filtered_stat$cv <- apply(junction_counts, 1, function(x) junction_counts_filtered_stat$mean/junction_counts_filtered_stat$sd)
write.xlsx(junction_counts_filtered_stat, "mydata_jc.xlsx")
####calculate sd/mean in excel, cause of memory
junction_counts_filtered_stat2 <- read.csv2("~/LMU/Binf/gobi/Blockteil/mydata_jc.csv", header= TRUE)

plot(x = junction_counts_filtered_stat2$mean  ,y = junction_counts_filtered_stat2$cv ,type = "p")


p4 <- qplot(data = junction_counts_filtered_stat2, mean, cv, xlab = "", ylab = "",
       geom_smooth(method = "auto", size = 1.5), theme_bw())

p4


print("Performing PCA...")

#Reduce dimensions with PCA
gene_counts_filtered <- t(gene_counts_filtered)
gene_counts_pca <- prcomp(gene_counts_filtered, center = TRUE, scale = FALSE)
plot (gene_counts_pca, type="l")


#junction_counts_filtered <- junction_counts[rowSums(junction_counts>5)>(0.2*ncol(junction_counts)),]
junction_counts_pca <- t(junction_counts)
junction_counts_pca <- prcomp(junction_counts_pca, center = TRUE, scale = FALSE)


#Take random subsample

#gc_pca <- gene_counts_pca$x
#save(gc_pca, junction_counts_filtered, phenotype_labels, file = paste(project_accession, ".RData", sep = ""))

#print("RData generated!")

filename <- paste(project_accession, "_python", ".csv", sep="")
write.csv(gene_counts_pca$x, file=paste("gc_", filename, sep=""), row.names = FALSE)
write.csv(junction_counts_pca$x, file=paste("jc_", filename, sep=""), row.names=FALSE)
write.csv(phenotype_labels$char, file=paste("labels_", filename, sep=""), row.names=FALSE)


temp <- data.frame(scale(gene_counts_pca$x[,1:3]))
plot3d(temp$PC1, temp$PC2, temp$PC3, col=phenotype_labels$char+1)
gc_cluster <- kmeans(gene_counts_pca$x[,1:3], length(unique(phenotype_labels)[,1]), nstart=20)
table(gc_cluster$cluster, phenotype_labels[,1])
