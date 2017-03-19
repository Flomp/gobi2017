library(rgl)
library(stats)
library(raster)
library(ggplot2)
library(devtools)
library(ggbiplot)
library(randomForest)

project_accession <- "SRP057500"
pheno_keyword <- "cancer type"

# Download and read ----
if(!dir.exists(project_accession)){
  dir.create(project_accession)
}
setwd(project_accession)
url<-"http://duffel.rail.bio/recount"

# Gene Counts
filename <- "counts_gene.tsv.gz"
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
                quiet = FALSE, mode = "w", cacheOK = TRUE)
}

gene_counts <- read.table(file =filename, sep = '\t', header = TRUE)

# Junction Counts 
filename <- "counts_jx.tsv.gz"
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
                quiet = FALSE, mode = "w", cacheOK = TRUE)
}
junction_counts <- read.table(file = filename, sep = '\t', header = TRUE)

# Phenotype
filename <- paste(project_accession, "tsv", sep=".")
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
                quiet = FALSE, mode = "w", cacheOK = TRUE)
}
phenotype <- read.table(file = filename, sep = '\t', header = TRUE)
phenotype <- phenotype[!phenotype$reads_downloaded==0,]
phenotype <- phenotype[!is.na(phenotype$auc),]

# Junction Positions
filename <- paste(project_accession, "junction_id_with_transcripts.bed.gz", sep=".")
if(!file.exists(filename)){
  download.file(paste(url, project_accession, filename, sep="/"), filename, "wget", 
                quiet = FALSE, mode = "w", cacheOK = TRUE)
}
junction_positions <-read.table(file = filename, sep = '\t', header = FALSE)

stopifnot(ncol(gene_counts)==ncol(junction_counts) && ncol(gene_counts) == nrow(phenotype)&&
            nrow(junction_counts)==nrow(junction_positions))

print("Download successfully completed!")

# Beautify tables ----

concat_row_names = paste (junction_positions[,1], junction_positions[,2], junction_positions[,3], sep="_")
rownames(junction_counts) <- concat_row_names

phenotype_labels <- data.frame(char=phenotype$characteristics)
rownames(phenotype_labels) <- phenotype$run
phenotype_labels$char <- tolower(phenotype_labels$char)
phenotype_labels$char <- gsub(paste(".*", pheno_keyword, ":", sep=""), "", phenotype_labels$char)
phenotype_labels$char <- gsub("(,|\\)$).*", "", phenotype_labels$char)

####### for binary classification: comment line with transform-function out (0: healthy, 1: cancer)
####### for 7 categories: comment line with ifelse-function out 
#phenotype_labels <- transform(phenotype_labels, char = as.integer(factor(char, unique(char))))-1
phenotype_labels$char <- ifelse(grepl('hc',phenotype_labels$char),phenotype_labels$char <- 0, phenotype_labels$char <- 1)

# Filter counts ----
FilterByLogFold<-function(dat, threshold){
  filtered <- dat[rowSums(dat)!=0,]
  filtered_healthy <- filtered[,phenotype_labels$char == 0]
  filtered_unhealthy <- filtered[,phenotype_labels$char == 1]
  
  healthy_means <- apply(filtered_healthy, 1, mean)
  unhealthy_means <-  apply(filtered_unhealthy, 1, mean)
  
  log_fold <- log(((unhealthy_means)+1)/((healthy_means)+1))
  
  means <- apply(filtered,1,mean) 
  sds <- apply(filtered,1,sd) 
  cv <- sds/means
  
  plot(log_fold,cv,xlab="Log Fold", ylab="Variationskoeffizient",main="Log Fold\nGene-Counts", cex=0.8)
  good <- which(log_fold>threshold | log_fold< -threshold)
  points(log_fold[good], cv[good], col="red", pch=19, cex=0.8)
  abline(v=threshold, lwd=2, lty=2)
  abline(v=-threshold, lwd=2, lty=2)
  
  return (filtered[good,])
}

FilterByVariance<-function(dat, threshold){
  filtered <- dat[rowSums(dat)!=0,]
  means <- apply(filtered,1,mean)
  ok <- means>threshold
  filtered <- filtered[ok,]

  means <- means[ok] 
  sds <- apply(filtered,1,sd)
  cv <- sqrt(sds/means)
  
  afit<-loess(cv~means) 
  resids<-afit$residuals
  plot(density(resids))
  good<-which(resids >= quantile(resids,0.90)) 
  
  #plots
  plot(cv~means, cex=0.8, xlab="Mittelwert", ylab = "Variationskoeffizient", main="Variationskoeffizient\nGene-Counts")
  lines(cbind(sort(means), predict(afit, newdata = sort(means))), 
        col="blue",lwd=3)
  points(means[good],cv[good],col="red",pch=19, cex=0.8) 
  return(filtered[good,])
}


FilterByExpression<-function(dat, threshold, percentage){
  filtered <- dat[rowSums(dat>threshold)>(percentage*ncol(dat)),]
  return(filtered)
}

FilterRandom<-function(dat, samplesize){
  filtered <- dat[sample(nrow(dat), samplesize), ]
}
#Parameters for the different variations of choosing good Data
print("Filtering features...")

#1. Try: Random sampling
gene_counts_filtered <- FilterRandom(gene_counts,2000)
junction_counts_filtered <- FilterRandom(junction_counts,2000)

#2.Try: Feature Selection by Expression 
#gene_counts_filtered <- FilterRandom(gene_counts,2000)
#junction_counts_filtered <- FilterByExpression(junction_counts, 5, 0.2)

#3.Try: Variance
#not done

#4.Try: Coefficient of variation vs. Mean
#gene_counts_filtered <- FilterByVariance(log(gene_counts+1),1)
#junction_counts_filtered <- FilterByVariance(log(junction_counts+1), 0.1)

#5.Try: LogFold
#gene_counts_filtered <- FilterByLogFold(gene_counts,2)
#junction_counts_filtered <- FilterByLogFold(junction_counts, 2)

# Write data for deep net ----
filename <- paste(project_accession, "_python", ".csv", sep="")
write.csv(t(gene_counts_filtered), file=paste("gc_", filename, sep=""), row.names = FALSE)
write.csv(t(junction_counts_filtered), file=paste("jc_", filename, sep=""), row.names=FALSE)
write.csv(phenotype_labels$char, file=paste("labels_", filename, sep=""), row.names=FALSE)


# Reduce dimensions with PCA ----
print("Performing PCA...")

gene_counts_filtered <- log(gene_counts_filtered+1)
junction_counts_filtered <- log(junction_counts_filtered+1)

gene_counts_filtered <- t(gene_counts_filtered)
gene_counts_pca <- prcomp(gene_counts_filtered, center = TRUE, scale = FALSE)

junction_counts_filtered <- t(junction_counts_filtered)
junction_counts_pca <- prcomp(junction_counts_filtered, center = TRUE, scale = FALSE)

eigen_gc <- gene_counts_pca$x
eigen_jc <- junction_counts_pca$x

plot(gene_counts_pca, type='l')
plot(junction_counts_pca, type ='l')
summary(junction_counts_pca)

# PCA plots ----
png("PCA_Gene.png")
gc_plot <- ggbiplot(gene_counts_pca, choices = 1:2, obs.scale = 1, var.scale = 1, groups = as.factor(phenotype_labels$char), ellipse = TRUE, 
                    circle = FALSE, arrow = 0.0,  labels = NULL,labels.size = 0, var.axes = FALSE)
gc_plot <- gc_plot + labs(color=("Patientengruppen"))
gc_plot <- gc_plot + ggtitle("PCA Component Plot (Gene)")
gc_plot <- gc_plot + theme(legend.direction = 'vertical', legend.position = 'right')
gc_plot <- gc_plot + theme(aspect.ratio=1.0)
print(gc_plot)
dev.off()

png("PCA_Junctions.png")
jc_plot <- ggbiplot(junction_counts_pca, choices = 1:2, obs.scale = 1, var.scale = 1, groups = as.factor(phenotype_labels$char), ellipse = TRUE, 
                    circle = FALSE, arrow = 0.0,  labels = NULL,labels.size = 0, var.axes = FALSE)
jc_plot <- jc_plot + labs(color=("Patientengruppen"))
jc_plot <- jc_plot + ggtitle("PCA Component Plot (Junctions)")
jc_plot <- jc_plot + theme(legend.direction = 'vertical', legend.position = 'right')
jc_plot <- jc_plot + theme(aspect.ratio=1.0)
print(jc_plot) 
dev.off()


### Random-Forest for Gene- & Junction Counts
gene_count_rf <- data.frame(outcome = as.character(phenotype_labels[,1]), data.matrix(eigen_gc))
junction_count_rf <- data.frame(outcome2 = as.character(phenotype_labels[,1]), data.matrix(eigen_jc))

# RF for Genes
outcome <- gene_count_rf$outcome
geneModel <- randomForest(outcome~.,data=gene_count_rf, importance = TRUE, ntree=1000, replace = TRUE, do.trace = TRUE, keep.forest=TRUE)
summary(geneModel)

# RF for Junctions
outcome2 <- junction_count_rf$outcome2
junctionModel <- randomForest(outcome2~.,data=junction_count_rf, importance = TRUE, ntree=1000, replace = TRUE, do.trace = TRUE, keep.forest=TRUE)
summary(junctionModel)

# Confusionmatrix and Error-Rate for all Phenotypes
print("Gen-Modell")
print(geneModel)
plot (geneModel)

print("Junction-Modell")
print(junctionModel)
plot(junctionModel)

#Plot for Random Forest
png(paste(project_accession, "_RF.png", sep = ""))
plot(1-geneModel$err.rate[,1], type="l", col="orange", main="Random Forrest", ylab="Accuracy", xlab="Anzahl Bäume", ylim=ylim=c(0.65,0.90), sub= "\n  Datenauswahl-Kriterium", cex.sub = 1, font.sub = 2)
lines(1-junctionModel$err.rate[,1], type="l", col="blue")
legend("bottomright", legend=c("Genes", "Junctions"), col=c("orange", "blue"), lty=1:1, cex=0.9, title="Input Datenart", bg='aliceblue')
dev.off()

print("Fi")
