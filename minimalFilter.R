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
  abline(v=2, lwd=2, lty=2)
  abline(v=-2, lwd=2, lty=2)
  
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

gene_counts_filtered <- FilterByVariance(log(gene_counts+1),1)
junction_counts_filtered <- FilterByVariance(log(junction_counts)+1, 0.1)

gene_counts_filtered <- FilterByLogFold(gene_counts,2)
junction_counts_filtered <- FilterByLogFold(junction_counts, 2)