library(gbm)
library(dplyr)

gene_count_gbm <- data.frame(outcome = as.numeric(phenotype_labels[,1]), data.matrix(gc_pca))
junction_count_gbm <- data.frame(outcome2 = as.numeric(phenotype_labels[,1]), t(data.matrix(junction_counts_random)))

## sample
n_gene = nrow(gene_count_gbm)
gene_trainIndex = sample(2:n_gene, size = round(0.7*n_gene), replace=FALSE)

gbm_gene_training_set <- gene_count_gbm[gene_trainIndex, ] 
gbm_gene_test_set <- gene_count_gbm[-gene_trainIndex, ] 

summary(gbm_gene_training_set)


model = gbm.fit(
  x = gene_count_gbm, y = gene_count_gbm$outcome, distribution = "bernoulli", n.trees = 5000,
  shrinkage = 0.01, interaction.depth = 3, n.minobsinnode = 10, nTrain = round (nrow(gene_count_gbm * 0.8)), 
  verbose = TRUE )





n_jxt = nrow(junction_count_gbm)
gbm_jxt_trainIndex = sample(2:n_jxt, size = round(0.7*n_jxt), replace=FALSE)

gbm_jxt_training_set <-  junction_count_gbm[jxt_trainIndex,] 
gbm_jxt_test_set <- junction_count_gbm[-jxt_trainIndex,] 


#gene_train_gbm1 <- gbm(outcome~.,data=gbm_gene_training_set, importance = TRUE, ntree=1000, replace = TRUE, do.trace = TRUE, keep.forest=TRUE))
summary(gene_train_gbm1)
plot.gbm(gene_train_gbm1)