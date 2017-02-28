
gene_count_glm <- data.frame(outcome = as.character(phenotype_labels[,1]), data.matrix(gc_pca))
junction_count_glm <- data.frame(outcome2 = as.character(phenotype_labels[,1]), t(data.matrix(junction_counts_random)))

n_gene = nrow(gene_count_glm)
gene_trainIndex = sample(2:n_gene, size = round(0.7*n_gene), replace=FALSE)

gene_training_set <- gene_count_glm[gene_trainIndex, ] 
gene_test_set <- gene_count_glm[-gene_trainIndex, ] 

n_jxt = nrow(junction_count_glm)
jxt_trainIndex = sample(2:n_jxt, size = round(0.7*n_jxt), replace=FALSE)

jxt_training_set <-  junction_count_glm[jxt_trainIndex,] 
jxt_test_set <- junction_count_glm[-jxt_trainIndex,] 

#logistic regression for more than binary outcome options

gene_train_glm1 <- glm(formula = outcome ~ ., family = binomial(link = "logit"), data = gene_training_set)
summary(gene_train_glm1)

gene_predict <- predict(gene_train_glm1,newdata=gene_test_set, type="response") 
gene_predict <- ifelse(gene_predict > 0.5,1,0)
misClasificError <- mean(gene_predict != gene_count_glm$outcome)
print(paste('Accuracy',1-misClasificError))


jxt_train_glm1 <- glm(formula = outcome2 ~ ., family = binomial(link = "logit"), data = jxt_training_set)
summary(jxt_train_glm1)

jxt_predict <- predict(jxt_train_glm1,newdata=jxt_test_set, type="response") 
jxt_predict <- ifelse(jxt_predict > 0.5,1,0)
misClasificError <- mean(jxt_predict != junction_count_glm$outcome2)
print(paste('Accuracy',1-misClasificError))
