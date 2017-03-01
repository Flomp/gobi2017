library(randomForest)
#library(rgl)
library(stats)

### Random-Forest gehoert zu Bagging (Datensaetze mit Bootstrap-Stichproben, weniger Overfitting, Varianzreduzierung, ohne Einteilung Training - Test)
# Formel, abhaengigkeit des Ergebnissen von allen Features
# Greedy-Algo
# am besten RF mehrfach ausfuehren

gene_count_rf <- data.frame(outcome = as.character(phenotype_labels[,1]), data.matrix(gc_pca))
#junction_count_rf <- data.frame(outcome2 = as.character(phenotype_labels[,1]), t(data.matrix(junction_counts_filtered)))
junction_count_rf <- data.frame(outcome2 = as.character(phenotype_labels[,1]), data.matrix(jxt_pca))

####### RF fuer Genes
outcome <- gene_count_rf$outcome
geneModel <- randomForest(outcome~.,data=gene_count_rf, importance = TRUE, ntree=1000, replace = TRUE, do.trace = TRUE, keep.forest=TRUE)

geneModel

summary(geneModel)
# Konfusionsmatrix mit OOB-Daten -> TP, TN, FP, TN
geneModel$obb.times
mean(geneModel$oob.times)

table(gene_count_rf$outcome, geneModel$predicted)

# OOB-Fehlerschaetzung statt Kreuzvalidierung
geneModel$err.rate

#Variablenwichtigkeit
imp <- varImpPlot(geneModel)

######### Junction RF
# Wald Junctions generieren
outcome2 <- junction_count_rf$outcome2
junctionModel <- randomForest(outcome2~.,data=junction_count_rf, importance = TRUE, ntree=1000, replace = TRUE, do.trace = TRUE, keep.forest=TRUE)

junctionModel

summary(junctionModel)

# Konfusionsmatrix mit OOB-Daten -> TP, TN, FP, TN
junctionModel$obb.times
mean(junctionModel$oob.times)

# OOB-Fehlerschaetzung statt Kreuzvalidierung
junctionModel$err.rate

#Variablenwichtigkeit
imp <- varImpPlot(junctionModel)

print("Accuracy Genes:")
#Accuarcy
1-geneModel$err.rate

print("Accuracy Junctions:")
#Accuarcy
1-junctionModel$err.rate

#########Plots fuer Random Forest
plot(1-geneModel$err.rate[,1], type="l", col="orange", main="Random Forrest", ylab="Accuracy", xlab="Anzahl Bäume")
lines(1-junctionModel$err.rate[,1], type="l", col="blue")
legend("bottomright", legend=c("Genes", "Junctions"), col=c("orange", "blue"), lty=1:2, cex=0.8)
