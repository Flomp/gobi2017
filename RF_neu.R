library(randomForest)
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

summary(geneModel)

mean(geneModel$oob.times)

# OOB-Fehlerschaetzung statt Kreuzvalidierung
geneModel$err.rate

#Variablenwichtigkeit
imp <- varImpPlot(geneModel)

######### Junction RF
# Wald Junctions generieren
outcome2 <- junction_count_rf$outcome2
junctionModel <- randomForest(outcome2~.,data=junction_count_rf, importance = TRUE, ntree=1000, replace = TRUE, do.trace = TRUE, keep.forest=TRUE)

summary(junctionModel)

mean(junctionModel$oob.times)

# OOB-Fehlerschaetzung statt Kreuzvalidierung
junctionModel$err.rate

#Variablenwichtigkeit
imp <- varImpPlot(junctionModel)

####Konfussionsmatrize aller Phenotypne
print("Gen-Modell")
print(geneModel)

print("Junction-Modell")
print(junctionModel)

######### Accuracies
print("Accuracy Genes:")
#Accuarcy Genes
print(1-geneModel$err.rate[1000,])

print("Accuracy Junctions:")
#Accuarcy Junctions
print(1-junctionModel$err.rate[1000,])

#########Plots fuer Random Forest
png(paste(project_accession, ".png", sep = ""))
plot(1-geneModel$err.rate[,1], type="l", col="orange", main="Random Forrest", ylab="Accuracy", xlab="Anzahl Bäume", ylim=range(1-geneModel$err.rate[,1],1-junctionModel$err.rate[,1]))
lines(1-junctionModel$err.rate[,1], type="l", col="blue")
legend("bottomright", legend=c("Genes", "Junctions"), col=c("orange", "blue"), lty=1:1, cex=0.9, title="Input Datenart", bg='aliceblue')
dev.off()
 