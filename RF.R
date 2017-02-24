library(randomForest)

### Random-Forest gehoert zu Bagging (Datensaetze mit Bootstrap-Stichproben, weniger Overfitting, Varianzreduzierung, ohne Einteilung Training - Test)
# Formel, abhaengigkeit des Ergebnissen von allen Features
# Greedy-Algo
# am besten RF mehrfach ausfuehren

gene_count_rf <- data.frame(outcome = as.character(phenotype_labels[,1]), data.matrix(gene_counts_random))

junction_count_rf <- data.frame(outcome = as.character(phenotype_labels[,1]), t(data.matrix(junction_counts_random)))
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

#geneModel$err.rate[100,1]
plot(geneModel$err.rate[,1], type="l")


#Vorhersage Datensatz
vorhersage <- predict(geneModel, newdata=gene_count_rf)
vorhersage_geneModel<- table(gene_count_rf$outcome,vorhersage)

(vorhersage_geneModel[1,1]+vorhersage_geneModel[2,2])/length(gene_count_rf$outcome) # Accuracy

#Variablenwichtigkeit
imp <- varImpPlot(geneModel)


# Wald Junctions generieren
outcome2 <- junction_count_rf$outcome
junctionModel <- randomForest(outcome2~.,data=junction_count_rf, importance = TRUE, ntree=1000, replace = TRUE, do.trace = TRUE, keep.forest=TRUE)

junctionModel

summary(junctionModel)

# Konfusionsmatrix mit OOB-Daten -> TP, TN, FP, TN
junctionModel$obb.times
mean(junctionModel$oob.times)

table(junction_count_rf$outcome, junctionModel$predicted)

# OOB-Fehlerschaetzung statt Kreuzvalidierung
junction_count_rf$err.rate

#geneModel$err.rate[100,1]
plot(junction_count_rf$err.rate[,1], type="l")


#Vorhersage Datensatz
vorhersage_junction <- predict(junctionModel, newdata=junction_count_rf)
vorhersage_junctionModel<- table(junction_count_rf$outcome,vorhersage_junction)

(vorhersage_junctionModel[1,1]+vorhersage_junctionModel[2,2])/length(junction_count_rf$outcome) # Accuracy

#Variablenwichtigkeit
imp <- varImpPlot(junctionModel)