library(h2o)
> h2o.init()
> setwd("/home/icb/lukas.simon/projects/GobiProjekt/")
> 
  > RunPredictions<-function(dat){
    >   load(dat)
    >   
      >   aframe<-data.frame(outcome=as.character(phenotype_labels[,1]),t(data.matrix(gene_counts_rando
                                                                                      > m)))
      >   
        >   train<-sample(nrow(aframe),round(nrow(aframe)*0.6))
        >   test<-setdiff(1:nrow(aframe),train)
        >   
          >   write.table(aframe[test,],file="~/test.tmp.txt",sep="\t",row.names=F,quote=F)
        >   write.table(aframe[train,],file="~/train.tmp.txt",sep="\t",row.names=F,quote=F)
        >   
          >   test<-h2o.importFile(path = normalizePath("~/test.tmp.txt"))
          >   train<-h2o.importFile(path = normalizePath("~/train.tmp.txt"))
          >   
            >   hidden<-c(200,200)
            >   genes.pred<-h2o.deeplearning(,1,training_frame = train,validation_frame = test,
                                             >                                epochs = 200, hidden=hidden)
            >   
              >   ####
              >   aframe<-data.frame(outcome=as.character(phenotype_labels[,1]),t(data.matrix(junction_counts_r
                                                                                              > andom)))
              >   
                >   train<-sample(nrow(aframe),round(nrow(aframe)*0.6))
                >   test<-setdiff(1:nrow(aframe),train)
                >   
                  >   write.table(aframe[test,],file="~/test.tmp.txt",sep="\t",row.names=F,quote=F)
                >   write.table(aframe[train,],file="~/train.tmp.txt",sep="\t",row.names=F,quote=F)
                >   
                  >   test<-h2o.importFile(path = normalizePath("~/test.tmp.txt"))
                  >   train<-h2o.importFile(path = normalizePath("~/train.tmp.txt"))
                  >   
                    >   hidden<-c(200,200)
                    >   junctions.pred<-h2o.deeplearning(,1,training_frame = train,validation_frame = test,
                                                         >                                    epochs = 200, hidden=hidden)
                    >   
                      >   return(list(genes.pred,junctions.pred))
                    > }
  > 
    > dat<-list.files("data",full.names = T)