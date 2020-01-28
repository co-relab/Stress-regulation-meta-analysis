
#Impostazione script
#devi prima creare tutte le variabili che non dipendono da altre
#o che sono "genitori" di altre variabili (come category)
#a questo punto crei un dataset<-con CBIND unisci tutte le colonne che hai creato
#dal dataset (di nome PINCOPALLO) applichi gli IFELSE per creare i "figli"
#delle variabili genitori... ad esempio PINCOPALLO$category===1 e via dicendo
#negli ifelse
#ATTENZIONE 
#puoi fare gli ifelse con il sample dentro è più preciso 
#e riproducibile
Type<-ifelse(category$category==1, 
             sample(c(1,2,3,4,5),length(category$category==1),replace = T), NA)

#DA CATEGORY ATTACCO  4 VARIABILI che cambiano per livello di category
set.seed(234)
category <-sample(c(1,2,3,4),200, replace=T)
category<-as.data.frame(category)
#creo la variabile figlia
Type_of_Sam<-ifelse(category$category==1, 
       c(1,2,3,4,5), NA)
Type_of_Sam
#attacco al dataset la variabile figlia
category$Type_of_Sam<-Type_of_Sam

type_of_environmet <-ifelse(category$category==2, 
                            c(1,2), NA)

category$type_of_environmet <-type_of_environmet 

type_of_exposure <-ifelse(category$category==2, 
                          c(1,2,3), NA)

category$type_of_exposure<-type_of_exposure
  

type_of_SocialSupport <-ifelse(category$category==4, 
                          c(1,2,3,4), NA) 
category$type_of_SocialSupport<-type_of_SocialSupport

head(category)
