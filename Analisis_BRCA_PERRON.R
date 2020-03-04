library(TCGAbiolinks)
library(Biobase)
library(SummarizedExperiment)

# Query miRNA -------------------------------------------------------------

query_mirna <- GDCquery(project = "TCGA-BRCA",
                              data.category = "Transcriptome Profiling",
                              data.type = "miRNA Expression Quantification",
                              sample.type = "Primary solid Tumor")

GDCdownload(query_mirna)
data_mirna <- GDCprepare(query_mirna, save = TRUE)

# Data frame --------------------------------------------------------------

####1: OBTENER PACIENTES QUE NO TIENEN MUTACIONES
uncommon_symbol <- read.csv("uncommon_symbol.csv")
correccion <- gsub(".", "-", colnames(uncommon_symbol), fixed = TRUE)
colnames(uncommon_symbol) <- correccion

datoschidos <- data_mirna
nombresmirna <- datoschidos$miRNA_ID

########1: QUEDARTE SOLO CON LOS READ COUNTS
mirnaPERRRON <- datoschidos[,c(grep("read_count",colnames(datoschidos)))]
mirnaPERRRON$miRNA<- nombresmirna

      #---Automatizando
mirnaPERRRON <- mirnaPERRRON[,c(length(colnames(mirnaPERRRON)),
                            1:length(colnames(mirnaPERRRON))-1)]

#########2: QUEDARTE SOLO CON LOS 12 DIGITOS DE TCGA
expresion <- gsub("read_count_", "",colnames(mirnaPERRRON))
expresion <- substr(expresion, 1, 12)

colnames(mirnaPERRRON) <-expresion

rownames(mirnaPERRRON) <- mirnaPERRRON$miRNA
mirnaPERRRON <- mirnaPERRRON[,-c(1)]


#########3: QUE PACIENTES DE DATASET DE MIRNAS TIENEN MUTACIONES
#Y CUALES NO usando uncommon

uncommon_filtrados <- mirnaPERRRON[,colnames(mirnaPERRRON)%in% colnames(uncommon_symbol)]
uncommon_filtrados

as.data.frame(uncommon_filtrados)
#########NOTA: HACER UN DATAFRAME PARA CADA GRUPO

#ID de uncommon y no uncommon (mutacione y sin mutaciones)
####UNCOMMON NO TIENEN MUTACIONES

id <- colnames(uncommon_symbol)
id <- id[-c(1,2,length(id))] #no mutación

id_mutation <- read.csv("mutacion.csv")
id_mutation <- as.character(id_mutation[,2])

################---Haciendo data.frame con mirnaPERRON
    #--- Sin mutaciones
matrix_nomutation <- mirnaPERRRON[,id] #no ideal pero sí salío lol

    #--- Con mutaciones
matrix_mutation <- mirnaPERRRON[,colnames(mirnaPERRRON) %in% id_mutation] #index

# Creando Control ---------------------------------------------------------

    ###---CREANDO CONTROL
query_control <- GDCquery(project = "TCGA-BRCA",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification",
                        sample.type = "Solid Tissue Normal")

GDCdownload(query_control)
data_control <- GDCprepare(query_control, save = TRUE)

rownames(data_control) <- data_control$miRNA_ID
   
   #ReadCounts
control <- data_control[,c(grep("read_count",colnames(data_control)))]


vacio_control <- gsub("read_count_", "",colnames(control))
vacio_control <- substr(vacio_control, 1, 12)

colnames(control) <-vacio_control
rownames(control) <- data_control$miRNA_ID


# Expresión diferencial ---------------------------------------------------

library(DESeq2)

grupoa <- "normal"
grupob <- "tumor"

numero_nomutation <- ncol(matrix_nomutation)
numero_mutation <- ncol(matrix_mutation)
numero_control <- ncol(control)

matrix_tumor <- cbind(matrix_nomutation,matrix_mutation) #Juntando todos los tumores
countdata <- matrix_tumor #Convirtiendo la matriz en la variable countdata
countdata <- as.matrix(countdata) #Regresando a matrix

(condition <- factor(c( rep(c(grupoa,grupob),
                            c(numero_nomutation,numero_mutation)))))

(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)

dds

dds$condition <- relevel(dds$condition, ref=grupoa)
dds <- DESeq(dds)

# Obtener los datos de Differential Expression
res <- results(dds)
table(res$padj<0.05)
## Ordenar por adjusted p-value
res <- res[order(res$padj), ]
## Unir con los datos de counts
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
head(resdata)
## Escribir resultados
write.csv(resdata, file="diffexpr-controlVStumor.csv")

                  ####NO MUTACION VS MUTACION (1)

grupomutation <- "mutation"
gruponomutation <- "no_mutation"

numero_grupo_mutation <- ncol(matrix_mutation)
numero_grupo_no_mutation <- ncol(matrix_nomutation)

matrix_nomutation_mutation <- cbind(matrix_nomutation,matrix_mutation)
countdata_nomutation_mutation <- matrix_nomutation_mutation 
countdata_nomutation_mutation <- as.matrix(countdata_nomutation_mutation) 


(condition_nomutation_mutation <- factor(c( rep(c(gruponomutation,grupomutation),
                            c(numero_grupo_no_mutation,numero_grupo_mutation)))))

(coldata_nomutation_mutation <- data.frame(row.names=colnames(countdata_nomutation_mutation), condition_nomutation_mutation))
dds_nomutation_mutation <- DESeqDataSetFromMatrix(countData=countdata_nomutation_mutation, 
                                       colData=coldata_nomutation_mutation, design=~condition_nomutation_mutation)

dds_nomutation_mutation

dds_nomutation_mutation$condition_nomutation_mutation <- relevel(dds_nomutation_mutation$condition_nomutation_mutation, 
                                                                 ref=gruponomutation)
dds_nomutation_mutation <- DESeq(dds_nomutation_mutation)

# Obtener los datos de Differential Expression
res_nomutation_mutation <- results(dds_nomutation_mutation)
table(res_nomutation_mutation$padj<0.05)
## Ordenar por adjusted p-value
res_nomutation_mutation <- res_nomutation_mutation[order(res_nomutation_mutation$padj), ]
## Unir con los datos de counts
resdata_nomutation_mutation <- merge(as.data.frame(res_nomutation_mutation), as.data.frame(counts(dds_nomutation_mutation, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_nomutation_mutation)[1] <- "Gene"
head(resdata_nomutation_mutation)
## Escribir resultados
write.csv(resdata_nomutation_mutation, file="diffexpr-mutationVSnomutation.csv")

###HASTA AQUÍ BIEN
#Aquí empieza el cagadero

                    ####NORMAL VS NO MUTACION (2)

gruponormal <- "normal"
gruponomutation2 <- "no_mutation"

no_normal <- ncol(control)
no_no_mutation <- ncol(matrix_nomutation)

#Dudoso
matrix_normal_nomutation <- cbind(control,matrix_nomutation)
countdata_normal_nomutation <- matrix_normal_nomutation
countdata_normal_nomutation <- as.matrix(countdata_normal_nomutation) 

(condition_normal_nomutation <- factor(c( rep(c(gruponormal,gruponomutation2),
                                     c(no_normal,no_no_mutation)))))
#Hay un paciente repetido
colnames(countdata_normal_nomutation) <- make.unique(colnames(countdata_normal_nomutation))

(coldata_normal_nomutation <- data.frame(row.names=colnames(countdata_normal_nomutation), 
                                         condition_normal_nomutation))

dds_normal_nomutation <- DESeqDataSetFromMatrix(countData=countdata_normal_nomutation,
                      colData=coldata_normal_nomutation, design=~condition_normal_nomutation)

dds_normal_nomutation

dds_normal_nomutation$condition_normal_nomutation <- relevel(dds_normal_nomutation$condition_normal_nomutation, ref=gruponormal)
dds_normal_nomutation <- DESeq(dds_normal_nomutation)

# Obtener los datos de Differential Expression
res_normal_nomutation <- results(dds_normal_nomutation)
table(res_normal_nomutation$padj<0.05)
## Ordenar por adjusted p-value
res_normal_nomutation <- res_normal_nomutation[order(res_normal_nomutation$padj), ]
## Unir con los datos de counts
resdata_normal_nomutation <- merge(as.data.frame(res_normal_nomutation), as.data.frame(counts(dds_normal_nomutation, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_normal_nomutation)[1] <- "Gene"
head(resdata_normal_nomutation)
## Escribir resultados
write.csv(resdata_normal_nomutation, file="diffexpr-controlVSnomutation.csv")


####NORMAL VS MUTACION (3)

grupo_normal <- "normal"
grupo_mutation <- "mutation"

numero_normal <- ncol(control)
no_mutation <- ncol(matrix_mutation)

matrix_normal_mutation <- cbind(control,matrix_mutation)
countdata_normal_mutation <- matrix_normal_mutation 
countdata_normal_mutation <- as.matrix(countdata_normal_mutation) 


(condition_normal_mutation <- factor(c(rep(c(grupo_normal,grupo_mutation),
                                                c(numero_normal,no_mutation)))))
#rownames duplicados
colnames(countdata_normal_mutation) <- make.unique(colnames(countdata_normal_mutation))

(coldata_normal_mutation <- data.frame(row.names=colnames(countdata_normal_mutation), condition_normal_mutation))
dds_normal_mutation <- DESeqDataSetFromMatrix(countData=countdata_normal_mutation, 
                                                  colData=coldata_normal_mutation, design=~condition_normal_mutation)

dds_normal_mutation
  
dds_normal_mutation$condition_normal_mutation <- relevel(dds_normal_mutation$condition_normal_mutation, 
                                                                 ref=grupo_normal)
dds_normal_mutation <- DESeq(dds_normal_mutation)

# Obtener los datos de Differential Expression
res_normal_mutation <- results(dds_normal_mutation)
table(res_normal_mutation$padj<0.05)
## Ordenar por adjusted p-value
res_normal_mutation <- res_normal_mutation[order(res_normal_mutation$padj), ]
## Unir con los datos de counts
resdata_normal_mutation <- merge(as.data.frame(res_normal_mutation), 
                                 as.data.frame(counts(dds_normal_mutation, normalized=TRUE)), 
                                 by="row.names", sort=FALSE)
names(resdata_normal_mutation)[1] <- "Gene"
head(resdata_normal_mutation)
## Escribir resultados
write.csv(resdata_normal_mutation, file="diffexpr-normalVSmutation.csv")

#Escribir resultados significativos  

#3
signif_normal_mutation <- resdata_normal_mutation[resdata_normal_mutation$padj< 0.05,]
  write.csv(signif_normal_mutation, file = "signif_normalVSmutation.csv")

#2
signif_normal_nomutation <- resdata_normal_nomutation[resdata_normal_nomutation$padj < 0.05,]
  write.csv(signif_normal_nomutation, file = "signif_normalVSnomutation.csv")

#1
signif_nomutation_mutation <- resdata_nomutation_mutation[resdata_nomutation_mutation$padj < 0.05,]
  write.csv(signif_nomutation_mutation, file = "signif_nomutationVSmutation.csv")
    
# miRNAs compartidos ------------------------------------------------------

####SE COMPARTEN MIRNAS ENTRE 3Y2? Y ENTRE 1 Y LAS OTRAS?
# 3 y 2 == normal_mutation VS normal_nomutation
# 1 y 2 == nomutation_mutation VS normal_nomutation
# 1 y 3 == nomutation_mutation VS normal_mutation

  #---Filtrando por P. Value
significativos_1 <- resdata_nomutation_mutation[resdata_nomutation_mutation$padj < 0.05,]
significativos_2 <- resdata_normal_nomutation[resdata_normal_nomutation$padj < 0.05,]
significativos_3 <- resdata_normal_mutation[resdata_normal_mutation$padj < 0.05,]

  #---Extraer compartidos

compartidos_3y2 <- intersect(as.character(significativos_3$Gene),
                                    as.character(significativos_2$Gene))

compartidos_1y2 <- intersect(as.character(significativos_1$Gene),
                             as.character(significativos_2$Gene))

compartidos_1y3 <- intersect(as.character(significativos_1$Gene),
                             as.character(significativos_3))

write.table(compartidos_3y2, file="compartidos_3y2.txt")
write.csv(compartidos_1y2, file="compartidos_1y2.csv")
write.csv(compartidos_1y3, file="compartidos_1y3.csv")

# Diagrama de Venn --------------------------------------------------------

####DIAGRAMA DE VENN DE LAS COMPARACIONES











