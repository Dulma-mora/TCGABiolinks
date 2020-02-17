library(TCGAbiolinks)
library(Biobase)
library(SummarizedExperiment)

#---

query_mirna <- GDCquery(project = "TCGA-BRCA",
                              data.category = "Transcriptome Profiling",
                              data.type = "miRNA Expression Quantification",
                              sample.type = "Primary solid Tumor")

GDCdownload(query_mirna)
data_mirna <- GDCprepare(query_no_mutation, save = TRUE)

#--
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
id <- colnames(uncommon_symbol)
id <- id[-c(1,2,length(id))] #no mutación

id_mutation <- read.csv("mutacion.csv")
id_mutation <- as.character(id_mutation[,2])

################---Haciendo data.frame con mirnaPERRON
    #--- Sin mutaciones
matrix_nomutation <- mirnaPERRRON[,id] #no ideal pero sí salío lol

    #--- Con mutaciones
matrix_mutation <- mirnaPERRRON[,colnames(mirnaPERRRON) %in% id_mutation] #index


    #---Creando control
query_control <- GDCquery(project = "TCGA-BRCA",
                        data.category = "Transcriptome Profiling",
                        data.type = "miRNA Expression Quantification",
                        sample.type = "Solid Tissue Normal")

GDCdownload(query_control)
data_control <- GDCprepare(query_no_mutation, save = TRUE)

rownames(data_control) <- data_control$miRNA_ID
   
   #ReadCounts
control <- data_control[,c(grep("read_count",colnames(data_control)))]


vacio_control <- gsub("read_count_", "",colnames(control))
vacio_control <- substr(vacio_control, 1, 12)

colnames(control) <-vacio_control
rownames(control) <- data_control$miRNA_ID

###---Expresión diferencial 
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

# Regularizar la transformacion de log en rlog para poder graficar en heatmap
#rld <- rlogTransformation(dds)
# Extraer la transformacion del s4class rld
#rld3 <- as.matrix(t(assay(rld)))

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

matrix_tumor <- cbind(matrix_nomutation,matrix_mutation)
countdata_nomutation_mutation <- matrix_tumor 
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
write.csv(resdata, file="diffexpr-mutationVSnomutation.csv")

                    ####NORMAL VS NO MUTACION (2)

gruponormal <- "normal"
gruponomutation2 <- "no_mutation"

no_normal <- ncol(control)
no_no_mutation <- ncol(matrix_nomutation)

matrix_normal_nomutation <- cbind(control,matrix_nomutation)
countdata_normal_nomutation <- matrix_normal_nomutation
countdata_normal_nomutation <- as.matrix(countdata_normal_nomutation) 

(condition_normal_nomutation <- factor(c( rep(c(gruponormal,gruponomutation2),
                                     c(no_normal,no_no_mutation)))))
#
(coldata_normal_nomutation <- data.frame(row.names=colnames(countdata_normal_nomutation), 
                                         condition_normal_nomutation))

dds_normal_nomutation <- DESeqDataSetFromMatrix(countData=countdata_normal_nomutation,
                      colData=coldata_normal_nomutation, design=~condition_normal_nomutation)

dds_mutation

dds_mutation$condition_mutation <- relevel(dds_mutation$condition_mutation, ref=grupomutation)
dds_mutation <- DESeq(dds_mutation)

# Obtener los datos de Differential Expression
res_mutation <- results(dds_mutation)
table(res_mutation$padj<0.05)
## Ordenar por adjusted p-value
res_mutation <- res_mutation[order(res_mutation$padj), ]
## Unir con los datos de counts
resdata_mutation <- merge(as.data.frame(res_mutation), as.data.frame(counts(dds_mutation, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata_mutation)[1] <- "Gene"
head(resdata_mutation)
## Escribir resultados
write.csv(resdata, file="diffexpr-mutationVSnomutation.csv")



####NORMAL VS MUTACION (3)

####SE COMPARTEN MIRNAS ENTRE 3Y2? Y ENTRE 1 Y LAS OTRAS?

####DIAGRAMA DE VENN DE LAS COMPARACIONES





