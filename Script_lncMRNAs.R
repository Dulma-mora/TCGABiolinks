library(TCGAbiolinks)
library(Biobase)
library(SummarizedExperiment)
library(dplyr)


# Query lncRNA | Esto no recuerdo que es honestly --------------------

query_lnc <- GDCquery(project = "TCGA-BRCA",
                        legacy = FALSE,
                        data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "HTSeq - Counts")

GDCdownload(query_lnc, method = "api")
data <- GDCprepare(query_lnc)

# CAGADERO-IGNORAR --------------------------------------------------------------

  #Obteniendo datos clinicos
data_clin <- GDCquery_clinic(project = "TCGA-BRCA", "Clinical")

  #Extrayendo Read Counts 
BRCAmatrix <- assay(data, "HTSeq - Counts")

# saveRDS(c(query_lnc, data, data_clin, BRCAmatrix), "data_lncBRCA.rds") 

# readRDS("data_lncBRCA.rds")       

write.csv(BRCAmatrix, "BRCAmatrix.csv") #Guardando para cuando la cague
BRCAmatrix <- read.csv("BRCAmatrix.csv") 

  #Quedandome con los primeros 12 digitos
correccion <- colnames(BRCAmatrix)
correccion <- substr(correccion, 1, 12)
  colnames(BRCAmatrix) <-correccion

  #Cambiando ensembl por synbol gene
library(org.Hs.eg.db)
  
#  cambio <- mapIds(org.Hs.eg.db,
                          #  keys=row.names(BRCAmatrix),
                          #  column="SYMBOL",
                          #  keytype="ENSEMBL",
                          #  multiVals="first")
#  cambio
  
#  symbol_gene <- merge(BRCAmatrix,cambio, by.x="row.names", by.y="row.names")
#  colnames(symbol_gene)[1] <- "Ensemble"
  
#  write.csv(symbol_gene, "BRCAmatrix_preparada.csv")



  
# Obteniendo lncRNA por su ensembl ID | Creando MB con BiomaRt ----------------------------

# Matriz de expresión de pacientes (BRCAmatrix)

matriz_chafa <- read.csv("BRCAmatrix.csv") 
  BRCAmatrix <- matriz_chafa[,-1]
  rownames(BRCAmatrix) <- matriz_chafa[,1]
  head(BRCAmatrix)
  
# Ahora se tiene que filtrar la matriz y tienen que quedar los
# transcritos que sean "lncRNA"
  
library("biomaRt")
library("XML") 

# Buscando servicios disponibles de BiomaRt
listMarts()

# Creando MART
ensembl <- useMart("ensembl") # MART
datasets <- listDatasets(ensembl)
head(datasets) 

# Seleccionando dataset (H. sapiens gene ensembl)
ensembl = useDataset("hsapiens_gene_ensembl", mart=ensembl)  


# Numerando atributos
atributes = listAttributes(ensembl)
atributes # Usaremos transcript_biotype |	Transcript type	| sequences

# Numerando filtros 
filters = listFilters(ensembl) 
filters # Usaremos ensembl_gene_id |Gene stable ID(s) [e.g. ENSG00000000003]

# Corroborando el filtro
filtro <- searchFilters(mart = ensembl, pattern = "ensembl.*id")
# Otra forma de buscar filtro, nos aparece inmediatamente 


# Preparando getMB

values_BRCAmatrix <- rownames(BRCAmatrix)
write.csv(values_BRCAmatrix, "ensembl_id.csv")

BM <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
      filters = 'ensembl_gene_id', 
      values = values_BRCAmatrix, 
      mart = ensembl, uniqueRows = TRUE)   # OMG FUNCIONO


# Haciendo pruebas

# Prueba 1. 9015 ENSG00000147761 - ENSG00000147761
valor_prueba1<- "ENSG00000147761"
prueba1 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
            filters = 'ensembl_gene_id', 
            values = valor_prueba1, 
            mart = ensembl, uniqueRows = TRUE) # COINCIDE!!


# Prueba 2 - 38328 | 38331 - ENSG00000249383 

valor_prueba2<- "ENSG00000249383"   # No están en el mismo lugar
prueba2 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
                 filters = 'ensembl_gene_id', 
                 values = valor_prueba2, 
                 mart = ensembl, uniqueRows = TRUE) #COINCIDE


# Prueba 3 - ENSG00000281398|56425 <- 56431

valor_prueba3<- "ENSG00000281398"   # No están en el mismo lugar
prueba3 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
                 filters = 'ensembl_gene_id', 
                 values = valor_prueba3, 
                 mart = ensembl, uniqueRows = TRUE) # COINCIDE

# Prueba 4 - 56493 | ENSG00000281920 <- 56499

valor_prueba4<- "ENSG00000281920"   # No están en el mismo lugar
prueba4 <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
                 filters = 'ensembl_gene_id', 
                 values = valor_prueba4, 
                 mart = ensembl, uniqueRows = TRUE)
#FUNCIONA



# Filtrando lncRNA de la matriz  ------------------------------------------

# Quedandome solo con lncRNA de BM

library(dplyr)

# Obteniendo lncRNA
lnc_BM <- filter(BM, gene_biotype=="lncRNA") # BM filtrado, solo lncRNA
# 14085 len

# Creando indice de compartidos entre la matriz y lncBM
lncRNA_index <- lnc_BM$ensembl_gene_id %in% rownames(BRCAmatrix)

  # --- MATRIZ CON SOLO lncRNA
# DATAFRAMEAFILTRAR[rows a filtrar , columnas a filtrar]

lncRNA_BRCA <- BRCAmatrix[lncRNA_index,]


#Confirmando que todos sean lncRNA

values_confirm <- rownames(lncRNA_BRCA)
BM_confirm <- getBM(attributes= c('gene_biotype', "ensembl_gene_id"), 
            filters = 'ensembl_gene_id', 
            values = values_confirm, 
            mart = ensembl, uniqueRows = TRUE)





















  