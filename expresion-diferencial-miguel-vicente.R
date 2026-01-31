# ===================================================

# Análisis de expresión diferencial funcional
# Fecha: "Thu Jan 29 16:41:47 2026"
# Autor: Miguel Angel Vicente Mesonero
# 
# ===================================================

# 0. Librerías
library(edgeR)

# 1. Cargo matriz de conteos y datos clínicos
mat <- read.delim("BRCA_exp_matrix.tsv",
                     header =T,
                     sep="\t",
                     fill =T,
                     stringsAsFactors = F)
clinical <- read.delim("clinical_info_TCGA-BRCA.tsv",
                       sep = "\t",
                       check.names = T,
                       header = TRUE,
                       fill=TRUE,
                       stringsAsFactors = F)

# 2. Exploro los datos
tibble::view(mat)
dim(mat) # Las dimensiones son 20500 genes (filas) x 1155 muestras (columnas)


colnames(clinical)  
dim(clinical) #114 variables clínicas

# 3. Normalizo datos (intersección, duplicados )

sum(colnames(mat) %in% clinical$sample) #1032/1054 muestras de mat están en clinical

  #Cojo solo las muestras de la matriz de expresión que tengo en clinical
snames<- intersect(colnames(mat), clinical$sample)

  #Vemos cuantos duplicados hay
sum(duplicated(clinical$sample)) 

  #Confirmamos que son duplicados
clinical [duplicated(clinical$sample),]

  #Me quedo con los que no son duplicados
clinical<- clinical [!duplicated(clinical$sample),]

  #Me quedo con los que están tanto en clínical como en matriz
rownames(clinical) <- clinical$sample
clinical <- clinical[snames,]
matriz <- mat[,snames]

  #Busco variables de er y edad
er <- grep("^er_", colnames(clinical))
colnames(clinical)[er]
colnames(clinical)[grep("age",colnames(clinical))]

#variables er_status_by_ihc & age_at_diagnosis

  #En clinical los que sean er positivos o negativos, es decir, quitamos [not evaluated]
clinical<-clinical[clinical$er_status_by_ihc == "Negative" |clinical$er_status_by_ihc == "Positive",]

  #En clinical: limpieza NAs e intersección con matriz
sum(is.na(clinical$age_at_diagnosis)) #Hay 0
matriz <- matriz[,rownames(clinical)]

  #Creo factor necesario para análisis
ER <- factor(clinical$er_status_by_ihc,levels = c("Positive","Negative"))
age <- clinical$age_at_diagnosis

design <- model.matrix()