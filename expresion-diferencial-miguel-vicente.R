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

# 3. Limpieza de duplicados y sincronización de datos 

sum(colnames(mat) %in% clinical$sample) #1032/1054 muestras de mat están en clinical

  #Cojo solo las muestras de la matriz de expresión que tengo en clinical
snames<- intersect(colnames(mat), clinical$sample)

  #Vemos cuantos duplicados hay
sum(duplicated(clinical$sample)) 

  #Confirmamos que son duplicados
clinical [duplicated(clinical$sample),]

  #Me quedo con los que no son duplicados
clinical<- clinical [!duplicated(clinical$sample),]

  #Sincronización de datos
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

  #En clinical: limpieza NAs e sincronización entre matrices
sum(is.na(clinical$age_at_diagnosis)) #Hay 0
matriz <- matriz[,rownames(clinical)]

  #Creo factor necesario para análisis
ER <- factor(clinical$er_status_by_ihc,levels = c("Positive","Negative"))
age <- clinical$age_at_diagnosis

  #Creamos matriz de diseño a partir de edad y ER
design <- model.matrix(~ER+age)

  #Creo la DGEList a partir de la matriz (información de los pacientes)
y <- DGEList(matriz)

# 4. Normalización de datos

y<- calcNormFactors(y)

  #Elimino aquellos con expresión baja
keep <- filterByExpr(y = y,design = design)
y <- y [keep,,keep.lib.sizes = TRUE]

y <-edgeR::estimateCommonDisp(y)

# 5. Modelado y Test 
 #Creación de modeloque explique expresión de cada gen, teniendo en cuenta la matriz de diseño, que se creó viendo las distribuciones de edad y ER
fit <- edgeR::glmQLFit(y,design)

  # Comparo diferencias de expresión en ER negativo contra positivo en el modelo creado previamente con glmQFit
qlf<- edgeR::glmQLFTest(fit,coef = "ERNegative")

  #Extraemos los genes que más explican diferencias entre ER+ y ER-
pvalues <-edgeR::topTags(qlf,n=Inf, adjust.method = "BH",sort.by = "PValue",p.value=1)
p_values <-as.data.frame (pvalues)
alpha <- 0.05
significativos <- p_values[p_values$FDR<alpha,]

