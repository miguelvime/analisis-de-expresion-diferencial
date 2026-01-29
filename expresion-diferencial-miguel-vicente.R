# ===================================================

# Análisis de expresión diferencial funcional
# Fecha: "Thu Jan 29 16:41:47 2026"
# Autor: Miguel Angel Vicente Mesonero
# 
# ===================================================

# 1. Cargo matriz de conteos y datos clínicos
mat <- read.delim("BRCA_exp_matrix.tsv",
                     header =T,
                     sep="\t",
                     fill =T,
                     stringsAsFactors = F)
library(data.table)
clinical <- fread("clinical_info_TCGA-BRCA.tsv")
  
  # read.delim("clinical_info_TCGA-BRCA.tsv",
  #                      sep = "\t",
  #                      check.names = T,
  #                      header = TRUE,
  #                      fill=TRUE,
  #                      stringsAsFactors = F)

#2. Exploro los datos
tibble::view(mat)
dim(mat)
###En el dataframe las filas corresponden a genes (20500-ALPL,ALPK1,etc-)
### Las columnas son muestras (1155)


#