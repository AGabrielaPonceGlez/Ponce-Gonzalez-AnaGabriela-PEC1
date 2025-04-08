# Intalación de librerias necesarias:
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("POMA")
library(POMA)
install.packages("ggtext", repos = "https://cloud.r-project.org")
library(ggtext)
library(magrittr)
library(tidyverse)


# PASOS PARA LA CREACIÓN DEL OBJETO:

# En primer lugar cargamos el dataset elegido. Este archivo contine datos y metadatos juntos:
features <- readr::read_csv("DataValues_S013.csv", col_names = TRUE)

# Vamos a realizar una exploración básica de los mismos para ver cómo están organizados:
head(features, 10)
str(features)

# Vemos que el dataset incluye 39 observaciones (una por paciente o muestra) y 696 variables (metabolitos y otras medidas clínicas).
# Los nombres de las columnas reflejan tanto la variable medida como el tiempo del muestreo (por ejemplo, "GLU_T0", "INS_T0", "bmi_T0"), lo que indica un seguimiento longitudinal de distintos parámetros.
- Los valores están expresados en formato numérico de tipo continuo, lo que es de esperarse al tratarse de este tipo de variables.
- El dataset incluye algunos valores nulos ("NA"), lo cual puede corresponderse con limitaciones técnicas o registros incompletos; pero que justifican la necesidad de hacer una limpieza previa al análisis exploratorio, y en este caso, una imputación de dichos valores.


# Después de haber visto como se organizan, las primeras variables, la presencia de valores nulos, el tipo de datos... vamos a proceder a limpiar los datos. Para ello, eliminamos la primera columna y la fila 26, ya que los metadatos no pueden tener valores nulos (NA):
features <- features[-1]
target <- features[-26,1:9]

# Luego, vamos a modificar el nombre de la observación "bypass" por "by_pass":
target$SURGERY[target$SURGERY == "by pass"] = "by_pass"

features <- features[-26,-c(1:9)]

# Como ya la base de datos está acomodada, podemos crear el objeto. Como comentamos anteriormente, hemos optado por el uso de POMA para crear la clase SummarizedExperiment (SE), ya que la función PomaCreateObject() del paquete POMA genera internamente un objeto de clase POMA que encapsula un SE
objeto <- PomaCreateObject(metadata = target, features = features)


# GUARDAR EL OBJETO EN EL FORMATO REQUERIDO:
save(objeto, file = "SummarizedExperiment_objeto.Rda")




# _________________________________________________________________________________________________________________________________


# ANÁLISIS EXPLORATORIO:

# Para poder hacer el análisis exploratorio, iniciamos con la imputación de los NA empleando el método kNN:
imputed <- objeto %>% 
  PomaImpute(method = "knn", zeros_as_na = TRUE, 
             remove_na = TRUE, cutoff = 20)


# Ahora, aplicamos la normalización por auto-escalado (z-score), que centra y escala cada variable para que tenga media cero y ds uno. De esta forma podemos comparar las variables sin que la escala afecte al resto del análisis:
normalized <- imputed %>% 
  PomaNorm(method = "auto_scaling")

# A continuación, vamos a generar gráficos de caja por muestra y variable para detectar outliers y ver el efecto de la normalización:

# Boxplot por muestras antes de normalizar:
PomaBoxplots(imputed, x = "samples") +
  ggtitle("No Normalizado") +
  theme(legend.position = "none")

# Boxplot por variables antes de normalizar (computacionalmente muy demandante):
PomaBoxplots(imputed, x = "features") +
  ggtitle("No Normalizado") +
  theme(legend.position = "none")

# Boxplot por muestras después de normalizar:
PomaBoxplots(normalized, x = "samples") +
  ggtitle("Normalizado") +
  theme(legend.position = "none")

# Boxplot por variables después de normalizar (computacionalmente muy demandante):
PomaBoxplots(normalized, x = "features") +
  ggtitle("Normalizado") +
  theme(legend.position = "none")

# Por otro lado, PomaDensity nos muestra la distribución de todas las características antes y después del proceso de normalización. La normalización debería homogeneizar la dispersión:
PomaDensity(imputed) +
  ggtitle("No normalizado") +
  theme(legend.position = "none")
PomaDensity(normalized) +
  ggtitle("Normalizado") +
  theme(legend.position = "none")

# Ahora eliminamos posibles muestras outlier con la función PomaOutliers. De esta forma, trabajamos con un conjunto más robusto para el análisis multivariante:
pre_processed <- PomaOutliers(normalized)$data


# Analisis de componentes principales:
poma_pca <- PomaPCA(pre_processed, load_length = 0)


# Ahora, representamos la varianza explicada por cada componente. Con esta primera se explican aprox el 60% de la varianza de los datos:
plot(poma_pca$eigenvalues_plot)

# Esta segunda muestra la distribucion de las muestras de acuerdo a su varianza:
plot(poma_pca$biplot)

# A continuación estudiamos las correlaciones entre metabolitos utilizando la matriz de correlaciones de Pearson:
poma_cor <- PomaCorr(pre_processed, method = "pearson", 
                     cluster = TRUE, sig_level = 0.6)

# Visualizamos las correlaciones entre las variables:
plot(poma_cor$corrplot)

# Ahora, con un análisis univariante identificamos los metabolitos diferencialmente expresados entre los grupos quirúrgicos. Empleamos el test de Mann-Whitney (no paramétrico) que es adecuado para datos sin distribución normal:
results <- pre_processed %>% 
  PomaUnivariate(method = "mann") %>% 
  dplyr::select(feature, fold_change, pvalue)

# Representamos los resultados mediante un volcano plot que muestra la relación entre la magnitud del cambio (log2 Fold Change) y la significancia estadística: (-log10 p-value):
results %>% 
  PomaVolcano(pval_cutoff = 0.05,
              log2fc_cutoff = NULL,
              labels = TRUE,
              x_label = "log2 (Fold Change)",
              y_label = "-log10 (P-value)")

# Finalmente, aplicamos un modelo lineal mediante el método LIMMA, mejorando la detección de diferencias significativas entre grupos al modelar mejor la varianza. Definimos el contraste entre los grupos "tubular" y "by_pass":
results <- pre_processed %>%
  PomaLimma(contrast = "tubular - by_pass") %>% 
  dplyr::select(feature, log2FC, pvalue)

# Representamos de nuevo los resultados con un volcano plot:
results %>% 
  PomaVolcano(pval_cutoff = 0.05,
              log2fc_cutoff = NULL,
              labels = TRUE,
              x_label = "log2 (Fold Change)",
              y_label = "-log10 (P-value)")
