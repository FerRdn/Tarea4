#########################
##### MICROARREGLOS #####
#########################

## 23 03 22
## Affymetrix Práctica

# http://datos.langebio.cinvestav.mx/~cei/cursos/R/

library (affy)
library (limma)
library (mouse4302.db)
library (mouse4302cdf)
library (pvclust)
library (vsn)

## Para trabajar solo en una carpeta
setwd ("C:/Users/monta/Documents/Tareas/Universidad/6to Semestre/Genómica Funcional/R Genómica/mir155")
getwd ()

## Cuadro de relación entre nombre de c/archivo y sus caract
pd <- read.table (choose.files (), header = T, as.is = T) # Archivo pdata.txt
# as.is permite que se lean como caracteres y no como factores
pd

## Archivos de imagen CEL
affyData <- ReadAffy (filenames = pd$filename)
# Mágicamente encuentra todos los archivos CEL
# en la carpeta en la que estoy trabajando
pData (affyData) <- pd # Esos archivos se guardan dentro del cuadro pd
sampleNames (affyData) <- pd$name # Cambia los nombres de identificación

affyData # :o


## Análisis

boxplot (affyData, col = rainbow (6))
# Boxplot de los seis microarreglos
# Con relleno de orgullosamente LGBT+
# Escoge específicamente seis colores para las 6 muestras

hist (affyData, col = rainbow (6))
# Lo mismo pero en histograma
# Es como la intensidad del color
# Y con ello su nivel de expresión
# Todas se ven en el mismo am lugar? jajaj

image (affyData [ , 3])
# Se supone que sale una imagen para esa muestra
# Pero como que no sale jaja
# Ya salió después de un ratillo
# Es puro bco y negro y supongo que si algo está mal
# Por cualquier error que exista
# Se vería de otro color y forma

heatmap (cor (exprs (affyData)), symm = T)
# Se ve la relación entre muestras
# Entre más oscuro, más parecidas?
# Es notorio que entre WTs y KOs son diferentes

corClust <- pvclust (exprs (affyData),
                     nboot = 1, method.dist = "correlation")
# Lo mismo pero con un dendrograma
# Requiere de pvclust y está haciento clusters
# de los microarreglos de expresión en affyData
# La distancia la obtiene de la correlación
# Es decir, entre más cercanos, más correlacionados
# nose q es nboot
plot (corClust)
# Los WTs se encuentran dentro de un mismo cluster
# Los KOs en otro, excepto por KO1
# No es anómalo pero sería mejor normalizar

# Para ver que, en efecto, no es anómalo
# Se hace un PCA para ver que se agrupe con los KOs
pca <- princomp (exprs (affyData))
# Objeto con análisis

plot (pca$loadings, main = "PCA",
      col = rainbow (6),  pch = 19, cex = 2)
# Plot con las muestras, título, color, forma y tamaño

text (pca$loadings, colnames (exprs (affyData)), pos = 3, cex = 0.8)
# Texto para cada objeto en la gráfica

# Sí se agrupa con los demás KOs


## Normalización
eset <- rma (affyData)
# Normalización de datos con
# RMA, Robust Microarray Average
# Especial para microarreglos

# Jaja no sale
# Sí sale pero si no sale
# Aplica la de los inges en sistemas

# No es el único, hay otros que se pueden ver con
normalize.AffyBatch.methods ()


## Comparación
## Pre-normalización vs Post-normalización

par (mfrow = c (1, 2))
# Para ver dos gráficas en una sola

boxplot (affyData, col = rainbow (6))
# Pre-normalización
boxplot (data.frame (exprs (eset)), col = rainbow (6))
# Post-normalización

par (mfrow = c (1, 1))
# Siempre regresar a una


par (mfrow = c (1,2))
# Dos gráficas en una
corClust <- pvclust (exprs (affyData),
                    nboot = 1, method.dist = "correlation")
plot (corClust, main = "Pre-normalización")
# El pre

corClustAfter <- pvclust (exprs (eset),
                          nboot = 1, method.dist = "correlation")
plot (corClustAfter, main = "Post-normalización")
# El post

par(mfrow = c (1, 1))
# Regresar a 1 gráfica


# Si bien, se normalizan los datos, no se quiere perder aquellos
# Datos que tuvieron una alta expresión
# Puede verse esto con:
meanSdPlot (exprs (eset))
# Gráfica del SD de los datos normalizados
# No se pierden los datos que tuvieron alta expresión
# jaja no sale
# Jaja ya sale


## Comportamiento de genes
# Como se eliminó gen 1428027_at
# Se puede ver su comportamiento en comparación
# Con los WTs, controles
boxplot (data.frame (exprs (eset)), col = "mediumpurple3");
lines (exprs (eset) ["1428027_at", ],
       lwd = 2, type = "b", col = "purple4")
# Boxplot de las expresiones con unas líneas que
# muestran la expresión específica del gen eliminado


## Guardar resultados
write.exprs (eset, file = "expr_normalizada.txt")
# Archivo de texto con datos normalizados

save (eset, file = "eset.Rdata")
# Archivo de R con datos normalizados


