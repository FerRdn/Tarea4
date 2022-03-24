###########################
##### MICROARREGLOS 2 #####
###########################

## 23 03 22
## Affymetrix Práctica 2
## Fer Rendón

## Librerías a usar

BiocManager::install ("Biobase")
library (Biobase)
library (limma)


## Cargar archivo con datos normalizados
exprs <- read.table (choose.files (), header = T, row.names = 1)
# Archivo expr_normalizada.txr
head (exprs)
# Muestra las primeros seis renglones de los microarrelgos normalizados

boxplot (exprs, col = rainbow (6))
# Boxplot de las expresiones de las seis muestras


## Calcular expresión diferencial
types <- factor (c ("KO", "KO", "KO", "WT", "WT", "WT"))
# Se describen las categorías de las seis muestras, KO y WT
types

design <- model.matrix (~ 0 + types)
# Crea matriz modelo con las categorías
# Me parece que también cuenta en cada renglón cuántas de esas hay
# Como solo hay una por renglón, y solo hay dos categorías
# Por eso pone 0s y 1s
# Ps es una matriz jaja

colnames (design) <- levels (types)
# A las columnas de esa matriz les asigna las categorías 
# creadas y guardadas en types

design
# ta daa

contMatrix <- makeContrasts (KO - WT, levels = design)
# Compara/contrasta los objetos en cada categoría
contMatrix
# Tiene dos niveles, KO y WT
# No termino de entender porqué da 1 y -1

fit <- lmFit (exprs, design)
# Ajuste lineal de los datos de expresión
# Como que ajusta cada muestra de exprs 
# de acuerdo a la matriz de design

fit2 <- contrasts.fit (fit, contMatrix)
# Estima contrastes deseados
# El ajuste que ya se había hecho se contrasta con
# El objeto que, vaya, contrasta los objetos en contMatrix

fit2 <- eBayes (fit2)
# Comprume las varianzas hacia valor global
# ... como¿

topTable (fit2, number = 20, sort.by = "p")
# Ya se tiene la lista de los genes diferencialmente expresados
# En base al P-value, se toman los 20 mejores
# Aquí sale 1428027, el gen eliminado


## Anotando los datos de expresión

library (mouse4302.db)

mouse4302 () # Datos del paquete

# Se deben obtener los probe IDs del objeto fit2
# Y se guardan en otro objeto
probes <- fit2$genes$ID

descriptions <- mget (probes, mouse4302GENENAME)
# Los probes los busca en los datos de Gene name y los guarda en un objeto

symbols <- mget (probes, mouse4302SYMBOL)
# Same pero ahora símbolos, nombres comunes o cortos

entrezids <- mget (probes, mouse4302ENTREZID)
# Same pero ahora entrezids, identificadores

# Todo lo obtenido se guarda en el objeto fit2
fit2$genes$EntrezID <- unlist (entrezids)
fit2$genes$Symbol <- unlist (symbols)
fit2$genes$Description <- unlist (descriptions)
# Unlist pasa anotaciones de texto a una columna

# Ahora cuando se quiere identificar a los 20 genes más
# "importantes" según su P-value, resulta más fácil
topTable (fit2, number = 20, sort.by = "p")


## Gráfica de volcán

volcanoplot (fit2, highlight = 10, names = fit2$genes$Symbol)
# Se hace el plot a partir de los datos en fit2
# Los 10 genes... ¿con un P-value más pequeño? son resaltados
# A esos 10 se les pone su nombre común (symbols)


# Cuadro con sondas elegidas
# Fold-change > 1.5 y P-value < 0.05
deTable <- topTable (fit2, number = nrow (fit2),
                     lfc = log2 (1.5), p.value = 0.05)
# FC siempre va a ser con log base 2

dim (deTable)

# Cuadro con datos del microarreglo pero ordenenados según logFC
fullTable <- topTable (fit2, number = nrow (fit2), sort.by = "logFC")
dim (fullTable)


## Guardar resultados
# Que no puedo guardar porque solo pude hacer la mitad
write.table (fullTable, file = "full_results.txt", row.names = F,
             sep = "\t", quote = F)


## Ejercicios
# Que tampoco puedo hacer porque no me muestra $genes
# Nada, no hay, no existe

# 1. ¿Cuántas sondas predicen como diferencialmente expresadas?
difexp <- topTable (fit2, number = nrow (fit2), p.value = 0.05)
# Es lo mismo de arriba pero solo con P-value < 0.05
dim (difexp)
# Son 417

# 2. ¿Cuántas decrementan y cuántas aumentan su expresión en el KO?


# 3. ¿Cuántos genes únicos hay en estas listas?
unique (fit2$genes)
# Me atrevería a decir que todos son únicos jaja



