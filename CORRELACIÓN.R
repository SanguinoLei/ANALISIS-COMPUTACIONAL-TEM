alibrary(qgraph)

#### CORRELACION DE LA VIA EGF, IGF, PDGF####
genes <- EGFR[, -1]

# Calcular la matriz de correlación de Pearson
correlation1_matrix <- cor(genes, method = "pearson")

qgraph(correlation1_matrix, layout="spring", legend=FALSE)

#### CORRELACION DE LA VIA WNT BETA CATENINA####
genes2 <- WNT[, -1]

# Calcular la matriz de correlación de Pearson
correlation2_matrix <- cor(genes2, method = "pearson")

qgraph(correlation2_matrix, layout="spring", legend=FALSE )

######CORELACION JAK STAT######

genesJAK<- JAK[, -1]

# Calcular la matriz de correlación de Pearson
correlationJAK_matrix <- cor(genesJAK, method = "pearson")
