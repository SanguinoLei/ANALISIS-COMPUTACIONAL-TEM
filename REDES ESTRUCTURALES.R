library(deSolve)
library(readxl)
library(igraph)

View(EGFR)

##trasponer matriz##

tEGFR<-t(EGFR)
tEGFR<-tEGFR[-1,]

##binarizacion###

binEGFR<-binarizeTimeSeries(tEGFR)
binEGFR
##reconstruir la red apartir de la binarizacion##


netEGFR<-reconstructNetwork(binEGFR$binarizedMeasurements,
                            method = "bestfit",
                            maxK = 1,
                            requiredDependencies = list(
                              "GRB2"= c("GRB2","EGFR"),
                              "SOS2"= c("SOS2","GRB2"),
                              "SOS1"= c("SOS1","GRB2"),
                              "MRAS...9"= c("MRAS...9","SOS2","SOS1"),
                              "MRAS...7"= c("MRAS...7","SOS2","SOS1"),
                              "KRAS"= c("KRAS","MRAS...7","MRAS...9"),
                              "RAF1"= c("RAF1","MRAS...7", "MRAS...9"),
                              "MMP2"= c("MMP2","JUN"),
                              "MMP7"= c("MMP7","JUN"),
                              "MMP11"= c("MMP11","JUN"),
                              "MMP24"= c("MMP24","JUN"),
                              "MMP28"= c("MMP28","JUN"),
                              "MMP3"= c("MMP3","JUN"),
                              "MMP9"= c("MMP9","JUN"),
                              "MMP10"= c("MMP10","JUN"),
                              "IL8"= c("IL8","JUN"),
                              "VEGFA"= c("VEGFA","JUN"),
                              "CCND1"= c("CCND1","JUN","MYC","ETS1","FOS"),
                              "CDK4"= c("CDK4","JUN","MYC","ETS1","FOS"),
                              "EGFR"= c("EGFR","EGF"),
                              "MMP2"= c("MMP2","MYC"),
                              "MMP7"= c("MMP7","MYC"),
                              "MMP11"= c("MMP11","MYC"),
                              "MMP24"= c("MMP24","MYC"),
                              "MMP28"= c("MMP28","MYC"),
                              "MMP3"= c("MMP3","MYC"),
                              "MMP9"= c("MMP9","MYC"),
                              "MMP10"= c("MMP10","MYC"),
                              "IL8"= c("IL8","MYC"),
                              "VEGFA"= c("VEGFA","MYC"),
                              "MMP2"= c("MMP2","FOS"),
                              "MMP7"= c("MMP7","FOS"),
                              "MMP11"= c("MMP11","FOS"),
                              "MMP24"= c("MMP24","FOS"),
                              "MMP28"= c("MMP28","FOS"),
                              "MMP3"= c("MMP3","FOS"),
                              "MMP9"= c("MMP9","FOS"),
                              "MMP10"= c("MMP10","FOS"),
                              "IL8"= c("IL8","FOS"),
                              "VEGFA"= c("VEGFA","FOS"),
                              "MMP2"= c("MMP2","ETS1"),
                              "MMP7"= c("MMP7","ETS1"),
                              "MMP11"= c("MMP11","ETS1"),
                              "MMP24"= c("MMP24","ETS1"),
                              "MMP28"= c("MMP28","ETS1"),
                              "MMP3"= c("MMP3","ETS1"),
                              "MMP9"= c("MMP9","ETS1"),
                              "MMP10"= c("MMP10","ETS1"),
                              "IL8"= c("IL8","ETS1"),
                              "VEGFA"= c("VEGFA","ETS1"),
                              "ITPKB"= c("ITPKB","PDGFA","PDGFB","PDGFC","PDGFD"),
                              "AKT3"= c("AKT3","PDGFA","PDGFB","PDGFC","PDGFD"),
                              "BAD"= c("BAD","ITPKB","AKT3"),
                              "MDM2"= c("MDM2","ITPKB","AKT3")))

                            
netEGFR


###graficar##

window()

plotNetworkWiring(netEGFR, edge.arrow.size=0.1, vertex.label.cex=0.7, vertex.color = "lightgrey", edge.width = 0.2, vertex.size=50, vertex.label.color = "black", rescale= FALSE,  xlim = c(-4, 4), ylim = c(-3.5, 3.5), margin = c(0.1, 0.1, 0.1, 0.1))


###booleana JAK-STAT microarreglos###

View(JAK)

tJAK<-t(JAK)
tJAK<-tJAK[-1,]

binJAK<-binarizeTimeSeries(tJAK)

netJAK <- reconstructNetwork(binJAK$binarizedMeasurements,
                             method="bestfit",
                             maxK=2,
                             requiredDependencies = list(
                               "STAT3"= c("STAT3","JAK1"),
                               "VEGFA"= c("VEGFA","STAT3"),
                               "PIM1"= c("PIM1","STAT3"),
                               "PIM2"= c("PIM2","STAT3")
                             ))

plotNetworkWiring(netJAK, edge.size= 2, edge.arrow.size=0.2, vertex.label.cex=0.4, vertex.color = "lightgrey", edge.width = 0.5, vertex.size=50, vertex.label.color = "black")
netJAK

###booleana WNT microarreglos###

View(WNT)

tWNT<-t(WNT)
tWNT<-tWNT[-1,]

binWNT<-binarizeTimeSeries(tWNT)

netWNT <- reconstructNetwork(binWNT$binarizedMeasurements,
                             method="bestfit",
                             maxK=2,
                             requiredDependencies = list(
                             "LRP6"= c("LRP6","WNT1"),
                             "LRP5"= c("LRP5","WNT1"),
                             "DVL2"= c("DVL2","LRP6"),
                             "DVL2"= c("DVL2","LRP5"),
                             "DVL3"= c("DVL3","LRP6"),
                             "DVL3"= c("DVL3","LRP5"),
                             "GBP3"= c("GBP3","LRP6"),
                             "GBP3"= c("GBP3","LRP5"),
                             "DVL1 /// LOC642469"= c("DVL1 /// LOC642469","LRP5"),
                             "DVL1 /// LOC642469"= c("DVL1 /// LOC642469","LRP6"),
                             "LEF1"= c("LEF1","AXIN1"),
                             "LEF1"= c("LEF1","AXIN2")
                             ))

plotNetworkWiring(netWNT, edge.arrow.size=0.1, vertex.label.cex=0.2, vertex.color = "lightgrey", edge.width = 0.3, vertex.size=60, vertex.label.color = "black", rescale= FALSE, xlim = c(-2.5, 2.5), ylim = c(-2.5, 2.5))

netWNT

