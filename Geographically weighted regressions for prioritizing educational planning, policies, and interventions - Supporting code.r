options(warn = -1) 
#install.packages("GWmodel")
library("GWmodel") 
#install.packages("gwrr")
library("gwrr")
#install.packages("sp")
library("sp")
#install.packages("gstat")
library("gstat")
#install.packages("ModelMap")
library("ModelMap")
#install.packages("RColorBrewer")
library("RColorBrewer")
#install.packages("ape")
library("ape")
#install.packages("usdm")
library("usdm")
#install.packages("rgdal")
library("rgdal")
#install.packages("tidyr")
library("tidyr")
#install.packages("MASS")
library("MASS")
#install.packages("Compositional")
library("Compositional")
#install.packages("dplyr")
library("dplyr")
#install.packages("spatialEco")
library("spatialEco")
#install.packages("broom")
library("broom")
install.packages("pastecs") 
library("pastecs")

setwd("C:\\Users\\g.vargas\\BOX\\IIEP_MyProjects\\MP_01000298_RND_SDA\\WorkFiles_Experts\\298-Issue-Papers\\298-Issue-Paper-GWR\\Replication files")

PolygonShape <- readOGR(dsn = "Data\\Colombia", layer = "GWR - Colombia", drop_unsupported_fields=FALSE, disambiguateFIDs=TRUE)
cat(" ")
cat("The final names of the variables that were imported are the following:")

colnames(PolygonShape@data)

ListVariables <- list("E_s11_to_1", "MurderRate", "ProxyGDP", "Poverty", "SISBEN1PC", "GC_indrura", "Vulnerabil", "Threat", "Lack_Respo", "TransEducP", "TransAlimE", "HS_Cober_7", "HS_Cober_1", "G_IGA_tota", "GC_discapi")

PolygonShapeTrimmed <- PolygonShape[c("Lat", "Long", "admin2Pcod", "E_s11_to_1", "SISBEN1PC", "MurderRate", "ProxyGDP", "Poverty", "GC_indrura", "Vulnerabil", "Threat", "Lack_Respo", "TransEducP", "TransAlimE", "HS_Cober_7", "HS_Cober_1", "G_IGA_tota", "GC_discapi")]
PolygonShapeTrimmed <- sp.na.omit(PolygonShapeTrimmed)

VariablesToSummarize <- PolygonShapeTrimmed[c("E_s11_to_1", "SISBEN1PC", "MurderRate", "ProxyGDP", "Poverty", "GC_indrura", "Vulnerabil", "Threat", "Lack_Respo", "TransEducP", "TransAlimE", "HS_Cober_7", "HS_Cober_1", "G_IGA_tota", "GC_discapi")]@data
stat.desc(VariablesToSummarize)

DeVar <- ListVariables[[1]] 
InDeVars <- unlist(ListVariables[2:length(ListVariables)], recursive = TRUE, use.names = FALSE)

model.sel <- model.selection.gwr(DeVar ,InDeVars, data = PolygonShapeTrimmed, 
                                 kernel = "bisquare", adaptive = TRUE, bw = 80) 
sorted.models <- model.sort.gwr(model.sel, numVars = length(InDeVars), 
                                ruler.vector = model.sel[[2]][,2]) 

model.list <- sorted.models[[1]]

model.view.gwr(DeVar, InDeVars, model.list = model.list)

plot(sorted.models[[2]][,2], col = "black", pch = 20, lty = 5, 
     main = "Alternative view of GWR model selection procedure", 
     ylab = "AICc", xlab = "Model number", type = "b")

attach(PolygonShapeTrimmed@data)
modelOLS <- lm(as.formula(sorted.models[[1]][[length(sorted.models[[1]])]][[1]]))
summary(modelOLS)
tidy_modelOLS <- tidy(modelOLS)
write.csv(tidy_modelOLS, "Tables\\ModelOLS.csv")

DependentVariable <- c(ListVariables[[1]])
# MoransIListVariables <- c(ListVariables[[1]])
# MoransIListVariables <- append(MoransIListVariables, "Lat")
# MoransIListVariables <- append(MoransIListVariables, "Long")
PolygonShapeData <- PolygonShapeTrimmed@data
PolygonShapeData.dists <- as.matrix(dist(cbind(PolygonShapeData$Long, PolygonShapeData$Lat)))
PolygonShapeData.dists.inv <- 1/PolygonShapeData.dists
diag(PolygonShapeData.dists.inv) <- 0
Moran.I(PolygonShapeData$E_s11_to_1, PolygonShapeData.dists.inv, na.rm=FALSE)

# gw.ss.bx <- gwss(PolygonShapeTrimmed, PolygonShapeTrimmed, vars = c("E_s11_to_1", "MurderRate", "ProxyGDP", "Poverty"), 
#                  kernel = "boxcar", adaptive = TRUE, bw = 48, quantile = TRUE)
gw.ss.bx <- gwss(PolygonShapeTrimmed, PolygonShapeTrimmed, vars = c(ListVariables), 
                 kernel = "boxcar", adaptive = TRUE, bw = 48, quantile = TRUE)

# gw.ss.bs <- gwss(PolygonShapeTrimmed, PolygonShapeTrimmed, vars = c("E_s11_to_1", "MurderRate", "ProxyGDP", "Poverty"), 
#                  kernel = "bisquare", adaptive = TRUE, bw = 48)
gw.ss.bs <- gwss(PolygonShapeTrimmed, PolygonShapeTrimmed, vars = c(ListVariables), 
                 kernel = "bisquare", adaptive = TRUE, bw = 48)

map.na = list("SpatialPolygonsRescale", layout.north.arrow(),
              offset = c(329000,261500), scale = 4000, col=1)
map.scale.1 = list("SpatialPolygonsRescale", layout.scale.bar(),
                   offset = c(326500,217000), scale = 5000, col=1, fill=c("transparent","blue"))
map.scale.2  = list("sp.text", c(326500,217900), "0", cex=0.9, col=1)
map.scale.3  = list("sp.text", c(331500,217900), "5km", cex=0.9, col=1)
map.layout <- list(map.na,map.scale.1,map.scale.2,map.scale.3)

mypalette.1 <- brewer.pal(9, "Reds") 
mypalette.2 <- brewer.pal(9, "Blues") 
mypalette.3 <- brewer.pal(9, "Greens")

LSD_DepVar <- paste(ListVariables[[1]],"LSD", sep="_")
IQR_DepVar <- paste(ListVariables[[1]],"IQR", sep="_")
CorrDepVarIndVar1 <- paste("Corr",paste(ListVariables[[1]], ListVariables[[2]], sep="."), sep="_")
CorrDepVarIndVar2 <- paste("Corr",paste(ListVariables[[1]], ListVariables[[3]], sep="."), sep="_")
SpeRhoDepVarIndVar1 <- paste("Spearman_rho", paste(ListVariables[[1]], ListVariables[[2]], sep="."), sep="_")
SpeRhoDepVarIndVar3 <- paste("Spearman_rho", paste(ListVariables[[1]], ListVariables[[4]], sep="."), sep="_")

# spplot(gw.ss.bx$SDF, "E_s11_to_1_LSD", key.space = "right", 
#        col.regions = mypalette.1, cuts = 7, 
#        main = "GW standard deviations for SABER11 (basic)", sp.layout = map.layout)
spplot(gw.ss.bx$SDF, LSD_DepVar, key.space = "right", 
       col.regions = mypalette.1, pretty=TRUE, cuts=8,
       main = "GW standard deviations for dependent variable (basic)", sp.layout = map.layout)

# spplot(gw.ss.bx$SDF, "E_s11_to_1_IQR", key.space = "right", 
#        col.regions = mypalette.1, cuts = 7, 
#        main = "GW inter-quartile ranges for SABER11 (robust)",
#        sp.layout = map.layout)
spplot(gw.ss.bx$SDF, IQR_DepVar, key.space = "right", 
       col.regions = mypalette.1, pretty=TRUE, cuts=8,
       main = "GW inter-quartile ranges for dependent variable (robust)",
       sp.layout = map.layout)



# spplot(gw.ss.bx$SDF, "Corr_E_s11_to_1.MurderRate", key.space = "right", 
#        col.regions = mypalette.2, at = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), 
#        main = "GW correlations: SABER11 and MurderRate (box-car kernel)", 
#        sp.layout = map.layout)
spplot(gw.ss.bx$SDF, CorrDepVarIndVar1, key.space = "right", 
       col.regions = mypalette.2, pretty=TRUE, cuts=8,
       main = "GW correlations: Dependent and first independent variable (box-car kernel)", 
       sp.layout = map.layout)



# spplot(gw.ss.bx$SDF, "Corr_E_s11_to_1.ProxyGDP", key.space = "right", 
#        col.regions = mypalette.2, at = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), 
#        main = "GW correlations: SABER11 and GDP (box-car kernel)", 
#        sp.layout = map.layout)
spplot(gw.ss.bx$SDF, CorrDepVarIndVar2, key.space = "right", 
       col.regions = mypalette.2, pretty=TRUE, cuts=8,
       main = "GW correlations: Dependent and second independent variable (box-car kernel)", 
       sp.layout = map.layout)




# spplot(gw.ss.bs$SDF, "Spearman_rho_E_s11_to_1.MurderRate",key.space = "right", 
#        col.regions = mypalette.3, at = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), 
#        main = "GW correlations: SABER11 and Murder rate (robust)", 
#        sp.layout = map.layout)
spplot(gw.ss.bs$SDF, SpeRhoDepVarIndVar1,key.space = "right", 
       col.regions = mypalette.3, pretty=TRUE, cuts=8, 
       main = "GW correlations: Dependent and first independent variable (robust)", 
       sp.layout = map.layout)



# spplot(gw.ss.bs$SDF, "Spearman_rho_E_s11_to_1.Poverty",key.space = "right", 
#        col.regions = mypalette.3, at = c(-0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1), 
#        main = "GW correlations: SABER11 and Poverty (robust)", 
#        sp.layout = map.layout)
spplot(gw.ss.bs$SDF, SpeRhoDepVarIndVar3,key.space = "right", 
       col.regions = mypalette.3, pretty=TRUE, cuts=8, 
       main = "GW correlations: Dependent and third independent variable (robust)", 
       sp.layout = map.layout)

#bw.gwr.1 <- bw.gwr(whatever_theoretical_model_is_selected, data = PolygonShapeTrimmed, approach = "AICc",
#                   kernel = "bisquare", adaptive = TRUE)
bw.gwr.1 <- bw.gwr(as.formula(sorted.models[[1]][[length(sorted.models[[1]])]][[1]]), data = PolygonShapeTrimmed, approach = "AICc",
                   kernel = "bisquare", adaptive = TRUE)
cat("")
cat("The optimal adaptative bandwidth is equal to ", bw.gwr.1)

#gwr.res <- gwr.basic(whatever_theoretical_model_is_selected, data = PolygonShapeTrimmed, bw = bw.gwr.1, 
#                     kernel = "bisquare", adaptive = TRUE, F123.test = TRUE)

gwr.res <- gwr.basic(as.formula(sorted.models[[1]][[length(sorted.models[[1]])]][[1]]), data = PolygonShapeTrimmed, bw = bw.gwr.1, 
                     kernel = "bisquare", adaptive = TRUE, F123.test = TRUE)
gwr.res

writeOGR(gwr.res$SDF, "Data\\Colombia", driver="ESRI Shapefile", layer='Results', overwrite_layer = TRUE)
names(gwr.res$SDF) 

#spplot(gwr.res$SDF, MurderRate, key.space = "right", 
#       col.regions = mypalette.2, pretty=TRUE, cuts=8,
#       main = "Basic GW regression coefficient estimates for murder rates",
#       sp.layout = map.layout)
spplot(gwr.res$SDF, ListVariables[[2]], key.space = "right", 
       col.regions = mypalette.2, pretty=TRUE, cuts=8,
       main = "Basic GW regression coefficient estimates for the first independent variable",
       sp.layout = map.layout)

spplot(gwr.res$SDF, ListVariables[[3]], key.space = "right", 
       col.regions = mypalette.2, pretty=TRUE, cuts=8,
       main = "Basic GW regression coefficient estimates for the second independent variable",
       sp.layout = map.layout)

spplot(gwr.res$SDF, ListVariables[[4]], key.space = "right", 
       col.regions = mypalette.2, pretty=TRUE, cuts=8,
       main = "Basic GW regression coefficient estimates for the third independent variable",
       sp.layout = map.layout)

#rgwr.res <- gwr.robust(as.formula(sorted.models[[1]][[length(sorted.models[[1]])]][[1]]), data = PolygonShapeTrimmed, bw = bw.gwr.1,
#                       kernel = "bisquare", adaptive = TRUE, F123.test = TRUE)
rgwr.res <- gwr.robust(as.formula(sorted.models[[1]][[length(sorted.models[[1]])]][[1]]), data = PolygonShapeTrimmed, bw = bw.gwr.1,
                       kernel = "bisquare", adaptive = TRUE, F123.test = TRUE)
rgwr.res

spplot(rgwr.res$SDF, ListVariables[[2]], key.space = "right", 
       col.regions = mypalette.1, pretty=TRUE, cuts=8,
       main = "Robust GW regression coefficient estimates for the first independent variable",
       sp.layout = map.layout)

spplot(rgwr.res$SDF, ListVariables[[3]], key.space = "right", 
       col.regions = mypalette.2, pretty=TRUE, cuts=8,
       main = "Robust GW regression coefficient estimates for the second independent variable",
       sp.layout = map.layout)

spplot(rgwr.res$SDF, ListVariables[[4]], key.space = "right", 
       col.regions = mypalette.3, pretty=TRUE, cuts=8,
       main = "Robust GW regression coefficient estimates for the third independent variable",
       sp.layout = map.layout)

writeOGR(rgwr.res$SDF, "Data\\Colombia", driver="ESRI Shapefile", layer='Results - Robust', overwrite_layer = TRUE)
