# Chapter 4

############################################
# Chapter 4: Estimated spatial effects
############################################

rm(list=ls())
load("data/workspaces/20220323_models_noexp_pred.RData")



# rezultati
# povprecje
# needs to be standardized!!

boundaryUE$randomeff.S.UE.small<- model.besag.S.UE.small.C.NA$summary.random$idUE[,"mean"]
boundaryUE$randomeff.HS.UE.small <- model.besag.HS.UE.small.C.NA$summary.random$idUE[,"mean"]

boundaryOB$randomeff.S.OB.small<- model.besag.S.OB.small.C.NA$summary.random$idOB[,"mean"]
boundaryOB$randomeff.HS.OB.small <- model.besag.HS.OB.small.C.NA$summary.random$idOB[,"mean"]


# sd

boundaryUE$randomeff.sd.S.UE.small <- model.besag.S.UE.small.C.NA$summary.random$idUE[,"sd"]
boundaryUE$randomeff.sd.HS.UE.small <- model.besag.HS.UE.small.C.NA$summary.random$idUE[,"sd"]

boundaryOB$randomeff.sd.S.OB.small <- model.besag.S.OB.small.C.NA$summary.random$idOB[,"sd"]
boundaryOB$randomeff.sd.HS.OB.small <- model.besag.HS.OB.small.C.NA$summary.random$idOB[,"sd"]


p_S_UE_mean <- ggplot(boundaryUE) + geom_sf(aes(fill = randomeff.S.UE.small)) +
  scale_fill_distiller(palette = "RdYlBu", guide = "colourbar", limits = c(-0.5, 0.5)) + theme_bw() + 
  # geom_text(aes(CEN_E, CEN_N, label = UE_ID)) + 
  ggtitle("Model GS - upravne enote") +
  theme(legend.position = "none")

p_HS_UE_mean <- ggplot(boundaryUE) + geom_sf(aes(fill = randomeff.HS.UE.small)) +
  scale_fill_distiller(palette = "RdYlBu", guide = "colourbar", limits = c(-0.5, 0.5), name = "Prostorski vpliv") + theme_bw() + 
  ggtitle("Model GHS - upravne enote")
#scale_fill_continuous(name = "Prostorski vpliv")

p_S_OB_mean <- ggplot(boundaryOB) + geom_sf(aes(fill = randomeff.S.OB.small)) +
  scale_fill_distiller(palette = "RdYlBu", guide = "colourbar", limits = c(-0.75, 0.75)) + theme_bw() + 
  # geom_text(aes(CEN_E, CEN_N, label = OB_ID)) + 
  ggtitle("Model GS - obèine") +
  theme(legend.position = "none")

p_HS_OB_mean <- ggplot(boundaryOB) + geom_sf(aes(fill = randomeff.HS.OB.small)) +
  scale_fill_distiller(palette = "RdYlBu", guide = "colourbar", limits = c(-0.75, 0.75),  name = "Prostorski vpliv") + theme_bw() + 
  ggtitle("Model GHS - obèine")




(p_S_UE_mean + p_HS_UE_mean) / (p_S_OB_mean + p_HS_OB_mean)




p_S_UE_sd <- ggplot(boundaryUE) + geom_sf(aes(fill = randomeff.sd.S.UE.small)) +
  scale_fill_distiller(palette = "Blues", guide = "colourbar", direction = 1, limits = c(0, 0.5)) + theme_bw() + 
  # geom_text(aes(CEN_E, CEN_N, label = UE_ID)) + 
  ggtitle("Model GS - upravne enote") +
  theme(legend.position = "none")

p_HS_UE_sd <- ggplot(boundaryUE) + geom_sf(aes(fill = randomeff.sd.HS.UE.small)) +
  scale_fill_distiller(palette = "Blues", guide = "colourbar", direction = 1, limits = c(0, 0.5), name = "sd") + theme_bw() + 
  ggtitle("Model GHS UE - upravne enote")

p_S_OB_sd <- ggplot(boundaryOB) + geom_sf(aes(fill = randomeff.sd.S.OB.small)) +
  scale_fill_distiller(palette = "Blues", guide = "colourbar", direction = 1, limits = c(0, 0.55)) + theme_bw() + 
  # geom_text(aes(CEN_E, CEN_N, label = OB_ID)) + 
  ggtitle("Model GS OB - obèine") +
  theme(legend.position = "none")

p_HS_OB_sd <- ggplot(boundaryOB) + geom_sf(aes(fill = randomeff.sd.HS.OB.small)) +
  scale_fill_distiller(palette = "Blues", guide = "colourbar", direction = 1, limits = c(0, 0.55), name = "sd") + theme_bw() + 
  ggtitle("Model GHS OB - obèine")




(p_S_UE_sd + p_HS_UE_sd) / (p_S_OB_sd + p_HS_OB_sd)

library('RColorBrewer')
library(INLA)
library(fields)
PlotSpdeField = function(mesh,sample.field1, param, zlim =NULL){
  proj <- inla.mesh.projector(mesh,dims=c(300, 300))
  field.proj1 <- inla.mesh.project(proj, sample.field1)
  
  poly1 <- sp::Polygon(sloveniaBoundary)
  firstPoly <- sp::Polygons(list(poly1), ID = "A")
  firstSpatialPoly <- sp::SpatialPolygons(list(firstPoly))
  ov <- sp::over(SpatialPoints(proj$lattice$loc), firstSpatialPoly)
  field.proj1[is.na(ov)] <- NA
  
  if (param == "povpr"){
    cols = rev(brewer.pal(11,'RdYlBu'))
    rf <- colorRampPalette(cols)   # make colors
    r <- rf(64)
    image.plot(list(x = proj$x, y=proj$y, z = (field.proj1)),col = r, zlim=zlim, xlim = c(370, 620), ylim = c(30, 195))#, ylim = c(0,1), xlim= c(0,1))
  }
  if (param == "sd") {
    cols = (brewer.pal(9,'Blues'))
    rf <- colorRampPalette(cols)   # make colors
    r <- rf(64)
    image.plot(list(x = proj$x, y=proj$y, z = (field.proj1)),col = r, zlim=zlim, xlim = c(370, 620), ylim = c(30, 195))#, ylim = c(0,1), xlim= c(0,1))
    
  }
  
}


par(mfrow = c(2,2))
PlotSpdeField(mesh,model.SPDE.S.small.NA$summary.random$fieldID$mean, param = "povpr", zlim =NULL)
# plot(mesh, add = T)
# points(phenoData.small$xpos, phenoData.small$ypos, col = rgb(0,0,0,0.7))
# plot(mapSlo, add=TRUE)
# 
points(sloveniaBoundary, type = "l")

PlotSpdeField(mesh,model.SPDE.HS.small.NA$summary.random$fieldID$mean, param = "povpr", zlim =NULL)
# plot(mesh, add = T)
# points(phenoData.small$xpos, phenoData.small$ypos, col = rgb(0,0,0,0.7))
# plot(mapSlo, add=TRUE)
# 
points(sloveniaBoundary, type = "l")

PlotSpdeField(mesh,model.SPDE.S.small.NA$summary.random$fieldID$sd, param = "sd", zlim =NULL)
# points(phenoData.small$xpos, phenoData.small$ypos, col = rgb(0,0,0,0.7))
points(sloveniaBoundary, type = "l")




PlotSpdeField(mesh,model.SPDE.HS.small.NA$summary.random$fieldID$sd, param = "sd", zlim =NULL)
# points(phenoData.small$xpos, phenoData.small$ypos, col = rgb(0,0,0,0.7))
points(sloveniaBoundary, type = "l")






par(mfrow=c(1,1))

PlotSpdeField(mesh,model.SPDE.HS.small.NA$summary.random$fieldID$mean, param = "povpr", zlim =NULL)
# plot(mesh, add = T)
# points(phenoData.small$xpos, phenoData.small$ypos, col = rgb(0,0,0,0.7))
# plot(mapSlo, add=TRUE)
# 
points(sloveniaBoundary, type = "l")



PlotSpdeField(mesh,model.SPDE.HS.small.NA$summary.random$fieldID$sd, param = "sd", zlim =NULL)
# points(phenoData.small$xpos, phenoData.small$ypos, col = rgb(0,0,0,0.7))
points(sloveniaBoundary, type = "l")

