# RUN MODELS (no exp)

# Spatial modelling in animal breeding
# Eva Lavrencic

library(INLA)
library(sf) #spatial
library(spdep) #spatial
library(rgdal) #spatial
library(tidyverse)
library(Matrix)
library(kableExtra)
library(gdata) # drop unused factor levels
library(patchwork) # combine plots
library(naniar) # missing values
library(ggpubr)
library(ggplot2)
# 
library(devtools)
library(RcppArmadillo)
find_rtools() # TRUE
has_devel()

############################################
# import data
############################################

rm(list = ls())
load("data/data_orig.RData")
load("data/data_prostor.RData")

cols_analysis <- c("idlok", "OB", "UE", "idziva", "idziv", "idosbocn", "datocen", "last004", 
                   "vpliv1", "vpliv2", "vpliv3", "last004scaled", "rowNumberAinv", "xpos", "ypos", "KO")

# phenoData <- phenoData[,cols_analysis]
# phenoData.small <- phenoData.small[,cols_analysis]


# phenoData.small$index <- as.numeric(row.names(phenoData.small))
phenoData.small$index = as.numeric(rownames(phenoData.small))

phenoData.small <- phenoData.small[order(phenoData.small$index), ]

head(phenoData.small)

mapUE <- mapUE[order(mapUE$UE_ID),]
mapOB <- mapOB[order(mapOB$OB_ID),]

boundaryUE <- boundaryUE[order(boundaryUE$UE_ID),]
boundaryOB <- boundaryOB[order(boundaryOB$OB_ID),]
############################################
# prepare graph of regions
############################################

# upravne enote

nb.mapUE <- poly2nb(mapUE) #class nb
nb2INLA("map.graphUE",nb.mapUE)
gUE <- inla.read.graph(filename = "map.graphUE")
# gUE
# n = 58123

# nnbs = number of neighbours
# nbs = lists of neighbours

#obcine
nb.mapOB <- poly2nb(mapOB)
nb2INLA("map.graphOB",nb.mapOB)
gOB <- inla.read.graph(filename = "map.graphOB")




############################################
# prepare region ids from graph (regions)
############################################

# prepare region ids
boundaryUE <- boundaryUE[order(boundaryUE$UE_ID),]
boundaryOB <- boundaryOB[order(boundaryOB$OB_ID),]

# grafID_UE <- data.frame("UE" = mapUE[,"UE_MID"], "idUE" = mapUE[,"UE_ID"])

grafID_UE <- data.frame("UE" = mapUE[,"UE_MID"], "idUE" = 1:nrow(mapUE))
colnames(grafID_UE) <- c("UE", "idUE")
grafID_OB <- data.frame("OB" = mapOB[,"OB_MID"], "idOB" = 1:nrow(mapOB))
# grafID_OB <- data.frame("OB" = mapOB[,"OB_MID"], "idOB" = mapOB[,"OB_ID"])
colnames(grafID_OB) <- c("OB", "idOB")
head(grafID_OB)
# merge region ids with phenoData
phenoData.small <- as.data.frame(merge(phenoData.small, grafID_UE, all.x=T))
phenoData.small <- as.data.frame(merge(phenoData.small, grafID_OB, all.x=T))

phenoData <- as.data.frame(merge(phenoData, grafID_UE, all.x=T))
phenoData <- as.data.frame(merge(phenoData, grafID_OB, all.x=T))

############################################
# define priors
############################################

hyperVarPed = list(theta = list(prior="pc.prec", param=c(sqrt(0.1),0.5)))
hyperVarIdlok.GH = list(theta = list(prior="pc.prec", param=c(sqrt(0.25),0.5)))
hyperVarIdlok.GHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))
hyperResVarG = list(theta = list(prior="pc.prec", param=c(sqrt(0.3),0.5)))
hyperResVarGHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))


# BESAG (scaled)
hyperVarBesag.S = list(theta = list(prior="pc.prec", param=c(sqrt(0.25),0.5)))
hyperVarBesag.HS = list(theta = list(prior="pc.prec", param=c(sqrt(0.10),0.5)))


# SPDE
hyperRange  = c(50, 0.8)
hyperVarSpdeS = c(sqrt(0.25), 0.5)
hyperVarSpdeHS = c(sqrt(0.10), 0.5)



############################################
# prepare data from last 2 years (NA)
############################################
phenoData.small <- phenoData.small[order(phenoData.small$index),]
phenoData.small$year = gsub("^.*/", "", phenoData.small$datocen)
phenoData.small$year = as.numeric(phenoData.small$year)

hist(phenoData.small$year) # take year 2018 and 2019
sum(phenoData.small$year == 2019 | phenoData.small$year == 2018) # 115 records
sum(is.na(phenoData.small$last004scaled)) # 50 missing phenotype data overall
# 5 missing phenotypes for last 3 years
sum(is.na(phenoData.small$last004scaled[phenoData.small$year == 2019 | phenoData.small$year == 2018 | phenoData.small$year == 2017]))

# years 2017, 2018, 2019

val_set <- phenoData.small$year == 2017 | phenoData.small$year == 2018 | phenoData.small$year == 2019
sum(val_set) # 246 records
year_pred <- which(val_set)

tmp = table(phenoData.small$idlok, val_set) # COLSUMS 3554, 246

# OZNACI PO KATEGORIJAH 2 IN 3
# naredi korelacijo po kategorijah
sum(tmp[,1] > 0 & tmp[,2] == 0) # 1625
sum(tmp[,1] == 0 & tmp[,2] > 0) # 48 niso v trening, so v napovednem
sum(tmp[,1] > 0 & tmp[,2] > 0) # 165 so v trening in napovednem




summary(phenoData.small$last004scaled)

phenoData.small$last004scaled[which.max(phenoData.small$last004scaled)] = NA

summary(phenoData.small$last004scaled)

############################################
# define model matrix and stack (for SPDE)
############################################

phenoData.small <- phenoData.small[order(phenoData.small$index),]

M = model.matrix(~1 + vpliv1 + vpliv2 + vpliv3, phenoData.small)
# Store in a stack
FactorNames =  tail(colnames(M), ncol(M)-1)
FactorStack = paste0( "list(", FactorNames, " = M[,", 2:(ncol(M)), "]),", collapse = "" )
FactorStack = ( gsub('.{1}$', '', FactorStack ) )
FactorStack= paste0("c(", FactorStack, ")")

# Make the formula
FormulaStack = paste("+", FactorNames, collapse = " ")  

formulaS.1 = 'last004scaled ~ intercept + f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + f(fieldID, model = spdeStatS) -1'

formulaHS.1 = 'last004scaled ~intercept + f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + f(idlok, model = "iid", hyper = hyperVarIdlok.GHS) + f(fieldID, model = spdeStatHS) -1'

formulaS = as.formula( paste0(formulaS.1, FormulaStack))
formulaHS = as.formula( paste0(formulaHS.1, FormulaStack))

# Make mesh and SPDE 
mesh = inla.mesh.2d(cbind(phenoData.small$xpos, phenoData.small$ypos), max.edge=c(10, 20), cutoff = 2.5, offset = 30)
A = inla.spde.make.A(mesh = mesh, loc = cbind(phenoData.small$xpos, phenoData.small$ypos) )
spdeStatS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatS$n.spde) 
spdeStatHS = inla.spde2.pcmatern(mesh = mesh, alpha = 2 ,  prior.range = hyperRange, prior.sigma = hyperVarSpdeS)
meshIndexHS = inla.spde.make.index(name = "fieldID", n.spde = spdeStatHS$n.spde) 

# Make stack
stack = inla.stack(data = list(last004scaled = phenoData.small$last004scaled),
                   A = list(A,1, 1),
                   effects = list(c(meshIndexS, list(intercept = 1)),
                                  c(list(rowNumberAinv = phenoData.small$rowNumberAinv),list(idlok = phenoData.small$idlok)),
                                  c(eval(parse(text = FactorStack ) ) )), tag = "est") 

stackHS = inla.stack(data = list(last004scaled = phenoData.small$last004scaled),
                     A = list(A,1, 1),
                     effects = list(c(meshIndexHS, list(intercept = 1)),
                                    c(list(rowNumberAinv = phenoData.small$rowNumberAinv),list(idlok = phenoData.small$idlok)),
                                    c(eval(parse(text = FactorStack ) ) ) ), tag = "est") 





############################################
# define other models
############################################

# base models
formulaB.G = last004scaled ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed)
# H MODEL (G + H: herd effect)
formulaB.H = last004scaled ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idlok, model = "iid", hyper = hyperVarIdlok.GH)

# S MODEL (G + field (regional effect)) BESAG
# UE = administrative regions
# OB = municipalities
formulaB.S.UE.C = last004scaled ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idUE, model = "besag", graph = gUE, scale.model = TRUE, hyper = hyperVarBesag.S)

formulaB.S.OB.C = last004scaled ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idOB, model = "besag", graph = gOB, scale.model = TRUE, hyper = hyperVarBesag.S)


# HS MODEL (G + location + field) BESAG
formulaB.HS.UE.C <- last004scaled ~ 1 + vpliv1 + vpliv2 + vpliv3 + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idlok, model = "iid", hyper = hyperVarIdlok.GHS) + 
  f(idUE, model = "besag", graph = gUE, scale.model = TRUE, hyper = hyperVarBesag.HS)

formulaB.HS.UE.C_full <- last004scaled ~ 1 + vpliv1 + vpliv2 + vpliv3 + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idlok, model = "iid", hyper = hyperVarIdlok.GHS) + 
  f(idUE, model = "besag", graph = gUE, scale.model = TRUE, hyper = hyperVarBesag.HS)


formulaB.HS.OB.C <- last004scaled ~ 1 + vpliv1 + vpliv2 + vpliv3 + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idlok, model = "iid", hyper = hyperVarIdlok.GHS) + 
  f(idOB, model = "besag", graph = gOB, scale.model = TRUE, hyper = hyperVarBesag.HS)



############################################
# run models
############################################
#phenoData.small <- phenoData.small[order(phenoData.small$rowNumberAinv),]
#phenoData.small <- phenoData.small[order(phenoData.small$idlok),]


# phenoData.small$index <- as.numeric(row.names(phenoData.small))
phenoData.small.full <- phenoData.small[order(phenoData.small$index), ]


model.G.small.NA.full <- inla(formulaB.G,  family = "gaussian", data = phenoData.small, 
                         control.compute = list(dic=T,cpo=F, openmp.strategy="huge", config=T), verbose=T)
model.H.small.NA.full <- inla(formulaB.H,  family = "gaussian", data = phenoData.small, 
                         control.compute = list(dic=T,cpo=F, openmp.strategy="huge", config=T), verbose=T)

model.besag.S.UE.small.C.NA.full <- inla(formulaB.S.UE.C,  family = "gaussian", data = phenoData.small, 
                                    control.compute = list(dic=T,cpo=F, openmp.strategy="huge", config=T), verbose=T)

model.besag.HS.UE.small.C.NA.full <- inla(formulaB.HS.UE.C,  family = "gaussian", data = phenoData.small, 
                                     control.compute = list(dic=T,cpo=F, openmp.strategy="huge", config=T), verbose=T)

model.besag.HS.UE.small.C.NA.full <- inla(formulaB.HS.UE.C_full,  family = "gaussian", data = phenoData.small, 
                                     control.compute = list(dic=T,cpo=F, openmp.strategy="huge", config=T), verbose=T)


model.besag.S.OB.small.C.NA.full <- inla(formulaB.S.OB.C,  family = "gaussian", data = phenoData.small, 
                                    control.compute = list(dic=T,cpo=F, openmp.strategy="huge", config=T), verbose=T)

model.besag.HS.OB.small.C.NA.full <- inla(formulaB.HS.OB.C,  family = "gaussian", data = phenoData.small, 
                                     control.compute = list(dic=T,cpo=F, openmp.strategy="huge", config=T), verbose=T)

# SPDE
model.SPDE.S.small.NA.full = inla(formula = formulaS, data = inla.stack.data(stack),
                             family = "normal", control.predictor =list(A=inla.stack.A(stack),compute = T),
                             control.family=list(list(hyper=hyperResVarGHS)),
                             control.compute = list(dic=T,cpo=F, config=T), verbose=T) 
model.SPDE.HS.small.NA.full = inla(formula = formulaHS, data = inla.stack.data(stackHS),
                               family = "normal", control.predictor =list(A=inla.stack.A(stackHS),compute = T),
                               control.family=list(list(hyper=hyperResVarGHS)),
                               control.compute = list(dic=T,cpo=F, config=T), verbose=T)



save.image("data/workspaces/20220404_models_noexp_full.RData")


