# RUN MODELS (exp)

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


load("data/workspaces/20220323_chapter3_data.RData")


############################################
# calculate distances between locations
############################################
Rcpp::sourceCpp("data/other/euclideandistance.cpp")
# DISTANCES==
rhoHS <- 27.5
rhoS <- 10
# Exponential distances
# small data
fun_distances <- function(data, rho, model = "S"){
  data <- data[order(data$idlok),] # *
  koordinate <- unique(paste(data$xpos, data$ypos))
  koordinate <- matrix(as.numeric(unlist(strsplit(koordinate, split=" "))), ncol=2, byrow = T)
  paired.distances <- fastPairDist(koordinate, koordinate)
  paired.distances.exp <- exp(-paired.distances/rho)
  Exp.inv <- solve(paired.distances.exp) # inverz matrike
  if (model == "S"){
    data$idExpS <- as.numeric(as.factor(data$idlok))
  }
  if (model == "HS"){
    data$idExpHS <- as.numeric(as.factor(data$idlok))
  }

  list(data = data, dist = paired.distances.exp, dist.inv = Exp.inv)
}

dist.small <- fun_distances(phenoData.small, rhoS, "S")
phenoData.small <- dist.small[[1]]
Exp.inv.small.S <- dist.small[[3]]

dist.small <- fun_distances(phenoData.small, rhoS, "HS")
phenoData.small <- dist.small[[1]]
Exp.inv.small.HS <- dist.small[[3]]


phenoData.small <- phenoData.small[order(phenoData.small$index), ]


############################################
# define priors
############################################
hyperVarPed = list(theta = list(prior="pc.prec", param=c(sqrt(0.1),0.5)))
hyperVarIdlok.GH = list(theta = list(prior="pc.prec", param=c(sqrt(0.25),0.5)))
hyperVarIdlok.GHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))
hyperResVarG = list(theta = list(prior="pc.prec", param=c(sqrt(0.3),0.5)))
hyperResVarGHS = list(theta = list(prior="pc.prec", param=c(sqrt(0.15),0.5)))

# Exp
hyperVarExp.S = list(theta = list(prior="pc.prec", param=c(sqrt(0.25),0.5)))
hyperVarExp.HS = list(theta = list(prior="pc.prec", param=c(sqrt(0.10),0.5)))



formulaB.S.Exp.prior.small = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idExpS, model = "generic0", Cmatrix = Exp.inv.small.S, hyper = hyperVarExp.S)
formulaB.HS.Exp.prior.small <- last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3 + 
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed) + 
  f(idlok, model = "iid", hyper = hyperVarIdlok.GHS) + 
  f(idExpHS, model = "generic0", Cmatrix = Exp.inv.small.HS, hyper = hyperVarExp.HS)




model.Exp.S.prior.small.NA <- inla(formulaB.S.Exp.prior.small,  family = "gaussian", 
                                   data = phenoData.small, control.compute = list(dic=T,cpo=F,config=T), verbose=T)
model.Exp.HS.prior.small.NA <- inla(formulaB.HS.Exp.prior.small,  family = "gaussian", 
                                    data = phenoData.small, control.compute = list(dic=T,cpo=F,config=T), verbose=T)







