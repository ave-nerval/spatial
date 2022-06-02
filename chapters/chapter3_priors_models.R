# Chapter 3 - change priors - models

rm(list=ls())
load("data/workspaces/20220323_chapter3_data.RData")


hyperVarPed01 = list(theta = list(prior="pc.prec", param=c(sqrt(0.1),0.5)))
hyperVarPed03 = list(theta = list(prior="pc.prec", param=c(sqrt(0.3),0.5)))
hyperVarPed05 = list(theta = list(prior="pc.prec", param=c(sqrt(0.5),0.5)))
hyperVarPed07 = list(theta = list(prior="pc.prec", param=c(sqrt(0.7),0.5)))
hyperVarPed09 = list(theta = list(prior="pc.prec", param=c(sqrt(0.9),0.5)))

formulaB.G = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv)

formulaB.G.01 = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed01)

formulaB.G.03 = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed03)

formulaB.G.05 = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed05)

formulaB.G.07 = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed07)

formulaB.G.09 = last004scaledNA ~ 1 + vpliv1 + vpliv2 + vpliv3  +  
  f(rowNumberAinv, model = "generic0", Cmatrix = Ainv,hyper = hyperVarPed09)

model.G.small <- inla(formulaB.G,  family = "gaussian", data = phenoData.small, 
                      control.compute = list(dic=T,cpo=F, openmp.strategy="huge", 
                                             config=T), verbose=T)

model.G.small.001 <- inla(formulaB.G.001,  family = "gaussian", data = phenoData.small, 
                          control.compute = list(dic=T,cpo=F, openmp.strategy="huge", 
                                                 config=T), verbose=T)

model.G.small.01 <- inla(formulaB.G.01,  family = "gaussian", data = phenoData.small, 
                         control.compute = list(dic=T,cpo=F, openmp.strategy="huge", 
                                                config=T), verbose=T)

model.G.small.03 <- inla(formulaB.G.03,  family = "gaussian", data = phenoData.small, 
                         control.compute = list(dic=T,cpo=F, openmp.strategy="huge", 
                                                config=T), verbose=T)

model.G.small.05 <- inla(formulaB.G.05,  family = "gaussian", data = phenoData.small, 
                         control.compute = list(dic=T,cpo=F, openmp.strategy="huge", 
                                                config=T), verbose=T)

model.G.small.07 <- inla(formulaB.G.07,  family = "gaussian", data = phenoData.small, 
                         control.compute = list(dic=T,cpo=F, openmp.strategy="huge", 
                                                config=T), verbose=T)

model.G.small.09 <- inla(formulaB.G.09,  family = "gaussian", data = phenoData.small, 
                         control.compute = list(dic=T,cpo=F, openmp.strategy="huge", 
                                                config=T), verbose=T)

