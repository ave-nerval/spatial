# Chapter 2

# rm(list=ls())
# load("data/workspaces/20220323_models_exp_pred.RData")
# load("data/workspaces/20220323_models_noexp_pred.RData")


############################################
# Chapter 2: DIC values and var comp
############################################

dic_small = c(model.G.small.NA$dic$dic,
              model.H.small.NA$dic$dic,
              model.besag.S.OB.small.C.NA$dic$dic,
              model.besag.S.UE.small.C.NA$dic$dic,
              model.besag.HS.OB.small.C.NA$dic$dic,
              model.besag.HS.UE.small.C.NA$dic$dic,
              model.SPDE.S.small.NA$dic$dic,
              model.SPDE.HS.small.NA$dic$dic,
              model.Exp.S.prior.small.NA$dic$dic,
              model.Exp.HS.prior.small.NA$dic$dic)

dic_names = c('model.G.small',
              'model.H.small',
              'model.besag.S.OB.small',
              'model.besag.S.UE.small',
              'model.besag.HS.OB.small',
              'model.besag.HS.UE.small',
              'model.SPDE.S.small',
              'model.SPDE.HS.small',
              'model.Exp.S.small',
              'model.Exp.HS.small')


dic = data.frame(name = dic_names, value = dic_small)

# write.csv2(dic, file='dic_values_pred_new_0323.csv')

### variance components


var_comp_model <- function(chosen_model, var_comp, region = "UE"){
  
  if(var_comp=="gen"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for rowNumberAinv`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp[c("quant0.025", "mean", "quant0.5" ,"quant0.975")]), 2)
  }
  
  else if(var_comp=="herd"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idlok`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp[c("quant0.025", "mean", "quant0.5" ,"quant0.975")]), 2)
  }
  
  else if(var_comp=="res"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for the Gaussian observations`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp[c("quant0.025", "mean", "quant0.5" ,"quant0.975")]), 2)
    
  }
  
  else if(var_comp == "exp"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idExp`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp[c("quant0.025", "mean", "quant0.5" ,"quant0.975")]), 2)
    
  }
  
  else if(var_comp=="spatial"){
    
    if(region=="UE"){
      var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idUE`)
      var_comp <- inla.zmarginal(var_comp)
      var_comp <- round(unlist(var_comp[c("quant0.025", "mean", "quant0.5" ,"quant0.975")]), 2)
      
    }
    
    else if(region=="OB"){
      var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idOB`)
      var_comp <- inla.zmarginal(var_comp)
      var_comp <- round(unlist(var_comp[c("quant0.025", "mean", "quant0.5" ,"quant0.975")]), 2)
    }
    
    
    
    
  }
  var_comp <- paste(var_comp, collapse= " ", sep = ", ")
  var_comp
  
}





var_comp_model_brezIZ <- function(chosen_model, var_comp, region = "UE"){
  
  if(var_comp=="gen"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for rowNumberAinv`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp["mean"]), 2)
  }
  
  else if(var_comp=="herd"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idlok`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp["mean"]), 2)
  }
  
  else if(var_comp=="res"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for the Gaussian observations`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp["mean"]), 2)
    
  }
  
  else if(var_comp == "exp"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idExp`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp["mean"]), 2)
    
  }
  
  else if(var_comp=="spatial"){
    
    if(region=="UE"){
      var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idUE`)
      var_comp <- inla.zmarginal(var_comp)
      var_comp <- round(unlist(var_comp["mean"]), 2)
      
    }
    
    else if(region=="OB"){
      var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for idOB`)
      var_comp <- inla.zmarginal(var_comp)
      var_comp <- round(unlist(var_comp["mean"]), 2)
    }
    
    
    
    
  }
  var_comp <- paste(var_comp, collapse= " ", sep = ", ")
  var_comp
  
}





G.small <- c(var_comp_model_brezIZ(chosen_model = model.G.small.NA, var_comp = "gen"),
             "/",
             "/",
             var_comp_model_brezIZ(chosen_model = model.G.small.NA, var_comp = "res"),
             round(model.G.small.NA$dic$dic,0))

GH.small <- c(var_comp_model_brezIZ(chosen_model = model.H.small.NA, var_comp = "gen"), 
              var_comp_model_brezIZ(chosen_model = model.H.small.NA, var_comp = "herd"),
              "/",
              var_comp_model_brezIZ(chosen_model = model.H.small.NA, var_comp = "res"),
              round(model.H.small.NA$dic$dic, 0))



GS.besag.OB.small.C <- c(var_comp_model_brezIZ(chosen_model = model.besag.S.OB.small.C.NA, var_comp = "gen"),
                         "/",
                         var_comp_model_brezIZ(chosen_model = model.besag.S.OB.small.C.NA, var_comp = "spatial", region = "OB"),
                         var_comp_model_brezIZ(chosen_model = model.besag.S.OB.small.C.NA, var_comp = "res"),
                         round(model.besag.S.OB.small.C.NA$dic$dic, 0))

GS.besag.UE.small.C <- c(var_comp_model_brezIZ(chosen_model = model.besag.S.UE.small.C.NA, var_comp = "gen"),
                         "/",
                         var_comp_model_brezIZ(chosen_model = model.besag.S.UE.small.C.NA, var_comp = "spatial"),
                         var_comp_model_brezIZ(chosen_model = model.besag.S.UE.small.C.NA, var_comp = "res"),
                         round(model.besag.S.UE.small.C.NA$dic$dic, 0))



GHS.besag.OB.small.C <- c(var_comp_model_brezIZ(chosen_model = model.besag.HS.OB.small.C.NA, var_comp = "gen"),
                          var_comp_model_brezIZ(chosen_model = model.besag.HS.OB.small.C.NA, var_comp = "herd"),
                          var_comp_model_brezIZ(chosen_model = model.besag.HS.OB.small.C.NA, var_comp = "spatial", region="OB"),
                          var_comp_model_brezIZ(chosen_model = model.besag.HS.OB.small.C.NA, var_comp = "res"),
                          round(model.besag.HS.OB.small.C.NA$dic$dic, 0))

GHS.besag.UE.small.C <- c(var_comp_model_brezIZ(chosen_model = model.besag.HS.UE.small.C.NA, var_comp = "gen"),
                          var_comp_model_brezIZ(chosen_model = model.besag.HS.UE.small.C.NA, var_comp = "herd"),
                          var_comp_model_brezIZ(chosen_model = model.besag.HS.UE.small.C.NA, var_comp = "spatial"),
                          var_comp_model_brezIZ(chosen_model = model.besag.HS.UE.small.C.NA, var_comp = "res"),
                          round(model.besag.HS.UE.small.C.NA$dic$dic, 0))

# Exp


GS.Exp.prior.small <- c(var_comp_model_brezIZ(chosen_model = model.Exp.S.prior.small.NA, var_comp = "gen"),
                        "/",
                        var_comp_model_brezIZ(chosen_model = model.Exp.S.prior.small.NA, var_comp = "exp"),
                        var_comp_model_brezIZ(chosen_model = model.Exp.S.prior.small.NA, var_comp = "res"),
                        round(model.Exp.S.prior.small.NA$dic$dic, 0))

GHS.Exp.prior.small <- c(var_comp_model_brezIZ(chosen_model = model.Exp.HS.prior.small.NA, var_comp = "gen"),
                         var_comp_model_brezIZ(chosen_model = model.Exp.HS.prior.small.NA, var_comp = "herd"),
                         var_comp_model_brezIZ(chosen_model = model.Exp.HS.prior.small.NA, var_comp = "exp"),
                         var_comp_model_brezIZ(chosen_model = model.Exp.HS.prior.small.NA, var_comp = "res"),
                         round(model.Exp.HS.prior.small.NA$dic$dic, 0))

# SPDE
GS.SPDE.small <- c(var_comp_model_brezIZ(chosen_model = model.SPDE.S.small.NA, var_comp = "gen"),
                   "/",
                   paste(round(c(model.SPDE.S.small.NA$summary.hyperpar[4,3]^2, model.SPDE.S.small.NA$summary.hyperpar[4,1]^2,
                                 model.SPDE.S.small.NA$summary.hyperpar[4,4]^2, model.SPDE.S.small.NA$summary.hyperpar[4,5]^2), 2), 
                         collapse= " "),
                   var_comp_model_brezIZ(chosen_model = model.SPDE.S.small.NA, var_comp = "res"),
                   round(model.SPDE.S.small.NA$dic$dic,0))

GHS.SPDE.small <- c(var_comp_model_brezIZ(chosen_model = model.SPDE.HS.small.NA, var_comp = "gen"),
                    var_comp_model_brezIZ(chosen_model = model.SPDE.HS.small.NA, var_comp = "herd"),
                    paste(round(c(model.SPDE.HS.small.NA$summary.hyperpar[5,3]^2, model.SPDE.HS.small.NA$summary.hyperpar[5,1]^2,
                                  model.SPDE.HS.small.NA$summary.hyperpar[5,4]^2, model.SPDE.HS.small.NA$summary.hyperpar[5,5]^2), 2), 
                          collapse= " "),
                    var_comp_model_brezIZ(chosen_model = model.SPDE.HS.small.NA, var_comp = "res"),
                    round(model.SPDE.HS.small.NA$dic$dic,0))



Modeli <- cbind(G.small, GS.besag.OB.small.C, GS.besag.UE.small.C, GS.SPDE.small, GS.Exp.prior.small, GH.small, GHS.besag.OB.small.C, GHS.besag.UE.small.C,
                GHS.SPDE.small, GHS.Exp.prior.small)


rownames(Modeli) <- c("genetic effect", "herd effect", "spatial effect", "residual", "DIC")



# write.csv2(Modeli, "modeli_var_comp_small_brezIZ_new.csv")



