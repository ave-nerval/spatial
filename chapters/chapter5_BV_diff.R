# Chapter 5

############################################
# Chapter 5: Difference in estimated genetic effects
############################################

rm(list=ls())
load("data/workspaces/20220323_models_noexp_pred.RData")




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



# merging effects


fun_merge_S_effect <- function(data, model, region = "UE"){
  # SPATIAL EFFECTS (GS + GHS model)
  
  if (region == 'UE'){
    spatial_effect <- cbind(model$summary.random$idUE[,c("ID","mean")])
    spatial_effect$mean = spatial_effect$mean/sqrt(as.numeric(var_comp_model_brezIZ(chosen_model = model, var_comp = "spatial", region = "UE")))
    
    names(spatial_effect) <- c("UE_ID", paste0('effect.S.', deparse(substitute(model))))
    
    ID_UE <- data.frame("UE" = boundaryUE$UE_MID, "UE_ID" = 1:nrow(boundaryUE))
    spatial_effect <- merge(ID_UE, spatial_effect, all.x=T)
  }
  #merge with phenodata, merge by UE
  #now we added an estimated spatial effect to each cow
  
  
  if (region == 'OB'){
    spatial_effect <- cbind(model$summary.random$idOB[,c("ID","mean")])
    spatial_effect$mean = spatial_effect$mean/sqrt(as.numeric(var_comp_model_brezIZ(chosen_model = model, var_comp = "spatial", region = "OB")))
    
    names(spatial_effect) <- c("OB_ID", paste0('effect.S.', deparse(substitute(model))))
    
    ID_OB <- data.frame("OB" = boundaryOB$OB_MID, "OB_ID" = 1:nrow(boundaryOB))
    spatial_effect <- merge(ID_OB, spatial_effect, all.x=T)
  }
  
  data <- as.data.frame(merge(data, spatial_effect, all.x=T))
  
  data
  
}


fun_merge_G_effect <- function(data, model){
  name_effect = paste0('effect.G.', deparse(substitute(model)))
  randomeff.Gen <- data.frame("rowNumberAinv" = model$summary.random$rowNumberAinv[,"ID"],
                              "Geff" = model$summary.random$rowNumberAinv[,"mean"])
  randomeff.Gen$Geff = randomeff.Gen$Geff/sqrt(as.numeric(var_comp_model_brezIZ(chosen_model = model, var_comp = "gen")))
  colnames(randomeff.Gen) = c("rowNumberAinv", name_effect)
  #merge together; genetic effect
  #merge by rowNumberAinv
  data <- as.data.frame(merge(data, randomeff.Gen, all.x=T))
  
  data
}


phenoData.small = fun_merge_S_effect(phenoData.small, model.besag.HS.OB.small.C.NA, region = "OB")
phenoData.small = fun_merge_G_effect(phenoData.small, model.G.small.NA)
phenoData.small = fun_merge_G_effect(phenoData.small, model.H.small.NA)
phenoData.small = fun_merge_G_effect(phenoData.small, model.besag.S.OB.small.C.NA)
phenoData.small = fun_merge_G_effect(phenoData.small, model.besag.HS.OB.small.C.NA)



phenoData.small = phenoData.small[order(phenoData.small$index),]

# differences in estimated BV
# small models
# G - GHS (OB)
BV_diff_G_GHS_OB.small  <- phenoData.small$effect.G.model.G.small[-year_pred] - phenoData.small$effect.G.model.besag.HS.OB.small.C[-year_pred]

# GH - GHS (OB)
BV_diff_GH_GHS_OB.small <- phenoData.small$effect.G.model.H.small[-year_pred] - phenoData.small$effect.G.model.besag.HS.OB.small.C[-year_pred]

# GS (OB) - GHS (OB)
BV_diff_GS_GHS_OB.small <- phenoData.small$effect.G.model.besag.S.OB.small.C[-year_pred]- phenoData.small$effect.G.model.besag.HS.OB.small.C[-year_pred]




df = data.frame(BV_diff_GH_GHS_OB.small,BV_diff_G_GHS_OB.small,BV_diff_GS_GHS_OB.small,"effect.S.model.besag.HS.OB.small.C" = phenoData.small$effect.S.model.besag.HS.OB.small.C[-year_pred])


y_labs = c("G - GHS", "GH - GHS", "GS - GHS")
names(y_labs) = c("BV_diff_G_GHS_OB.small", "BV_diff_GH_GHS_OB.small", "BV_diff_GS_GHS_OB.small")

p <- df %>%
  pivot_longer(cols = c("BV_diff_GH_GHS_OB.small", "BV_diff_G_GHS_OB.small", "BV_diff_GS_GHS_OB.small")) %>%
  
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin=cut_width(effect.S.model.besag.HS.OB.small.C, width=0.5) ) %>%
  
  # plot
  ggplot( aes(x=bin, y=value) ) +
  geom_boxplot(fill="#69b3a2") +
  theme_bw() +
  facet_grid(name~., labeller = labeller(name = y_labs)) +
  ylab("Razlika v ocenjenih genetskih vrednostih") +
  xlab("Ocenjeni prostorski vpliv") +
  theme(text = element_text(size = 16), axis.text = element_text(size = 8))

p



phenoData.small = fun_merge_S_effect(phenoData.small, model.besag.HS.UE.small.C.NA, region = "UE")
# phenoData.small = fun_merge_G_effect(phenoData.small, model.G.small.NA)
# phenoData.small = fun_merge_G_effect(phenoData.small, model.H.small.NA)
phenoData.small = fun_merge_G_effect(phenoData.small, model.besag.S.UE.small.C.NA)
phenoData.small = fun_merge_G_effect(phenoData.small, model.besag.HS.UE.small.C.NA)



# differences in estimated BV
# small models
# G - GHS (UE)
BV_diff_G_GHS_UE.small  <- phenoData.small$effect.G.model.G.small[-year_pred] - phenoData.small$effect.G.model.besag.HS.UE.small.C[-year_pred]

# GH - GHS (UE)
BV_diff_GH_GHS_UE.small <- phenoData.small$effect.G.model.H.small[-year_pred] - phenoData.small$effect.G.model.besag.HS.UE.small.C[-year_pred]

# GS (UE) - GHS (UE)
BV_diff_GS_GHS_UE.small <- phenoData.small$effect.G.model.besag.S.UE.small.C[-year_pred]- phenoData.small$effect.G.model.besag.HS.UE.small.C[-year_pred]




df = data.frame(BV_diff_GH_GHS_UE.small,BV_diff_G_GHS_UE.small,BV_diff_GS_GHS_UE.small,"effect.S.model.besag.HS.UE.small.C" = phenoData.small$effect.S.model.besag.HS.UE.small.C[-year_pred])


y_labs = c("G - GHS", "GH - GHS", "GS - GHS")
names(y_labs) = c("BV_diff_GH_GHS_UE.small", "BV_diff_G_GHS_UE.small", "BV_diff_GS_GHS_UE.small")

p <- df %>%
  pivot_longer(cols = c("BV_diff_GH_GHS_UE.small", "BV_diff_G_GHS_UE.small", "BV_diff_GS_GHS_UE.small")) %>%
  
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin=cut_width(effect.S.model.besag.HS.UE.small.C, width=0.5) ) %>%
  
  # plot
  ggplot( aes(x=bin, y=value) ) +
  geom_boxplot(fill="#69b3a2") +
  theme_bw() +
  facet_grid(name~., labeller = labeller(name = y_labs)) +
  ylab("Razlika v ocenjenih genetskih vrednostih") +
  xlab("Ocenjeni prostorski vpliv") +
  theme(text = element_text(size = 16), axis.text = element_text(size = 10))   

p

