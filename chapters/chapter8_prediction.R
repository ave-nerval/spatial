# Chapter 6

############################################
# Chapter 6: Phenotype prediction
############################################

rm(list=ls())
load("data/workspaces/20220323_models_noexp_pred.RData")
load("data/workspaces/20220404_models_noexp_full.RData")




library(INLA)
library(tidyverse)
library(psych)

phenoData.small <- phenoData.small[order(phenoData.small$index), ]

# PART A: CORRELATIONS (PHENOTYPE)

# correlations
# here you need original values!!
# 2 are NAs in original data!!!!

cor(phenoData.small$last004scaled[year_pred], model.G.small.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
cor(phenoData.small$last004scaled[year_pred], model.H.small.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
cor(phenoData.small$last004scaled[year_pred], model.besag.S.UE.small.C.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
cor(phenoData.small$last004scaled[year_pred], model.besag.S.OB.small.C.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
cor(phenoData.small$last004scaled[year_pred], model.besag.HS.UE.small.C.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
cor(phenoData.small$last004scaled[year_pred], model.besag.HS.OB.small.C.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
cor(phenoData.small$last004scaled[year_pred], model.SPDE.S.small.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
cor(phenoData.small$last004scaled[year_pred], model.SPDE.HS.small.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
# cor(phenoData.small$last004scaled[year_pred], model.Exp.S.prior.small.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")
# cor(phenoData.small$last004scaled[year_pred], model.Exp.HS.prior.small.NA$summary.fitted.values["mean"][,1][year_pred], use = "complete.obs")


df_wide = data.frame("last004scaled" = phenoData.small$last004scaled[year_pred],
                     "G" = model.G.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GS_UE" = model.besag.S.UE.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GS_OB" = model.besag.S.OB.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GS_spde" = model.SPDE.S.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GH" = model.H.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GHS_UE" = model.besag.HS.UE.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GHS_OB" = model.besag.HS.OB.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GHS_spde" = model.SPDE.HS.small.NA$summary.fitted.values["mean"][,1][year_pred])



# pair plots including correlations

pairs.panels(df_wide, 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             digits = 3
)

phenoData.small <- phenoData.small[order(phenoData.small$index),]

############################################
# adding groups
############################################

tmp = table(phenoData.small$idlok, val_set) # COLSUMS 3554, 246
# tmp rownames = idlok

# OZNACI PO KATEGORIJAH 2 IN 3
# naredi korelacijo po kategorijah
# to se navezuje na lokacije!!!
sum(tmp[,1] > 0 & tmp[,2] == 0) # 1625
sum(tmp[,1] == 0 & tmp[,2] > 0) # 48 niso v trening, so v napovednem
sum(tmp[,1] > 0 & tmp[,2] > 0) # 165 so v trening in napovednem
# sum = 1838 = locations


idlok_only_train = as.numeric(names(which(tmp[,1] == 0 & tmp[,2] > 0)))
idlok_both_train_val = as.numeric(names(which(tmp[,1] > 0 & tmp[,2] > 0)))



df_wide = data.frame("last004scaled" = phenoData.small$last004scaled[year_pred],
                     "G" = model.G.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GS_UE" = model.besag.S.UE.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GS_OB" = model.besag.S.OB.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GS_spde" = model.SPDE.S.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     # "GS_exp" = model.Exp.S.prior.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GH" = model.H.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GHS_UE" = model.besag.HS.UE.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GHS_OB" = model.besag.HS.OB.small.C.NA$summary.fitted.values["mean"][,1][year_pred],
                     "GHS_spde" = model.SPDE.HS.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     # "GHS_exp" = model.Exp.HS.prior.small.NA$summary.fitted.values["mean"][,1][year_pred],
                     "idlok" = phenoData.small$idlok[year_pred],
                     "UE" = phenoData.small$UE[year_pred],
                     "idUE" = phenoData.small$idUE[year_pred],
                     "OB" = phenoData.small$OB[year_pred],
                     "idOB" = phenoData.small$idOB[year_pred])


# to so pa opazovanja!!!
df_wide$group = NA
df_wide$group[df_wide$idlok %in% idlok_only_train] = "train"
df_wide$group[df_wide$idlok %in% idlok_both_train_val] = "train_val"


df_wide_only_train = df_wide[df_wide$idlok %in% idlok_only_train,] # 59
df_wide_both_train_val = df_wide[df_wide$idlok %in% idlok_both_train_val,] # 187

drops <- c("idlok", "group", "UE")

# pearson correlation
# all data (whole df_wide) (same as above)
pairs.panels(df_wide[,!(names(df_wide) %in% drops)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             digits = 3,
             pch = 21,
             bg=c("royalblue4","lightcyan2")[factor(df_wide$group)]
)


# only train
pairs.panels(df_wide[df_wide$group == "train",][,!(names(df_wide) %in% drops)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             digits = 3,
             pch = 20
             #bg=c("royalblue4","lightcyan2")[factor(df_wide$group)]
)



pairs.panels(df_wide[df_wide$group == "train_val",][,!(names(df_wide) %in% drops)], 
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE, # show correlation ellipses
             digits = 3,
             pch = 20
             #bg=c("royalblue4","lightcyan2")[factor(df_wide$group)]
)




# PART B: diff phenotype vs spatial (GHS)


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


df_wide <- fun_merge_S_effect(df_wide, model.besag.HS.OB.small.C.NA.full, region = "OB")



df_long <- pivot_longer(df_wide, cols = c("G", "GS_UE", "GS_OB", "GS_spde", "GH",
                                          "GHS_UE", "GHS_OB", "GHS_spde"), names_to = "model")





# GS (OB) - GHS (OB)
df_wide$pheno_diff_GH_GHS <- df_wide$GH- df_wide$GHS_OB
df_wide$pheno_diff_G_GHS <- df_wide$G- df_wide$GHS_OB
df_wide$pheno_diff_GS_GHS <- df_wide$GS_spde- df_wide$GHS_OB


#dfGH = data.frame("diff" = pheno_diff_GH_GHS, "s_eff" = df_wide$effect.S.model.besag.HS.OB.small.C.NA.full)
#dfG = data.frame("diff" = pheno_diff_G_GHS, "s_eff" = df_wide$effect.S.model.besag.HS.OB.small.C.NA.full)
#dfGS = data.frame("diff" = pheno_diff_GS_GHS, "s_eff" = df_wide$effect.S.model.besag.HS.OB.small.C.NA.full)


# y_labs = c("G - GHS", "GH - GHS", "GS - GHS")
# names(y_labs) = c("BV_diff_GH_GHS_OB.small", "BV_diff_G_GHS_OB.small", "BV_diff_GS_GHS_OB.small")

pG <- dfG %>%
  # pivot_longer(cols = c("BV_diff_GH_GHS_OB.small", "BV_diff_G_GHS_OB.small", "BV_diff_GS_GHS_OB.small")) %>%
  
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin=cut_width(s_eff, width=0.5) ) %>%
  
  # plot
  ggplot( aes(x=bin, y=diff) ) +
  geom_boxplot(fill="#69b3a2") +
  theme_bw() +
  # facet_grid(name~., labeller = labeller(name = y_labs)) +
  ylab("Razlika v ocenjenih genetskih vrednostih") +
  xlab("Ocenjeni prostorski vpliv")

pG



y_labs = c("G - GHS", "GH - GHS", "GS - GHS")
names(y_labs) = c("pheno_diff_G_GHS", "pheno_diff_GH_GHS", "pheno_diff_GS_GHS")


p <- df_wide %>%
  pivot_longer(cols = c("pheno_diff_G_GHS", "pheno_diff_GH_GHS", "pheno_diff_GS_GHS")) %>%
  
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin=cut_width(effect.S.model.besag.HS.OB.small.C.NA.full, width=0.5) ) %>%
  
  # plot
  ggplot( aes(x=bin, y=value) ) +
  geom_boxplot(fill="#69b3a2") +
  theme_bw() +
  #facet_grid(name~.)+
  facet_grid(name~., labeller = labeller(name = y_labs)) +
  ylab("Razlika v ocenjenem fenotipu") +
  xlab("Ocenjeni prostorski vpliv") +
  theme(text = element_text(size = 16), axis.text = element_text(size = 8))

p



ptrain <- df_wide[df_wide$group == 'train',] %>%
  pivot_longer(cols = c("pheno_diff_G_GHS", "pheno_diff_GH_GHS", "pheno_diff_GS_GHS")) %>%
  
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin=cut_width(effect.S.model.besag.HS.OB.small.C.NA.full, width=0.5) ) %>%
  
  # plot
  ggplot( aes(x=bin, y=value) ) +
  geom_boxplot(fill="#69b3a2") +
  theme_bw() +
  #facet_grid(name~.)+
  facet_grid(name~., labeller = labeller(name = y_labs)) +
  ylab("Razlika v ocenjenem fenotipu") +
  xlab("Ocenjeni prostorski vpliv") +
  theme(text = element_text(size = 16), axis.text = element_text(size = 8)) +
  ggtitle('Lokacije ni v uèni množici')


ptrain



ptrainval <- df_wide[df_wide$group == 'train_val',] %>%
  pivot_longer(cols = c("pheno_diff_G_GHS", "pheno_diff_GH_GHS", "pheno_diff_GS_GHS")) %>%
  
  # Add a new column called 'bin': cut the initial 'carat' in bins
  mutate( bin=cut_width(effect.S.model.besag.HS.OB.small.C.NA.full, width=0.5) ) %>%
  
  # plot
  ggplot( aes(x=bin, y=value) ) +
  geom_boxplot(fill="#69b3a2") +
  theme_bw() +
  #facet_grid(name~.)+
  facet_grid(name~., labeller = labeller(name = y_labs)) +
  ylab("Razlika v ocenjenem fenotipu") +
  xlab("Ocenjeni prostorski vpliv") +
  theme(text = element_text(size = 16), axis.text = element_text(size = 8)) +
  ggtitle('Lokacija v uèni množici')

ptrainval





