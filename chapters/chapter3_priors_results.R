# Chapter 3

# rm(list=ls())
# load("data/workspaces/20220323_chapter3_models.RData")

############################################
# Chapter 3: Vpliv apriorne porazdelitve na aposteriorno oceno variance
############################################

var_comp_model <- function(chosen_model, var_comp){
  
  if(var_comp=="gen"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for rowNumberAinv`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp[c("quant0.025", "quant0.25","mean", "quant0.5", "quant0.75", "quant0.975")]), 2)
    var_comp <- c(var_comp, "gen")
  }
  
  else if(var_comp=="res"){
    var_comp <- inla.tmarginal(function(x) 1/x,chosen_model$marginals.hyperpar$`Precision for the Gaussian observations`)
    var_comp <- inla.zmarginal(var_comp)
    var_comp <- round(unlist(var_comp[c("quant0.025", "quant0.25","mean", "quant0.5", "quant0.75", "quant0.975")]), 2)
    var_comp <- c(var_comp, "res")
    
  }
  
  c(deparse(substitute(chosen_model)), var_comp)
  
}


df_gen <- rbind(#var_comp_model(model.G.small, var_comp = "gen"),
  #var_comp_model(model.G.small, var_comp = "res"),
  var_comp_model(model.G.small, var_comp = "gen"),
  var_comp_model(model.G.small, var_comp = "res"),
  var_comp_model(model.G.small.01, var_comp = "gen"),
  var_comp_model(model.G.small.01, var_comp = "res"),
  var_comp_model(model.G.small.03, var_comp = "gen"),
  var_comp_model(model.G.small.03, var_comp = "res"),
  var_comp_model(model.G.small.05, var_comp = "gen"),
  var_comp_model(model.G.small.05, var_comp = "res"),
  var_comp_model(model.G.small.07, var_comp = "gen"),
  var_comp_model(model.G.small.07, var_comp = "res"),
  var_comp_model(model.G.small.09, var_comp = "gen"),
  var_comp_model(model.G.small.09, var_comp = "res"))

df_gen <- as.data.frame(df_gen)
df_gen$prior <- c(rep(c(0, 0.01, 0.03, 0.05, 0.07, 0.09), each = 2))
colnames(df_gen) <- c("model", "lower", "Q1", "mean", "median", "Q3", "upper", "variance", "prior")
df_gen$lower <- as.numeric(df_gen$lower)
df_gen$mean <- as.numeric(df_gen$mean)
df_gen$median <- as.numeric(df_gen$median)
df_gen$upper <- as.numeric(df_gen$upper)
df_gen$Q1 <- as.numeric(df_gen$Q1)
df_gen$Q3 <- as.numeric(df_gen$Q3)

ggplot(df_gen, aes(x = prior, y = mean, color = variance)) +
  geom_point(size = 3, shape = 1) +
  geom_pointrange(aes(ymax = upper, ymin = lower, alpha = 0.1), linetype=2, show_guide=FALSE) +
  geom_pointrange(aes(ymax = Q3, ymin = Q1)) +
  theme_bw() +
  ylab("Aposteriorna ocena variance") +
  xlab("u") +
  scale_color_discrete(name = "Komponenta variance", labels = c("Genetska", "Ostanek")) +
  #scale_x_continuous(breaks = round(seq(0, 1, by = 0.1),1)) 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10))


