# Chapter 1

# rm(list=ls())
# load("data/workspaces/20220323_models_noexp_pred.RData")

############################################
# Chapter 1: posterior precision/variance
############################################

# # Plot posterior genomic precision 
par(mfrow = c(2,2))
#genetic variance
plot(inla.tmarginal(function(x) 1/x,model.G.small.NA$marginals.hyperpar$`Precision for rowNumberAinv`), type = "l",xlim = c(0,0.7), main = "Genetska varianca", ylim = c(0,20), xlab = "Varianca", ylab = "Gostota porazdelitve")
lines(inla.tmarginal(function(x) 1/x, model.H.small.NA$marginals.hyperpar$`Precision for rowNumberAinv` ), col = 2)
lines(inla.tmarginal(function(x) 1/x, model.besag.S.OB.small.C.NA$marginals.hyperpar$`Precision for rowNumberAinv` ), col = 3)
lines(inla.tmarginal(function(x) 1/x, model.besag.S.UE.small.C.NA$marginals.hyperpar$`Precision for rowNumberAinv` ), col = 3, lty = 2)
lines(inla.tmarginal(function(x) 1/x, model.SPDE.S.small.NA$marginals.hyperpar$`Precision for rowNumberAinv` ), col = 3, lty = 3)

lines(inla.tmarginal(function(x) 1/x, model.besag.HS.OB.small.C.NA$marginals.hyperpar$`Precision for rowNumberAinv` ), col = 4)
lines(inla.tmarginal(function(x) 1/x, model.besag.HS.UE.small.C.NA$marginals.hyperpar$`Precision for rowNumberAinv` ), col = 4, lty = 2)
lines(inla.tmarginal(function(x) 1/x, model.SPDE.HS.small.NA$marginals.hyperpar$`Precision for rowNumberAinv` ), col = 4, lty = 3)

# legend("right", col = 1:4, legend = c("G", "GH", "GS", "GHS"), lty  =1)

#residual variance
plot(inla.tmarginal(function(x) 1/x,model.G.small.NA$marginals.hyperpar$`Precision for the Gaussian observations`), type = "l",xlim = c(0.4,0.95), main = "Ostanek",ylim = c(0,20), xlab = "Varianca", ylab = "Gostota porazdelitve")
lines(inla.tmarginal(function(x) 1/x, model.H.small.NA$marginals.hyperpar$`Precision for the Gaussian observations` ), col = 2)
lines(inla.tmarginal(function(x) 1/x, model.besag.S.OB.small.C.NA$marginals.hyperpar$`Precision for the Gaussian observations` ), col = 3)
lines(inla.tmarginal(function(x) 1/x, model.besag.S.UE.small.C.NA$marginals.hyperpar$`Precision for the Gaussian observations` ), col = 3, lty = 2)
lines(inla.tmarginal(function(x) 1/x, model.SPDE.S.small.NA$marginals.hyperpar$`Precision for the Gaussian observations` ), col = 3, lty = 3)

lines(inla.tmarginal(function(x) 1/x, model.besag.HS.OB.small.C.NA$marginals.hyperpar$`Precision for the Gaussian observations` ), col = 4)
lines(inla.tmarginal(function(x) 1/x, model.besag.HS.UE.small.C.NA$marginals.hyperpar$`Precision for the Gaussian observations` ), col = 4, lty = 2)
lines(inla.tmarginal(function(x) 1/x, model.SPDE.HS.small.NA$marginals.hyperpar$`Precision for the Gaussian observations` ), col = 4, lty = 3)

# legend("right", col = 1:3, legend = c("G", "GH", "GS", "GHS"), lty  =1)

#herd effect variance
plot(inla.tmarginal(function(x) 1/x,model.H.small.NA$marginals.hyperpar$`Precision for idlok`), type = "l",xlim = c(0.1,0.5), main = "Varianca vpliva èrede",ylim = c(0,20), xlab = "Varianca", ylab = "Gostota porazdelitve", col = 2)
lines(inla.tmarginal(function(x) 1/x, model.besag.HS.OB.small.C.NA$marginals.hyperpar$`Precision for idlok` ), col = 4)
lines(inla.tmarginal(function(x) 1/x, model.besag.HS.UE.small.C.NA$marginals.hyperpar$`Precision for idlok` ), col = 4, lty = 2)
lines(inla.tmarginal(function(x) 1/x, model.SPDE.HS.small.NA$marginals.hyperpar$`Precision for idlok` ), col = 4, lty = 3)
# legend("right", col = 3, legend = c("BESAG (UE)", "BESAG (OB)", "SPDE"), lty  =1:3)

# legend("left", col = 1:2, legend = c("GH", "GHS"), lty  =1)

#spatial variance
spdeHS =inla.spde2.result(inla =model.SPDE.HS.small.NA, name = "fieldID",spde = spdeStatHS)
spdeS =inla.spde2.result(inla =model.SPDE.S.small.NA, name = "fieldID",spde = spdeStatS)


plot(inla.tmarginal(function(x) 1/x,model.besag.S.OB.small.C.NA$marginals.hyperpar$`Precision for idOB`), col = 3, type = "l",xlim = c(0,0.35), main = "Prostorska varianca",ylim = c(0,20), xlab = "Varianca", ylab = "Gostota porazdelitve")
lines(inla.tmarginal(function(x) 1/x, model.besag.S.UE.small.C.NA$marginals.hyperpar$`Precision for idUE` ), col = 3, lty = 2)
# lines(inla.tmarginal(function(x) 1/x, model.SPDE.S.small.NA$marginals.hyperpar$`Precision for idUE` ), col = 1)

lines(inla.tmarginal(function(x) 1/x, model.besag.HS.OB.small.C.NA$marginals.hyperpar$`Precision for idOB` ), col = 4)
lines(inla.tmarginal(function(x) 1/x, model.besag.HS.UE.small.C.NA$marginals.hyperpar$`Precision for idUE` ), col = 4, lty = 2)
# lines(inla.tmarginal(function(x) 1/x, model.SPDE.HS.small.NA$marginals.hyperpar$`Precision for idUE` ), col = 2)
lines(spdeS$marginals.variance.nominal$variance.nominal.1, type = "l",xlim = c(0,0.3), main = "Variance for field", col = 3, lty = 3)
lines(spdeHS$marginals.variance.nominal$variance.nominal.1, type = "l",xlim = c(0,0.3), main = "Variance for field", col = 4, lty = 3)





# par(mfrow = c(1,1))
# plot(inla.tmarginal(function(x) 1/x,model.besag.S.OB.small.C.NA$marginals.hyperpar$`Precision for idOB`), col = 3, type = "l",xlim = c(0,0.4), main = "Prostorska varianca",ylim = c(0,20), xlab = "Varianca", ylab = "Gostota porazdelitve")

# legend("right", col = 1:2, legend = c("GS", "GHS"), lty  =1)
# legend("right", col = 1:4, legend = c("G", "GH", "GS", "GHS"), lty  =1, title = "Model")
# legend("right", col = 4, title = "Prostorski model", legend = c("BESAG (UE)", "BESAG (OB)", "SPDE"), lty  =1:3)




