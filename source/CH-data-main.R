#################
## CH-data-main.R
#################
##
## Descriptuion:  This R script generates performs data analysis under 
##                model "nModel".

# Check and install packages if needed
packages <- c("spsUtil", "Matrix")

for (pkg in packages) {
  if (!(pkg %in% rownames(installed.packages()))) {
    install.packages(pkg)
  }
}

# Load libraries
library(spsUtil)
library(Matrix)

# Load functions
source("./source/CH-fun.R")

# Receive number of model as input
nModel <- commandArgs(trailingOnly = TRUE)
nModel <- as.numeric(nModel)

#########################
## Folders for storage ##
#########################

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime,
                  "nMod", nModel, sep = "_")

## current directory
currDir <- getwd()
setwd(currDir)

## folders
resPath <- paste("./results", sep = "")
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
mBICPath <- paste0(resPath, "/mBICs/")
cAICPath <- paste0(resPath, "/cAICs/")
cBICPath <- paste0(resPath, "/cBICs/")
psi_path <- paste0(resPath, "/all_psi/")
resPath <- paste0(paste(resPath, currTime, sep = "/"), "/")
figPath <- paste0(resPath, "figures/")


ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(figPath)), dir.create(file.path(figPath)), FALSE)
ifelse(!dir.exists(file.path(mBICPath)), dir.create(file.path(mBICPath)), FALSE)
ifelse(!dir.exists(file.path(cAICPath)), dir.create(file.path(cAICPath)), FALSE)
ifelse(!dir.exists(file.path(cBICPath)), dir.create(file.path(cBICPath)), FALSE)
ifelse(!dir.exists(file.path(psi_path)), dir.create(file.path(psi_path)), FALSE)

# load variant data
variant_dataset <- read.csv(file = "./data/VariantData/dataset_cancer_cell.csv", sep = ',')


################################################
# data preprocessing
################################################
proc.data <- get.data(variant_dataset)
tps <- proc.data$tps
B0 <- proc.data$B0
R0 <- proc.data$R0
Y0 <- proc.data$Y0
Z0 <- proc.data$Z0
X0 <- proc.data$X0
ijk_sampled <- proc.data$ijk_sampled
ijktyr <- proc.data$ijktyr
IasJKt <- proc.data$IasJKt
as <- proc.data$as
################################################

## starting value for psi
psi_c <- rep(0, 14)

## candidate models
allComb <- as.matrix(expand.grid(rep(list(c(0, 1)), 4)))
colnames(allComb) <- c("0a", "0s", "ra", "rs")
allMods <- cbind("01" = 1, allComb[,c(1,2)], 
                 "r1" = 1, allComb[,c(3,4)])
allMods <- allMods[nrow(allMods):1,]

## model to fit
model <- allMods[nModel,]
covNames <- c("01", "0a", "0s", 
                  "r1", "ra", "rs")[as.logical(model)]

X0 <- X0[,which(colnames(X0) %in% covNames)]

psi_c <- head(psi_c, length(psi_c) - (6 - length(covNames)))

cat(paste0("Fitting model:\t", paste(covNames, collapse = "_"), "\n"))

print(head(X0))

EM.res <- get.EM(psi = psi_c,
                 X = X0,
                 Y = Y0,
                 Z = Z0,
                 R = R0,
                 B = B0,
                 control = list(verbose = TRUE,
                                maxIT = 10000,
                                maxTime = 23, # in hours
                                xtol = 1e-6, # 1e-8,
                                ftol = 1e-6, # 1e-8,
                                local = FALSE),
                 BFGS_ctrl = list(maxit = 10000,
                                  factr = 0, # 1e5
                                  lmm = 250, # 250
                                  trace = 1,
                                  REPORT = 1))

## store results
save.image(paste0(resPath,"/EM_res_data_nMod", 
                  sprintf("%02d", as.integer(nModel)),".RData"))

## optimal conditional moments
Ev_c <- EM.res$EVuc$Eu
Eu_c <- matrix(data = as.numeric(B0 %*% EM.res$EVuc$Eu),
               ncol = 2,
               byrow = TRUE)

## optimal beta parameters
beta_opt <- head(EM.res$x, ncol(X0))
names(beta_opt) <- colnames(X0)

## time grid for prediction
tps_pred <- seq(min(tps), max(tps), length.out = 100)

## information criteria
cAIC <- get.cAIC(psi = EM.res$x,
                 X = X0,
                 Y = Y0,
                 Z = Z0 %*% B0,
                 R = R0,
                 B = B0,
                 Euc = EM.res$EVuc$Eu)

cBIC <- get.cBIC(psi = EM.res$x,
                 X = X0,
                 Y = Y0,
                 Z = Z0 %*% B0,
                 R = R0,
                 B = B0,
                 Euc = EM.res$EVuc$Eu)

mBIC <- get.mBIC(psi = EM.res$x,
                 X = X0,
                 Y = Y0,
                 Z = Z0 %*% B0,
                 R = R0,
                 B = B0)

## store results
save.image(paste0(resPath,"/EM_res_data_wAIC_nMod", 
                  sprintf("%02d", as.integer(nModel)),".RData"))

names(cAIC) <- "cAIC"
names(cBIC) <- "cBIC"
names(mBIC) <- "mBIC"
write.table(cAIC,
            paste0(resPath, "/cAIC_nMod", sprintf("%02d", as.integer(nModel)), ".txt"), 
            sep = "\t", row.names = TRUE, col.names = FALSE)
write.table(cAIC,
            paste0(cAICPath, "/cAIC_nMod", sprintf("%02d", as.integer(nModel)), ".txt"), 
            sep = "\t", row.names = TRUE, col.names = FALSE)

write.table(cBIC,
            paste0(resPath, "/cBIC_nMod", sprintf("%02d", as.integer(nModel)), ".txt"), 
            sep = "\t", row.names = TRUE, col.names = FALSE)
write.table(cBIC,
            paste0(cBICPath, "/cBIC_nMod", sprintf("%02d", as.integer(nModel)), ".txt"), 
            sep = "\t", row.names = TRUE, col.names = FALSE)

write.table(mBIC,
            paste0(resPath, "/mBIC_nMod", sprintf("%02d", as.integer(nModel)), ".txt"),
            sep = "\t", row.names = TRUE, col.names = FALSE)
write.table(mBIC,
            paste0(mBICPath, "/mBIC_nMod", sprintf("%02d", as.integer(nModel)), ".txt"),
            sep = "\t", row.names = TRUE, col.names = FALSE)

write.table(c(cAIC, cBIC, mBIC),
            paste0(resPath, "/ICs_nMod", sprintf("%02d", as.integer(nModel)), ".txt"), 
            sep = "\t", row.names = TRUE, col.names = FALSE)


##################################################
## Generate results at clone level (i,j,k)
##################################################

## list all clones
ijk_wNames <- IasJKt[, c("i", "j", "k", 
                          "gene", "cDNA_change", "indiv_ID")]
ijk_wNames <- unique(ijk_wNames)

## sort clones by estimated growth rate
ordClones <- order(sapply(1:nrow(ijk_sampled), function(r){
  ijk_r <- ijk_sampled[r,]
  x_r <- rep(c(1, as[ijk_r$k, c("a", "s")]), times = 2)
  names(x_r) <-  c("01", "0a", "0s", 
                  "r1", "ra", "rs")
  x_r <- x_r[colnames(X0)]
  wr <- grep(pattern = "r", x = names(beta_opt), fixed = TRUE)
  as.numeric(x_r[wr] %*% beta_opt[wr] + Eu_c[r, 2])
}), decreasing = TRUE)


ordRatesClones <- sort(sapply(1:nrow(ijk_sampled), function(r){
  ijk_r <- ijk_sampled[r,]
  x_r <- rep(c(1, as[ijk_r$k, c("a", "s")]), times = 2)
  names(x_r) <-  c("01", "0a", "0s", 
                  "r1", "ra", "rs")
  x_r <- x_r[colnames(X0)]
  wr <- grep(pattern = "r", x = names(beta_opt), fixed = TRUE)
  return(as.numeric(x_r[wr] %*% beta_opt[wr] + Eu_c[r, 2]))
}), decreasing = TRUE)


## sort clones by estimated initial VAF level
x0Clones <- sapply(1:nrow(ijk_sampled), function(r){
  ijk_r <- ijk_sampled[r,]
  x_r <- rep(c(1, as[ijk_r$k, c("a", "s")]), times = 2)
  names(x_r) <-  c("01", "0a", "0s", 
                  "r1", "ra", "rs")
  x_r <- x_r[colnames(X0)]
  w0 <- grep(pattern = "0", x = names(beta_opt), fixed = TRUE)
  return(as.numeric(x_r[w0] %*% beta_opt[w0] + Eu_c[r, 1]))
})

## list all clones, sorted by estimated growth rates, including rates
ijk_wNames_wEst <- cbind(ijk_wNames[ordClones,],
                         x0 = get.logistic(x0Clones)[ordClones],
                         r = ordRatesClones,
                         u0 = Eu_c[ordClones,1],
                         u1 = Eu_c[ordClones,2])
rownames(ijk_wNames_wEst) <- NULL

## export results
write.csv(ijk_wNames_wEst, paste0(resPath, "/ijk_estimates.csv"))

pdf(file = paste0(figPath, "/ijk_plots.pdf"), width = 10, height = 7)
par(mar = c(2,2,1,.5), mfrow = c(5,5))
for (l in 1:length(ordClones)) {
  
  r <- ordClones[l]
  ijk_r <- ijk_sampled[r,]
  
  X0_r <- cbind(1, t(replicate(length(tps_pred), as[ijk_r$k, c("a", "s")])),
                tps_pred, t(replicate(length(tps_pred), as[ijk_r$k, c("a", "s")]))*tps_pred)
  colnames(X0_r) <- c("01", "0a", "0s", 
                  "r1", "ra", "rs")
  X0_r <- X0_r[, names(beta_opt)]

  ## random-effects prediction matrix
  Z0_r <- cbind(1, tps_pred)
  
  plot(x = tps[ijktyr[which(ijktyr$i == ijk_r$i &
                                 ijktyr$j == ijk_r$j &
                                 ijktyr$k == ijk_r$k), c("t")]],
       y = ijktyr[which(ijktyr$i == ijk_r$i &
                          ijktyr$j == ijk_r$j &
                          ijktyr$k == ijk_r$k), c("VAF_01")],
       xlim = range(tps),
       ylim = c(0,1),
       pch = 20,
       cex = 2,
       cex.axis = 1,
       cex.lab = 1,
       main = paste0(ijk_wNames_wEst[l, 4:6], collapse = "_"),
       xlab = "time",
       ylab = "VAF")
  
  lines(x = tps_pred,
        y = get.logistic(X0_r %*% beta_opt + Z0_r %*% Eu_c[r,]),
        lwd = 3, col = "blue")
}
dev.off()


####################################
## Generate results at level (i,j) 
####################################

L2 <- build_L2_sparse(ijk_sampled = ijk_sampled, 
                      p = 2)

ij_sampled <- unique(ijk_sampled[, c("i", "j")])

Euij_c <- matrix(data = as.numeric(L2 %*% Ev_c),
                 ncol = 2,
                 byrow = TRUE)

## estimated parameters
parms_js <- t(sapply(1:nrow(ij_sampled), function(r){
  ij_r <- ij_sampled[r,]
  ijks_r <- ijk_sampled[which(ijk_sampled$i == ij_r$i &
                                ijk_sampled$j == ij_r$j),]
  
  theta_ijs <- sapply(1:nrow(ijks_r), function(s){
    ijk_r <- ijks_r[s,]
    
    X_ijk_r <- bdiag(lapply(1:2, function(l){t(c(1, as[ijk_r$k, c("a", "s")]))}))
    colnames(X_ijk_r) <- c("01", "0a", "0s", 
                           "r1", "ra", "rs")
    
    X_ijk_r <- X_ijk_r[,names(beta_opt)]
    
    as.numeric(X_ijk_r %*% beta_opt)
    
  })
  
  mean_theta_ijs <- apply(theta_ijs, 1, mean)
  return(mean_theta_ijs + Euij_c[r,])
}))


ordcDNAchs <- order(parms_js[,2], decreasing = TRUE)

ij_wNames <- IasJKt[, c("i", "j", 
                          "gene", "cDNA_change")]
ij_wNames <- unique(ij_wNames)

ij_wNames_wEst <- cbind(ij_wNames[ordcDNAchs,],
                         x0 = get.logistic(parms_js[ordcDNAchs,1]),
                         r = parms_js[ordcDNAchs,2],
                         u0 = Euij_c[ordcDNAchs,1],
                         u1 = Euij_c[ordcDNAchs,2])
rownames(ij_wNames_wEst) <- NULL

## export results
write.csv(ij_wNames_wEst, paste0(resPath, "/ij_estimates.csv"))

pdf(file = paste0(figPath, "/ij_plots.pdf"), width = 10, height = 7)
par(mar = c(2,2,1,.5), mfrow = c(5,5))
for (l in 1:length(ordcDNAchs)){
  r <- ordcDNAchs[l]
  ij_r <- ij_sampled[r,]
  ijks_r <- ijk_sampled[which(ijk_sampled$i == ij_r$i &
                                ijk_sampled$j == ij_r$j),]
  
  theta_ijs <- sapply(1:nrow(ijks_r), function(s){
    ijk_r <- ijks_r[s,]
    
    X_ijk_r <- bdiag(lapply(1:2, function(l){t(c(1, as[ijk_r$k, c("a", "s")]))}))
    colnames(X_ijk_r) <- c("01", "0a", "0s", 
                           "r1", "ra", "rs")
    
    X_ijk_r <- X_ijk_r[,names(beta_opt)]
    
    as.numeric(X_ijk_r %*% beta_opt)
    
  })
  
  mean_theta_ijs <- apply(theta_ijs, 1, mean)
  
  plot(x = tps[ijktyr[which(ijktyr$i == ij_r$i &
                                 ijktyr$j == ij_r$j), c("t")]],
       y = ijktyr[which(ijktyr$i == ij_r$i &
                          ijktyr$j == ij_r$j), c("VAF_01")],
       xlim = range(tps),
       ylim = c(0,1),
       pch = 20,
       cex = 2,
       cex.axis = 1,
       cex.lab = 1,
       main = paste0(ij_wNames_wEst[l, 3:4], collapse = "_"),
       xlab = "time",
       ylab = "VAF")
  
  lines(x = tps_pred,
        y = get.logistic(cbind(1, tps_pred) %*% (mean_theta_ijs + Euij_c[r,])),
        lwd = 2,
        col = "blue")
  
}
dev.off()

#################################
## Generate results at level (i) 
#################################

L1 <- build_L1_sparse(ijk_sampled = ijk_sampled, 
                      p = 2)

i_sampled <- unique(ijk_sampled[, c("i")])

Eui_c <- matrix(data = as.numeric(L1 %*% Ev_c),
                ncol = 2,
                byrow = TRUE)

parms_is <- t(sapply(1:length(i_sampled), function(r){
  i_r <- i_sampled[r]
  ijks_r <- ijk_sampled[which(ijk_sampled$i == i_r),]
  
  theta_is <- sapply(1:nrow(ijks_r), function(s){
    ijk_r <- ijks_r[s,]
    
    X_ijk_r <- bdiag(lapply(1:2, function(l){t(c(1, as[ijk_r$k, c("a", "s")]))}))
    colnames(X_ijk_r) <- c("01", "0a", "0s", 
                           "r1", "ra", "rs")
    
    X_ijk_r <- X_ijk_r[,names(beta_opt)]
    
    as.numeric(X_ijk_r %*% beta_opt)
    
  })
  
  mean_theta_is <- apply(theta_is, 1, mean)
  
  return(mean_theta_is + Eui_c[r,])
}))


ordGenes <- order(parms_is[,2], decreasing = TRUE)
i_wNames <- IasJKt[, c("i", "gene")]
i_wNames <- unique(i_wNames)
i_wNames_wEst <- cbind(i_wNames[ordGenes,],
                         x0 = get.logistic(parms_is[ordGenes,1]),
                         r = parms_is[ordGenes,2],
                         u0 = Eui_c[ordGenes,1],
                         u1 = Eui_c[ordGenes,2])
rownames(i_wNames_wEst) <- NULL

## export results
write.csv(i_wNames_wEst, paste0(resPath, "/i_estimates.csv"))

pdf(file = paste0(figPath, "/i_plots.pdf"), width = 10, height = 7)
par(mar = c(2,2,1,.5), mfrow = c(5,5))
for (l in 1:length(ordGenes)) {
  r <- ordGenes[l]
  i_r <- i_sampled[r]
  ijks_r <- ijk_sampled[which(ijk_sampled$i == i_r),]
  
  theta_is <- sapply(1:nrow(ijks_r), function(s){
    ijk_r <- ijks_r[s,]
    
    X_ijk_r <- bdiag(lapply(1:2, function(l){t(c(1, as[ijk_r$k, c("a", "s")]))}))
    colnames(X_ijk_r) <- c("01", "0a", "0s", 
                           "r1", "ra", "rs")
    
    X_ijk_r <- X_ijk_r[,names(beta_opt)]
    
    as.numeric(X_ijk_r %*% beta_opt)
    
  })
  
  mean_theta_is <- apply(theta_is, 1, mean)

  plot(x = tps[ijktyr[which(ijktyr$i == i_r), c("t")]],
       y = ijktyr[which(ijktyr$i == i_r), c("VAF_01")],
       xlim = range(tps),
       ylim = c(0,1),
       pch = 20,
       cex = 2,
       cex.axis = 1,
       cex.lab = 1,
       main = i_wNames_wEst[l,2],
       xlab = "time",
       ylab = "VAF")
  
  lines(x = tps_pred,
        y = get.logistic(cbind(1, tps_pred) %*% (mean_theta_is + Eui_c[r,])),
        lwd = 3, col = "blue")
}
dev.off()

## export results
write.csv(EM.res$x, paste0(resPath,"/psi_estimates_nmod", 
                           sprintf("%02d", as.integer(nModel)),".csv"))

write.csv(EM.res$x, paste0(psi_path,"/psi_estimates_nmod", 
                           sprintf("%02d", as.integer(nModel)),".csv"))

