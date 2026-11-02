#################
## CH-sim-main.R
#################
##
## Descriptuion:  Simulation studies for validation.
##                This script performs a single simulation.

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

nSim <- commandArgs(trailingOnly = TRUE)
nSim <- as.numeric(nSim) # unlist(strsplit(nSim, split = "-"))[1]

#########################
## Folders for storage ##
#########################

currTime <- Sys.time()
currTime <- gsub(" ", "_", currTime)
currTime <- gsub("-", "", currTime)
currTime <- gsub(":", "", currTime)
currTime <- paste(currTime,
                  "nSim", nSim, sep = "_")

## current directory
currDir <- getwd()
setwd(currDir)

## folders
resPath <- paste("./results", sep = "")
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
resPath <- paste0(paste(resPath, currTime, sep = "/"), "/")
figPath <- paste0(resPath, "figures/")

ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
ifelse(!dir.exists(file.path(figPath)), dir.create(file.path(figPath)), FALSE)

simResPath <- paste("./results/sim/", sep = "")
simFigPath <- paste0(simResPath, "/figures/")
psi_path <- paste0(simResPath, "/all_psi/")
Ev_path <- paste0(simResPath, "/all_Ev/")

Eu_ijk_path <- paste0(simResPath, "/all_Euijk/")
Eu_ij_path <- paste0(simResPath, "/all_Euij/")
Eu_i_path <- paste0(simResPath, "/all_Eui/")

ifelse(!dir.exists(file.path(simResPath)), dir.create(file.path(simResPath)), FALSE)
ifelse(!dir.exists(file.path(simFigPath)), dir.create(file.path(simFigPath)), FALSE)
ifelse(!dir.exists(file.path(psi_path)), dir.create(file.path(psi_path)), FALSE)
ifelse(!dir.exists(file.path(Ev_path)), dir.create(file.path(Ev_path)), FALSE)
ifelse(!dir.exists(file.path(Eu_ijk_path)), dir.create(file.path(Eu_ijk_path)), FALSE)

ifelse(!dir.exists(file.path(Eu_ij_path)), dir.create(file.path(Eu_ij_path)), FALSE)
ifelse(!dir.exists(file.path(Eu_i_path)), dir.create(file.path(Eu_i_path)), FALSE)

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
ijkt <- proc.data$ijkt
IasJKt <- proc.data$IasJKt
as <- proc.data$as
################################################

## true vector parameter
psi.true <- c(-0.7669597469819380286537580104777589440345764160156250000000000000000000000000000000000000000000000000,
              0.0077029788264883101642221951976807758910581469535827636718750000000000000000000000000000000000000000,
              0.0273227815909087012880007705462048761546611785888671875000000000000000000000000000000000000000000000,
              0.0108505257266413007027905379686671949457377195358276367187500000000000000000000000000000000000000000,
              -0.0000647175321981660947850273357673245300247799605131149291992187500000000000000000000000000000000000,
              -0.0001396207920001559909650018509452706894080620259046554565429687500000000000000000000000000000000000,
              0.1189805137526380035550843672353948932141065597534179687500000000000000000000000000000000000000000000,
              -0.3761189859179580041725898809090722352266311645507812500000000000000000000000000000000000000000000000,
              2.0060115003182699133787991740973666310310363769531250000000000000000000000000000000000000000000000000,
              -10.9841614428155995852876003482379019260406494140625000000000000000000000000000000000000000000000000000,
              -2.6054123218335898393149818730307742953300476074218750000000000000000000000000000000000000000000000000,
              -12.8128441110841997385705326450988650321960449218750000000000000000000000000000000000000000000000000000,
              0.0323251819583681987424483850190881639719009399414062500000000000000000000000000000000000000000000000,
              -10.3351694224250003628640115493908524513244628906250000000000000000000000000000000000000000000000000000)
b.true <- psi.true[1:6]
s.true <- psi.true[7:8]
tau.true <- tail(psi.true, 6)

Dg <- diag(exp(tau.true[1:2]))
Dc <- diag(exp(tau.true[3:4]))
Dp <- diag(exp(tau.true[5:6]))

## simulate true random-effects
set.seed(123)
U <- NULL
U_ij <- NULL

u_i <- matrix(data = rnorm(n = 2*length(unique(ijk_sampled$i)),
                           mean = 0,
                           sd = sqrt(diag(Dg))),
              ncol = 2,
              byrow = TRUE)

for (i in unique(ijk_sampled$i)) {
  ijk_sampled_i <- ijk_sampled[which(ijk_sampled$i == i),]
  
  u_ij <- t(matrix(data = rnorm(n = 2*length(unique(ijk_sampled_i$j)),
                                mean = u_i[i,],
                                sd = sqrt(diag(Dc))),
                   nrow = 2,
                   ncol = length(unique(ijk_sampled_i$j))))
  U_ij <- rbind(U_ij, cbind(i = i, j = unique(ijk_sampled_i$j), u_ij))
  
  for (j in unique(ijk_sampled_i$j)) {
    ijk_sampled_ij <- ijk_sampled_i[which(ijk_sampled_i$j == j),]
    
    u_ijk <- t(matrix(data = rnorm(n = 2*length(unique(ijk_sampled_ij$k)),
                                   mean = u_ij[which(unique(ijk_sampled_i$j) == j),],
                                   sd = sqrt(diag(Dp))),
                      nrow = 2,
                      ncol = length(unique(ijk_sampled_ij$k))))
    
    U <- rbind(U, cbind(ijk_sampled_ij, u_ijk))
    
  }
  
}

colnames(U) <- c("i", "j", "k", "u1", "u2")
colnames(U_ij) <- c("i", "j", "u1", "u2")

U_i <- cbind(i = I, u_i)
colnames(U_i) <- c("i", "u1", "u2")

u.true <- c(t(U[, c("u1", "u2")]))


## generate synthetic data
Y <- as.numeric(X0 %*% b.true + Z0 %*% u.true)
set.seed(nSim)
Sigma.true <- exp(s.true[1] + s.true[2]*R0)
noise <- rnorm(n = length(Y),
               mean = rep(0, length(Y)),
               sd = sqrt(Sigma.true))
Y0 <- Y + noise
rm(Y)
ijkty <- cbind(ijkt, Y0)

## start value for vector parameter
psi_c <- rep(0, 14)

## candidate models
allComb <- as.matrix(expand.grid(rep(list(c(0, 1)), 4)))
colnames(allComb) <- c("0a", "0s", "ra", "rs")
allMods <- cbind("01" = 1, allComb[,c(1,2)], 
                 "r1" = 1, allComb[,c(3,4)])
allMods <- allMods[nrow(allMods):1,]

## true model
nModel <- 1
model <- allMods[nModel,]
covNames <- c("01", "0a", "0s", 
                  "r1", "ra", "rs")[as.logical(model)]
X0 <- X0[,which(colnames(X0) %in% covNames)]

## parameters of true model
psi_c <- head(psi_c, length(psi_c) - (6 - length(covNames)))

cat(paste0("Fitting model:\t", paste(covNames, collapse = "_"), "\n"))

print(head(X0))

## Run E-M inference algorithm
EM.res <- get.EM(psi = psi_c,
                 X = X0,
                 Y = Y0,
                 Z = Z0,
                 R = R0,
                 B = B0,
                 control = list(verbose = TRUE,
                                maxIT = 10000,
                                maxTime = 23, # in hours
                                xtol = 1e-5, # 1e-8,
                                ftol = 1e-6, # 1e-8,
                                local = FALSE),
                 BFGS_ctrl = list(maxit = 10000,
                                  factr = 0, # 1e5
                                  lmm = 250, # 250
                                  trace = 1,
                                  REPORT = 1))

## export results
save.image(paste0(resPath,"/EM_res_data_nSim", sprintf("%03d", as.integer(nSim)),".RData"))

write.csv(EM.res$x, paste0(psi_path, "/psi_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))
write.csv(EM.res$x, paste0(resPath, "/psi_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))

pdf(file = paste0(figPath, "/psi_est.pdf"), width = 5, height = 4)
par(mar = c(2,2,1,.5), mfrow = c(1,1))
plot(psi.true, EM.res$x)
lines(range(c(psi.true, EM.res$x)),
      range(c(psi.true, EM.res$x)),
      col = "red")
dev.off()

Evc <- EM.res$EVuc$Eu


write.csv(Evc, paste0(Ev_path, "/Ev_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))
write.csv(Evc, paste0(resPath, "/Ev_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))

write.csv(as.numeric(B0 %*% Evc), paste0(Eu_ijk_path, "/Eu_ijk_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))
write.csv(as.numeric(B0 %*% Evc), paste0(resPath, "/Eu_ijk_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))

pdf(file = paste0(figPath, "/u_ijk.pdf"), width = 5, height = 4)
par(mar = c(2,2,1,.5), mfrow = c(1,1))
plot(x = c(t(U[,c("u1", "u2")])), 
     y = as.numeric(B0 %*% Evc),
     xlim = range(c(t(U[,c("u1", "u2")])), 
                  as.numeric(B0 %*% Evc)),
     ylim = range(c(t(U[,c("u1", "u2")])), 
                  as.numeric(B0 %*% Evc)))
lines(range(c(t(U[,c("u1", "u2")])), 
            as.numeric(B0 %*% Evc)),
      range(c(t(U[,c("u1", "u2")])), 
            as.numeric(B0 %*% Evc)),
      col = "red")
dev.off()

L2 <- build_L2_sparse(ijk_sampled = ijk_sampled, p = 2)
L1 <- build_L1_sparse(ijk_sampled = ijk_sampled, p = 2)

pdf(file = paste0(figPath, "/u_ij.pdf"), width = 5, height = 4)
par(mar = c(2,2,1,.5), mfrow = c(1,1))
plot(c(t(U_ij[,c("u1", "u2")])), 
     L2 %*% Evc,
     xlim = range(c(t(U_ij[,c("u1", "u2")])), as.numeric(L2 %*% Evc)),
     ylim = range(c(t(U_ij[,c("u1", "u2")])), as.numeric(L2 %*% Evc)))
lines(range(c(t(U_ij[,c("u1", "u2")])), as.numeric(L2 %*% Evc)),
      range(c(t(U_ij[,c("u1", "u2")])), as.numeric(L2 %*% Evc)),
      col = "red")
dev.off()

pdf(file = paste0(figPath, "/u_i.pdf"), width = 5, height = 4)
par(mar = c(2,2,1,.5), mfrow = c(1,1))
plot(c(t(U_i[,c("u1", "u2")])), 
     L1 %*% Evc,
     xlim = range(c(t(u_i)), as.numeric(L1 %*% Evc)),
     ylim = range(c(t(u_i)), as.numeric(L1 %*% Evc)))
lines(range(c(t(u_i)), as.numeric(L1 %*% Evc)),
      range(c(t(u_i)), as.numeric(L1 %*% Evc)),
      col = "red")
dev.off()

write.csv(as.numeric(L2 %*% Evc), paste0(Eu_ij_path, "/Euij_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))
write.csv(as.numeric(L2 %*% Evc), paste0(resPath, "/Euij_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))

write.csv(as.numeric(L1 %*% Evc), paste0(Eu_i_path, "/Eui_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))
write.csv(as.numeric(L1 %*% Evc), paste0(resPath, "/Eui_estimates_sim", sprintf("%03d", as.integer(nSim)), ".csv"))



