#################
## CH-sim-main.R
#################
##
## Descriptuion:  Aggregate all results from simulations.
##

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

#########################
## Folders for storage ##
#########################

## current directory
currDir <- getwd()
setwd(currDir)

## folders
resPath <- paste("./results", sep = "")
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)
figPath <- paste0(resPath, "figures/")

simResPath <- paste("./results/sim/", sep = "")
simFigPath <- paste0(simResPath, "/figures/")
psi_path <- paste0(simResPath, "/all_psi/")
Ev_path <- paste0(simResPath, "/all_Ev/")

Eu_ijk_path <- paste0(simResPath, "/all_Euijk/")
Eu_ij_path <- paste0(simResPath, "/all_Euij/")
Eu_i_path <- paste0(simResPath, "/all_Eui/")

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

# Get all psi estimates
files <- list.files(path = psi_path, 
                    pattern = "\\.csv$", 
                    full.names = TRUE)

# Read all files into a list of data frames
data_list <- lapply(files, read.csv)



# Optional: combine into one data frame
all_psi <- sapply(data_list, function(X){X[,2]})


# Get all estimated Eu_ijk
files <- list.files(path = Eu_ijk_path, 
                    pattern = "\\.csv$", 
                    full.names = TRUE)

# Read all files into a list of data frames
data_list <- lapply(files, read.csv)

all_u <- sapply(data_list, function(X){X[,2]})

## prediction plots
nMonths <- 1:10 # months for prediction
pdf(file = paste0(simFigPath, "/pred_accuracy.pdf"), width = 12, height = 4)
par(mar = c(2,2.5,.5,.5), mfrow = c(2,5))
for (nMonth in nMonths) {
  # tps_pred <- max(tps) + 365/12*nMonth
  tps_pred <- max(tps) + 12*nMonth
  
  # ijk_sampled <- unique(ijkt[, c("i", "j", "k")])
  ijkt_pred <- cbind(ijk_sampled[rep(1:nrow(ijk_sampled), each = length(tps_pred)),],
                     t = rep(1:length(tps_pred), times = nrow(ijk_sampled)))
  
  X0_pred <- cbind(1, as[ijkt_pred$k, c("a", "s")],
                   tps_pred[ijkt_pred$t], as[ijkt_pred$k, c("a", "s")]*tps_pred[ijkt_pred$t])
  
  Z0_pred <- bdiag(lapply(1:nrow(ijk_sampled), function(r){
    cbind(1, tps_pred[ijkt_pred[which(ijkt_pred$i == ijk_sampled[r,"i"] & 
                                        ijkt_pred$j == ijk_sampled[r,"j"] &
                                        ijkt_pred$k == ijk_sampled[r,"k"]), "t"]])
  }))
  
  
  Y_true <- as.numeric(X0_pred %*% b.true + Z0_pred %*% u.true)
  
  
  
  b.est <- head(all_psi[,1], 6)
  u.est <- all_u[,1]
  
  Y_est <- as.numeric(X0_pred %*% b.est + Z0_pred %*% u.est)
  
  plot(Y_true, 
       Y_est,
       xlim = range(Y_true, 
                    Y_est),
       ylim = range(Y_true, 
                    Y_est),
       cex.axis = 1.5,
       cex.lab = 1.5)
  lines(range(Y_true, 
              Y_est),
        range(Y_true, 
              Y_est),
        col = "red",
        lwd = 3)
  text(x = .05*min(Y_true, 
                    Y_est),
       y = .8*max(Y_true, 
                    Y_est),
       labels = paste0(nMonth, " years"),
       font = 2,
       cex = 1.5)
}
dev.off()

## UQ on vector parameter
pdf(file = paste0(simFigPath, "/psi_est.pdf"), width = 5, height = 4)
par(mar = c(4.5,5.5,1,1))
plot(replicate(ncol(all_psi), psi.true), 
     all_psi,
     xlim = range(psi.true, all_psi),
     ylim = range(psi.true, all_psi),
     xlab = expression(psi),
     ylab = expression(hat(psi)),
     cex.axis = 2,
     cex.lab = 2,
     pch = 21,
     cex = 2,
     lwd = 2)
lines(range(psi.true, all_psi),
      range(psi.true, all_psi),
      lwd = 4,
      col = "red")
dev.off()

## UQ on conditional expectations
pdf(file = paste0(simFigPath, "/Eu_est.pdf"), width = 5, height = 4)
par(mar = c(4.5,5.5,1,1))
plot(u.true,
     all_u[,1],
     cex.axis = 2,
     cex.lab = 2,
     xlab = expression(paste(E, "[", u, ";",psi[true], "]")),
     ylab = expression(paste(E, "[", u, ";",hat(psi), "]")),
     pch = 21,
     cex = 2,
     lwd = 2)
lines(range(u.true,
     u.est),
     range(u.true,
     u.est),
     col = "red",
     lwd = 3)
dev.off()

## export results
save.image(paste0(simResPath,"/sim-aggregated.RData"))

