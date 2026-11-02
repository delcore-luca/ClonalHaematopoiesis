######################
## CH-data-mAVG-main.R
######################
##
## Descriptuion:  This R script generates the results of the averaged model,
##                using the model specific results obtaibed from CH-data-main.R
##

# Load functions
source("./source/CH-fun.R")

# Set folders path
folder_path <- "./results/mBICs/"
figPath <- "./results/mAVGfigures/"
resPath <- "./results/mAVGres/"
ifelse(!dir.exists(file.path(figPath)), dir.create(file.path(figPath)), FALSE)
ifelse(!dir.exists(file.path(resPath)), dir.create(file.path(resPath)), FALSE)

# Get all .txt files
files <- list.files(path = folder_path, pattern = "\\.txt$", full.names = TRUE)
t(t(files))

# Read all files into a list
data_list <- lapply(files, read.table, header = FALSE, sep = "", stringsAsFactors = FALSE)


mBICs <- as.numeric(sapply(data_list, function(l){unlist(l)[2]}))
sprintf("%.100f",mBICs)


w <- exp(-.5*(mBICs - min(mBICs)))/sum(exp(-.5*(mBICs - min(mBICs))))

pdf(file = paste0(figPath, "/ws.pdf"), width = 7, height = 3)
par(mar = c(3,3,.5,0), mfrow = c(1,1))
# Create barplot without x-axis labels
bp <- barplot(
  w[order(w, decreasing = T)],
  border = NA,               # remove bar borders
  names.arg = FALSE,         # suppress default labels
  col = "grey",
  cex.axis = 1.5
)

# Add rotated labels at 45 degrees
text(
  x = bp,                    # bar midpoints
  y = par("usr")[3] - 0.005,  # slightly below x-axis
  labels = paste0("M", 1:length(w))[order(w, decreasing = T)],
  srt = 45,                  # rotation angle
  adj = 1,                   # adjust text position
  cex = 1.5,
  xpd = TRUE                  # allow text outside plot region
)
dev.off()

# Get all psi estimates
psi_path <- "./results/all_psi/"
psi_files <- list.files(path = psi_path, pattern = "\\.csv$", full.names = TRUE)
data_list <- lapply(psi_files, read.csv, header = TRUE, sep = ",", stringsAsFactors = FALSE)
all.psi.est <- lapply(data_list, function(X){X[,2]})
names(all.psi.est) <- paste0("M", 1:length(all.psi.est))

## candidate models
allComb <- as.matrix(expand.grid(rep(list(c(0, 1)), 4)))
colnames(allComb) <- c("0a", "0s", "ra", "rs")
allMods <- cbind("01" = 1, allComb[,c(1,2)], 
                 "r1" = 1, allComb[,c(3,4)])
allMods <- allMods[nrow(allMods):1,]

allMods_wpsi <- cbind(allMods, matrix(data = 1, nrow = 16, ncol = 8))
all_psi <- matrix(data = 0, nrow = 16, ncol = 14)

## calculate averaged model parameters
for (m in 1:16) {
  model <- allMods_wpsi[m,]
  all_psi[m, as.logical(model)] <- all.psi.est[[m]]
}

psi_avg <- colSums(all_psi*w)
all_psi_wavg <- rbind(all_psi,
                      psi_avg)

rownames(all_psi_wavg) <- c(1:16, "m")

sink(paste0(resPath, "/all_psi_wavg_table.txt"))
print(xtable::xtable(formatC(all_psi_wavg, format = "e", digits = 2)), 
      include.rownames = TRUE)
sink()

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
B0_lvls <- proc.data$B0_lvls
################################################

# averaged-model conditional moments
EVuc_avg <- get.EVu(psi = psi_avg,
                    X = X0,
                    Z = Z0 %*% B0,
                    Y = Y0,
                    R = R0,
                    B_lvls = B0_lvls)

Ev_c <- EVuc_avg$Eu

Eu_c <- matrix(data = as.numeric(B0 %*% EVuc_avg$Eu),
               ncol = 2,
               byrow = TRUE)

beta_opt <- head(psi_avg, 6)
names(beta_opt) <- c("01", "0a", "0s", 
                     "r1", "ra", "rs")

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

sink(paste0(resPath, "/ijk_table.txt"))
print(xtable::xtable(cbind(cbind(ijk_wNames_wEst[1:50, 4:6],
                                 formatC(as.matrix(ijk_wNames_wEst[1:50, 7:8]), format = "e", digits = 3)),
                           cbind(ijk_wNames_wEst[51:100, 4:6],
                                 formatC(as.matrix(ijk_wNames_wEst[51:100, 7:8]), format = "e", digits = 3)))), 
      include.rownames = FALSE)
sink()

## prediction time grid (in months)
tps_pred <- seq(0, max(tps) + 120, length.out = 500)

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
  
  Z0_r <- cbind(1, tps_pred)
  
  plot(x = tps[ijktyr[which(ijktyr$i == ijk_r$i &
                              ijktyr$j == ijk_r$j &
                              ijktyr$k == ijk_r$k), c("t")]],
       y = ijktyr[which(ijktyr$i == ijk_r$i &
                          ijktyr$j == ijk_r$j &
                          ijktyr$k == ijk_r$k), c("VAF_01")],
       xlim = range(tps_pred),
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

L2 <- build_L2_sparse(ijk_sampled = ijk_sampled, p = 2)
ij_sampled <- unique(ijk_sampled[, c("i", "j")])
Euij_c <- matrix(data = as.numeric(L2 %*% Ev_c),
                 ncol = 2,
                 byrow = TRUE)

Z0_r <- cbind(1, tps_pred)

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
ij_wNames <- IasJKt[, c("i", "j", "gene", "cDNA_change")]
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
       xlim = range(tps_pred),
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

L1 <- build_L1_sparse(ijk_sampled = ijk_sampled, p = 2)
i_sampled <- unique(ijk_sampled[, c("i")])

Eui_c <- matrix(data = as.numeric(L1 %*% Ev_c),
                ncol = 2,
                byrow = TRUE)

Z0_r <- cbind(1, tps_pred)

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
       xlim = range(tps_pred),
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

sink(paste0(resPath, "/i_table.txt"))
print(xtable::xtable(cbind(i_wNames_wEst[,2],
                           formatC(as.matrix(i_wNames_wEst[,c(3:4)]), format = "e", digits = 3))), 
      include.rownames = FALSE)
sink()

############################################################
## results for participants with at least 5 clones recorded
############################################################

ijk_wNames_wEst_g5clones <- ijk_wNames_wEst[which(ijk_wNames_wEst$k %in% which(table(ijk_sampled$k) >= 5)),]

pdf(file = paste0(figPath, "/ijks_kwg5clones_sorted.pdf"), width = 15, height = 10)
layout(matrix(1:30, nrow = 6, ncol = 5, byrow = FALSE),
       heights = c(2, 2, 1, 2, 2, 1))
for (k in unique(ijk_wNames_wEst_g5clones$k)) {
  ijk_wNames_wEst_g5clones_k <- ijk_wNames_wEst_g5clones[ijk_wNames_wEst_g5clones$k == k,]
  Y_k <- sapply(1:nrow(ijk_wNames_wEst_g5clones_k), function(r){
    ijk_r <- ijk_wNames_wEst_g5clones_k[r,]
    
    return(ijk_r[,"x0"]/(ijk_r[,"x0"] + (1 - ijk_r[,"x0"])*exp(-ijk_r[,"r"]*tps_pred)))
    
  })
  
  XY_k <- Reduce(rbind, lapply(1:nrow(ijk_wNames_wEst_g5clones_k), function(r){
    ijk_r <- ijk_wNames_wEst_g5clones_k[r,]
    
    xy_r <- cbind(x = tps[ijktyr[which(ijktyr$i == ijk_r$i &
                                         ijktyr$j == ijk_r$j &
                                         ijktyr$k == ijk_r$k), c("t")]],
                  y = ijktyr[which(ijktyr$i == ijk_r$i &
                                     ijktyr$j == ijk_r$j &
                                     ijktyr$k == ijk_r$k), c("VAF_01")],
                  r = r)
    
    return(xy_r)
    
  }))
  
  par(mar = c(2,2.5,1,.5))
  plot(x = XY_k[,"x"], 
       y = XY_k[,"y"], 
       xlim = range(tps_pred),
       ylim = range(c(XY_k[,"y"], Y_k)),
       pch = 20,
       cex = 2,
       lty = 1, 
       lwd = 2,
       cex.axis = 1.5,
       cex.lab = 1.5,
       col = XY_k[,"r"],
       main = paste0("ID-", unique(ijk_wNames_wEst_g5clones_k$indiv_ID),
                     "-a", unique(variant_dataset[which(variant_dataset$indiv_ID == unique(ijk_wNames_wEst_g5clones_k$indiv_ID)), c("age")]),
                     "-s", ifelse(unique(variant_dataset[which(variant_dataset$indiv_ID == unique(ijk_wNames_wEst_g5clones_k$indiv_ID)), c("sex")]) == "Male", "M", "F"))
  )
  matplot(tps_pred, Y_k, 
          type = 'l', 
          lty = 1, 
          lwd = 3, 
          col = 1:ncol(Y_k), 
          add = TRUE)
  abline(v = max(tps), col = "grey", lwd = 4, lty = 2)
  
  
  mat <- Y_k/rowSums(Y_k)
  
  colnames(mat) <- apply(ijk_wNames_wEst_g5clones_k[,4:6], 1, paste, collapse = "_")
  
  # Number of rows and columns
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  # Compute cumulative sums for stacking
  mat_cum <- t(apply(mat, 1, cumsum))
  mat_cum <- cbind(0, mat_cum)  # add 0 for bottom of first area
  
  # Colors
  cols <- 1:nc
  
  # Create empty plot
  par(mar = c(2,2.5,1,.5))
  plot(tps_pred, 
       tps_pred,
       cex.axis = 1.5,
       cex.lab = 1.5,
       type = "n",
       ylim = range(mat_cum),
       xlab = "Time", 
       ylab = "Proportion")
  
  # Draw polygons for each category
  for (j in 1:nc) {
    polygon(
      x = c(tps_pred, rev(tps_pred)),
      y = c(mat_cum[, j], rev(mat_cum[, j+1])),
      col = cols[j],
      border = NA
    )
  }
  abline(v = max(tps), col = "grey", lwd = 4, lty = 2)
  
  par(mar = c(.5,.5,.5,.5))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), type = "n",
       axes = FALSE, xlab = "", ylab = "", ann = FALSE, bty = "n")
  legend(x = "topleft",
         legend = colnames(mat),
         col = 1:ncol(mat),
         lwd = 10,
         bty = "n")
}
dev.off()


############################################################
## results for participants with at least 3 clones recorded
############################################################

ijk_wNames_wEst_g3clones <- ijk_wNames_wEst[which(ijk_wNames_wEst$k %in% which(table(ijk_sampled$k) >= 3)),]

pdf(file = paste0(figPath, "/ijks_kwg3clones_sorted.pdf"), width = 15, height = 10)
layout(matrix(1:30, nrow = 6, ncol = 5, byrow = FALSE),
       heights = c(2, 2, 1, 2, 2, 1))
for (k in unique(ijk_wNames_wEst_g3clones$k)) {
  ijk_wNames_wEst_g3clones_k <- ijk_wNames_wEst_g3clones[ijk_wNames_wEst_g3clones$k == k,]
  Y_k <- sapply(1:nrow(ijk_wNames_wEst_g3clones_k), function(r){
    ijk_r <- ijk_wNames_wEst_g3clones_k[r,]
    
    return(ijk_r[,"x0"]/(ijk_r[,"x0"] + (1 - ijk_r[,"x0"])*exp(-ijk_r[,"r"]*tps_pred)))
    
  })
  
  XY_k <- Reduce(rbind, lapply(1:nrow(ijk_wNames_wEst_g3clones_k), function(r){
    ijk_r <- ijk_wNames_wEst_g3clones_k[r,]
    
    xy_r <- cbind(x = tps[ijktyr[which(ijktyr$i == ijk_r$i &
                                         ijktyr$j == ijk_r$j &
                                         ijktyr$k == ijk_r$k), c("t")]],
                  y = ijktyr[which(ijktyr$i == ijk_r$i &
                                     ijktyr$j == ijk_r$j &
                                     ijktyr$k == ijk_r$k), c("VAF_01")],
                  r = r)
    
    return(xy_r)
    
  }))
  
  par(mar = c(2,2.5,1,.5))
  plot(x = XY_k[,"x"], 
       y = XY_k[,"y"], 
       xlim = range(tps_pred),
       ylim = range(c(XY_k[,"y"], Y_k)),
       pch = 20,
       cex = 2,
       lty = 1, 
       lwd = 2,
       cex.axis = 1.5,
       cex.lab = 1.5,
       col = XY_k[,"r"],
       main = paste0("ID-", unique(ijk_wNames_wEst_g3clones_k$indiv_ID),
                     "-a", unique(variant_dataset[which(variant_dataset$indiv_ID == unique(ijk_wNames_wEst_g3clones_k$indiv_ID)), c("age")]),
                     "-s", ifelse(unique(variant_dataset[which(variant_dataset$indiv_ID == unique(ijk_wNames_wEst_g3clones_k$indiv_ID)), c("sex")]) == "Male", "M", "F"))
  )
  matplot(tps_pred, Y_k, 
          type = 'l', 
          lty = 1, 
          lwd = 3, 
          col = 1:ncol(Y_k), 
          add = TRUE)
  abline(v = max(tps), col = "grey", lwd = 4, lty = 2)
  
  
  mat <- Y_k/rowSums(Y_k)
  
  colnames(mat) <- apply(ijk_wNames_wEst_g3clones_k[,4:6], 1, paste, collapse = "_")
  
  # Number of rows and columns
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  # Compute cumulative sums for stacking
  mat_cum <- t(apply(mat, 1, cumsum))
  mat_cum <- cbind(0, mat_cum)  # add 0 for bottom of first area
  
  # Colors
  cols <- 1:nc
  
  # Create empty plot
  par(mar = c(2,2.5,1,.5))
  plot(tps_pred, 
       tps_pred,
       cex.axis = 1.5,
       cex.lab = 1.5,
       type = "n", 
       # ylim = c(0,1),
       # xlim = c(0, max(tps_pred)*1.5),
       ylim = range(mat_cum),
       xlab = "Time", 
       ylab = "Proportion")
  
  # Draw polygons for each category
  for (j in 1:nc) {
    polygon(
      x = c(tps_pred, rev(tps_pred)),
      y = c(mat_cum[, j], rev(mat_cum[, j+1])),
      col = cols[j],
      border = NA
    )
  }
  abline(v = max(tps), col = "grey", lwd = 4, lty = 2)
  
  par(mar = c(.5,.5,.5,.5))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), type = "n",
       axes = FALSE, xlab = "", ylab = "", ann = FALSE, bty = "n")
  legend(x = "topleft",
         legend = colnames(mat),
         col = 1:ncol(mat),
         lwd = 10,
         bty = "n")
}
dev.off()

############################################################
## results for participants with at least 2 clones recorded
############################################################

ijk_wNames_wEst_g2clones <- ijk_wNames_wEst[which(ijk_wNames_wEst$k %in% which(table(ijk_sampled$k) >= 2)),]

pdf(file = paste0(figPath, "/ijks_kwg2clones_sorted.pdf"), width = 15, height = 10)
layout(matrix(1:30, nrow = 6, ncol = 5, byrow = FALSE),
       heights = c(2, 2, 1, 2, 2, 1))
for (k in unique(ijk_wNames_wEst_g2clones$k)) {
  ijk_wNames_wEst_g2clones_k <- ijk_wNames_wEst_g2clones[ijk_wNames_wEst_g2clones$k == k,]
  Y_k <- sapply(1:nrow(ijk_wNames_wEst_g2clones_k), function(r){
    ijk_r <- ijk_wNames_wEst_g2clones_k[r,]
    
    return(ijk_r[,"x0"]/(ijk_r[,"x0"] + (1 - ijk_r[,"x0"])*exp(-ijk_r[,"r"]*tps_pred)))
    
  })
  
  XY_k <- Reduce(rbind, lapply(1:nrow(ijk_wNames_wEst_g2clones_k), function(r){
    ijk_r <- ijk_wNames_wEst_g2clones_k[r,]
    
    xy_r <- cbind(x = tps[ijktyr[which(ijktyr$i == ijk_r$i &
                                         ijktyr$j == ijk_r$j &
                                         ijktyr$k == ijk_r$k), c("t")]],
                  y = ijktyr[which(ijktyr$i == ijk_r$i &
                                     ijktyr$j == ijk_r$j &
                                     ijktyr$k == ijk_r$k), c("VAF_01")],
                  r = r)
    
    return(xy_r)
    
  }))
  
  par(mar = c(2,2.5,1,.5))
  plot(x = XY_k[,"x"], 
       y = XY_k[,"y"], 
       xlim = range(tps_pred),
       ylim = range(c(XY_k[,"y"], Y_k)),
       pch = 20,
       cex = 2,
       lty = 1, 
       lwd = 2,
       cex.axis = 1.5,
       cex.lab = 1.5,
       col = XY_k[,"r"],
       main = paste0("ID-", unique(ijk_wNames_wEst_g2clones_k$indiv_ID),
                     "-a", unique(variant_dataset[which(variant_dataset$indiv_ID == unique(ijk_wNames_wEst_g2clones_k$indiv_ID)), c("age")]),
                     "-s", ifelse(unique(variant_dataset[which(variant_dataset$indiv_ID == unique(ijk_wNames_wEst_g2clones_k$indiv_ID)), c("sex")]) == "Male", "M", "F"))
  )
  matplot(tps_pred, Y_k, 
          type = 'l', 
          lty = 1, 
          lwd = 3, 
          col = 1:ncol(Y_k), 
          add = TRUE)
  abline(v = max(tps), col = "grey", lwd = 4, lty = 2)
  
  
  mat <- Y_k/rowSums(Y_k)
  
  colnames(mat) <- apply(ijk_wNames_wEst_g2clones_k[,4:6], 1, paste, collapse = "_")
  
  # Number of rows and columns
  nr <- nrow(mat)
  nc <- ncol(mat)
  
  # Compute cumulative sums for stacking
  mat_cum <- t(apply(mat, 1, cumsum))
  mat_cum <- cbind(0, mat_cum)  # add 0 for bottom of first area
  
  # Colors
  cols <- 1:nc
  
  # Create empty plot
  par(mar = c(2,2.5,1,.5))
  plot(tps_pred, 
       tps_pred,
       cex.axis = 1.5,
       cex.lab = 1.5,
       type = "n",
       ylim = range(mat_cum),
       xlab = "Time", 
       ylab = "Proportion")
  
  # Draw polygons for each category
  for (j in 1:nc) {
    polygon(
      x = c(tps_pred, rev(tps_pred)),
      y = c(mat_cum[, j], rev(mat_cum[, j+1])),
      col = cols[j],
      border = NA
    )
  }
  abline(v = max(tps), col = "grey", lwd = 4, lty = 2)
  
  par(mar = c(.5,.5,.5,.5))
  plot(NA, xlim = c(0, 1), ylim = c(0, 1), type = "n",
       axes = FALSE, xlab = "", ylab = "", ann = FALSE, bty = "n")
  legend(x = "topleft",
         legend = colnames(mat),
         col = 1:ncol(mat),
         lwd = 10,
         bty = "n")
}
dev.off()


######################
## Generate IC tables
######################

# Get all .txt files
cAICfiles <- list.files(path = "./results/cAICs/", pattern = "\\.txt$", full.names = TRUE)

# Read all files into a list
cAIC_data_list <- lapply(cAICfiles, read.table, header = FALSE, sep = "", stringsAsFactors = FALSE)
cAICs <- as.numeric(sapply(cAIC_data_list, function(l){unlist(l)[2]}))

# Get all .txt files
cBICfiles <- list.files(path = "./results/cBICs/", pattern = "\\.txt$", full.names = TRUE)

# Read all files into a list
cBIC_data_list <- lapply(cBICfiles, read.table, header = FALSE, sep = "", stringsAsFactors = FALSE)
cBICs <- as.numeric(sapply(cBIC_data_list, function(l){unlist(l)[2]}))

sink(paste0(resPath, "/ICs_table.txt"))
print(xtable::xtable(cbind(
  # 1:nrow(allMods),
  # ifelse(allMods, "YES", "NO"),
  allMods,
  cAICs, cBICs, mBICs), digits = 3))
sink()


########################
## Confidence Intervals
########################

X0_full <- X0

## Calculate Hessian for each candidate model
mll.H_allMods <- sapply(1:(nrow(all_psi_wavg) - 1), function(m){
  
  model <- allMods[m,]
  covNames <- c("01", "0a", "0s", 
                "r1", "ra", "rs")[as.logical(model)]
  
  X0 <- X0_full[,which(colnames(X0_full) %in% covNames)]
  
  psi_m <- all_psi_wavg[m,]
  psi_m <- psi_m[c(which(colnames(X0_full) %in% covNames), 7:14)]
  
  mll.H.m <- get.mll.hess(psi = psi_m,
                          X = X0,
                          Y = Y0,
                          Z = Z0 %*% B0,
                          R = R0,
                          B = B0)
  
  mll.H.m_full <- matrix(data = 0, nrow = 14, ncol = 14)
  
  mll.H.m_full[c(which(colnames(X0_full) %in% covNames), 7:14),
               c(which(colnames(X0_full) %in% covNames), 7:14)] <- as.matrix(mll.H.m)
  
  
  return(as.matrix(mll.H.m_full))
}, simplify = "array")

## Calculate Covariance matrix of the estimates for each candidate model
mll.Cov_allMods <- sapply(1:dim(mll.H_allMods)[3], function(k){
  
  model <- allMods[k,]
  covNames <- c("01", "0a", "0s", 
                "r1", "ra", "rs")[as.logical(model)]
  I_hat <- -mll.H_allMods[,,k]
  I_hat <- I_hat[c(which(colnames(X0_full) %in% covNames), 7:14), 
                 c(which(colnames(X0_full) %in% covNames), 7:14)]
  
  # Check positive definiteness
  eigenvalues <- eigen(I_hat,
                       symmetric = TRUE,
                       only.values = TRUE)$values
  
  if (min(eigenvalues) <= 0) {
    I_hat <- nearPD(I_hat)$mat
    # I_hat <- nearest_spd(I_hat)
    # warning("Observed information is not positive definite.")
  }
  
  # Covariance matrix and standard errors
  vcov_hat <- solve(I_hat)
  # vcov_hat <- pinv(I_hat)
  
  vcov_hat_full <- matrix(data = 0, nrow = 14, ncol = 14)
  vcov_hat_full[c(which(colnames(X0_full) %in% covNames), 7:14),
                c(which(colnames(X0_full) %in% covNames), 7:14)] <- as.matrix(vcov_hat)
  
  return(as.matrix(vcov_hat_full))
}, simplify = "array")


## Calculate Wald confidence intervals
CIs <- t(format(Reduce(cbind, lapply(1:17, function(k){
  if(k <= 16){
    ## Use covariance of the k-th candidate model
    Cov_psi_avg <- mll.Cov_allMods[,,k]
    psi_hat <- all_psi_wavg[k,]
  }else{
    ## Derive the covariance of the averaged model
    Cov_psi_avg <- w[1]*mll.Cov_allMods[,,1] + w[1]*(all_psi_wavg[1,] - all_psi_wavg["m",])%*%t(all_psi_wavg[1,] - all_psi_wavg["m",])
    for (j in 2:dim(mll.Cov_allMods)[3]) {
      Cov_psi_avg <- Cov_psi_avg + w[j]*mll.Cov_allMods[,,j] + w[j]*(all_psi_wavg[j,] - all_psi_wavg["m",])%*%t(all_psi_wavg[j,] - all_psi_wavg["m",])
    }
    psi_hat <- all_psi_wavg["m",]
  }
  
  isposdef(Cov_psi_avg)
  vcov_hat <- Cov_psi_avg
  se_hat <- sqrt(diag(vcov_hat))
  
  # 95% Wald confidence intervals
  ci_hat <- cbind(Estimate = psi_hat,
                  SE = se_hat,
                  Lower95 = psi_hat - 1.96 * se_hat,
                  Upper95 = psi_hat + 1.96 * se_hat)
  
  rownames(ci_hat) <- c("$\\beta_{10}$",
                        "$\\beta_{a0}$",
                        "$\\beta_{s0}$",
                        "$\\beta_{11}$",
                        "$\\beta_{a1}$",
                        "$\\beta_{s1}$",
                        "$\\sigma_0$",
                        "$\\sigma_1$",
                        "$\\log\\tau_{g1}$",
                        "$\\log\\tau_{g2}$",
                        "$\\log\\tau_{c1}$",
                        "$\\log\\tau_{c2}$",
                        "$\\log\\tau_{p1}$",
                        "$\\log\\tau_{p2}$")
  
  if(k <= 16){
    ci_hat[which(!as.logical(allMods[k,])),] <- NA
  }
  
  return(ci_hat[,c(3,4)])
})), scientific = TRUE, digits = 3))


## export results
sink(paste0(resPath, "/ConfidenceIntervals_table.txt"))
print(xtable::xtable(CIs),
      sanitize.text.function = identity)
sink()


pdf(file = paste0(figPath, "/allModelsCIs.pdf"), width = 7, height = 7)
par(mar = c(1,2,.7,.5))
lay <- matrix(c(matrix(1:18, nrow = 9, ncol = 2, byrow = FALSE),
                c(19, 20)), ncol = 2, byrow = TRUE)
layout(lay, heights = c(rep(1, 9), 0.5))
for (k in 1:16) {
  Cov_psi_avg <-mll.Cov_allMods[,,k]
  
  isposdef(Cov_psi_avg)
  
  vcov_hat <- Cov_psi_avg
  se_hat <- sqrt(diag(vcov_hat))
  
  psi_hat <- all_psi_wavg[k,]
  
  # 95% Wald confidence intervals
  ci_hat <- cbind(Estimate = psi_hat,
                  SE = se_hat,
                  Lower95 = psi_hat - 1.96*se_hat,
                  Upper95 = psi_hat + 1.96*se_hat)
  
  rownames(ci_hat) <- c("beta01",
                        "beta0age",
                        "beta0sex",
                        "r1",
                        "rage",
                        "rsex",
                        "s0",
                        "s1",
                        "taugene1",
                        "taugene2",
                        "taucDNA1",
                        "taucDNA2",
                        "tauparticipant1",
                        "tauparticipant2")
  
  
  ci_hat[which(!as.logical(allMods[k,])),] <- NA
  
  matplot(ci_hat[,c(1,3,4)], 
          type = 'l',
          col = c("black", "grey", "grey"),
          lty = c(1,2,2),
          lwd = 2,
          ylim = c(-43.65932, 16.08409),
          bty = "n",
          xaxt = 'n')
  title(main = k, adj = .05)
}

Cov_psi_avg <- w[1]*mll.Cov_allMods[,,1] + w[1]*(all_psi_wavg[1,] - all_psi_wavg["m",])%*%t(all_psi_wavg[1,] - all_psi_wavg["m",])
for (k in 2:dim(mll.Cov_allMods)[3]) {
  Cov_psi_avg <- Cov_psi_avg + w[k]*mll.Cov_allMods[,,k] + w[k]*(all_psi_wavg[k,] - all_psi_wavg["m",])%*%t(all_psi_wavg[k,] - all_psi_wavg["m",])
}

vcov_hat <- Cov_psi_avg
se_hat <- sqrt(diag(vcov_hat))
psi_hat <- all_psi_wavg["m",]

# 95% Wald confidence intervals
ci_hat <- cbind(Estimate = psi_hat,
                SE = se_hat,
                Lower95 = psi_hat - 1.96 * se_hat,
                Upper95 = psi_hat + 1.96 * se_hat)

matplot(ci_hat[,c(1,3,4)], 
        type = 'l',
        col = c("black", "grey", "grey"),
        lty = c(1,2,2),
        lwd = 2,
        ylim = c(-43.65932, 16.08409),
        bty = "n",
        xaxt = 'n')
title(main = "avg", adj = .05)

plot.new()
par(mfg = c(10, 1))
plot.new()
plot.window(xlim = c(1, 14), ylim = c(0, 1))
axis(1, at = 1:14, labels = FALSE, line = -3)
text(x = 1:14, y = par("usr")[3] + 1.5,
     labels = expression(beta[10], beta[a0], beta[s0],
                         beta[11], beta[a1], beta[s1],
                         sigma[0], sigma[1],
                         log(tau[g1]), log(tau[g2]),
                         log(tau[c1]), log(tau[c2]),
                         log(tau[p1]), log(tau[p2])),
     srt = 45, adj = 1,
     xpd = TRUE, cex = 1)

plot.new()
plot.window(xlim = c(1, 14), ylim = c(0, 1))
axis(1, at = 1:14, labels = FALSE, line = -3)
text(x = 1:14, y = par("usr")[3] + 1.5,
     labels = expression(beta[10], beta[a0], beta[s0],
                         beta[11], beta[a1], beta[s1],
                         sigma[0], sigma[1],
                         log(tau[g1]), log(tau[g2]),
                         log(tau[c1]), log(tau[c2]),
                         log(tau[p1]), log(tau[p2])),
     srt = 45, adj = 1,
     xpd = TRUE, cex = 1)
dev.off()

## export results
save.image(paste0(resPath,"/avgModelRes.RData"))
