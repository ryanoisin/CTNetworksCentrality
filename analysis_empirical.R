# This files contians supplementary materials for Ryan & Hamaker (2020)
# 25 Feb 2020, o.ryan@uu.nl

# ----------------------------------------------------------------------------------
# -------------------- Load functions and pacakges ---------------------------------
# ----------------------------------------------------------------------------------

library(ctsem)
library(tibble)
library(qgraph)
library(RColorBrewer)

source("functions/fun_plots.R") # functions to aid plotting
source("functions/EI_VAR.R") # Expected Influence from VAR model
source("functions/IE.R") # CT indirect effect
source("functions/DE.R") # CT direct effect
source("functions/TE.R") # CT total effect
source("functions/ct_centrality.R") # 


# ----------------------------------------------------------------------------------
# -------------------------- Load data and process ---------------------------------
# ----------------------------------------------------------------------------------

# ---------- Load Raw Data -------------
# data available for download from https://osf.io/c6xt4/download #
rawdata <- read.csv("empirical_data/ESMdata.csv",header=TRUE, stringsAsFactors = FALSE)

# Extract time variable
t1 <- as.POSIXct(paste(rawdata$date,rawdata$resptime_s),format="%d/%m/%y %H:%M:%S")
time <- as.numeric(difftime(t1,t1[1], units="hours"))

# Select required variables and scale
data_esm <- rawdata[,c("se_selfdoub","phy_tired","mood_irritat","pat_restl")]
data_esm <- apply(data_esm, 2, scale)

# create ID variable
id <- rep(1,nrow(data_esm))

# Create long form dataset for ctsem
data_long <- cbind(id, time, data_esm)


# ----------------------------------------------------------------------------------
# ------------------------------- Model Estimation ---------------------------------
# ----------------------------------------------------------------------------------

set.seed(123)
# make variable labels
lab <- c("S","F","I","R")
p <- 4

# Specify CT-VAR model
ctmodel <- ctModel(type='stanct',
                   manifestNames= c("se_selfdoub","phy_tired","mood_irritat","pat_restl"),
                   latentNames= paste0(lab,"_l"),
                   LAMBDA = diag(nrow=p),
                   DRIFT = "auto",
                   MANIFESTMEANS = matrix(data=0, nrow=p, ncol=1),
                   MANIFESTVAR=diag(0,p),
                   CINT = "auto",
                   DIFFUSION ="auto")
ctfit <- ctStanFit(data_long, ctmodel, optimize=TRUE, cores=2)


# ---------- DT-VAR Model Estimation -------------
dtmodel <- ctModel(type='standt',
                   manifestNames= c("se_selfdoub","phy_tired","mood_irritat","pat_restl"),
                   latentNames=  paste0(lab,"_l"),
                   LAMBDA = diag(nrow=p),
                   DRIFT = "auto",
                   MANIFESTMEANS = matrix(data=0, nrow=p, ncol=1),
                   MANIFESTVAR=diag(0,p), 
                   CINT = "auto",
                   DIFFUSION ="auto")
dtfit <- ctStanFit(data_long, dtmodel, optimize=TRUE, cores=2)


# Get parameters from output

ctres <- getdrift(summary(ctfit))
dtres <- getdrift(summary(dtfit), mode = "DT")

# Save
output <- tibble::lst(time, ctmodel, ctfit, dtmodel, dtfit,
            ctres, dtres)
saveRDS(output, file = "estimates.RDS")

# Save session information for reproducibility
sink("SessionInfo.txt")
sessionInfo()
sink()

# ----------------------------------------------------------------------------------
# ---------------------------- Analysis and Figures --------------------------------
# ----------------------------------------------------------------------------------

# Load if necessary
res <- readRDS("estimates.RDS")
ctres <- res$ctres ; dtres <- res$dtres

# extract parameter matrices
drift <- ctres$drift
phi <- dtres$phi

#----------------------------------------------------
#---------- Figure 7: Histogram of TI ---------------
#----------------------------------------------------

# Make time-intervals between subsequent occassions
dtime <- time[-1] - time[-length(time)]
median(dtime)
dtime_quant <- dtime[dtime <= quantile(dtime, c(0.975))]

# for pdf scaling
sc <- .7

pdf("figures/hist_ti.pdf", 
    height = sc * 8, width = sc * 16)
hist(dtime_quant, 
     #main = "TI Distribution (97.5 percentile)",
     main = NULL,
     xlab = "Time-Interval", col = "#666666", 
     cex.lab  =  1.75, cex.main = 2.75, cex.axis = 1.5)
abline(v=median(dtime),  col = "#E6AB02", lty = 2, lwd = 3)
dev.off()

#----------------------------------------------------
#---------- Figure 8: CT and DT networks ------------
#----------------------------------------------------

# DT network
pdf("figures/CTnet_emp.pdf",width = 8, height = 8)
netplot(drift, labels = lab)
dev.off()

# CT network
pdf("figures/DTnet_emp.pdf",width = 8, height = 8)
netplot(phi, labels = lab, shape = "square")
dev.off()

#----------------------------------------------------
#---------- Figure 9: Phi(DT) plots -----------------
#----------------------------------------------------

# set range of Delta t
dts = seq(0,5,.1)

# get phi across these intervals
phidt.ar <- get_phidt_ar(drift, dts)

# define subsets of parameters
s1 <- which(drift == diag(drift), arr.ind = TRUE) # Auto-regressive effects
s2 <- rbind(c(1,2), c(1,3), c(1,4), c(2,3), c(2,4), c(3,4)) # Cross-lags
s3 <- rbind(c(2,1), c(3,1), c(3,2), c(4,1), c(4,2), c(4,3)) # Cross-lags

slist <- list(s1,s2,s3)

# helper part to define labels
l1 <- n_to_lab(s1)
l2 <- n_to_lab(s2)
l3 <- n_to_lab(s3)

lablist <- list(l1,l2,l3)

# Start loop of making figures
for(i in 1:3){
  smat <- slist[[i]]
  
   pdf(paste0("figures/phidt_emp_", i,"_lab.pdf"), 
       height = sc * 8, width = sc * 11)
  
   # set up layout
  lmat <- matrix(c(1,2),1,2, byrow = TRUE)
  lo <- layout(mat = lmat, 
               heights = c(1),
               widths = c(1,.33))
  
  if(i ==1){ ylim = c(-0.5,1) } else{ ylim = c(-0.2,0.3)  }
  
  # Make the plot
  par(mar = c(5.3,4.1,4.1,.5))
  phidt_plot_subset(input = phidt.ar, dts = dts, cols = cols, 
                    smat = smat,
                    xlines = FALSE,
                    ylim = ylim)
  
  # Make legend
  par(mar = c(.1,.1,.1,.1))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
  
  lmat <- lablist[[i]]
  legtext <- as.vector(unlist(apply(lmat, 1, function(row){
    bquote(phi[.(paste0(row[1], row[2]))](Delta*t))
  })))
  
  legend("center",
         as.expression(legtext),
         lty = smat[, 1], col = cols[smat] , lwd = lwd,
         cex = cex.leg)
  
  
  dev.off()
} # end of for loop


#----------------------------------------------------
#---------- Figure 10: CT Centrality plots ----------
#----------------------------------------------------

ctcent_ar <- array(0, dim = c(length(dts),p,3))

# ct centrality
for(k in 1:length(dts)){
  temp <- ct_centrality(drift, dt=dts[k])
  ctcent_ar[k,,1] <- temp$TotalEffectCentrality
  ctcent_ar[k,,2] <- temp$IndirectEffectCentrality
}

# dt centrality
ei <- EI_VAR(phi)
between <- centrality(t(phi))$Betweenness

# Set up plotting parameters
cex.axis = cex.leg = 2
cex.lab  = 2.25
labT <- paste0(lab, "(t)")


pdf("figures/cent_emp_total.pdf",width = sc*19, height = sc*8)
cent_compare(ctmat =  ctcent_ar[,,1] , dtcent = ei$step2, ctype = NULL,
             dts = dts, p = p, cols = cols, labels = labT, xat2 = c(-0.4,0,0.4), 
             xlim = c(-.5,.5),
             cex.axis = cex.axis, cex.leg = cex.leg, cex.main = cex.main, cex.lab = cex.lab, 
             lwd = lwd,
             linex= linex, liney = liney,  cex.cline = 2, cex.point = lwd)
dev.off()

pdf("figures/cent_emp_indirect.pdf",width = sc*19, height = sc*8)
cent_compare(ctmat =  ctcent_ar[,,2] , dtcent = between , ctype = NULL,
             dts = dts, p = p, cols = cols, labels = labT, xat2 = c(0,1,2),
             cex.axis = cex.axis, cex.leg = cex.leg, cex.main = cex.main, cex.lab = cex.lab, 
             lwd = lwd,
             linex= linex, liney = liney,  cex.cline = 2, cex.point = lwd)
dev.off()



#----------------------------------------------------
#---------- Table: Parameter Estimates---------------
#----------------------------------------------------

res <- readRDS("estimates.RDS")
ctres <- res$ctres ; dtres <- res$dtres

# round to three places for table
ctres$drift <- round(ctres$drift,3) 
ctres$lower <- round(ctres$lower,3)
ctres$upper <- round(ctres$upper,3)

cmat <- matrix(as.character(ctres$drift),4,4)

for(i in 1:4){
  cmat[i,] <- paste0(cmat[i,]," (",ctres$lower[i,],", ", 
                     ctres$upper[i,],")")
}
lab <- c("SD","Ti","Ir","Rl")
dimnames(cmat) = list(lab,lab)
cmat
xtable::xtable(cmat)



dtres$phi <- round(dtres$phi,3) 
dmat <- matrix(as.character(dtres$phi),4,4)

for(i in 1:4){
  dmat[i,] <- paste0(dmat[i,]," (",round(dtres[[3]][i,],3),", ", 
                     round(dtres[[4]][i,],3),")")
}
lab <- c("SD","Ti","Ir","Rl")
dimnames(dmat) = list(lab,lab)
dmat
xtable::xtable(dmat)

