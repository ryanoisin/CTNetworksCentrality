# This files contians supplementary materials for Ryan & Hamaker (2020)
# 25 Feb 2020, o.ryan@uu.nl

# -------------------------------------------------------------------------------------
# ------------------ Source functions and packages ------------------------------------
# -------------------------------------------------------------------------------------

library(qgraph)
library(RColorBrewer)
library(expm)
library(ggplot2)
library(combinat)
library(networktools)
library(animation)
require(xtable)

source("functions/fun_plots.R") # functions to aid plotting
source("functions/EI_VAR.R") # Expected Influence from VAR model
source("functions/IE.R") # CT indirect effect
source("functions/DE.R") # CT direct effect
source("functions/TE.R") # CT total effect
source("functions/ct_centrality.R") # CT centrality measures
source("functions/simPress.R")

figwd <- paste0(getwd(),"/figures/")

# -------------------------------------------------------------------------------------
# ----------------------- Preliminary set-up ------------------------------------------
# -------------------------------------------------------------------------------------

# Example Drift Matrix
p <- 4
drift <- matrix(c(-6,  1.25,     0,  5,
                  5.5, -2.5,     0,   0,
                  0,  5.9, -6,   0 ,
                  0, -7.3,     2.5, -6),4,4,byrow=T)

# Plotting parameters
cols <- brewer.pal(p, "Dark2")
cols <- matrix(cols, p , p ,byrow=T)

lwd = 4
cex.axis = cex.leg = 1.5
cex.main = 2.75
cex.lab  = 1.75
linex = 3.1
liney = 2.75
font.main = 1

# for pdf scaling
sc <- .7

# -------------------------------------------------------------------------------------
# ------------------------ Figure 1: DT Networks  -------------------------------------
# -------------------------------------------------------------------------------------

##### Get CT and DT networks #####

# Set up range of dt
steps = 0.01
dt.min = 0
dt.max = 2
p <- dim(drift)[1]
dts <- seq(dt.min,dt.max,steps)

# Get phi at a range of dts
phidt.ar <- get_phidt_ar(drift, dts)

# Choose two specific Delta t values to plot
dt_1 <- 51
dt_2 <- 101

phi1 <- phidt.ar[,,dt_1]
phi2 <- phidt.ar[,,dt_2]

# -------- Figure 1 ----------------

# DT(2) network
 pdf(paste0(figwd,"DTnet2_col.pdf"),width = 8, height = 8)
  netplot(phi2, asize = 10,  edge.labels = T, greyscale = F,
        maximum = 0.5, shape = "square",
        labels = c(expression(Y[1]),
                   expression(Y[2]),
                   expression(Y[3]),
                   expression(Y[4])))
 dev.off()

# DT(1) network
 pdf(paste0(figwd,"DTnet1_col.pdf"),width = 8, height = 8)
  netplot(phi1,asize = 10,  edge.labels = T, greyscale = F,
          maximum = 0.5, shape ="square",
          labels = c(expression(Y[1]),
                     expression(Y[2]),
                     expression(Y[3]),
                     expression(Y[4])))
 dev.off()


# -------------------------------------------------------------------------------------
# -------------------------- Figure 2: CT Network -------------------------------------
# -------------------------------------------------------------------------------------
  
# CT network
   pdf(paste0(figwd,"CTnet_col.pdf"),width = 8, height = 8)
  netplot(drift, asize = 10, edge.labels = T, greyscale = F, shape = "circle",
          labels = c(expression(Y[1]),
                     expression(Y[2]),
                     expression(Y[3]),
                     expression(Y[4])))
 dev.off()
 
  
  
  
 # -------------------------------------------------------------------------------------
 # ------------ In-text Section 2 - Path-effects and centrality  -----------------------
 # -------------------------------------------------------------------------------------
 
 # For Delta t  = 2
 # Expected Influence
  eip2 <- EI_VAR(phi2)
 # Betweenness
  bp2 <- qgraph::centrality(t(phi2))$Betweenness

  # For Delta t  = 1
  # Expected Influence
  eip1 <- EI_VAR(phi1)
  # Betweenness
  bp1 <- qgraph::centrality(t(phi1))$Betweenness
  
  # Make a table
  tab <- round(cbind(eip2$step2, eip2$step1, bp2, eip1$step2, eip1$step1, bp1),3)

  # ------ Table 1 --------------
  # xtable(tab, digits = 3)
  
  # -------------------------------------------------------------------------------------
  # -------------------------- Figure 3: Phi dt Plot ------------------------------------
  # -------------------------------------------------------------------------------------
  
  # Phi-dt plot of all parameters in three parts
  # First make selections of parameters to plot seperately
   s1 <- which(drift == diag(drift), arr.ind = TRUE) # Auto-regressive effects
   s2 <- rbind(c(1,2), c(1,4), c(2,1), c(3,2), c(4,2), c(4,3)) # Cross-lags, drift =/= 0
   s3 <- rbind(c(1,3), c(2,3), c(2,4), c(3,1), c(3,4), c(4,1)) # Cross-lags, drift = 0
   slist <- list(s1,s2,s3)
   

   
   # Start loop of making figures
    for(i in 1:3){
      smat <- slist[[i]]
   
   pdf(paste0(figwd,"phidt_threeparts_vert", i,".pdf"), 
       height = sc * 8, width = sc * 11)
   
   lmat <- matrix(c(1,2),1,2, byrow = TRUE)
   lo <- layout(mat = lmat, 
                heights = c(1),
                widths = c(1,.33))

   # Make the plot
    par(mar = c(5.3,4.1,4.1,.5))
    phidt_plot_subset(input = phidt.ar, dts = dts, lwd = 4, cols = cols, 
                     smat = smat,
                     xl = 1.5, yl = 1.06,
                     titletext = "",
                     liney = liney,
                     linex = linex,
                     cex.axis = cex.axis,
                     cex.main  = cex.main,
                     cex.lab = cex.lab,
                     cex.leg = cex.leg)
   # Make legend
   par(mar = c(.1,.1,.1,.1))
   plot.new()
   plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
   legtext <- as.vector(unlist(apply(smat, 1, function(row){
     bquote(phi[.(paste0(row[1], row[2]))](Delta*t))
   })))
   
   legend("center",
          as.expression(legtext),
          lty = smat[, 1], col = cols[smat] , lwd = lwd,
          cex = cex.leg)
   

   dev.off()
    } # end of for loop
   
# -------------------------------------------------------------------------------------
# ----------------- Supp Materials: Animated Figure -----------------------------------
# -------------------------------------------------------------------------------------

# Scaling tuning parameter for animation
sca <- 1.5

lags <- c(seq(0.5, 2, .01), seq(0.01, 0.49, .01))
   
par(mfrow = c(1, 2))
   
ani.options(
  ani.dev = 'pdf',
  ani.type = 'pdf',
  ani.width = sca * 13.6,
  ani.height = sca * 6.8,
  interval = .08,
  nmax = 100
)

sall <- as.matrix(expand.grid(1:4,1:4))
sall <- sall[!sall[,1] == sall[,2],]


colleg <-
  c(bquote(Phi[.(paste0(".", 1))]),  bquote(Phi[.(paste0(".", 2))]),
    bquote(Phi[.(paste0(".", 3))]),  bquote(Phi[.(paste0(".", 4))]))
typeleg <-
  c(bquote(Phi[.(paste0(1, "."))]),  bquote(Phi[.(paste0(2, "."))]),
    bquote(Phi[.(paste0(3, "."))]),  bquote(Phi[.(paste0(4, "."))]))


# Create PDF animated figure
saveLatex({
  par(mar = c(5, 6, 5, 6),mfrow=c(1,2))
     for(j in 1:length(lags)){
       phidt_plot_subset(input = phidt.ar, dts = dts, lwd = 4, 
                         cols = cols, 
                         smat = sall,
                         xl = 1.5, yl = 1.06,
                         titletext = "",
                         liney = liney,
                         linex = linex,
                         cex.axis = cex.axis,
                         cex.main  = cex.main,
                         cex.lab = cex.lab,
                         cex.leg = cex.leg,
                         legtf = FALSE, xlines = FALSE,
                         ylim = c(-.6,.6))
       abline(v = lags[j])
       
       legend("topright",
              as.expression(colleg),
              lty = 1, col = cols[1,1:4] , lwd = lwd,
              cex = cex.leg)
       
       legend("bottomright",
              as.expression(typeleg),
              lty = 1:4, col = "grey" , lwd = lwd,
              cex = cex.leg)
       
       phitmp <- expm(drift*lags[j])
       diag(phitmp) <- 0
       # par(mar = c(50, 50, 50, 50))
       netplot(phitmp, asize = 10, edge.labels = F, greyscale = F,
               maximum = 0.5, shape = "square",
               labels = c(expression(Y[1]),
                          expression(Y[2]),
                          expression(Y[3]),
                          expression(Y[4])))
     } 
   }, ani.basename = "BM", ani.opts = "controls,loop,width=0.8\\textwidth",
   latex.filename = paste0(figwd,"netanim_ex.tex"),
   img.name=paste0(figwd,"netanim"))




# -------------------------------------------------------------------------------------
# -------------------------- Figure 4: DT Centrality ----------------------------------
# -------------------------------------------------------------------------------------
  
# First calculate centrality across a range of delta t
  DT_total <- DT_direct <- DT_between <- matrix(0,length(dts),4)

# Loop through phi values at different dt
 for(i in 1:length(dts)){
   matrix <- phidt.ar[,,i]
   diag(matrix) <- 0
   eitemp  <- EI_VAR(matrix)
   DT_total[i,] <- eitemp$step2
   DT_direct[i,] <-eitemp$step1
   DT_between[i,] <- centrality(t(matrix))$Betweenness
 }
 
# set up plotting parameters
   cols <- brewer.pal(p, "Dark2")
   cols <- matrix(cols, p , p ,byrow=T)

# ------------ Plot Figure 4 panel a -----------------
   
pdf(paste0(figwd,"DTcent_deltat_vert1.pdf"), 
         height = sc * 8, width = sc * 11)
     
     lmat <- matrix(c(1,2),1,2, byrow = TRUE)
     lo <- layout(mat = lmat, 
                  heights = c(1),
                  widths = c(1,.33))   

plot.new()
  plot.window(xlim=c(0, max(dts)), ylim=c(-.6,max(DT_total)))
  abline(h = 0)
for(c in 1:4) lines(y = DT_total[,c], x = dts, type = "l", col = cols[c,c], lwd = 4)
  axis(1, cex.axis = cex.axis)
  axis(2, cex.axis = cex.axis, at = c(-.5,0,.5,1))
    title(xlab = expression(paste("Time Interval (", Delta, "t)")), 
        main = "", cex.main  = cex.main,
        cex.lab = cex.lab,  font.main = font.main)
    title(ylab = expression(paste("Centrality value")),
          line = liney,
          cex.lab = cex.lab)
    abline(v = 0.5, lty = 2); abline(v = 1, lty = 2) 
  
    par(mar = c(.1,.1,.1,.1))
    
plot.new()
  plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
  legtext <- list(bquote(Y[1](t)), bquote(Y[2](t)),
                    bquote(Y[3](t)), bquote(Y[4](t)))
    
  legend("center",
           as.expression(legtext),
           lty = 1, col = cols[1,1:4] , lwd = lwd,
           cex = cex.leg)
    
    
  dev.off()
    
    
# ----------- Panel b ---------------

pdf(paste0(figwd,"DTcent_deltat_vert2.pdf"), 
        height = sc * 8, width = sc * 11)
    lmat <- matrix(c(1,2),1,2, byrow = TRUE)
    lo <- layout(mat = lmat, 
                 heights = c(1),
                 widths = c(1,.33))   
    
plot.new()
  plot.window(xlim=c(0, max(dts)), ylim=c(-.6,max(DT_total)))
    abline(h = 0)
    for(c in 1:4) lines(y = DT_direct[,c], x = dts, type = "l", col = cols[c,c], lwd = 4)
    axis(1, cex.axis = cex.axis)
    axis(2, cex.axis = cex.axis, at = c(-.5,0,.5,1))
    title(xlab = expression(paste("Time Interval (", Delta, "t)")), 
          main = "", cex.main  = cex.main,
          cex.lab = cex.lab,  font.main = font.main)
    title(ylab = expression(paste("Centrality value")),
          line = liney,
          cex.lab = cex.lab)
    abline(v = 0.5, lty = 2); abline(v = 1, lty = 2) 
    
    par(mar = c(.1,.1,.1,.1))
    plot.new()
    plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
    legtext <- list(bquote(Y[1](t)), bquote(Y[2](t)),
                    bquote(Y[3](t)), bquote(Y[4](t)))
    
    legend("center",
           as.expression(legtext),
           lty = 1, col = cols[1,1:4] , lwd = lwd,
           cex = cex.leg)
dev.off()

# ----------- Panel c ---------------
    
pdf(paste0(figwd,"DTcent_deltat_vert3.pdf"), 
        height = sc * 8, width = sc * 11)
  lmat <- matrix(c(1,2),1,2, byrow = TRUE)
  lo <- layout(mat = lmat, 
                 heights = c(1),
                 widths = c(1,.33))   
  
  # select some dt values to plot betweenness for  
  sel1 <- sapply(seq(0,2,0.25),function(vec) which(dts == vec))
  DT_between_sel <- DT_between[sel1,]
  
  # set seed for jitter
  set.seed(1)
  
  plot.new()
    plot.window(xlim=c(0, max(dts)), ylim=c(-0.25,3.25))
    
    for(c in 1:4) lines(y = jitter(DT_between_sel[,c], amount = 0.02), dts[sel1], type = "b", col = cols[c,c],
                        lwd = 4, pch = 20, cex = 5)  
    axis(1, cex.axis = cex.axis)
    axis(2, seq(0,3,by = 0.5), cex.axis = cex.axis)
    title(xlab =expression(paste("Time Interval (", Delta, "t)")), 
          main = "", cex.main  = cex.main,
          cex.lab = cex.lab,  font.main = font.main)
    title(ylab = "Centrality Value", line = liney, cex.lab = cex.lab)
    abline(v = 0.5, lty = 2); abline(v = 1, lty = 2) 
    
    par(mar = c(.1,.1,.1,.1))
    plot.new()
    plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
    legtext <- list(bquote(Y[1](t)), bquote(Y[2](t)),
                    bquote(Y[3](t)), bquote(Y[4](t)))
    
    legend("center",
           as.expression(legtext),
           lty = 1, col = cols[1,1:4] , lwd = lwd,
           cex = cex.leg)
    dev.off()
    
 
# -------------------------------------------------------------------------------------
# --------------- Figure 5: CT Interventions Trajectories -----------------------------
# -------------------------------------------------------------------------------------

# plotting parameter
cex.lab = 1.5

# Figure 5 a ------- Pulse Intervention
  pdf(paste0(figwd,"InterventionTrajectories_impulse.pdf"),height = sc*(8), width = sc*(8))
  trajplot(drift,
           delta.t = 1,
           lwd = lwd,
           ymin = -1,
           xmin = -0.05,
           xaxis = seq(0,1,length = 5),
           yaxis = seq(-1,1,length = 5),
           leg = F,
           cols = NULL,
           type = "TE",
           cex.lab = cex.lab
  )
  dev.off()
  
 # Figure 5 b ------- Press Intervention
  
  pdf(paste0(figwd,"InterventionTrajectories_press.pdf"),height = sc*(8), width = sc*(8+2))
  pressplot(drift, M = 2, press = 1, start = NULL, # arguments to pass to simPress
                        xaxis = seq(0,1,length = 5),
                        yaxis = seq(-1,1,length = 5),
                        delta.t = 1,
                        lwd = 3,
                        lwdp = 2,
                        ymin = -1,
                        xmin = -0.05,
                        cex.lab  = cex.lab,
                        cols = NULL,
                        leg = T
  )
  dev.off()

  # Figure 5 c -------  Direct Effect Intervention
pdf(paste0(figwd,"InterventionTrajectories_DE.pdf"),height = sc*(8), width = sc*(8))
  trajplot(drift, ymin = -1.1, 
           xaxis = seq(0,1,length = 5),
           yaxis = c(-1,-.5,0,.5,1),
           lwd = lwd,
           leg = F,
           cols = NULL,
           type = "DE",
           cex.lab = cex.lab
  )
  dev.off()
  
  
  # Figure 5 d -------  Direct, Indirect and Total effects
  
  pdf(paste0(figwd,"InterventionTrajectories_IE.pdf"),height = sc*(8), width = sc*(8+2))
  IEplot(drift = drift,
         x = 2,
         y = 4,
         m = c(1,3),
         x.int = 1,
         delta.t = 1,
         cols = c("#66A61E", "#E6AB02", "#A6761D"),
         xmin = -0.05,
         ymin = -1,
         ymax = 1,
         xaxis = seq(0,1,length = 5),
         yaxis = seq(-1,1,length = 5),
         lwd = 3,
         cex.lab = cex.lab)
  dev.off()
 
# -------------------------------------------------------------------------------------
# --------------------- Figure 6: CT Centrality ---------------------------------------
# -------------------------------------------------------------------------------------
  

# Get CT centrality over time
ctcent_ar <- array(0, dim = c(length(dts), 4, 2))

for(i in 1:length(dts)){
  temp <- ct_centrality(drift, dt=dts[i])
  ctcent_ar[i,,1] <- temp$TotalEffectCentrality
  ctcent_ar[i,,2] <- temp$IndirectEffectCentrality
}


# ----- Make Figure -------

  # Loop over centrality measures
  for(j in 1:2){
  pdf(paste0(figwd,"CTcent_deltat_vert",j,".pdf"), 
      height = sc * 8, width = sc * 11)
  
  # Plot layout set-up
  lmat <- matrix(c(1,2),1,2, byrow = TRUE)
  lo <- layout(mat = lmat, 
               heights = c(1),
               widths = c(1,.33))   

  # Main plot
  plot.new()
  plot.window(xlim=c(0, 2), ylim=c(-1.2,1.2))
  abline(h = 0)
  for(c in 1:4) lines(y = ctcent_ar[,c,j], x=dts, type = "l", col = cols[c,c], lwd = 4)
  axis(1, seq(0,2,0.5), cex.axis = 1.25)
  axis(2, c(-1,0,1), cex.axis = 1.25)
  title(xlab = expression(paste("Time Interval (", Delta, "t)")), 
        main = "", cex.main  = cex.main,
        cex.lab = cex.lab,  font.main = font.main)
  title(ylab = expression(paste("Centrality value")),
        line = liney,
        cex.lab = cex.lab)
  
  # Plot Legend
  par(mar = c(.1,.1,.1,.1))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
  legtext <- list(bquote(Y[1](t)), bquote(Y[2](t)),
                  bquote(Y[3](t)), bquote(Y[4](t)))
  
  legend("center",
         as.expression(legtext),
         lty = 1, col = cols[1,1:4] , lwd = lwd,
         cex = cex.leg)
  
  
  dev.off()
  }
  
