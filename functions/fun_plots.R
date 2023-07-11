# o.ryan@uu.nl; December 2019

# This file contains functions to aid in plotting DT and CT networks and parameters

# ----------------------------------------------------------------------
# ---------- Plot matrix in network form  --------------------
# ----------------------------------------------------------------------

netplot <- function(mat, greyscale = FALSE, maximum = .5, asize = 6, edge.labels = TRUE,
                    edge.label.cex = 2, fade = FALSE, shape = "circle",
                    labels = TRUE,
                    vsize = 20,
                    esize = 12){
  
  layout <- rbind(c(0,1), 
                  c(1,1), 
                  c(1,0),
                  c(0,0))
  
  if(isTRUE(labels)){ labels = c("X1", "X2", "X3", "X4") }
  
  m_lty <- matrix(1, 4, 4)
  m_lty[mat<0] <- 2 
  
  m_col <- matrix("blue",4,4)
  m_col[mat > 0 ] <- "firebrick2"
  if(greyscale){
    qgraph::qgraph(t(mat), 
                   layout = layout,
                   directed = T, 
                   edge.color = "darkgrey",
                   edge.labels = edge.labels,
                   edge.label.cex = edge.label.cex,
                   edge.label.color = "darkgrey",
                   # curved = FALSE, 
                   lty = t(m_lty), 
                   vsize = vsize, 
                   esize = esize,
                   asize= asize,
                   # color = cols,
                   mar = c(8, 10, 8, 8), maximum=maximum,
                   fade = fade,
                   shape = shape,
                   maximum = maximum,
                   labels = labels)
  } else{
    qgraph::qgraph(t(mat),
                   edge.color = t(m_col),
                   layout = layout,
                   directed = T, 
                   edge.labels = edge.labels,
                   edge.label.cex = edge.label.cex,
                   # curved = FALSE, 
                   lty = t(m_lty), 
                   vsize = 20, 
                   esize = 12,
                   asize= asize,
                   # color = cols,
                   mar = c(8, 10, 8, 8), maximum=maximum,
                   fade = fade,
                   shape = shape,
                   maximum = maximum,
                   labels = labels)
  }
}



# ----------------------------------------------------------------------
# ---------- Calculate Phi(Dt) for a range of Dt  ----------------------
# ----------------------------------------------------------------------

get_phidt_ar <- function(drift,dts){
  p <- ncol(drift)
  out <- vapply(
    lapply(
      dts,
      FUN = function(inter)
        as.matrix(expm::expm(drift * inter))
    ), 
    identity, matrix(0, p, p))
  return(out)}


# ----------------------------------------------------------------------
# ---------- Plot Phi as function of Dt  -------------------------------
# ----------------------------------------------------------------------


phidt_plot_subset <- function(input, dts, lwd = 4, cols = cols, smat,
                              xl = 1.5, yl = 1.06, liney = 2.5, linex = 3.1,
                              cex.axis = 1.5,
                              cex.main  = 2.75,
                              cex.lab = 1.75,
                              titletext = "",
                              font.main = 1,
                              vline = c(0.5,1),
                              cex.leg = 1.5,
                              legtf = FALSE,
                              maintf = FALSE,
                              xlines = TRUE,
                              ylim = NULL){
  # needs phidt_array as input type
  
  if(is.null(ylim)){  ylim = c(min(input),max(input)) }
 
  plot.new()
  plot.window(xlim = range(dts), ylim = ylim)
  axis(1, cex.axis = cex.axis)
  axis(2, cex.axis = cex.axis)
  title(ylab = "Parameter Value", line=liney, cex.lab = cex.lab)
  title(xlab = expression(paste("Time Interval (", Delta, "t)")), line=linex, cex.lab = cex.lab)
  title(main = titletext, cex.main = cex.main, font.main = 1)
  
  if(xlines == TRUE){
  abline(v = vline[1], lty = 2); abline(v = vline[2], lty = 2) 
  }
  for(i in 1:nrow(smat)){
    lines(y=input[smat[i,1],smat[i,2],],x=dts,lty=smat[i,1],col=cols[smat[i,1],smat[i,2]],lwd=lwd)
  }
  abline(h =0)
  
  if(isTRUE(legtf)){
    # Make a legend
    legtext <- as.vector(unlist(apply(smat, 1, function(row){
      bquote(Phi[.(paste0(row[1], row[2]))])
    })))
    
    
    legend(x = xl, y = yl,
           as.expression(legtext),
           lty = smat[, 1], col = cols[smat] , lwd = lwd,
           cex = cex.leg)
  }
}

# ----------------------------------------------------------------------
# ---------- Plot TE and DE Intervention Trajectories ------------------
# ----------------------------------------------------------------------

trajplot <- function(drift, 
                     x = 2, 
                     m = c(1,3), 
                     delta.t = 1, 
                     m.int = 0, 
                     x.int = 1,
                     lwd = 3,
                     lwdp = 2,
                     ymin = -0.75,
                     xmin = -0.05,
                     cex.axis = 1.5,
                     cex.leg = 1.5,
                     cex.main = 2.75,
                     cex.lab  = 1.75,
                     xaxis = seq(0,1,length = 5),
                     yaxis = seq(-0.5,1,length = 4),
                     leg = TRUE,
                     type = "TE",
                     cols = NULL,
                     return_values = FALSE
){
  
  # ------ Preliminary set up -----------
  
  n <- nrow(drift)
  
  if(is.null(cols)){
    cols <- RColorBrewer::brewer.pal(4, "Dark2")
  }
  
  # Assign mediators
  M <- m
  S <- diag(nrow(drift))
  S[M,M] <- 0 
  
  # Set vector of starting values
  start <- rep(0,n)
  start[x] <- x.int
  
  # Get CT matrix under indirect effects intervention (CT)
  if(type == "TE"){
    drift_tilde <- drift
  } else{
    drift_tilde <- S%*%drift%*%S
  }

  
  # Set up the time vector
  dts2 <- seq(0,delta.t,.001)
  
  # Set up storage vectors
  y_ct <- matrix(0,length(dts2),n+1)
  
  
  # ------ Generate the trajectories  -----------
  
  # Get Trajectories for CT and total (DT2) effects
  
  for(i in 1:length(dts2)){
    # Total Effect
    y_ct[i,1:n] <- expm(drift_tilde*dts2[i])%*%start
    y_ct[i,n+1] <- dts2[i]
  }
  
  
  # -------------------------------------------
  # ------ Plot trajectories ------------------
  # -------------------------------------------
  
  # Tune allows us to see overlapping lines more clearly
  tune <- .01
  
  
  # Set up layout
 if(!leg){
    layout(mat = 1, 
           heights = 1,
           widths = 1)
    par(mar = c(5.1,4.1,4.1,.5))
    main1 <- ""
  } else if(leg){
    lmat <- matrix(c(1,2),1,2, byrow = TRUE)
    lo <- layout(mat = lmat, 
                 heights = c(1,1),
                 widths = c(1,.25))
    par(mar = c(5.1,4.1,4.1,.5))
    main2 <- ""
  }
  
  
  
  # TE
  plot.new()
  plot.window(xlim=c(xmin, delta.t), ylim=c(ymin,1))
  abline(h = 0)
  
  if(type == "DE"){
    # move mediators around a little for better visualization
    y_ct <- y_ct ; y_ct[,M[1]] <- y_ct[,M[1]] + tune ; y_ct[,M[2]] <- y_ct[,M[2]] - tune 
  }
  
  for(c in 1:n) lines(y = y_ct[,c], x = y_ct[,n+1], type = "l", col = cols[c], lwd = lwd)
  axis(1, at = xaxis,  cex.axis = cex.axis)
  axis(2, at = yaxis, cex.axis = cex.axis)
  
  points(y = y_ct[1,x], x = y_ct[1, n+1], pch =
           5, cex = 2, col = cols[x], lwd = lwdp)
    
  if(type == "DE"){
    for(j in 1:length(m)){
    # start intervention points
    points(y = y_ct[1,m[j]], x = 0, pch =
             18, cex = 3, col = cols[m[j]], lwd = lwdp)
    # end intervention points
    points(y = y_ct[nrow(y_ct),m[j]], x = y_ct[nrow(y_ct),n+1], pch =
             18, cex = 3, col = cols[m[j]], lwd = lwdp)
    }
  }
  
    title(xlab = "Time (t)", ylab = "Variable Values", main = main1,
          font.main = 1, cex.main = cex.main, cex.lab = cex.lab)
  
  if(leg){
    par(mar=c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim=c(0, 1), ylim = c(0, 1))
      legend("left",
             c(expression(Y[1](t)),
               expression(Y[2](t)),
               expression(Y[3](t)),
               expression(Y[4](t))), 
             pch = 19, col=cols, cex = cex.leg, yjust = 0, xjust = 0)
    }
  if(return_values){
  return(y_ct)
  }
}

# ----------------------------------------------------------------------
# ------- Plot Trajectories under Press Interventions ------------------
# ----------------------------------------------------------------------


pressplot <- function(drift, M, press, start = NULL, # arguments to pass to simPress
                      xaxis = seq(0,1,length = 5),
                      yaxis = seq(-0.5,1,length = 4),
                      delta.t = 1,
                      lwd = 3,
                      lwdp = 2,
                      ymin = -1,
                      xmin = -0.05,
                      cex.axis = 1.5,
                      cex.leg = 1.5,
                      cex.main = 2.75,
                      cex.lab  = 1.75,
                      cols = NULL,
                      leg = TRUE
){
  
  n <- nrow(drift)
  
  if(is.null(cols)){
    cols <- RColorBrewer::brewer.pal(4, "Dark2")
  }
  
  # Set up the time vector
  dts2 <- seq(0,delta.t,.001)
  
  # # Set up storage vectors
  y_ct <- matrix(0,length(dts2),n+1)
  
  # Get Trajectories for CT and total (DT2) effects
  
  for(i in 1:length(dts2)){
    # Total Effect
    y_ct[i,1:n] <- simPress(drift,dt=dts2[i],M,press, start)
    y_ct[i,n+1] <- dts2[i]
  }
  
  # Set up layout
  
  if(!leg){
    layout(mat = 1, 
           heights = 1,
           widths = 1)
    par(mar = c(5.1,4.1,4.1,.5))
    main1 <- ""
  } else if(leg){
    lmat <- matrix(c(1,2),1,2, byrow = TRUE)
    lo <- layout(mat = lmat, 
                 heights = c(1,1),
                 widths = c(1,.25))
    par(mar = c(5.1,4.1,4.1,.5))
    main2 <- ""
  }
  
  # Make plot
  
  plot.new()
  plot.window(xlim=c(xmin, delta.t), ylim=c(ymin,1))
  abline(h = 0)
  for(c in 1:4) lines(y = y_ct[,c], x = y_ct[,n+1], type = "l", col = cols[c], lwd = lwd)
  axis(1, at = xaxis,  cex.axis = cex.axis)
  axis(2, at = yaxis, cex.axis = cex.axis)
  
  # start intervention points
  points(y = y_ct[1,M], x = 0, pch =
           18, cex = 3, col = cols[M], lwd = lwdp)
  # end intervention points
  points(y = y_ct[nrow(y_ct),M], x = y_ct[nrow(y_ct),n+1], pch =
           18, cex = 3, col = cols[M], lwd = lwdp)
  
  title(xlab = "Time (t)", ylab = "Variable Values", main = main2,
        font.main = font.main, cex.main = cex.main, cex.lab = cex.lab)
  
  if(leg){
    par(mar=c(0, 0, 0, 0))
    plot.new()
    plot.window(xlim=c(0, 1), ylim = c(0, 1))
    legend("left",
           c(expression(Y[1](t)),
             expression(Y[2](t)),
             expression(Y[3](t)),
             expression(Y[4](t))), 
           pch = 19, col=cols, cex = cex.leg, yjust = 0, xjust = 0)
  }
  
}


# ----------------------------------------------------------------------
# ------- Plot Effects as Function of Time-Interval --------------------
# ----------------------------------------------------------------------

# Function to extract the TE, IE and DE

trajIE <- function(drift,
                   x = 2,
                   y = 4,
                   m = c(1,3),
                   x.int = 1,
                   delta.t = 1){
  n <- nrow(drift)
  start <- rep(0,n)
  start[x] <- x.int
  
  # Set up the time vector
  dts2 <- seq(0,delta.t,.001)
  
  out_TE <- out_DE <- matrix(0,nrow = length(dts2), ncol = n)
  for(i in 1:nrow(out_TE)){
    out_TE[i,] <- TE(drift,IV= 1:n, DV = 1:n,dt =dts2[i])%*%start
    out_DE[i,] <- DE(drift,IV= 1:n, DV = 1:n, M =m, dt =dts2[i])%*%start
  }
  
  out <- cbind(out_TE[,y], out_DE[,y], out_TE[,y] - out_DE[,y], dts2 )
  colnames(out) <- c("TE","DE","IE","dt")
  out
}

# Function to plot the effects over time
IEplot <- function(drift,
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
                   cex.axis = 1.5,
                   cex.leg = 1.5,
                   cex.main = 2.75,
                   cex.lab  = 1.75,
                   lwd = 3){
  
  out <- trajIE(drift,
                x = 2,
                y = 4,
                m = c(1,3),
                x.int = 1,
                delta.t = 1)
  
  
  lmat <- matrix(c(1,2),1,2, byrow = TRUE)
  lo <- layout(mat = lmat, 
               heights = c(1,1),
               widths = c(1,.25))
  par(mar = c(5.1,4.1,4.1,.5))
  main2 <- ""
  
  plot.new()
  plot.window(xlim=c(xmin, delta.t), ylim=c(ymin,1))
  abline(h = 0)
  
  lines(y = out[,1], x = out[,4], type = "l", col = cols[1], lwd = lwd)
  lines(y = out[,2], x = out[,4], type = "l", col = cols[2], lwd = lwd)
  lines(y = out[,3], x = out[,4], type = "l", col = cols[3], lwd = lwd)
  axis(1, at = xaxis,  cex.axis = cex.axis)
  axis(2, at = yaxis, cex.axis = cex.axis)
  
  title(xlab = "Time (t)", ylab = "Effect Size", main = main2,
        font.main = 1, cex.main = cex.main, cex.lab = cex.lab)
  
  par(mar=c(0, 0, 0, 0))
  plot.new()
  plot.window(xlim=c(0, 1), ylim = c(0, 1))
  legend("left",
         c(expression(paste("TE(",Delta,"t)",sep = "")),
           expression(paste("DE(",Delta,"t)",sep = "")),
           expression(paste("IE(",Delta,"t)",sep = ""))), 
         pch = 19, col=cols, cex = cex.leg, yjust = 0, xjust = 0)
  
}

# ----------------------------------------------------------------------
# ------- Extract Drift matrix parameters from ctsem model summary -----
# ----------------------------------------------------------------------

getdrift <- function(sumobj, mode = "CT", byrow = TRUE){
  ind <- which(sumobj$parmatrices$matrix=="DRIFT")
  p = sqrt((ind[length(ind)] - ind[1])+1)
  drift_df <-sumobj$parmatrices[ind,]
  drift <- matrix(drift_df[,"Mean"],p,p, byrow = byrow)
  lower <- matrix(drift_df[,"2.5%"],p,p, byrow = byrow)
  upper <- matrix(drift_df[,"97.5%"],p,p, byrow = byrow)
  drift_sig <- drift
  drift_sig[!(sign(lower) == sign(upper))] <- 0
  
  if(mode == "CT"){
    out <- tibble::lst(drift_df,drift,lower,upper,drift_sig)
  }
  if(mode == "DT"){
    out <- list(phi_df = drift_df,phi = drift,lower,upper,phi_sig = drift_sig)
  }
  return(out)
}

# ----------------------------------------------------------------------
# --------------------- Plot Centrality Comparisons --------------------
# ----------------------------------------------------------------------

cent_compare <- function(ctmat, dtcent, dts, p, ctype, cols, labels = NULL, xat2 = NULL,
                         cex.axis = 1.5, cex.leg = 1.5, cex.main = 2.75, cex.lab = 1.75, lwd = 4,
                         linex= 3.1, liney = 2.5,
                         cex.cline = 1, cex.point = 5,
                         xlim2 = NULL){
  
  if(is.null(labels)){ labels <- LETTERS[1:p] }
  
  lwd = 4
  linex = 3.1
  liney = 2.5
  font.main = 1
  
  lmat <- matrix(c(1,2,3),1,3, byrow = TRUE)
  lo <- layout(mat = lmat, 
               heights = c(1),
               widths = c(1,1,.33))   
  par(mar= c(5.1,4.5, 4.1, 2.1))
  
  mm <- max(abs(range(ctmat)))
  ax <- pretty(seq(-mm,mm, length = 3))
  
  plot.new()
  plot.window(xlim=range(dts), ylim = range(ax))
  abline(h = 0)
  for(c in 1:p) lines(y = ctmat[,c], x=dts, type = "l", col = cols[c,c], lwd = 4)
  axis(1, cex.axis = cex.axis)
  axis(2, ax, cex.axis = cex.axis)
  title(xlab = expression(paste("Time Interval (", Delta, "t)")), 
        main = "", cex.main  = cex.main,
        cex.lab = cex.lab,  font.main = font.main)
  title(ylab = expression(paste("Centrality value")),
        line = liney,
        cex.lab = cex.lab)
  title(main = ctype,
        cex.main = cex.lab)
  
  # Plot DT centrality
  if(is.null(xlim2)){
    xlim2 <- range(dtcent)
    if(min(xlim2) < -1 | max(xlim2) > 1){
      xlim2 <- xlim2
      xlim2[1] <- xlim2[1] -.1
      xlim2[2] <- xlim2[2] +.1
    } else { xlim2 = range(ax) }
    
    xlim2[1] <- xlim2[1] -.1
    xlim2[2] <- xlim2[2] +.1
  }
  plot.new()
  plot.window(xlim=xlim2, ylim = c(1,p))
  lines(y = 1:p, x = dtcent, cex = cex.cline)
  points(y = 1:p, x = dtcent, col = cols[1,1:p], cex = cex.point, pch = 20)
  if(!is.null(xat2)){axis(1, at = xat2, cex.axis = cex.axis, cex.lab = cex.lab)} else {
    axis(1, cex.axis = cex.axis, cex.lab = cex.lab)}
  
  abline(v=0,lty = 2, cex.axis = cex.axis, cex.lab = cex.lab)
  axis(2,at = 1:p, labels = labels, cex.axis = cex.axis, cex.lab = cex.lab)
  
  title(xlab = "Centrality Value", 
        cex.lab = cex.lab)
  title(ylab = "Variable", 
        cex.lab = cex.lab, line = liney)
  
  
  par(mar = c(.1,.1,.1,.1))
  plot.new()
  plot.window(xlim=c(-1, 1), ylim = c(-1, 1))
  
  legend("center",
         labels,
         lty = 1, col = cols[1,1:p] , lwd = lwd,
         cex = cex.leg)
  
  
}

# ----------------------------------------------------------------------
# ------------------------- Create Legend Labels -----------------------
# ----------------------------------------------------------------------
n_to_lab <- function(s1){
  matrix(dplyr::recode(s1,
                       "1" = "S",
                       "2" = "F",
                       "3" = "I",
                       "4" = "R"),nrow(s1),2)
}

