# Function to calculate the Expected Influence measures for a DT-VAR model

EI_VAR <- function(phi){
  # One-step expected influence based on qgraph (no auto-regressive effects)
  ei1 <- colSums(phi) - diag(phi)
  # Two-step expected influence
  # Similar to networktools implementation but without auto-regressive effects
  # te2 <- apply(phi%*%phi,2,sum) - diag(phi%*%phi)
  
  # Essentially - total effects at lag 2 (excluding effects on itself)
  # does include the autoregressive parameter, but does not include paths
  # Y_1 -> Y_2 -> Y_1
  # does include
  # Y_1 -> Y_1 -> Y_2
  p <- ncol(phi)
  te2 <- rep(NA,p)
  for(i in 1:p){
    tmp <- rep(NA,p)
    for(j in 1:p){
      tmp[j] <- sum(phi[j,i] *phi[-i,j])
    }
    te2[i] <- sum(tmp)
  }
    # phi[1,1]*(phi[2,1] + phi[3,1]) +
    # phi[2,1]*(phi[3,2]+ phi[2,2]) +
    # phi[3,1]*(phi[2,3] + phi[3,3])
    # 
    #   phi[1,2]*(phi[1,1] + phi[3,1]) +
    #   phi[2,2]*(phi[1,2]+ phi[3,2]) +
    #   phi[3,2]*(phi[1,3] + phi[3,3])
    # 
  ei2 <- te2 + ei1
  
  if(is.null(dimnames(phi))){ names(ei1) <- names(ei2) <- paste0("node",1:nrow(phi))
     }else{ names(ei1) <- names(ei2) <- colnames(phi) }
     
  list(step1 = ei1, step2 = ei2)
  
}

