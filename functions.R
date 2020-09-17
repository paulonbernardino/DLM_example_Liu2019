###### Functions for DLM example from Liu et al. 2019 translated to R
### September 17th 2020
### Paulo Negri Bernardino

## Forward filtering function
forwardFilteringM <- function(M){
  Y <- M$Y
  X <- M$X
  rseas <- M$rseas
  delta <- M$delta
  Prior <- M$prior
  period <- 365.25/16
  deltrend <- delta[1]; delregn <- delta[2]; delseas <- delta[3]; delvar <- delta[4]
  Ftrend <- c(1,0)
  ntrend <- length(Ftrend)
  Gtrend <- matrix(c(1,0,1,1), nrow=2)
  itrend <- c(0:(ntrend-1)) 
  nregn <- ncol(X)
  Fregn <- rep(0, nregn)
  Gregn <- matrix(0, ncol = nregn, nrow = nregn); diag(Gregn)<-1
  iregn <- c((ntrend):(ntrend+nregn-1)) 
  pseas <- length(rseas)
  nseas <- pseas*2
  iseas <- c((ntrend+nregn):(ntrend+nregn+nseas-1))
  Fseas <- rep(c(1,0), pseas)
  
  Gseas <- matrix(0, ncol = nseas, nrow = nseas)
  c<-c(); s<-c()
  for (j in 1:pseas){
    c[j] <- cos(2*pi*rseas[j]/period)
    s[j] <- sin(2*pi*rseas[j]/period)
  }
  Gseas <- matrix(c(c[1], -s[1], 0, 0,
                    s[1], c[1], 0, 0,
                    0, 0, c[2], -s[2],
                    0, 0, s[2], c[2]), nrow = 4, ncol = 4)
  
  Fi <- c(Ftrend, Fregn, Fseas)
  G <- bdiag(Gtrend, Gregn, Gseas); G <- as.matrix(G)
  m <- M$prior[[1]]; C <- M$prior[[2]]; S <- M$prior[[3]]; nu <- M$prior[[4]]
  
  Ti <- length(Y)
  sm <- rep(0, length(m))
  sCi <- rep(0,nrow(C))
  sC <- list()
  for (i in 1:nrow(C)){
    sC[[i]] <- sCi
  }
  sS <- c(0)
  snu <- c(0)
  slik <- c(0)
  
  for (t in 1:Ti){
    a <- G %*% m
    R1 <- G %*% C
    R <- R1 %*% t(G) 
    # add 1 to indices because it's R, not python
    R[(itrend+1), (itrend+1)] <- R[(itrend+1), (itrend+1)] / deltrend
    R[(iregn+1), (iregn+1)] <- R[(iregn+1), (iregn+1)] / delregn
    R[(iseas+1), (iseas+1)] <- R[(iseas+1), (iseas+1)] / delseas
    nu <- delvar*nu
    Fi[(iregn+1)] <- X[t,]
    
    A1 <- R %*% Fi
    Q <- (t(Fi) %*% A1) + S
    A <- c(A1)/Q
    f <- t(Fi) %*% a
    y <- Y[t]
    
    if (!is.na(y)){
      e <- y-f
      ac <- (nu + e^2/Q) / (nu+1)
      rQ <- sqrt(Q)
      mlik <- dt(e/rQ, nu) /rQ # prob. dist. function
      m <- a + A*e
      C <- as.numeric(ac) * (R - (A %*% t(A) * as.numeric(Q)) )
      nu <- nu+1
      S <- ac*S
    }else{
      m <- a
      C <- R
      if (t<(Ti-1)){
        X[t+1,1] <- f
      }
      mlik <- NA
    }
    
    # assigning results
    sm <- cbind(sm, m)
    for (i in 1:nrow(C)){
      sC[[i]] <- cbind(sC[[i]], C[,i])
    }
    snu <- c(snu, nu)
    sS <- c(sS, S)
    slik <- c(slik, mlik)
  }
  return(list(sm, sC, snu, slik))
}

Index_low <- function(N, date0, percentile){
  percentile <- percentile
  interval <- 16
  dd <- seq(date0, as.Date("2015-12-03"), by = 16)
  dd_num <- as.numeric(dd)
  idq <- !is.na(tsN)
  tt1 <- as.numeric(date0):as.numeric(as.Date("2015-11-17")) # denser time series
  tt1_num <- tt1[-length(tt1)]
  N_idq <- N[idq]
  time_idq <- dd_num[idq]
  
  lm_itp <- approxfun(x=time_idq, y=N_idq, method="linear") 
  nn_itp <- lm_itp(tt1_num) # estimating NDVI for denser ts
  
  yday <- as.numeric( strftime( as.Date(tt1_num), format = "%j" )) # day of year
  
  ndvi_mean<-c()
  for (i in 1:365){
    avg <- mean(nn_itp[yday==i])
    ndvi_mean <- c(ndvi_mean, avg)
  }
  ndvi_sd<-c()
  for (i in 1:365){
    sd <- sd(nn_itp[yday==i])
    ndvi_sd <- c(ndvi_sd, sd)
  }
  
  if (length(ndvi_mean)==365){ndvi_mean <- c(ndvi_mean, ndvi_mean[length(ndvi_mean)])}
  ndvi_sd <- c(ndvi_sd, ndvi_sd[length(ndvi_sd)])
  
  tt2 <- seq(dd[1], dd[length(dd)], by="days")
  tt2 <- tt2[-length(tt2)]
  tt2_num <- as.numeric(tt2)
  yday2 <- as.numeric( strftime( as.Date(tt2_num), format = "%j" )) # day of year
  
  nv <- qnorm(1-(1-percentile)/2) #percentage points function
  
  lowboundary <- c()
  for (i in 1:length(tt2)){
    lb <- ndvi_mean[(yday2[i])] - nv * ndvi_sd[(yday2[i])]
    lowboundary <- c(lowboundary, lb)
  }
  
  index_low <- c()
  for (i in 1:length(dd)){
    if (!is.na(N[i])){
      if(N[i] < lowboundary[tt2_num==dd_num[i]]){
        index_low <- c(index_low, i)
      }
    }
  } 
  return(index_low)
}