####### Functions

iwpQ <- function(v,x) {
  n  <- nrow(x)
  z  <- matrix(rnorm(v*n),v,n) %*% chol(x)
  return(chol2inv(chol(crossprod(z))))
}

###### Data aug

augmentData <- function(Y, l, tau, d, lambda, gamma, delta) {
  
  # Y       N+1 dim data matrix with Date as first column
  # l       lag length
  # tau     controls prior on own 1st lags
  # d       decay for higher lags
  # lambda  prior for the constant
  # gamma   sum of coefficients unit roots
  # delta   co-integration prior
  
  N   <- ncol(Y)-1  # Number of variables
  
  # Data matrix
  X  <- matrix(1, nrow=nrow(Y), ncol=1)
  for (i in 1:l) {
    X <- cbind(X, lag(Y[,-1], i)) 
  }
  
  X    <- as.matrix(X[-(1:l),])
  Y    <- as.matrix(Y[-(1:l),-1])
  T    <- length(Y[,1])
  
  # mean of the data
  mu   <- colMeans(Y)
  dmu  <- diag(mu)
  
  sd <- matrix(0,1,N)
  # Individual AR(2) regressions to get s1, s2
  for (j in 1:N) {
    Xa <- X[,c(1, seq(1+j,1+N*l,N))]
    XX <- chol2inv(chol(crossprod(Xa)))
    e  <- Y[,j] - Xa %*% XX %*% crossprod(Xa, Y[,j])
    sd[1,j] <- sqrt(sum(e^2)/(T-2))
  }
  ds <- diag(c(sd))
  
  ######## Dummy observation priors
  
  # specify dummy observations for first lag
  yd1 <- ds/tau
  xd1 <- cbind(0, ds/tau, matrix(0, N, N*(l-1)))
  
  # specify dummies for higher lags
  yd2 <- matrix(0, N*(l-1), N)
  xd2 <- matrix(0, N*(l-1), 1+N*l)
  for (j in 1:(l-1))
  {
    xd2[((j-1)*N+1):(j*N), (2+j*N):(1+(j+1)*N)] <- ds*((j+1)^d)/tau
  }
  
  # specify priors for the constants
  yd3 <- matrix(0, 1, N)
  xd3 <- cbind(matrix(lambda, 1, 1), matrix(0, 1, N*l))
  
  # specify dummy variables for covariance matrix
  yd4 <- ds
  xd4 <- matrix(0, N, N*l+1)
  
  # specify sum of coefficient dummies
  yd5 <- gamma*dmu
  xd5 <- cbind(0, matrix(rep(dmu,l), N, N*l))
  
  # specify common stochastic trend dummies
  yd6 <- matrix(delta*mu, 1, N)
  xd6 <- delta*cbind(1, matrix(mu, 1, N*l))
  
  # Data and all dummy observations
  Ystar <- rbind(Y, yd1, yd2, yd3, yd4, yd5, yd6)
  Xstar <- rbind(X, xd1, xd2, xd3, xd4, xd5, xd6)
  
  return(list(Ystar,Xstar))
  
}

##################

Gibbs_estimate <- function(Y, X, reps, burn, ce, nf){
  
  N <- ncol(Y)
  T <- nrow(Y)
  l <- (ncol(X)-1)/N
  m <- N*(l+2)+2
  
  # compute posterior mean
  bhat   <- chol2inv(chol(crossprod(X))) %*% crossprod(X,Y)
  # compute initial value of sigma
  e      <- Y - X %*% bhat
  sigma  <- crossprod(e)/T
  
  # vectorise bhat, calculate (X'X)^-1
  M    <- cbind(c(bhat))
  iXs  <- chol2inv(chol(crossprod(X)))
  
  # Store Gibbs samples of forecasts
  fstore <- matrix(0,reps-burn,nf*N)
  bstore <- matrix(0,reps-burn,N*(N*l+1))
  sstore <- matrix(0,reps-burn,N)
  
  # Gibbs loop
  pb <- progress_bar$new(total=reps, clear=TRUE)
  
  for (i in 1:reps) { 
    
    pb$tick()
    
    # Draw VAR coeffs
    V     <- kronecker(sigma, iXs)
    b     <- M + t(chol(V)) %*% rnorm(N*(N*l+1))
    
    # Draw sigma
    bs    <- matrix(b, nrow=N*l+1, ncol=N, byrow=FALSE)
    e     <- Y - X%*%bs
    sca   <- chol2inv(chol(crossprod(e)))
    sigma <- iwpQ(T, sca)
    
    if (i>burn) { 
      if(ce > 0) {
        bstore[i-burn,] <- c(bs)
        sstore[i-burn,] <- diag(sigma)
      }
      if(nf > 0) {
        # forecast GDP growth and inflation for nf periods
        yhat       <- matrix(0,nf+l,N)
        yhat[1:l,] <- Y[(T-l+1-m):(T-m),]
        for (j in (l+1):(nf+l)) {
          ldata    <- c(1, c(t(yhat[(j-1):(j-l),])))
          yhat[j,] <- ldata %*% bs + t(rnorm(N))%*%chol(sigma)
        }
        fstore[i-burn,] <- c(yhat[-(1:l),])
      }
    }
  }
  
  return(list(bstore, sstore, fstore))
  
}

############################

fan_chart <- function(Y, out, controls, nb){
  
  # Dimensions
  N   <- ncol(Y)-1 # No variables
  nf  <- ncol(out)/N # No forecast periods
  
  # Variable names 
  nms <- colnames(Y)[-1]
  
  # Dates
  fdt <- seq.Date(tail(Y$Date,1), by="quarter", length.out=nf+1)
  
  # Colours for lines/centres
  col_tail <- "grey95"
  if (N <= 8) {
    centre <- brewer.pal(n = ifelse(N<3, 3, N), name = "Dark2")  
  } else {centre <- rep(c("green", "red", "purple"), round(N/3))}
  
  # Quantiles
  qv  <- c(.05,.2,.35,.65,.8,.95)
  lq  <- length(qv)
  
  # Colors for fans sections
  col <- NULL
  for (k in 1:N) {
    CRP <- c(col_tail,centre[k],col_tail)
    col <- c(col, colorRampPalette(CRP)(lq+1)[2:lq])
  }
  
  ##### Stack-band version #########################################
  
  mm <- out %>% 
    apply(2, quantile, probs=qv) %>% 
    t() %>% 
    as_tibble() %>% 
    mutate(Variable = rep(nms, each=nf), F = rep(1:nf, N)) %>% 
    pivot_longer(cols = -c(Variable, F)) %>% 
    unite(names, c(Variable, name)) %>% 
    pivot_wider(names_from = names, values_from = value) %>% 
    select(-F) %>% 
    rbind(rep(c(as.matrix(tail(Y[,-1],1))), each=lq), .) %>%
    mutate(Date = fdt) 
  
  mma <- mm %>% 
    select(-ends_with("_95%")) %>% 
    rename_all(~ c(paste(rep(nms, each=(lq-1)), letters[1:(lq-1)], sep="_"), "Date")) %>%
    pivot_longer(cols=-Date, names_to = "Band", values_to = "LV")
  
  mmb <- mm %>% 
    select(-ends_with("_5%")) %>% 
    rename_all(~ c(paste(rep(nms, each=(lq-1)), letters[1:(lq-1)], sep="_"), "Date")) %>%
    pivot_longer(cols=-Date, names_to = "Band", values_to = "UV")
  
  jmm <- left_join(mma, mmb) %>% 
    mutate(name = gsub("_.*", "", Band))
  
  hist <- tail(Y, nb) %>% 
    pivot_longer(cols=-Date)
  
  bind_rows(hist, jmm) %>% 
    ggplot() + 
    geom_rect(aes(xmin=max(Date[!is.na(value)], na.rm=TRUE),
                  xmax=max(Date, na.rm=TRUE), 
                  ymin=-Inf, ymax=Inf), 
              fill=col_tail, alpha=.1) +
    geom_ribbon(aes(x=Date, ymin=LV, ymax=UV, fill=Band), alpha=.8) +
    geom_line(aes(x=Date, y=value, color=name)) + 
    scale_fill_manual(values=col) +
    scale_color_manual(values=centre) +
    scale_x_date(expand = c(0,0)) +
    theme_minimal() +
    facet_grid(~ name) + 
    theme(legend.position="none") +
    labs(title="Dummy variable BVAR model", 
         subtitle=controls, 
         y="Percentage points", x="")
  
}

###############

coeff_plot <- function(Y, l, b, s, si, controls){
  ###### Plot estimated param densities #####
  # Retrieve constants/AR1/AR2/sigma
  
  N   <- ncol(Y)-1
  nms <- colnames(Y)[-1]
  
  #Sigma
  sc <- s %>%
    as_tibble(.name_repair = ~ paste0("sigma^2 -", nms)) %>%
    mutate(iter = 1:nrow(s), w = "Variance" ) %>% 
    pivot_longer(cols=-c(w,iter))
  # Constant
  pc <- b[,((0:(N-1))*(N*l+1)+1)] %>% 
    as_tibble(.name_repair = ~ paste0("alpha -", nms)) %>%
    mutate(iter = 1:nrow(b), w = "Constant" ) %>% 
    pivot_longer(cols=-c(w,iter))
  # AR1
  ac <- b[,((0:(N-1))*(N*l+1)+1)+(1:N)] %>%
    as_tibble(.name_repair = ~ paste("rho[1] -", nms)) %>%
    mutate(iter = 1:nrow(b), w = "AR1" ) %>% 
    pivot_longer(cols=-c(w,iter))
  
  csequences <- bind_rows(sc, pc, ac)
  
  if (l > 1) {
    #AR2
    dc <- b[,((0:(N-1))*(N*l+1)+1)+2*(1:N)] %>%
      as_tibble(.name_repair = ~ paste0("rho[2] -", nms)) %>%
      mutate(iter = 1:nrow(b), w = "AR2" ) %>% 
      pivot_longer(cols=-c(w,iter))
    csequences <- bind_rows(csequences, dc)
  }
  
  p <- list()
  # All as histograms/densities
  p[[1]] <- csequences %>%
    ggplot() + 
    geom_histogram(aes(x=value, y=after_stat(density), fill=w), alpha=0.33, bins=150) +
    geom_density(aes(x=value, colour=w)) +
    facet_wrap(~ name, labeller=label_parsed, ncol=N, scales="free_y") +
    theme_minimal() +
    theme(legend.position="none") + 
    labs(title="Coefficients: Histograms and kernel-smoothed estimated densities", 
         subtitle=controls, x="", y="")
  
  # Plot selected sequences
  if(si > 0) {
    p[[2]] <- csequences %>%
      filter(iter < si) %>% 
      ggplot() + 
      geom_line(aes(x=iter, y=value, colour=name), show.legend=FALSE) +
      facet_wrap(~name, labeller=label_parsed, ncol=N, scales="free") +
      theme_minimal() + 
      labs(title=paste("Sequences from the Gibbs Sampler, first",si,"iterations"), 
           subtitle=controls, x="", y="")
  }
  
  return(p)
}