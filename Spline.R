#library(splines2)
library(tidyverse)

ex1 <- TRUE
ex1 <- FALSE

if(ex1) {
  X     <- c(-0.23, -0.1, -0.05, -0.025, 0.001, 0.02, 0.05, 0.11, 0.22)
  Y     <- c( 24, 21, 18.75, 17.25, 15.75, 15, 14.25, 12.95, 12.5)
  xaxis <- seq(-0.4, 0.4, 0.01)
} else {
  X     <- c(27., 32., 37., 42., 47., 52., 57., 62., 67., 72., 77., 82., 87., 92.)
  Y     <- c(3.89, 2.45, 2.49, 3.81, 6.34, 10.49, 15.94, 26.91, 40.74, 59.55, 86.02, 145.42, 172.15, 230.80)/1000.
  xaxis <- seq(20, 97, 1)
}

D <- tibble(X=X, Y=Y)

ggplot(D) +
  geom_point(aes(x=X, y=log(Y)), shape=23, fill="green") +
  theme_light()

sf <- splinefun(X, y=Y, method = "natural")
ss <- spline(X, y=Y, method = "natural")

###################################
# a binary value for whether we want natural or clamped
###################################

cubic_spline = function(X,Y,nc) {

  # set nc == 0 for clamped at zero, anything else for natural

  # Define differences
  n  <- length(X)-1
  dy <- Y[-1] - Y[-(n+1)]
  h  <- X[-1] - X[-(n+1)]
  
  # Cubic spline: Natural and clamped use H and u
  
  H1 <- diag(c(h[2:(n-1)], 0))
  H  <- diag(2*(h[-n] + h[-1])) + cbind(H1[,(n-1)], H1[,-(n-1)]) + rbind(H1[(n-1),], H1[-(n-1),])
  u  <- matrix(-6*(dy[-n]/h[-n] - dy[-1]/h[-1]), n-1, 1)
  
  if(nc == 0) {     # Clamped flat (zero 1st deriv at boundary)
    H[1,1]     <- (3/2)*h[1] + 2*h[2]
    H[n-1,n-1] <- (3/2)*h[n] + 2*h[n-1]
    u[1]       <- u[1]     - 3*dy[1]/h[1]
    u[n-1]     <- u[n-1]   + 3*dy[n-1]/h[n-1]
    mc         <- solve(H,u)
    m          <- rbind(3*dy[1]/(h[1]^2)-mc[1]/2,
                        mc,
                       -3*dy[n-1]/(h[n-1]^2)-mc[n-1]/2)
    } else {        # Natural spline (zero 2nd deriv at boundary)
    m          <- rbind(0, solve(H,u), 0)
  }

  # To make calculation easier shift m by 1.
  me <- m[-1]
  m1 <- m[-length(m)]

  # Calculate coefficients
  return(tibble("a" = Y[1:n],
                "b" = dy/h - h*(2*m1 + me)/6,
                "c" = m1/2,
                "d" = (me-m1)/(6*h)))
}

########################

# Clamped
C0 <- cubic_spline(X, Y, 0L)
a  <- C0$a
b  <- C0$b
c  <- C0$c
d  <- C0$d

interp0 <- tibble(x=xaxis, s=1L)
for (i in 1:(length(X)-1)) {
  interp0 <- mutate(interp0 , s = if_else((x >= X[i] & x < X[i+1]), i, s))
  }

interp0 <- interp0 %>%
  mutate(f = a[s] + b[s]*(x-X[s]) + c[s]*(x-X[s])^2 + d[s]*(x-X[s])^3) %>%
  mutate(f = if_else(x < X[1], Y[1], 
                     if_else(x >= tail(X,1), tail(Y,1), f))) %>%
  fill(f)

# Natural
C1 <- cubic_spline(X, Y, 1L)
a  <- C1$a
b  <- C1$b
c  <- C1$c
d  <- C1$d

interp1 <- tibble(x=xaxis, s=1L)
for (i in 1:(length(X)-1)) {
  interp1 <- mutate(interp1 , s = if_else((x >= X[i] & x < X[i+1]), i, s))
  }

interp1 <- interp1 %>%
  mutate(s = if_else(x >= tail(X,2)[1], (length(X)-1L), s)) %>%
  mutate(f = a[s] + b[s]*(x-X[s]) + c[s]*(x-X[s])^2 + d[s]*(x-X[s])^3) %>%
  mutate(f = if_else(x < X[1],  a[s] + b[s]*(x-X[s]), f))  %>%
  mutate(f = if_else(x >= X[length(X)], f, f))

#############################################################################

ggplot(D) +
  geom_point(aes(x=X, y=log(Y)), shape=23, fill="green", size=3) +
  geom_point(data=interp1, aes(x=x, y=log(sf(x))), shape=21, fill="pink", size=2.2) +
  geom_line(data=interp0,  aes(x=x, y=log(f)), color="red", linetype=1, linewidth=0.9) +
  geom_line(data=interp1,  aes(x=x, y=log(f)), color="blue", linetype=2, linewidth=0.9) +
  theme_light()
