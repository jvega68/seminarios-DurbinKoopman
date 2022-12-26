nilo <- scan("data/nile.dat", skip =1)

llm <- function(y, a,P,sigma2_eps,sigma2_eta){
# Estimación del modelo lineal local vía recursiones 
# del filtro de Kalman.
  TT <- length(y)
  att  <- Ptt <- K <- F_ <- v <- alfahat <- r <- L <- V <- N <- NULL
  for(i in 1:TT){
    F_[i] <- P[i] + sigma2_eps
    K[i] <- P[i]/F_[i]
    v[i] <- y[i] - a[i]
    att[i] <- a[i] + K[i]*v[i]
    a[i+1] <- a[i] + K[i]*v[i]
    Ptt[i] <- P[i]*(1 - K[i])
    P[i+1] <- P[i]*(1 - K[i]) + sigma2_eta
  }
  
  
# Parte para el suavizamiento del estado
  
r[TT] <- N[TT] <- 0
for(i in TT:1){
  L[i] <- 1- K[i]
  if(i>=2) r[i-1] <- v[i]/F_[i] + L[i]*r[i]
  alfahat[i] <- a[i] + P[i]*ifelse(i==1,NA,r[i-1])
  V[i] <- P[i] -P[i]^2*ifelse(i==1,NA,N[i-1])
  if(i>=2) N[i-1] = 1/F_[i] + L[i]^2*N[i]
}

return(data.frame(  y = c(y,NA),
                    att = c(att,NA),
                      a = a,
                    Ptt = c(Ptt,NA),
                      P = P,
                      K = c(K, NA), 
                      alfahat = c(NA,alfahat),
                      L = c(NA, L),
                      N = c(NA, N),
                      V = c(NA, V),
                      r = c(NA, r)
                    ))

  
}



a <- llm(nilo, a = 0, P = 10e7, sigma2_eps = 15099, sigma2_eta = 1469.1)

par(mfrow = c(1,2))
plot(nilo)
a$a[1] <- NA
lines(a$a, lwd = 2)
lines(a$a - 1.65*sqrt(a$P), col = "red", lty=3, lwd = 2)
lines(a$a + 1.65*sqrt(a$P), col = "red", lty=3, lwd = 2)
lines(a$alfahat, lwd = 2, col = "blue")

plot(a$y - a$a, type = "l")
abline(h=0)

y <- datasets::Nile
struct <- StructTS(y, type = "level")
if (struct$code != 0) stop("optimizer did not converge")
print(struct$coef)


lines(struct$fitted, col  = "red", lwd = 5)
