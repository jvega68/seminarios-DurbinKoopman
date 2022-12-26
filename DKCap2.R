llm <- function(y,a,P,sigma2_eps,sigma2_eta){
  # Estimación del modelo lineal local vía recursiones 
  # del filtro de Kalman.
  # Tambièn estima el suavizamiento del estado. 
  
  TT <- length(y)
  att  <- Ptt <- K <- F_ <- v <- alfahat <- r <- L <- V <- N <- epshat <- etahat <- Varepshat <- Varetahat <- u <- D <- Vs<-NULL
  
  # ìndices de casos faltantes
  m <- which(is.na(y))
  
  # Recursiones:
  for(i in 1:TT){
    F_[i] <- P[i] + sigma2_eps
    K[i] <- ifelse(i %in% m, 0,P[i]/F_[i])
    v[i] <- y[i] - a[i]
    att[i] <- a[i] + ifelse(i %in% m, 0, K[i]*v[i])
    Ptt[i] <- P[i]*(1 - K[i])
    a[i+1] <- a[i] + ifelse(i %in% m, 0, K[i]*v[i])
    P[i+1] <- P[i]*(1 - K[i]) + sigma2_eta  
  }
  
  # Parte para el suavizamiento del estado. Esta parte es recursiva inversa
  
  r[TT] <- N[TT] <- 0 # valores iniciales
  for(i in TT:1){
    L[i] <- 1- K[i]
    if(i>=2) r[i-1] <- v[i]/F_[i] + L[i]*ifelse(is.na(r[i]),0,r[i])
    if(i>=2) N[i-1] = 1/F_[i] + L[i]^2*N[i]
    alfahat[i] <- a[i] + P[i]*ifelse(i == 1, NA, r[i-1])
    V[i] <- abs(P[i] -P[i]^2*ifelse(i == 1, NA, N[i-1]))
  }
  
  
  # Suavizamiento de las perturbaciones
  
  for(i in TT:1){
    u[i] <- v[i]/F_[i] - K[i]*r[i]
    epshat[i] <- sigma2_eps*u[i]
    D[i] <- 1/F_[i] + K[i]^2*N[i]
    Varepshat[i] <- sigma2_eps - sigma2_eps^2*D[i]
    etahat[i] <- sigma2_eta*r[i]
    Varetahat[i] <- sigma2_eta - sigma2_eta^2*N[i]
  }
  
  return(data.frame(  y = c(y,NA),
                      att = c(att,NA),
                      a = a,
                      v = c(v,NA),
                      Ptt = c(Ptt,NA),
                      P = P,
                      K = c(K, NA), 
                      alfahat = c(NA,alfahat),
                      L = c(NA, L),
                      N = c(NA, N),
                      V = c(NA, V),
                      r = c(NA, r),
                      F_ = c(NA,F_),
                      u = c(NA,u),
                      epshat = c(NA,epshat),
                      Varepshat = c(NA,Varepshat),
                      etahat = c(NA,etahat),
                      Varetahat = c(NA,Varetahat)
  ))
  
  
}

y <- datasets::Nile
a <- llm(y, a = 0, P = 10e7, sigma2_eps = 15099, sigma2_eta = 1469.1)

RecKalman <-function(datos){
  modelo <- StructTS(datos, type="level")
  n <- length(modelo$data)
  Ft <- NULL
  Pt <- NULL
  at <- NULL
  Kt <- NULL
  v2t <-NULL
  vt <- modelo$data-modelo$fitted
  Pt[1] <- modelo$model0$P
  at[1] <- modelo$data[1]
  s2eps <- modelo$coef[2]
  s2eta <- modelo$coef[1]
  for(i in 1:n){
    v2t[i] <- modelo$data[i] - at[i]
    Ft[i] <- Pt[i]+s2eps
    Kt[i] <- Pt[i]/Ft[i]
    at[i+1] <- at[i]+ Kt[i]*v2t[i]
    Pt[i+1] <- Pt[i]*(1-Kt[i]) + s2eta
  }
  #Recursiones para suavizamiento:
  rt <- NULL
  Nt <- NULL
  Dt <- NULL
  Lt <- 1-Kt
  rt[n] <- 0
  Nt[n] <- 0
  alfahatt <- NULL
  var.alfat <- NULL
  for(i in n:2){
    rt[i-1] <- v2t[i]/Ft[i]+Lt[i]*rt[i]
    Nt[i-1] <- 1/Ft[i]+Lt[i]^2*Nt[i]
    alfahatt[i] <- at[i]+Pt[i]*rt[i-1]
    var.alfat[i] <- Pt[i]-Pt[i]^2*Nt[i-1]
  }
  r0 <- vt[1]/Ft[1]+Lt[1]*rt[1]
  N0 <- 1/Ft[1]+Lt[1]^2*Nt[1]
  alfahatt[1] <- at[1]+Pt[1]*r0
  var.alfat[1] <- Pt[1]-Pt[1]^2*N0
  return(list(rec = cbind(v2t=v2t,
                          Ft=Ft,
                          Kt=Kt,
                          at=at[-1],
                          Pt=Pt[1:n],
                          fitted=modelo$fitted,
                          vt=vt,
                          Lt=Lt,
                          Nt=Nt,
                          rt=rt,
                          alfahat=alfahatt,var.alfat=var.alfat),
              errores = cbind(Dt=1/Ft+Kt^2*Nt,epshat=s2eps*vt/Ft-Kt*rt,
                              vareps=s2eps-s2eps^2*(1/Ft+Kt^2*Nt),etahat=s2eta*rt,
                              vareta=s2eta-s2eta^2*Nt),
              init = c(r0=r0,N0=N0,sigma2eps= s2eps,sigma2eta=s2eta)))
  
}