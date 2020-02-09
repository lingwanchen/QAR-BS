#########  required functions for the proposed method  ##########

## Survival function for univariate censoring setting
Surv.uni<-function(Y, delta) {
  pts<-sort(unique(Y))
  npts<-length(pts)
  n<-length(Y)
  obs.mat<-matrix(Y,npts,n,byrow=T)
  delta.mat<-matrix(delta,npts,n,byrow=T)
  Dt.mat<-ifelse(obs.mat==pts & delta.mat == 1 , 1, 0)  # npts x n
  Yt.mat<-ifelse(obs.mat >= pts,1,0)  # npts x n
  lambda<-rowSums(Dt.mat)/rowSums(Yt.mat) # npts-vector
  Sxx<-cumprod((1-lambda))  # Survival for sort.time
  Sx<-pmax(Sxx[colSums(Yt.mat)],1e-20)
  return(list(St=sort(Sx[delta==1],decreasing=T), Yt=sort(Y[delta==1])))
}

est.surv.fun<-function(Y, pts, St) {
  npts<-length(pts)
  n<-length(Y)
  obs.mat<-matrix(Y,npts,n,byrow=T)
  Yt.mat<-ifelse(obs.mat >= pts,1,0)  
  prob<-c(1,St)
  Sx<-(prob[colSums(Yt.mat)+1])
  return(Sx)
}

## Function for estimating beta
PF.method<-function(Y, tr.Y, delta, Zmat, tauseq, M=10^7, Surv.C){
  Gx<-est.surv.fun(Y, Surv.C$Yt, Surv.C$St) # est of survival of censoring                 
  wi<-(delta==1)/Gx               # weight
  Yn<-c(tr.Y*wi,M,M)                 # new Yn for PF's method
  hat.beta.tau<-c()
  for( i in 1: length(tauseq))
  {
    tau<-tauseq[i]
    Xmat <- rbind(Zmat*wi, apply(-Zmat*wi,2, sum), apply(2*Zmat*tau, 2, sum))
    fit<-rq(Yn~Xmat-1)$coef   # est.beta by PF's method
    hat.beta.tau<-cbind(hat.beta.tau, fit)  # each row contains fit from all tau. 
  }
  return(list(hat.b=hat.beta.tau, Gx=Gx))
}

## Function: qor to copula
X.func<-function(Zmat,gamma,u,v){
  y<-c(exp(Zmat%*%gamma))
  copula<- (u+v)/2 + (1-sqrt( (y-1)^2*(u-v)^2+2*(y-1)*(u+v-2*u*v)+1 ))/(2*(y-1))
  copula[y==1]<-u*v
  return(copula)
}

dX.func<-function(y, u, v){
  diff.X<- -(y-1)^(-2)/2 + ( (u-v)^2 + 2*(u+v-2*u*v)*(y-1)^(-1) + (y-1)^(-2) )^(-1/2)*((u+v-2*u*v)*(y-1)^(-2) + (y-1)^(-3))/2
  diff.X[y==1]<-0
  return(diff.X)
}

creat.list<-function(Y1.d1,Y2.d2, Zmat, tauseq, invw){ 
  list.obj<-c() 
  for(i in 1:length(tauseq)){
    list.obj[[i]]<-cbind(tauseq[i], Y1.d1[,i], Y2.d2[,i], Zmat, invw)
  }
  return(list.obj)
}

## Estimating Equation for gamma
est.func<-function(gamma, para.mat){
  p<-dim(para.mat)[2]-4
  Zmat<-para.mat[ , c(4:(3+p))] 
  u<-v<-para.mat[1, 1]  
  Gx12<-para.mat[, (4+p)]
  ind.12<-para.mat[, 2]*para.mat[, 3]*1
  est.copula<-X.func(Zmat,gamma,u,v)
  estimate.func<-colMeans(Zmat*(ind.12/Gx12 - est.copula))
  return(estimate.func)
}

d.est.func.<-function(gamma, para.mat){
  p<-dim(para.mat)[2]-4
  Zmat<-para.mat[ , c(4:(3+p))] 
  u<-v<-para.mat[1, 1]  
  Gx12<-para.mat[, (4+p)]
  ind.12<-para.mat[, 2]*para.mat[, 3]*1
  est.copula<-X.func(Zmat,gamma,u,v)
  y<-c(exp(Zmat%*%gamma))
  -t(Zmat)%*%(Zmat*y*dX.func(est.copula,u,v))/n
  
  estimate.func<-colMeans(Zmat*(ind.12/Gx12 - est.copula))
  return(estimate.func)
}

## Function for estimating gamma 
opt.gamma.func<-function(para.mat){
  p<-dim(para.mat)[2]-4
  output<-multiroot(est.func, start=rep(0.05,p), #start=rep(0.05,p), 
                    parms=para.mat, rtol = 1e-12, atol = 1e-12)
  result<-list(hat.gamma=output$root)
  return(result)
}

### try 
# creat.list2<-function(Y1.d1,Y2.d2, Zmat, tauseq, invw, de1 ,de2){ 
#   list.obj<-c() 
#   for(i in 1:length(tauseq)){
#     list.obj[[i]]<-cbind(tauseq[i], Y1.d1[,i], Y2.d2[,i], Zmat, invw, de1, de2)
#   }
#   return(list.obj)
# }
# ## Estimating Equation for gamma
# est.func2<-function(gamma, para.mat){
#   p<-dim(para.mat)[2]-6
#   Zmat<-para.mat[ , c(4:(3+p))] 
#   u<-v<-para.mat[1, 1]  
#   Gx12<-para.mat[, (4+p)]
#   ind.12<-para.mat[, 2]*para.mat[, 3]*1
#   de.12<-para.mat[, 8]*para.mat[, 9]*1
#   est.copula<-X.func(Zmat,gamma,u,v)
#   estimate.func<-colMeans(Zmat*(ind.12 - est.copula)*de.12/Gx12)
#   return(estimate.func)
# }
# 
# ## Function for estimating gamma 
# opt.gamma.func2<-function(para.mat){
#   output<-multiroot(est.func2, start=c(1,0,2), 
#                     parms=para.mat, rtol = 1e-12, atol = 1e-12)
#   result<-list(hat.gamma=output$root)
#   return(result)
# }


###############################
## Influence functions for beta
# sub-functions
Gij.mat<-function(Y, delta, case, tseq){
  n<-length(Y)
  pts<-sort(unique(Y))
  npts<-length(pts)
  
  obs.mat<-matrix(Y,npts,n,byrow=T)
  delta.mat<-matrix(delta,npts,n,byrow=T)
  
  dN.mat<-ifelse(obs.mat==pts & delta.mat == case, 1, 0)  # npts x n
  yt.mat<-ifelse(obs.mat >= pts, 1, 0)  # npts x n
  yt.vec<-rowMeans(yt.mat)  # npts-vector
  lambda<-rowSums(dN.mat)/rowSums(yt.mat) # npts-vector
  Sxx<-cumprod((1-lambda))  # Survival for sort.unique.time
  dM.mat<-dN.mat-yt.mat*lambda
  gij<-dM.mat/yt.vec  # npts x n 
  
  # Calculate the Gij for special time.seq #
  ntseq<-length(tseq)
  tseq.mat<-matrix(tseq,npts,ntseq,byrow=T)
  ytseq.pts.mat<-ifelse(tseq.mat >= pts,1,0)
  posi<-colSums(ytseq.pts.mat)
  Gx.ytseq<-pmax(Sxx[posi],1e-20)  # ntseq-vector Survival function
  G.vec<-c()
  for( i in 1:ntseq)
  {
    G.vec<-cbind(G.vec, colSums(gij*ytseq.pts.mat[,i]*Gx.ytseq[i]))
    #    G.vec<-cbind(G.vec, colSums(matrix(gij[1:n.yt.pts[i],],,n.yt.pts[i]))) 
  }
  return(G.vec) # row: unique.time; col: subject
}

Est.Sigma.gi<-function (fit, Y, tr.Y, delta, Gx, Zmat, case, tau){ 
  ri<-Zmat%*%fit-tr.Y
  w<-Zmat*c(ri>=0)*c(delta==1)/Gx^2 # n x p matrix
  ww<-(ri>=0)*(delta==1)/Gx
  eta1<-Zmat*c(ww-tau)
  
  gij<-Gij.mat(Y, delta, case, Y)  # Case=0 : censoring state
  eta2<-c()
  for ( p in 1:dim(Zmat)[2])
  {
    eta2<-cbind(eta2, colMeans(gij*w[,p]))
  }
  eta<-eta1-eta2
  return(eta)
}

## Generalized Estimating Equation for beta
#Sn<-function(Y, Zmat, fit, Gx, tau){
#  ri<-Zmat%*%fit-Y
#  w<-(ri>=0)/Gx
#  Sn.b<-colMeans(Zmat*c(w-tau))*sqrt(n)
#  return(Sn.b)
#}

## induced smoothing est for GEE of beta (score func of beta)
Un.b<-function(fit, tr.Y, Zmat, H, wi, tau){
  zwz<-diag(Zmat%*%H%*%t(Zmat))
  ri<-Zmat%*%fit-tr.Y
  Un<-colMeans(Zmat*c(pnorm(ri/sqrt(zwz))*wi-tau))
  return(Un)
}

## induced smoothing est for derivative of score func of beta
An.b<-function(fit, tr.Y, Zmat, H, wi, ZZ, pr, pc){
  zwz<-diag(Zmat%*%H%*%t(Zmat))
  ri<-Zmat%*%fit-tr.Y
  An.mat0<-c(wi*dnorm(ri/sqrt(zwz))/sqrt(zwz))*ZZ
  An<-matrix(colMeans(An.mat0), pr, pc)
  return(list(An=An, An.mat=An.mat0))
}

#### Influence functions for beta ######
induced.std.beta.gij<-function(fit, Y, tr.Y, delta, Gx, Zmat, ZZ, tau){
  wi<-(delta==1)/Gx 
  n<-length(Y)
  p<-dim(Zmat)[2]
  H<-diag(1/n, p)
  Un<-Un.b(fit, tr.Y, Zmat, H, wi, tau)
  A0<-An.b(fit, tr.Y, Zmat, H, wi, ZZ, p, p)$An
  
  An.v<-c(A0)    
  An.n<-ifelse(An.v==0,1e-10,An.v)
  A0<-matrix(An.n, p, p)
  
  inf.0<-Est.Sigma.gi(fit, Y, tr.Y, delta, Gx, Zmat, 0, tau) # never change
  Sigma.0<-var(inf.0)
  b0<-fit
  
  Conv<-FALSE
  iter<-1
  
  while(Conv!=TRUE)
  {
    b.n<-b0-solve(A0)%*%matrix(Un) # update beta
    #b.n<-multiroot(Un.b, b0, tr.Y=tr.Y, Zmat=Z, H=H, wi=wi, tau=tau)$root 
    if ( max(abs(b.n)) >= 30 ){ Conv=TRUE }
    else{
      An<-An.b(b.n, tr.Y, Zmat, H, wi, ZZ, p, p)$An # update An
      An.v<-c(An)    
      An.n<-ifelse(An.v==0,1e-10,An.v)
      An<-matrix(An.n, p, p)
      H<-solve(An)%*%Sigma.0%*%solve(An)/n # update Hn
      Conv=((max(abs(b0-b.n))<=1e-10)*(sum(abs(An/A0-1))<=1e-15)==1)
      if ( iter >= 250 ){ Conv=TRUE }
      iter<-iter+1
      b0<-b.n
      A0<-An
      Un<-Un.b(b0, tr.Y, Zmat, H, wi, tau) #updata Un
    }
  }  
  An.arr<-An.b(b.n, tr.Y, Zmat, H, wi, ZZ, p, p)
  An<-An.arr$An
  An.mat<-An.arr$An.mat
  
  A0.arr<-An.b(fit, tr.Y, Zmat, H, wi, ZZ, p, p)
  A0<-A0.arr$An
  A0.mat<-A0.arr$An.mat
  
  return(list(An=An, An.mat=An.mat, H=H, beta=b.n, A0=A0, A0.mat=A0.mat, infy=inf.0)) # H is covariance of beta
}

##############################

## induced smoothing est for partial derivative of score func of gamma
Pn.b<-function(fit.j, tr.Yj, Zj, Hj, deltaj, ZZj3, invw, in.deltab.jc, pr, pc){
  An.mat.j<-An.b(fit.j, tr.Yj, Zj, Hj, wi=deltaj, ZZj3, pr, pc)$An.mat
  Pn.j<-matrix(colMeans(An.mat.j*in.deltab.jc/invw), pr, pc)
  return(Pn.j)
}

Est.gij<-function(w1.psi, Zmat, Y.c, delta.c, invw, case){ 
  w<-Zmat*w1.psi/invw # n x p matrix
  gij<-Gij.mat(Y.c, delta.c, case, Y.c)  # Case=1 for censoring status in delta.c
  Eta<-c()
  for ( p in 1:dim(Zmat)[2])
  {
    Eta<-cbind(Eta, colMeans(gij*w[,p]))
  }
  return(Eta)
}

## est of derivative of score func of gamma
# Jn<-function(Zmat, ZZ, gamma, u, v){
#   p<-dim(Zmat)[2]
#   y<-c(exp(Zmat%*%gamma))
#   diff.x<-dX.func(y, u, v)
#   Jn.r<-matrix(colMeans(ZZ*diff.x*y), p, p)
#   return(Jn.r)
# }

## approximate normal term 
Psi.n<-function(Zmat, w1.psi, gamma, u, v, P1.b, P2.b, A1.b, A2.b, b1.inf, b2.inf){
  w2.psi<-X.func(Zmat,gamma,u,v)
  psi1<-Zmat*(w1.psi-w2.psi)
  psi2<-t(P1.b%*%solve(A1.b)%*%t(b1.inf))
  psi3<-t(P2.b%*%solve(A2.b)%*%t(b2.inf))
  Psi<-psi1-psi2-psi3
  return(Psi)
}

Psi.n2<-function(Zmat, w1.psi, wc.psi, gamma, u, v, P1.b, P2.b, A1.b, A2.b, b1.inf, b2.inf){
  w2.psi<-X.func(Zmat,gamma,u,v)
  psi1<-Zmat*(w1.psi-w2.psi)
  psi2<-t(P1.b%*%solve(A1.b)%*%t(b1.inf))
  psi3<-t(P2.b%*%solve(A2.b)%*%t(b2.inf))
  Psi<-psi1-psi2-psi3-wc.psi
  return(Psi)
}
