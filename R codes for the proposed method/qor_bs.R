## qor_bs is the main function for using the proposed method 

qor_bs<-function(Y1, Y2, delta1, delta2, Zmat, tauseq)
{
  ## Y1: time to the first event 
  ## Y2: time to the second event 
  ## delta1: observed indicator function for the first event
  ## delta2: observed indicator function for the second event
  ## Zmat: the design matrix of covariates, including the intercept
  ## tauseq: the quantile of interest
  
  n<-length(Y1)
  m<-length(tauseq)
  p<-dim(Zmat)[2]
  delta.c<-1-delta1*delta2
  time.c<-pmax(Y1,Y2)
  sum(delta.c)/n  # overall censoring rate= 0.315
  
  ZZ<-c()
  for(k in 1: p) { ZZ<-cbind(ZZ, Zmat[,k]*Zmat) } # nx(p^2) matrix 
  
  Surv.C<-Surv.uni(time.c, delta.c)
  invw<-est.surv.fun(time.c, Surv.C$Yt, Surv.C$St)
  wi12<-delta1*delta2/invw
  
  Surv.c1<-Surv.uni(Y1, 1-delta1)
  Surv.c2<-Surv.uni(Y2, 1-delta2)
  
  ## marginal coefficient estimates: \beta_1 \beta_2
  fit.Y1<-PF.method(Y1, log(Y1), delta1, Zmat, tauseq, M=10^7, Surv.c1)
  fit.Y2<-PF.method(Y2, log(Y2), delta2, Zmat, tauseq, M=10^7, Surv.c2)
  
  Gx1<-fit.Y1$Gx
  Gx2<-fit.Y2$Gx
  
  Q1.u<-Zmat%*%fit.Y1$hat.b
  Q2.v<-Zmat%*%fit.Y2$hat.b
  
  Ind.Y1.delta1<-((log(Y1)-Q1.u)<=0)*delta1
  Ind.Y2.delta2<-((log(Y2)-Q2.v)<=0)*delta2
  
  para.mat.all<-creat.list(Ind.Y1.delta1, Ind.Y2.delta2, Zmat, tauseq, invw)
  
  ## estimates for the association coefficient: \gamma
  (fit.r.seq<-matrix(unlist(lapply(para.mat.all, opt.gamma.func)),p))  
  
  sd.b1.mat<-diag(0,p,m)
  sd.b2.mat<-diag(0,p,m)
  sd.r.mat<-diag(0,p,m)
  for (j in 1: m)
  {
    tau.j<-tauseq[j]
    b1.est<-fit.Y1$hat.b[,j]
    b2.est<-fit.Y2$hat.b[,j]
    
    gamma.est<-fit.r.seq[, j]
    
    is.Y1<-induced.std.beta.gij(b1.est, Y1, log(Y1), delta1, Gx1, Zmat, ZZ, tau.j)
    is.Y2<-induced.std.beta.gij(b2.est, Y2, log(Y2), delta2, Gx2, Zmat, ZZ, tau.j)
    
    A01.b<-is.Y1$A0
    A02.b<-is.Y2$A0
    
    A01.mat.b<-is.Y1$A0.mat # for P1.b
    A02.mat.b<-is.Y2$A0.mat # for P2.b
    
    b1.inf<-is.Y1$infy
    b2.inf<-is.Y2$infy
    
    H1.b<-is.Y1$H        # for P1.b
    H2.b<-is.Y2$H        # for P2.b
    
    sd.b1.mat[,j]<-sqrt(diag(H1.b)) # sd for beta1
    sd.b2.mat[,j]<-sqrt(diag(H2.b)) # sd for beta2
    
    P1.b<-Pn.b(b1.est, log(Y1), Zmat, H1.b, delta1, ZZ, invw, Ind.Y2.delta2[,j], p, p)
    P2.b<-Pn.b(b2.est, log(Y2), Zmat, H2.b, delta2, ZZ, invw, Ind.Y1.delta1[,j], p, p)
    
    Jn.gamma<-jacobian(est.func, x=gamma.est,  para.mat=para.mat.all[[j]])
    
    w1.psi<-(Ind.Y1.delta1*Ind.Y2.delta2)[,j]/invw
    wc.psi<-Est.gij(w1.psi, Zmat, time.c, delta.c, invw, 1)
    
    Psi.b<-Psi.n2(Zmat, w1.psi, wc.psi, gamma.est, tau.j, tau.j, P1.b, P2.b, A01.b, A02.b, b1.inf, b2.inf)
    inf.mat<-t(Psi.b)%*%(Psi.b)/n
    cov.gamma<-solve(Jn.gamma)%*%inf.mat%*%solve(Jn.gamma)
    sd.r.mat[,j]<-sqrt(diag(cov.gamma)/n)
  }
  
  
  ### estimates of coefficients and standard deviations in marginal and association models.
  fit.Y1$hat.b
  fit.Y2$hat.b
  fit.r.seq
  colnames(fit.Y1$hat.b)<-colnames(fit.Y2$hat.b)<-colnames(fit.r.seq)<-paste0("tau_",tauseq)
  colnames(sd.b1.mat)<-colnames(sd.b2.mat)<-colnames(sd.r.mat)<-paste0("tau_",tauseq)
  
  
  return(list(est.beta1=fit.Y1$hat.b, est.beta2=fit.Y2$hat.b, est.gamma=  fit.r.seq,
              est.sd.beta1=sd.b1.mat, est.sd.beta2=sd.b2.mat, est.sd.gamma=sd.r.mat))
  
  
}
