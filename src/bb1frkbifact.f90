! Subroutines lbb1derivs cbb1derivs are in file bb1facts.f90
! Subroutines lfrk2derivs is in file frk12fact-lpdf.f90
! Subroutines lbb1pdf lfrkpdf cbb1cdf are in file cop12factlik.f90

! bi-factor model with BB1 linking copulas for common latent
!                 and Frank linking copula for group latent
! inputs
!   npar = 3*dvar = number of parameters
!   th = vector of parameters; 
!   mgrp = number of groups 
!   n = sample size
!   dvar = sum(grsize) = dimension of data
!   grsize = vector with group sizes
!   udata = nxd matrix of uniform scores 
!   nq = number of quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes 
! outputs 
!   nllk = negative log-likelihood, 
!   grad = gradient of nllk, 
!   hess = hessian of nllk
subroutine strbb1frk2(npar,th,mgrp,n,dvar,grsize,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,nq,ip,jp,ind
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision lpdf(2*dvar),der1(npar),der2(dvar+npar)
  double precision ccdf(dvar),cder1(2*dvar),cder2(npar),der1bb1(2),der2bb1(3),cder1bb1(2),cder2bb1(3)
  double precision nllk,liki,llk,lk,lk2,grad(npar),hess(npar,npar)
  integer i,iq,iq2,jg,mj,mj2,ind1,ind2,grsize(mgrp),ibb1(2)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision der1u(dvar), der2u(dvar),derumix(dvar), dermxj(npar)
  double precision hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

  ! npar = 3*dvar; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)    
    liki = 0.d0; grdi = 0.d0; hssi = 0.d0; 
    !i0 ---> iq
    do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; hssaux = 0.d0 
      ind = 0; intj = 0.d0;  grdj = 0.d0; der2j = 0.d0; dermxj = 0.d0 
      !j ---> jg
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        !ij ---> iq2
        do iq2 =1,nq   
          lk = 1.d0; lk2 = 1.d0 
          do mj = ind1,ind2
            !C_{ij|V_0}
            call cbb1derivs(uvec(mj),xl(iq),th((/mj,mj+dvar/)),ccdf(mj),cder1bb1,cder2bb1)                                         
            if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
            if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999  
            cder1((/mj,mj+dvar/)) = cder1bb1
            cder2((/mj,mj+2*dvar,mj+dvar/)) = cder2bb1
            !c_{i,V_j;V_0} 
            call lfrk2derivs(ccdf(mj),xl(iq2),th(2*dvar+mj),lpdf(dvar+mj),der1(2*dvar+mj),der2(3*dvar+mj),&
                 der1u(mj),der2u(mj),derumix(mj))   
            !c_{ij,V_0}
            call lbb1derivs(uvec(mj),xl(iq),th((/mj,mj+dvar/)),lpdf(mj),der1bb1,der2bb1)                                           
            der1((/mj,mj+dvar/)) = der1bb1
            der2((/mj,mj+2*dvar,mj+dvar/)) = der2bb1  
            lk = lk*exp(lpdf(mj))
            lk2 = lk2*exp(lpdf(dvar+mj))
          end do
          llk = lk*lk2
          intj(jg) = intj(jg) + wl(iq2)*llk !intj value of the jth integral
          intjj(ind1:ind2) = intj(jg)  !intjj: repeat intj(j) grsize[j] times
          intjj((dvar+ind1):(dvar+ind2)) = intj(jg)
          intjj((2*dvar+ind1):(2*dvar+ind2)) = intj(jg)
          do mj = ind1,ind2
            ibb1 = (/mj,mj+dvar/) 
            !1st o.d. wrt th.bb1,dl.bb1
            grdj(ibb1) = grdj(ibb1) + wl(iq2)*llk*(der1u(mj)*cder1(ibb1)+der1(ibb1))                                        
            !1st o.d. wrt th.frk
            grdj(2*dvar+mj) = grdj(2*dvar+mj) + wl(iq2)*llk*der1(2*dvar+mj)
            !2nd o.d. wrt th.bb1,dl.bb1
            der2j(ibb1) = der2j(ibb1)  + wl(iq2)*llk*(der2u(mj)+der1u(mj)*der1u(mj))*cder1(ibb1)*cder1(ibb1)                  
            der2j(ibb1) = der2j(ibb1)  + wl(iq2)*llk*(der1u(mj)*cder2(ibb1)+der2(ibb1)+der1(ibb1)*der1(ibb1))       
            der2j(ibb1) = der2j(ibb1)  + wl(iq2)*llk*2*der1(ibb1)*der1u(mj)*cder1(ibb1)                            
            !2nd o.d. wrt th.frk
            der2j(2*dvar+mj) = der2j(2*dvar+mj)+ wl(iq2)*llk*(der2(3*dvar+mj)+der1(2*dvar+mj)*der1(2*dvar+mj)) 
            !mixed d. wrt th.bb1 and dl.bb1
            dermxj(ibb1) = dermxj(ibb1) + wl(iq2)*llk &
                   *(der1(ibb1)*der1(2*dvar+mj)+(derumix(mj)+der1u(mj)*der1(2*dvar+mj))*cder1(ibb1))  
            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + wl(iq2)*llk &
                   *(der2(mj+2*dvar)+der1(mj)*der1(mj+dvar)+der1u(mj)*der1(dvar+mj)*cder1(mj))
            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + wl(iq2)*llk &
                   *(der1u(mj)*der1(mj)*cder1(mj+dvar)+der1u(mj)*cder2(2*dvar+mj))
            dermxj(2*dvar+mj) = dermxj(2*dvar+mj) + wl(iq2)*llk &
                   *((der2u(mj)+der1u(mj)*der1u(mj))*cder1(mj)*cder1(dvar+mj))
            do mj2 = ind1,(mj-1)  !2nd order der. of the jth int. wrt th(mj), th(mj2)
              hssaux(mj,mj2)  = hssaux(mj,mj2) + wl(iq2)*llk*(der1u(mj)*cder1(mj)+der1(mj))*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj+dvar,mj2) = hssaux(mj+dvar,mj2) + wl(iq2)*llk &
                    *(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj,mj2+dvar) = hssaux(mj,mj2+dvar) + wl(iq2)*llk &
                    *(der1u(mj)*cder1(mj)+der1(mj))*(der1u(mj2)*cder1(mj2+dvar) &
                    +der1(mj2+dvar))
              hssaux(mj+dvar,mj2+dvar) = hssaux(mj+dvar,mj2+dvar) + wl(iq2)*llk &
                    *(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))*(der1u(mj2)*cder1(mj2+dvar)+der1(mj2+dvar))
              hssaux(mj+2*dvar,mj2) = hssaux(mj+2*dvar,mj2) + wl(iq2)*llk*der1(2*dvar+mj)*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(mj2+2*dvar,mj) = hssaux(mj2+2*dvar,mj) + wl(iq2)*llk*der1(2*dvar+mj2)*(der1u(mj)*cder1(mj)+der1(mj))
              hssaux(mj+2*dvar,mj2+dvar) = hssaux(mj+2*dvar,mj2+dvar) + wl(iq2)*llk &
                    *der1(2*dvar+mj)*(der1u(mj2)*cder1(mj2+dvar)+der1(mj2+dvar))
              hssaux(mj2+2*dvar,mj+dvar) = hssaux(mj2+2*dvar,mj+dvar) + wl(iq2)*llk &
                    *der1(2*dvar+mj2)*(der1u(mj)*cder1(mj+dvar)+der1(mj+dvar))
              hssaux(mj+2*dvar,mj2+2*dvar) = hssaux(mj+2*dvar,mj2+2*dvar) + wl(iq2)*llk*der1(2*dvar+mj)*der1(2*dvar+mj2)
              hssaux(mj2,mj)  = hssaux(mj,mj2)
              hssaux(mj2,mj+dvar) = hssaux(mj+dvar,mj2)
              hssaux(mj2+dvar,mj) = hssaux(mj,mj2+dvar) 
              hssaux(mj2+dvar,mj+dvar) = hssaux(mj+dvar,mj2+dvar)
              hssaux(mj2,mj+2*dvar) = hssaux(mj+2*dvar,mj2)
              hssaux(mj,mj2+2*dvar) = hssaux(mj2+2*dvar,mj) 
              hssaux(mj2+dvar,mj+2*dvar) = hssaux(mj+2*dvar,mj2+dvar) 
              hssaux(mj+dvar,mj2+2*dvar) = hssaux(mj2+2*dvar,mj+dvar)  
              hssaux(mj2+2*dvar,mj+2*dvar) = hssaux(mj+2*dvar,mj2+2*dvar)  
            end do                    
          end do
        end do          
      end do

      int0 = product(intj) !product of j inner integrals  
      grd0 = int0*(grdj/intjj)!gradient of the product
      !m1 ---> ip; m2 ---> jp
      do ip=2,npar !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      ind = 0
      !j ---> jg
      do jg=1,mgrp    !2nd order derivatives within groups 
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
        hss0(ind1:ind2,dvar+ind1:dvar+ind2)=int0*hssaux(ind1:ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,ind1:ind2)=int0*hssaux(dvar+ind1:dvar+ind2,ind1:ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,ind1:ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(ind1:ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
        hss0(2*dvar+ind1:2*dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(2*dvar+ind1:2*dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,2*dvar+ind1:2*dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,2*dvar+ind1:2*dvar+ind2)/intj(jg)
      end do

      do ip=1,npar
        hss0(ip,ip)  = int0*der2j(ip)/intjj(ip)
      end do
      do ip=1,dvar
        hss0(ip,ip+dvar)  = int0*dermxj(ip+2*dvar)/intjj(ip)
        hss0(ip+2*dvar,ip)  = int0*dermxj(ip)/intjj(ip)
        hss0(ip+2*dvar,ip+dvar) = int0*dermxj(ip+dvar)/intjj(ip)
        hss0(ip+dvar,ip)  = hss0(ip,ip+dvar)
        hss0(ip,ip+2*dvar)  = hss0(ip+2*dvar,ip)
        hss0(ip+dvar,ip+2*dvar) = hss0(ip+2*dvar,ip+dvar)
      end do
      liki = liki + wl(iq)*int0   !integrating over V0
      grdi = grdi + wl(iq)*grd0
      hssi = hssi + wl(iq)*hss0
    end do
    nllk = nllk - log(liki)       !updating loglikelihood
    grad = grad - grdi/liki       !updating gradient
    do ip=1,npar                  !updating hessian
      do jp=1,ip
        hessi(ip,jp)=hssi(ip,jp)/liki-grdi(ip)*grdi(jp)/liki/liki 
        hessi(jp,ip)=hessi(ip,jp)
      end do
    end do 
    hess = hess - hessi 
  end do
  return
  end

! bi-factor model with BB1 linking copulas for common latent
!                 and Frank linking copula for group latent
! This function computes the likelihood only, without derivatives
! inputs
!   npar = 3*dvar = number of parameters
!   th = vector of parameters; 
!   mgrp = number of groups 
!   n = sample size
!   dvar = sum(grsize) = dimension of data
!   grsize = vector with group sizes
!   udata = nxd matrix of uniform scores 
!   nq = number of quadrature points
!   wl = vector of quadrature weights
!   xl = vector of quadrature nodes 
! output 
!   nllk = negative log-likelihood
subroutine strbb1frk2nllk(npar,th,mgrp,n,dvar,grsize,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,mgrp,dvar,n,nq,ind
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision lpdf(2*dvar), ccdf(dvar)
  double precision nllk,liki,llk,lk,lk2
  integer i,iq,iq2,jg,mj,ind1,ind2,grsize(mgrp)
  double precision intj(mgrp), int0

  ! npar = 3*dvar; dvar = sum(grsize)
  nllk=0.d0 
  do i =1,n 
    uvec = udata(i,:) 
    liki = 0.d0
    do iq =1,nq
      int0 = 1.d0; intj = 0.d0; ind = 0;
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        do iq2 =1,nq   
          lk = 1.d0; lk2 = 1.d0        
          do mj = ind1,ind2
            !C_{ij|V_0}
            call cbb1cdf(uvec(mj),xl(iq),th((/mj,mj+dvar/)),ccdf(mj))                                         
            if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
            if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999               
            !c_{i,V_j;V_0} 
            call lfrkpdf(ccdf(mj),xl(iq2),th(2*dvar+mj),lpdf(dvar+mj))   
            !c_{ij,V_0}
            call lbb1pdf(uvec(mj),xl(iq),th((/mj,mj+dvar/)),lpdf(mj))  
            lk = lk*exp(lpdf(mj)) 
            lk2 = lk2*exp(lpdf(dvar+mj)) 
          end do
          llk = lk*lk2
          intj(jg) = intj(jg) + wl(iq2)*llk !intj value of the jth integral
        end do 
      end do
      int0 = product(intj)  !product of j inner integrals  
      liki = liki + wl(iq)*int0  !integrating over V0
    end do
    nllk = nllk - log(liki) !updating loglikelihood
  end do
  return
  end

