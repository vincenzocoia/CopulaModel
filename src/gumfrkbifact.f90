! Subroutines lgum1derivs cgumderivs are in file gum12fact.f90
! Subroutine lfrk2derivs is in file frk12fact-lpdf.f90
! Subroutines lgumpdf lfrkpdf cgumcdf are in file cop12factlik.f90

! bi-factor model with Gumbel linking copulas for common latent
!                 and Frank linking copula for group latent
! inputs
!   npar = 2*dvar = number of parameters
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
subroutine strgumfrk2(npar,th,mgrp,n,dvar,grsize,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,nq,ip,jp,ind,jj,jj2
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision lpdf(npar),der1(npar),der2(npar)
  double precision ccdf(dvar),cder1(dvar),cder2(dvar)
  double precision nllk,liki,llk,lk,lk2,grad(npar),hess(npar,npar)
  integer i,iq,iq2,jg,mj,mj2,ind1,ind2,grsize(mgrp)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision der1u(dvar), der2u(dvar),derumix(dvar), dermxj(dvar)
  double precision hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

  ! npar = 2*dvar; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)  
    liki = 0.d0; grdi = 0.d0; hssi = 0.d0; 
    do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; hssaux = 0.d0 
      ind = 0; intj = 0.d0;  grdj = 0.d0; der2j = 0.d0; dermxj = 0.d0 
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        do iq2 =1,nq   
          lk = 1.d0; lk2 = 1.d0  
          do mj = ind1,ind2  ! within group index
            call cgumderivs(uvec(mj),xl(iq),th(mj),ccdf(mj),cder1(mj),cder2(mj))  !C_{ij|V_0} 
            if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
            if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999  
            call lfrk2derivs(ccdf(mj),xl(iq2),th(dvar+mj),lpdf(dvar+mj),  &
             der1(dvar+mj),der2(dvar+mj),der1u(mj),der2u(mj),derumix(mj)) !c_{i,V_j;V_0} 
            call lgum1derivs(uvec(mj),xl(iq),th(mj),lpdf(mj),der1(mj),der2(mj)) !c_{ij,V_0}  
            lk = lk*exp(lpdf(mj)) 
            lk2 = lk2*exp(lpdf(dvar+mj)) 
          end do
          llk = lk*lk2
          intj(jg) = intj(jg) + wl(iq2)*llk  !intj value of the jth integral
          intjj(ind1:ind2) = intj(jg) !intjj: repeat intj(j) grsize[j] times
          intjj((dvar+ind1):(dvar+ind2)) = intj(jg)
          do mj = ind1,ind2
            jj=dvar+mj  ! index for parameter for latent for subgroup
            grdj(mj) = grdj(mj) + wl(iq2)*llk*(der1u(mj)*cder1(mj)+der1(mj))
            grdj(jj) = grdj(jj) + wl(iq2)*llk*der1(jj) !grdj: der. of the jth int. wrt th(mj)
            der2j(mj) = der2j(mj) + wl(iq2)*llk*(der2u(mj)+der1u(mj)*der1u(mj))*cder1(mj)*cder1(mj) !der2j: second order der. wrt th(mj)
            der2j(mj) = der2j(mj) + wl(iq2)*llk*(der1u(mj)*cder2(mj)+der2(mj)+der1(mj)*der1(mj)) !mistake1: added der1(mj)*der1(mj)
            der2j(mj) = der2j(mj) + wl(iq2)*llk*2*der1(mj)*der1u(mj)*cder1(mj) !mistake2: added one extra term 
            der2j(jj)= der2j(jj)+ wl(iq2)*llk*(der2(jj)+der1(jj)*der1(jj))
            dermxj(mj) = dermxj(mj) + wl(iq2)*llk*(der1(mj)*der1(jj)+(derumix(mj)+der1u(mj)*der1(jj))*cder1(mj))
            do mj2 = ind1,(mj-1)  !2nd order der. of the jth int. wrt th(mj), th(mj2)
              jj2=dvar+mj2  ! 2nd index for parameter for latent for subgroup
              hssaux(mj,mj2) = hssaux(mj,mj2) + wl(iq2)*llk*(der1u(mj)*cder1(mj)+der1(mj))*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(jj,mj2) = hssaux(jj,mj2) + wl(iq2)*llk*der1(jj)*(der1u(mj2)*cder1(mj2)+der1(mj2))
              hssaux(jj2,mj) = hssaux(jj2,mj) + wl(iq2)*llk*der1(dvar+mj2)*(der1u(mj)*cder1(mj)+der1(mj))
              hssaux(jj,jj2) = hssaux(jj,jj2) + wl(iq2)*llk*der1(jj)*der1(jj2)
              hssaux(mj2,mj) = hssaux(mj,mj2)
              hssaux(mj2,jj) = hssaux(jj,mj2)
              hssaux(mj,jj2) = hssaux(jj2,mj)
              hssaux(jj2,jj) = hssaux(jj,jj2)
            end do   
          end do  
        end do
      end do
      int0 = product(intj)     !product of j inner integrals 
      grd0 = int0*(grdj/intjj) !gradient of the product
      do ip=2,npar  !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      ind = 0
      do jg=1,mgrp  !2nd order derivatives within groups   
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,ind1:ind2)=int0*hssaux(dvar+ind1:dvar+ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,dvar+ind1:dvar+ind2)=int0*hssaux(ind1:ind2,dvar+ind1:dvar+ind2)/intj(jg)
      end do
      do ip=1,npar
        hss0(ip,ip) = int0*der2j(ip)/intjj(ip)
      end do
      do ip=1,dvar
        hss0(ip+dvar,ip) = int0*dermxj(ip)/intjj(ip)
        hss0(ip,ip+dvar) = hss0(ip+dvar,ip)
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

! bi-factor model with Gumbel linking copulas for common latent
!                 and Frank linking copula for group latent
! This function computes the likelihood only, without derivatives
! inputs
!   npar = 2*dvar = number of parameters
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
subroutine strgumfrk2nllk(npar,th,mgrp,n,dvar,grsize,udata,nq,wl,xl,nllk)
  implicit none
  integer npar,mgrp,dvar,n,nq,ind
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision lpdf(npar),ccdf(dvar)
  double precision nllk,liki,llk,lk,lk2
  integer i,iq,iq2,jg,mj,ind1,ind2,grsize(mgrp)
  double precision intj(mgrp), int0
  
  ! npar = 2*dvar; dvar = sum(grsize)
  nllk=0.d0; 
  do i =1,n 
    uvec = udata(i,:)   
    liki = 0.d0; 
    do iq =1,nq
      int0 = 1.d0; ind = 0; intj = 0.d0;  
      do jg =1,mgrp   !jth group                         
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        do iq2 =1,nq   
          lk = 1.d0; lk2 = 1.d0   
          do mj = ind1,ind2 ! within group index
            call cgumcdf(uvec(mj),xl(iq),th(mj),ccdf(mj)) !C_{ij|V_0}
            if(ccdf(mj) < 0.00001) ccdf(mj) = 0.00001
            if(ccdf(mj) > 0.99999) ccdf(mj) = 0.99999    
            call lfrkpdf(ccdf(mj),xl(iq2),th(dvar+mj),lpdf(dvar+mj)) !c_{i,V_j;V_0}  
            call lgumpdf(uvec(mj),xl(iq),th(mj),lpdf(mj)) !c_{ij,V_0}  
            lk = lk*exp(lpdf(mj))
            lk2 = lk2*exp(lpdf(dvar+mj))
          end do
          llk = lk*lk2
          intj(jg) = intj(jg) + wl(iq2)*llk !intj value of the jth integral
        end do   
      end do
      int0 = product(intj)   !product of j inner integrals  
      liki = liki + wl(iq)*int0   !integrating over V0
    end do
    nllk = nllk - log(liki)       !updating loglikelihood
  end do
  return
  end

