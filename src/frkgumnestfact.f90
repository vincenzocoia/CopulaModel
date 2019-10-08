! Code written by Pavel Krupskii and revised by HJ
! The subroutine lgum1derivs is in file gum12fact.f90
! The subroutine lfrkderivs is in file frk12fact-lpdf.f90

! nested-factor model with Frank linking copulas for group to common latent
!         and Gumbel linking copulas from observed to group latent.
! inputs
!   npar = mgrp+dvar = number of parameters
!   th = vector of parameters; 
!        parameters for links of group latent to common latent go first
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
subroutine strfrkgum1(npar,th,mgrp,n,dvar,grsize,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,nq,ip,jp, ind
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision lpdf(npar),der1(npar),der2(npar),llk
  double precision nllk,liki,lk,grad(npar),hess(npar,npar)
  integer i,iq,iq2,jg,jg2,mj,mj2,nj,nj2,ind1,ind2,grsize(mgrp)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

  ! npar = dvar+mgrp; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)  
    liki = 0.d0; grdi = 0.d0; hssi = 0.d0; 
    do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; hssaux = 0.d0 
      ind = 0; intj = 0.d0;  grdj = 0.d0; der2j = 0.d0 
      do jg =1,mgrp   !jth group                         
        !jg2=jg+dvar
        jg2=jg
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        do iq2 =1,nq   
          lk = 1.d0      
          do mj = ind1,ind2 ! index for within subgroup
            nj=mj+mgrp
            call lgum1derivs(uvec(mj),xl(iq2),th(nj),lpdf(nj),der1(nj),der2(nj))   !c_{ij,V_j}   
            lk = lk*exp(lpdf(nj))
          end do
          call lfrk1derivs(xl(iq2),xl(iq),th(jg2),lpdf(jg2),der1(jg2),der2(jg2))   !c_{V_j,V_0}  
          llk = exp(lpdf(jg2))*lk  !llk: values of the jth integrand
          intj(jg) = intj(jg) + wl(iq2)*llk !intj value of the jth integral
          !intjj(ind1:ind2) = intj(jg) !intjj: repeat intj(j) grsize[j] times
          intjj(ind1+mgrp:ind2+mgrp) = intj(jg) !intjj: repeat intj(j) grsize[j] times
          do mj = ind1,ind2
            nj=mj+mgrp
            grdj(nj) = grdj(nj) + wl(iq2)*llk*der1(nj) !grdj: der. of the jth int. wrt th(mj)
            der2j(nj)= der2j(nj)+ wl(iq2)*llk*(der2(nj)+der1(nj)*der1(nj)) !der2j: second order der. wrt th(mj)
            hssaux(jg2,nj) = hssaux(jg2,nj) + wl(iq2)*llk*der1(nj)*der1(jg2)
            hssaux(nj,jg2) = hssaux(jg2,nj)
            if (mj < ind2) then
              do mj2 = (mj+1),ind2  !2nd order der. of the jth int. wrt th(mj), th(mj2)
                nj2=mj2+mgrp
                hssaux(nj,nj2) = hssaux(nj,nj2) + wl(iq2)*llk*der1(nj)*der1(nj2)
                hssaux(nj2,nj) = hssaux(nj,nj2)
              end do 
            endif         
          end do
          grdj(jg2) = grdj(jg2) + wl(iq2)*llk*der1(jg2) 
          der2j(jg2)= der2j(jg2)+ wl(iq2)*llk*(der2(jg2)+der1(jg2)*der1(jg2))
        end do
      end do
      !intjj(dvar+1:dvar+mgrp)=intj
      intjj(1:mgrp)=intj
      int0 = product(intj)   !product of j inner integrals  
      grd0 = int0*(grdj/intjj) !gradient of the product
      do ip=2,npar !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      !ind = 0
      ind = mgrp
      do jg=1,mgrp   !2nd order derivatives within groups 
        !jg2=jg+dvar
        jg2=jg
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,jg2)=int0*hssaux(ind1:ind2,jg2)/intj(jg)
        hss0(jg2,ind1:ind2)=int0*hssaux(jg2,ind1:ind2)/intj(jg)
      end do
      do ip=1,npar
        hss0(ip,ip) = int0*der2j(ip)/intjj(ip)
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

