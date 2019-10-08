! Code written by Pavel Krupskii and revised by HJ
! The subroutine lfrk1derivs is in file frk12fact-lpdf.f90
! The subroutine lbb1derivs is in file bb1facts.f90

! nested-factor model with Frank linking copulas for group to common latent
!         and BB1 linking copulas from observed to group latent.
! inputs
!   npar = mgrp+2*dvar = number of parameters
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
subroutine strfrkbb1(npar,th,mgrp,n,dvar,grsize,udata,nq,wl,xl,nllk,grad,hess)
  implicit none
  integer npar,mgrp,dvar,n,nq,ip,jp, ind
  double precision th(npar),udata(n,dvar),wl(nq),xl(nq),uvec(dvar)
  double precision lpdf(dvar+mgrp),der1(npar),der2(npar+dvar),llk, der1bb1(2), der2bb1(3)
  double precision nllk,liki,lk,grad(npar),hess(npar,npar)
  integer i,iq,iq2,jg,jg2,mj,mj2,nj,nj2,ind1,ind2,grsize(mgrp),ibb1(2)
  double precision intj(mgrp), intjj(npar), grdj(npar), der2j(npar), int0, grd0(npar), grdi(npar)
  double precision dermxj(dvar), hss0(npar,npar), hssi(npar,npar), hessi(npar,npar), hssaux(npar,npar)

  ! npar = 2*dvar+mgrp; dvar = sum(grsize)
  nllk=0.d0; grad = 0.d0; hess=0.d0 
  do i =1,n 
    uvec = udata(i,:)    
    liki = 0.d0; grdi = 0.d0; hssi = 0.d0; 
    do iq =1,nq
      int0 = 1.d0; grd0 = 1.d0; dermxj = 0.d0; hssaux = 0.d0 
      ind = 0; intj = 0.d0;  grdj = 0.d0; der2j = 0.d0  
      do jg =1,mgrp   !jth group
        !jg2=jg+dvar                         
        jg2=jg
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg);      
        do iq2 =1,nq   
          lk = 1.d0 
          do mj = ind1,ind2
            nj=mj+mgrp
            call lbb1derivs(uvec(mj),xl(iq2),th((/nj,nj+dvar/)),lpdf(nj),der1bb1,der2bb1)   !c_{ij,V_j}   
            der1((/nj,nj+dvar/)) = der1bb1
            der2((/nj,nj+2*dvar,nj+dvar/)) = der2bb1
            lk = lk*exp(lpdf(nj)) 
            end do
          call lfrk1derivs(xl(iq2),xl(iq),th(jg2),lpdf(jg2),der1(jg2),der2(jg2))    !c_{V_j,V_0}  
          llk = exp(lpdf(jg2))*lk !llk: values of the jth integrand
          intj(jg) = intj(jg) + wl(iq2)*llk !intj value of the jth integral
          intjj(ind1+mgrp:ind2+mgrp) = intj(jg) !intjj: repeat intj(j) grsize[j] times
          do mj = ind1,ind2
            nj=mj+mgrp
            ibb1 = (/nj,nj+dvar/)
            grdj(ibb1) = grdj(ibb1) + wl(iq2)*llk*der1(ibb1)  !grdj: der. of the jth int. wrt th(mj)
            der2j(ibb1)= der2j(ibb1)+ wl(iq2)*llk*(der2(ibb1)+der1(ibb1)*der1(ibb1)) !der2j: second order der. wrt th(mj)^2 (th(m+mj)^2)
            dermxj(mj) = dermxj(mj) + wl(iq2)*llk*(der2(nj+2*dvar)+der1(nj)*der1(nj+dvar)) !dermxj: mixed der. wrt th(mj)=th and th(m+mj)=dl
            hssaux(jg2,ibb1) = hssaux(jg2,ibb1) + wl(iq2)*llk*der1(ibb1)*der1(jg2)
            hssaux(ibb1,jg2) = hssaux(jg2,ibb1)
            if (mj < ind2) then
              do mj2 = (mj+1),ind2  !2nd order der. of the jth int. wrt th(mj), th(mj2)
                nj2=mj2+mgrp                                  
                hssaux(nj,nj2) = hssaux(nj,nj2) + wl(iq2)*llk*der1(nj)*der1(nj2) !th,th
                hssaux(nj+dvar,nj2+dvar) = hssaux(nj+dvar,nj2+dvar) + wl(iq2)*llk*der1(nj+dvar)*der1(nj2+dvar)  !dl,dl
                hssaux(nj+dvar,nj2) = hssaux(nj+dvar,nj2) + wl(iq2)*llk*der1(nj+dvar)*der1(nj2) !dl,th
                hssaux(nj,nj2+dvar) = hssaux(nj,nj2+dvar) + wl(iq2)*llk*der1(nj)*der1(nj2+dvar) !th,dl
                hssaux(nj2,nj) = hssaux(nj,nj2)
                hssaux(nj2+dvar,nj+dvar) = hssaux(nj+dvar,nj2+dvar)
                hssaux(nj2,nj+dvar) = hssaux(nj+dvar,nj2)
                hssaux(nj2+dvar,nj) = hssaux(nj,nj2+dvar)
              end do 
            endif         
          end do
          grdj(jg2) = grdj(jg2) + wl(iq2)*llk*der1(jg2) 
          der2j(jg2)= der2j(jg2)+ wl(iq2)*llk*(der2(jg2)+der1(jg2)*der1(jg2))
        end do
      end do
      intjj(1+dvar+mgrp:dvar+dvar+mgrp)=intjj(1+mgrp:mgrp+dvar)
      intjj(1:mgrp)=intj
      int0 = product(intj) !product of j inner integrals  
      grd0 = int0*(grdj/intjj) !gradient of the product
      do ip=2,npar !hessian of the product
        do jp=1,ip-1
          hss0(ip,jp) = int0*grdj(ip)*grdj(jp)/intjj(ip)/intjj(jp)
          hss0(jp,ip) = hss0(ip,jp)
        end do
      end do
      !ind = 0
      ind = mgrp
      do jg=1,mgrp  !2nd order derivatives within groups  
        !jg2=jg+dvar
        jg2=jg                           
        ind1 = 1 + ind;  ind2 = grsize(jg) + ind;  ind = ind + grsize(jg)
        hss0(ind1:ind2,ind1:ind2)=int0*hssaux(ind1:ind2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)=int0*hssaux(dvar+ind1:dvar+ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,ind1:ind2)=int0*hssaux(dvar+ind1:dvar+ind2,ind1:ind2)/intj(jg)
        hss0(ind1:ind2,dvar+ind1:dvar+ind2)=int0*hssaux(ind1:ind2,dvar+ind1:dvar+ind2)/intj(jg)
        hss0(ind1:ind2,jg2)=int0*hssaux(ind1:ind2,jg2)/intj(jg)
        hss0(jg2,ind1:ind2)=int0*hssaux(jg2,ind1:ind2)/intj(jg)
        hss0(dvar+ind1:dvar+ind2,jg2)=int0*hssaux(dvar+ind1:dvar+ind2,jg2)/intj(jg)
        hss0(jg2,dvar+ind1:dvar+ind2)=int0*hssaux(jg2,dvar+ind1:dvar+ind2)/intj(jg)
      end do

      do ip=1,npar
        hss0(ip,ip) = int0*der2j(ip)/intjj(ip)
      end do
      do ip=1+mgrp,dvar+mgrp
        hss0(ip,ip+dvar) = int0*dermxj(ip-mgrp)/intjj(ip)
        hss0(ip+dvar,ip) = hss0(ip,ip+dvar)
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
