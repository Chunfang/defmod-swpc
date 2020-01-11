! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module local

  use elems
  use esh3d
  implicit none
  real(8),allocatable :: ipoint(:,:),weight(:)

contains

  ! Form element [K]
  subroutine FormElK(ecoords,estress,E,nu,visc,expn,dt,k,strng)
    implicit none
    real(8) :: ecoords(npel,dmn),estress(nip,cdmn),E,nu,visc,expn,dt,          &
       k(eldof,eldof),alpha,D(cdmn,cdmn),S(cdmn,cdmn),dN(dmn,npel),detj,       &
       B(cdmn,eldof),betad(cdmn,cdmn)
    integer :: i
    character(2) :: strng
    k=f0
    call DMat(D,E,nu)
    if (strng=="Kv") then
       alpha=f1
       call MatInv(D,S)
    end if
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       if (strng=="Kv") then
          call Matbetad(betad,estress(i,:),visc,expn)
          call MatInv(S+alpha*dt*betad,D)
       end if
       k=k+matmul(transpose(B),matmul(D,B))*weight(i)*detj
    end do
  end subroutine FormElK

  ! Form element [K] for RVE (visco = false, poro = false)
  subroutine FormElKRVE(ecoords,E,nu,k,ellip)
    implicit none
    real(8) :: ecoords(npel,dmn),E,nu,k(eldof,eldof),D(cdmn,cdmn),dN(dmn,npel),&
       detj,B(cdmn,eldof),ellip(:,:),Dstr(cdmn,cdmn)
    integer :: i
    k=f0
    call DMat(D,E,nu)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       call FormDstr(ipoint(i,:),ecoords,E,nu,ellip,Dstr)
       D=matmul(D,Dstr)
       k=k+matmul(transpose(B),matmul(D,B))*weight(i)*detj
    end do
  end subroutine FormElKRVE

  ! Rescale element [Kv] for new dt
  subroutine RscElKv(ecoords,estress,E,nu,visc,expn,dt,k,ddt)
    implicit none
    real(8) :: ecoords(npel,dmn),estress(nip,cdmn),E,nu,visc,expn,dt,          &
       k(eldof,eldof),alpha,D(cdmn,cdmn),S(cdmn,cdmn),dN(dmn,npel),detj,       &
       B(cdmn,eldof),betad(cdmn,cdmn),ddt
    integer :: i
    k=f0
    call DMat(D,E,nu)
    alpha=f1
    call MatInv(D,S)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       call Matbetad(betad,estress(i,:),visc,expn)
       call MatInv(S+alpha*dt*betad,D)
       k=k-matmul(transpose(B),matmul(D,B))*weight(i)*detj
       call MatInv(S+alpha*(dt+ddt)*betad,D)
       k=k+matmul(transpose(B),matmul(D,B))*weight(i)*detj
    end do
  end subroutine RscElKv

  ! Form element [Kp]
  subroutine FormElKp(ecoords,estress,E,nu,visc,expn,H,Bc,fi,Kf,theta,scale,   &
    dt,k,strng)
    implicit none
    real(8) :: ecoords(npel,dmn),estress(nip,cdmn),E,nu,visc,expn,H,Bc,fi,Kf,  &
       theta,scale,dt,D(cdmn,cdmn),alpha,V(cdmn,cdmn),Q(dmn,dmn),s,m(cdmn,1),  &
       pN,G,Kb,Ksinv,dN(dmn,npel),detj,B(cdmn,eldof),N(1,npel),betad(cdmn,cdmn)
    real(8),target :: k(eldof+eldofp,eldof+eldofp)
    real(8),pointer :: kl(:,:)
    integer :: i
    character(2) :: strng
    k=f0; Q=f0
    call DMat(D,E,nu)
    if (strng=="Kv") then
       alpha=f1
       call MatInv(D,V)
    end if
    do i=1,dmn
       Q(i,i)=H
    end do
    s=scale
    if (dmn==2) m(:,1)=(/f1,f1,f0/)
    if (dmn==3) m(:,1)=(/f1,f1,f1,f0,f0,f0/)
    pN=f1/dble(npel)
    G=E/(f2*(f1+nu))
    Kb=E/(f3*(f1-f2*nu))
    Ksinv=(f1-Bc)/Kb
    if (strng=="Kp" .or. strng=="Kv") then
       do i=1,nip
          call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
          call BMat(dN,B)
          call ShapeFunc(N(1,:),ipoint(i,:))
          kl=>k(:eldof,:eldof)
          if (strng=="Kv") then
             call Matbetad(betad,estress(i,:),visc,expn)
             call MatInv(V+alpha*dt*betad,D)
          end if
          kl=kl+matmul(transpose(B),matmul(D,B))*weight(i)*detj
          kl=>k(eldof+1:,eldof+1:)
          kl=kl+matmul(transpose(dN),matmul(Q,dN))*weight(i)*detj*theta*dt*s*s
          kl=kl+matmul(transpose(N-pN),N-pN)*weight(i)*detj*s*s*0.5d0/G
          kl=kl+matmul(transpose(N),N)*((Bc-fi)*Ksinv+fi/Kf)*weight(i)*detj*s*s
          kl=>k(:eldof,eldof+1:)
          kl=kl+matmul(transpose(B),matmul(m,N))*Bc*weight(i)*detj*s
          kl=>k(eldof+1:,:eldof)
          kl=kl-transpose(matmul(transpose(B),matmul(m,N))*Bc*weight(i)*detj)*s
       end do
    end if
    if (strng=="Kc") then
       do i=1,nip
          call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
          kl=>k(eldof+1:,eldof+1:)
          kl=kl+matmul(transpose(dN),matmul(Q,dN))*weight(i)*detj*s*s
       end do
    end if
  end subroutine FormElKp

  ! Dynamic poroelastic stiffness
  subroutine AddElHSinvHt(ecoords,E,nu,Bc,fi,Kf,k)
    implicit none
    real(8) :: ecoords(npel,dmn),Bc,Kf,dN(dmn,npel),detj,B(cdmn,eldof),SinvB2, &
       m(cdmn,1),N(1,npel),E,nu,fi,Kb
    real(8),target :: k(eldof,eldof)
    integer :: i
    Kb=E/(f3*(f1-f2*nu))
    SinvB2=f1/((Bc-fi)*(f1-Bc)/Kb+fi/Kf)*Bc**2
    if (dmn==2) m(:,1)=(/f1,f1,f0/)
    if (dmn==3) m(:,1)=(/f1,f1,f1,f0,f0,f0/)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       call ShapeFunc(N(1,:),ipoint(i,:))
       k=k+matmul(transpose(B),matmul(matmul(matmul(m,N),transpose(matmul(m,   &
          N))),B))*SinvB2*weight(i)*detj
    end do
  end subroutine AddELHSinvHt

  ! Kp with tensor valued fluid mobility Q(dmn,dmn)
  subroutine FormElKPerm(ecoords,estress,E,nu,visc,expn,Q,Bc,fi,Kf,theta,      &
    scale,dt,k,strng)
    implicit none
    real(8) :: ecoords(npel,dmn),estress(nip,cdmn),E,nu,visc,expn,Bc,fi,Kf,    &
       theta,scale,dt,D(cdmn,cdmn),alpha,V(cdmn,cdmn),Q(dmn,dmn),s,m(cdmn,1),  &
       pN,G,Kb,Ksinv,dN(dmn,npel),detj,B(cdmn,eldof),N(1,npel),                &
       betad(cdmn,cdmn)
    real(8),target :: k(eldof+eldofp,eldof+eldofp)
    real(8),pointer :: kl(:,:)
    integer :: i
    character(2) :: strng
    k=f0
    call DMat(D,E,nu)
    if (strng=="Kv") then
       alpha=f1
       call MatInv(D,V)
    end if
    s=scale
    if (dmn==2) m(:,1)=(/f1,f1,f0/)
    if (dmn==3) m(:,1)=(/f1,f1,f1,f0,f0,f0/)
    pN=f1/dble(npel)
    G=E/(f2*(f1+nu))
    Kb=E/(f3*(f1-f2*nu))
    Ksinv=(f1-Bc)/Kb
    if (strng=="Kp" .or. strng=="Kv") then
       do i=1,nip
          call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
          call BMat(dN,B)
          call ShapeFunc(N(1,:),ipoint(i,:))
          kl=>k(:eldof,:eldof)
          if (strng=="Kv") then
             call Matbetad(betad,estress(i,:),visc,expn)
             call MatInv(V+alpha*dt*betad,D)
          end if
          kl=kl+matmul(transpose(B),matmul(D,B))*weight(i)*detj
          kl=>k(eldof+1:,eldof+1:)
          kl=kl+matmul(transpose(dN),matmul(Q,dN))*weight(i)*detj*theta*dt*s*s
          kl=kl+matmul(transpose(N-pN),N-pN)*weight(i)*detj*s*s*0.5d0/G
          kl=kl+matmul(transpose(N),N)*((Bc-fi)*Ksinv+fi/Kf)*weight(i)*detj*s*s
          kl=>k(:eldof,eldof+1:)
          kl=kl+matmul(transpose(B),matmul(m,N))*Bc*weight(i)*detj*s
          kl=>k(eldof+1:,:eldof)
          kl=kl-transpose(matmul(transpose(B),matmul(m,N))*Bc*weight(i)*detj)*s
       end do
    end if
    if (strng=="Kc") then
       do i=1,nip
          call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
          kl=>k(eldof+1:,eldof+1:)
          kl=kl+matmul(transpose(dN),matmul(Q,dN))*weight(i)*detj*s*s
       end do
    end if
  end subroutine FormElKPerm

  ! Rescale element [Kp] for new dt
  subroutine RscElKp(ecoords,estress,E,nu,visc,expn,H,theta,scale,dt,k,ddt,    &
    strng)
    implicit none
    real(8) :: ecoords(npel,dmn),estress(nip,cdmn),E,nu,visc,expn,H,theta,     &
       scale,dt,ddt,D(cdmn,cdmn),alpha,V(cdmn,cdmn),Q(dmn,dmn),s,dN(dmn,npel), &
       detj,B(cdmn,eldof),N(1,npel),betad(cdmn,cdmn)
    real(8),target :: k(eldof+eldofp,eldof+eldofp)
    real(8),pointer :: kl(:,:)
    integer :: i
    character(2) :: strng
    k=f0; Q=f0
    call DMat(D,E,nu)
    if (strng=="Kv") then
       alpha=f1
       call MatInv(D,V)
    end if
    do i=1,dmn
       Q(i,i)=H
    end do
    s=scale
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       call ShapeFunc(N(1,:),ipoint(i,:))
       if (strng=="Kv") then
          kl=>k(:eldof,:eldof)
          call Matbetad(betad,estress(i,:),visc,expn)
          call MatInv(V+alpha*dt*betad,D)
          kl=kl-matmul(transpose(B),matmul(D,B))*weight(i)*detj
          call MatInv(V+alpha*(dt+ddt)*betad,D)
          kl=kl+matmul(transpose(B),matmul(D,B))*weight(i)*detj
       end if
       kl=>k(eldof+1:,eldof+1:)
       kl=kl+matmul(transpose(dN),matmul(Q,dN))*weight(i)*detj*theta*ddt*s*s
    end do
  end subroutine RscElKp

  ! Form element Hs (-ve)
  subroutine FormElHs(ecoords,epr,E,nu,scale,f)
    implicit none
    real(8) :: ecoords(npel,dmn),epr(npel),E,nu,scale,f(eldofp),detj,N(npel),  &
       pN,pepr,G,prip
    integer :: i
    f=f0
    pN=f1/dble(npel)
    pepr=sum(epr)/dble(npel)
    G=E/(f2*(f1+nu))
    do i=1,nip
       call FormdetJ(ipoint(i,:),ecoords,detj)
       call ShapeFunc(N,ipoint(i,:))
       prip=dot_product(N,epr)
       f=f-(N-pN)*(prip-pepr)*weight(i)*detj*scale*0.5d0/G
    end do
  end subroutine FormElHs

  ! Form element [M]
  subroutine FormElM(ecoords,rho,m)
    implicit none
    real(8) :: ecoords(npel,dmn),rho,m(eldof,eldof),ms(npel,npel),N(1,npel),   &
       detj,val
    integer :: i
    m=f0; ms=f0
    do i=1,nip
       call FormdetJ(ipoint(i,:),ecoords,detj)
       call ShapeFunc(N(1,:),ipoint(i,:))
       ms=ms+rho*matmul(transpose(N),N)*weight(i)*detj
    end do
    do i=1,dmn
       m(i::dmn,i::dmn)=ms
    end do
    do i=1,eldof
       val=sum(m(i,:)); m(i,:)=f0; m(i,i)=val ! Diagonalize [M]
    end do
  end subroutine FormElM

  ! Form element damping (Abs) matrix [C]
  subroutine FormElAbsC(enodes,ecoords,eside,dir,E,nu,rho,c)
    implicit none
    real(8) :: ecoords(npel,dmn),E,nu,rho,c(eldof,eldof),Vp,Vs,area
    integer :: enodes(npel),eside,dir,i,j,j1,j2,snodes(nps)
    c=f0
    Vs=sqrt(E/(f2*(f1+nu)*rho))
    Vp=sqrt(((f1-nu)*E)/((f1-f2*nu)*(f1+nu)*rho))
    call EdgeAreaNodes(enodes,ecoords,eside,area,snodes)
    do i=1,nps
       do j=1,npel
          if (snodes(i)==enodes(j)) then
             do j1=1,dmn
                j2=dmn*j-dmn+j1
                if (j1==dir) c(j2,j2)=area*Vp*rho/dble(nps)
                if (j1/=dir) c(j2,j2)=area*Vs*rho/dble(nps)
             end do
          end if
       end do
    end do
  end subroutine FormElAbsC

  subroutine FormElAbsC1(enodes,ecoords,eside,matrot,E,nu,rho,c)
    implicit none
    real(8) :: ecoords(npel,dmn),E,nu,rho,c(eldof,eldof),Vp,Vs,area,           &
       matc(dmn,dmn),matrot(dmn,dmn)
    integer :: enodes(npel),eside,i,j,j1,j2,snodes(nps)
    c=f0
    Vs=sqrt(E/(f2*(f1+nu)*rho))
    Vp=sqrt(((f1-nu)*E)/((f1-f2*nu)*(f1+nu)*rho))
    call EdgeAreaNodes(enodes,ecoords,eside,area,snodes)
    do i=1,nps
       do j=1,npel
          if (snodes(i)==enodes(j)) then
             matc=f0
             j2=dmn*j-dmn
             ! Viscous matrix in facet coordinate
             do j1=1,dmn
                if (j1==dmn) matc(j1,j1)=area*Vp*rho/dble(nps)
                if (j1/=dmn) matc(j1,j1)=area*Vs*rho/dble(nps)
             end do
             ! Viscous matrix in global coordinate
             matc=matmul(matmul(matrot,matc),transpose(matrot))
             c(j2+1:j2+dmn,j2+1:j2+dmn)=matc
          end if
       end do
    end do
  end subroutine FormElAbsC1

  ! Form element index
  subroutine FormElIndx(enodes,indx)
    implicit none
    integer :: enodes(npel),indx(eldof),j,j1
    do j=1,npel
       do j1=1,dmn
          indx(dmn*j-j1+1)=dmn*enodes(j)-j1+1
       end do
    end do
  end subroutine FormElIndx

  ! Form element indexp
  subroutine FormElIndxp(enodes,indxp)
    implicit none
    integer :: enodes(npel),indxp(eldof+eldofp),j,j1
    do j=1,npel
       do j1=1,dmn
          indxp(dmn*j-j1+1)=(dmn+1)*enodes(j)-j1
       end do
       indxp(npel*dmn+j)=(dmn+1)*enodes(j)
    end do
  end subroutine FormElIndxp

  ! Calculate element stress
  subroutine CalcElStress(ecoords,edisp,E,nu,estress)
    implicit none
    real(8) :: ecoords(npel,dmn),edisp(eldof),E,nu,estrain(nip,cdmn),          &
       estress(nip,cdmn),D(cdmn,cdmn),dN(dmn,npel),detj,B(cdmn,eldof)
    integer :: i
    call DMat(D,E,nu)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       estrain(i,:)=matmul(B,edisp); estress(i,:)=matmul(D,estrain(i,:))
    end do
  end subroutine CalcElStress

  ! Calculate element strain
  subroutine CalcElStrain(ecoords,edisp,estrain)
    implicit none
    real(8) :: ecoords(npel,dmn),edisp(eldof),estrain(nip,cdmn),dN(dmn,npel),  &
       detj,B(cdmn,eldof)
    integer :: i
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       estrain(i,:)=matmul(B,edisp)
    end do
  end subroutine CalcElStrain

  ! Calculate element stress for RVE (visco = false, poro = false)
  subroutine CalcElStressRVE(ecoords,edisp,E,nu,estress,ellip)
    implicit none
    real(8) :: ecoords(npel,dmn),edisp(eldof),E,nu,estrain(nip,cdmn),          &
       estress(nip,cdmn),D(cdmn,cdmn),dN(dmn,npel),detj,B(cdmn,eldof),         &
       ellip(:,:),Dstr(cdmn,cdmn)
    integer :: i
    call DMat(D,E,nu)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       call FormDstr(ipoint(i,:),ecoords,E,nu,ellip,Dstr)
       D=matmul(D,Dstr)
       estrain(i,:)=matmul(B,edisp); estress(i,:)=matmul(D,estrain(i,:))
    end do
  end subroutine CalcElStressRVE

  ! Reform element RHS
  subroutine ReformElRHS(ecoords,estress,E,nu,visc,expn,dt,f)
    implicit none
    real(8) :: ecoords(npel,dmn),estress(nip,cdmn),E,nu,visc,expn,dt,f(eldof), &
       alpha,D(cdmn,cdmn),S(cdmn,cdmn),dN(dmn,npel),detj,B(cdmn,eldof),        &
       beta(cdmn),betad(cdmn,cdmn)
    integer :: i
    call DMat(D,E,nu); call MatInv(D,S)
    alpha=f1;f=f0
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       call Matbetad(betad,estress(i,:),visc,expn)
       call MatInv(S+alpha*dt*betad,D)
       call Matbeta(beta,estress(i,:),visc,expn)
       f=f+matmul(transpose(B),matmul(D,dt*beta))*weight(i)*detj
    end do
  end subroutine ReformElRHS

  ! Calculate element VStress
  subroutine CalcElVStress(ecoords,edisp,estress_inout,E,nu,visc,expn,dt)
    implicit none
    real(8) :: ecoords(npel,dmn),edisp(eldof),estress_inout(nip,cdmn),E,nu,    &
       visc,expn,dt,estress_in(nip,cdmn),estress_update(nip,cdmn),alpha,       &
       D(cdmn,cdmn),S(cdmn,cdmn),dN(dmn,npel),detj,B(cdmn,eldof),beta(cdmn),   &
       betad(cdmn,cdmn)
    integer :: i
    alpha=f1; estress_in=estress_inout
    call DMat(D,E,nu); call MatInv(D,S)
    do i=1,nip
       call FormdNdetJ(ipoint(i,:),ecoords,dN,detj)
       call BMat(dN,B)
       call Matbetad(betad,estress_in(i,:),visc,expn)
       call MatInv(S+alpha*dt*betad,D)
       call Matbeta(beta,estress_in(i,:),visc,expn)
       estress_update(i,:)=matmul(D,(matmul(B,edisp)-dt*beta))
    end do
    estress_inout=estress_in+estress_update
  end subroutine CalcElVStress

  ! Computes the strain displacement matrix 'B'
  subroutine BMat(dN,B)
    implicit none
    real(8) :: dN(dmn,npel),B(cdmn,eldof)
    integer :: j
    B=f0
    select case(dmn)
    case(2)
       do j=1,npel
          B(1,2*j-1)=dN(1,j); B(1,2*j)=f0
          B(2,2*j-1)=f0     ; B(2,2*j)=dN(2,j)
          B(3,2*j-1)=dN(2,j); B(3,2*j)=dN(1,j)
       end do
    case(3)
       do j=1,npel
          B(1,3*j-2)=dN(1,j); B(1,3*j-1)=f0     ; B(1,3*j)=f0
          B(2,3*j-2)=f0     ; B(2,3*j-1)=dN(2,j); B(2,3*j)=f0
          B(3,3*j-2)=f0     ; B(3,3*j-1)=f0     ; B(3,3*j)=dN(3,j)
          B(4,3*j-2)=dN(2,j); B(4,3*j-1)=dN(1,j); B(4,3*j)=f0
          B(5,3*j-2)=f0     ; B(5,3*j-1)=dN(3,j); B(5,3*j)=dN(2,j)
          B(6,3*j-2)=dN(3,j); B(6,3*j-1)=f0     ; B(6,3*j)=dN(1,j)
       end do
    end select
  end subroutine BMat

  ! (ipcoord, ecoords, ellip) -> Dstr(cdmn,cdmn)
  subroutine FormDstr(ipcoord,ecoords,Em,vm,ellip,Dstr)
    ! ipcoord(3), Gauss points of FE;
    ! ellip(nellip,11), nellip inclusions in one element, center, semi-axes,
    ! rotation angels, Eh, vh;
    ! Dstr(6,6), Strain=>elastic strain at Gauss points;
    ! Dstr(i,:,:)=[I+sum_e(Dei(Cm-dCeSe)^-1dCe)]^-1, i inp, e ellip;
    ! S2(6,6), Rank 2 first Eshelby tensors of each inclusion;
    ! D2(6,6), Rank 2 second Eshelby tensors of each inclusion and gauss point;
    ! D2e(6,6), elastic stran=>eigen strain (Cm-dCeSe)^-1dCe of each inclusion;
    ! Cm(6,6), Rank 2 host stiffness matrix;
    ! Ch(6,6), Rank 2 stiffness matrix of each inclusion;
    ! dC(6,6), Cm-Ch
    implicit none
    real(8) :: ipcoord(dmn),ecoords(npel,dmn),Dstr(cdmn,cdmn),dN(dmn,npel),    &
       jacob(dmn,dmn),invj(dmn,dmn),xobs(3),Em,vm,DInv(6,6),Cm(6,6),Ch(6,6),   &
       dC(6,6),PIvec(3),fderphi(3),tderpsi(3,3,3),S2(6,6),D4(3,3,3,3),D2(6,6), &
       D2e(6,6),R(3,3),Rb(3,3),R_init(3,3),Rb_init(3,3),R2(3,3),R2b(3,3),      &
       D2R(6,6),ang(3),a(3),ellip(:,:),Eh,vh
    integer :: i,nellip
    nellip=size(ellip,1)
    call ShapeFuncd(dN,ipcoord) ! Form dN/de, dN/dn, (dN/ds)
    jacob=matmul(dN,ecoords)
    call MatInv(jacob,invj)
    xobs=matmul(invj,ipcoord) ! (eps,eta,sig) -> (x,y,z)
    xobs=xobs*sqrt((ipcoord(1)**2+ipcoord(2)**2+ipcoord(3)**2)/(xobs(1)**2+    &
       xobs(2)**2+xobs(3)**2)) ! Normalize
    DInv=f0
    do i=1,cdmn
       DInv(i,i)=f1
    end do
    call CMat(Em,vm,Cm)
    do i=1,nellip
       a=ellip(i,4:6)
       call AxesSort(a,R_init,Rb_init)
       ! Rotation matrices w.r.t the ellipsoid
       ang=ellip(i,7:9)
       call Ang2Mat(ang,R,f1)
       call Ang2Mat(ang,Rb,-f1)
       R2=matmul(R_init,Rb)  ! Gauss=>ellipsoid
       R2b=matmul(R,Rb_init) ! Ellipsoid=>Gauss
       call EshS2(vm,a,S2,PIvec)
       Eh=ellip(i,10); vh=ellip(i,11)
       call CMat(Eh,vh,Ch); dC=Cm-Ch
       call MatInv(Cm-matmul(dC,S2),D2e)
       D2e=matmul(D2e,dC)
       xobs=xobs-ellip(i,:3) ! Relative coordinate
       xobs=matmul(R2,xobs)  ! Element=>ellipsoid
       if (xobs(1)**2/a(1)**2+xobs(2)**2/a(2)**2+xobs(3)**2/a(3)**2<=f1) then
          D2=S2
       else
          call EshD4(vm,a,xobs,D4,fderphi,tderpsi)
          call T4T2(D4,D2)
       end if
       D2=matmul(D2,D2e)
       call T2Rot(D2,R2b,D2R) ! Rotate to Cartesian before add
       DInv=DInv+D2R
    end do
    call MatInv(DInv,Dstr)
    ! Inspect
    !call MatDet(Dstr(:3,:3),tmp)
    !print'(7(F0.6X))',xobs,DStr(1,1),DStr(2,2),DStr(3,3),tmp
  end subroutine FormDstr

  ! Stage a1>=a2>=a3 and retruen rotation matrces
  subroutine AxesSort(a,R,Rb)
    implicit none
    integer :: i,j
    real(8) :: a(3),R1(3,3),R2(3,3),R3(3,3),R(3,3),Rb(3,3),exh(3,3),tmp,ang(3)
    exh=f0; R=f0
    do i=1,2
       R(i,i)=f1
       do j=2,3
          if (a(i)<a(j)) then
             exh(i,j)=f1
             tmp=a(i); a(i)=a(j); a(j)=tmp
          end if
       end do
    end do
    ang=pi/f2*(/f0,f0,exh(1,2)/)
    call Ang2Mat(ang,R1,f1)
    ang=pi/f2*(/f0,exh(1,3),f0/)
    call Ang2Mat(ang,R2,f1)
    ang=pi/f2*(/exh(2,3),f0,f0/)
    call Ang2Mat(ang,R3,f1)
    R=matmul(R3,matmul(R2,R1))
    Rb=transpose(R)
  end subroutine AxesSort

  ! Computes dN/dx(yz) and the determinant of Jacobian 'detj'
  subroutine FormdNdetJ(ipcoord,ecoords,dN,detj)
    implicit none
    real(8) :: ipcoord(dmn),detj,ecoords(npel,dmn),dN(dmn,npel),               &
       jacob(dmn,dmn),invj(dmn,dmn)
    call ShapeFuncd(dN,ipcoord) ! Form dN/de, dN/dn, (dN/ds)
    jacob=matmul(dN,ecoords)
    call MatDet(jacob,detj)
    call MatInv(jacob,invj)
    dN=matmul(invj,dN) ! Form dN/dx, dN/dy, (dN/dz)
  end subroutine FormdNdetJ

  ! Computes only the determinant of Jacobian 'detj' for use in Hs
  subroutine FormdetJ(ipcoord,ecoords,detj)
    implicit none
    real(8) :: ipcoord(dmn),detj,ecoords(npel,dmn),dN(dmn,npel),               &
       jacob(dmn,dmn)
    call ShapeFuncd(dN,ipcoord) ! Form dN/de, dN/dn, (dN/ds)
    jacob=matmul(dN,ecoords)
    call MatDet(jacob,detj)
  end subroutine FormdetJ

  ! Computes D, the matrix of elastic properties
  subroutine DMat(D,E,nu)
    implicit none
    real(8) :: D(:,:),E,nu
    if (size(D,1)==3) call DMat2d(D,E,nu)
    if (size(D,1)==6) call DMat3d(D,E,nu)
  end subroutine DMat

  ! Computes D, the matrix of elastic properties
  subroutine DMat2d(D,E,nu)
    implicit none
    real(8) :: D(3,3),E,nu,c
    c=E/((f1+nu)*(f1-f2*nu))
    D=c*reshape((/(f1-nu),nu,f0,nu,(f1-nu),f0,f0,f0,(f1-f2*nu)/f2/),           &
       shape=(/3,3/))
  end subroutine DMat2d

  ! Computes D, the matrix of elastic properties
  subroutine DMat3d(D,E,nu)
    implicit none
    real(8) :: D(6,6),E,nu,c
    c=E/((f1+nu)*(f1-f2*nu))
    D=c*reshape((/(f1-nu),nu,nu,f0,f0,f0,nu,(f1-nu),nu,f0,f0,f0,nu,nu,(f1-nu), &
       f0,f0,f0,f0,f0,f0,(f1-f2*nu)/f2,f0,f0,f0,f0,f0,f0,(f1-f2*nu)/f2,f0,f0,  &
       f0,f0,f0,f0,(f1-f2*nu)/f2/),shape=(/6,6/))
  end subroutine DMat3d

  ! Computes \beta(stress), viscoelastic strain rate
  subroutine Matbeta(beta,estress,visc,expn)
    implicit none
    real(8) :: beta(:),estress(:),visc,expn
    if (size(estress,1)==3) call Matbeta2d(beta,estress,visc,expn)
    if (size(estress,1)==6) call Matbeta3d(beta,estress,visc,expn)
  end subroutine Matbeta

  ! Computes \beta(stress), viscoelastic strain rate
  subroutine Matbeta2d(beta,estress,visc,expn)
    implicit none
    real(8) :: beta(3),estress(3),kappa,visc,expn,cMat(3,3),s1,s2,s3
    s1=estress(1); s2=estress(2); s3=estress(3)
    kappa=sqrt(((s1-s2)/f2)**2+s3**2)
    cMat=reshape((/f1,-f1,f0,-f1,f1,f0,f0,f0,f4/),shape=(/3,3/))
    beta=((kappa**(expn-f1))/(f4*visc))*matmul(cMat,estress)
  end subroutine Matbeta2d

  ! Computes \beta(stress), viscoelastic strain rate
  subroutine Matbeta3d(beta,estress,visc,expn)
    implicit none
    real(8) :: beta(6),estress(6),kappa,visc,expn,cMat(6,6),s1,s2,s3,s4,s5,s6
    s1=estress(1); s2=estress(2); s3=estress(3)
    s4=estress(4); s5=estress(5); s6=estress(6)
    kappa=sqrt(((s1-s2)**2+(s2-s3)**2+(s1-s3)**2)/f6+s4**2+s5**2+s6**2)
    cMat=reshape((/                                                            &
        f4/f3,-f2/f3,-f2/f3, f0 , f0 , f0 ,                                    &
       -f2/f3, f4/f3,-f2/f3, f0 , f0 , f0 ,                                    &
       -f2/f3,-f2/f3, f4/f3, f0 , f0 , f0 ,                                    &
          f0 ,   f0 ,   f0 , f4 , f0 , f0 ,                                    &
          f0 ,   f0 ,   f0 , f0 , f4 , f0 ,                                    &
          f0 ,   f0 ,   f0 , f0 , f0 , f4 /),shape=(/6,6/))
    beta=((kappa**(expn-f1))/(f4*visc))*matmul(cMat,estress)
  end subroutine Matbeta3d

  ! Computes \beta', Jacobian Matrix (differentiate \beta w.r.t. components of
  ! stress)
  subroutine Matbetad(betad,estress,visc,expn)
    implicit none
    real(8) :: betad(:,:),estress(:),visc,expn
    if (size(estress,1)==3) call Matbetad2d(betad,estress,visc,expn)
    if (size(estress,1)==6) call Matbetad3d(betad,estress,visc,expn)
  end subroutine Matbetad

  ! Computes \beta', Jacobian Matrix
  subroutine Matbetad2d(betad,estress,visc,expn)
    implicit none
    real(8) :: betad(3,3),estress(3),kappa,c1,c2,c3,visc,expn,s1,s2,s3
    s1=estress(1); s2=estress(2); s3=estress(3)
    kappa=sqrt(((s1-s2)/f2)**2+s3**2)
    if (kappa==f0) then
       betad=f0; return
    end if
    c1=f1+(expn-f1)*((s1-s2)/(f2*kappa))**2
    c2=f1+(expn-f1)*(s3/kappa)**2
    c3=(expn-f1)*(s1*s3-s2*s3)/(kappa)**2
    betad=((kappa**(expn-f1))/(f4*visc))*reshape((/c1,-c1,c3,-c1,c1,-c3,c3,    &
       -c3,f4*c2/),shape=(/3,3/))
  end subroutine Matbetad2d

  ! Computes \beta', Jacobian Matrix
  subroutine Matbetad3d(betad,estress,visc,expn)
    implicit none
    real(8) :: betad(6,6),estress(6),kappa,visc,expn,c,Sx,Sy,Sz,T1,T2,T3,s1,   &
       s2,s3,s4,s5,s6
    s1=estress(1); s2=estress(2); s3=estress(3)
    s4=estress(4); s5=estress(5); s6=estress(6)
    kappa=sqrt(((s1-s2)**2+(s2-s3)**2+(s1-s3)**2)/f6+s4**2+s5**2+s6**2)
    if (kappa==f0) then
       betad=f0; return
    end if
    c=sqrt(expn-f1)
    Sx=c*(f2*s1-s2-s3)/(f3*kappa)
    Sy=c*(f2*s2-s3-s1)/(f3*kappa)
    Sz=c*(f2*s3-s1-s2)/(f3*kappa)
    T1=c*f2*s4/kappa; T2=c*f2*s5/kappa; T3=c*f2*s6/kappa
    betad=((kappa**(expn-f1))/(f4*visc))*reshape((/                            &
          f4/f3+Sx**2,-f2/f3+Sx*Sy,-f2/f3+Sx*Sz,  Sx*T1 ,  Sx*T2 ,  Sx*T3 ,    &
         -f2/f3+Sx*Sy, f4/f3+Sy**2,-f2/f3+Sy*Sz,  Sy*T1 ,  Sy*T2 ,  Sy*T3 ,    &
         -f2/f3+Sx*Sz,-f2/f3+Sy*Sz, f4/f3+Sz**2,  Sz*T1 ,  Sz*T2 ,  Sz*T3 ,    &
             Sx*T1   ,    Sy*T1   ,    Sz*T1   ,f4+T1**2,  T1*T2 ,  T1*T3 ,    &
             Sx*T2   ,    Sy*T2   ,    Sz*T2   ,  T2*T1 ,f4+T2**2,  T2*T3 ,    &
             Sx*T3   ,    Sy*T3   ,    Sz*T3   ,  T3*T1 ,  T3*T2 ,f4+T3**2/),  &
          shape=(/6,6/))
  end subroutine Matbetad3d

end module local
