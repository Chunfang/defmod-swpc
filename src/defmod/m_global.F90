! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module global

#include <petscversion.h>

  use local
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
  implicit none
#include "petscdef.h"
#else
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
#endif
  ! Global variables
  integer :: nnds,nels,nmts,nceqs,nfrcs,ntrcs,nabcs,p,frq,dsp,dsp_hyb,lm_str,  &
     bod_frc,steps,tstep,steps_dyn,tstep_dyn,nfnd,hyb,frq_dyn,nobs,n_log,      &
     n_log_dyn,n_log_wave,ih,igf,rsf,n_log_slip,nceqs_ncf,init,frq_wave,       &
     frq_slip,n_lmnd,lmnd0,nobs_loc,nfnd_loc,ngp,ngp_loc,fvin
  real(8) :: alpha,beta,t,dt,t_dyn,dt_dyn,t_lim,val,wt,scale,rslip,t_sta,t_abs,&
     t_hyb,v_bg,vtol,trunc
  integer,allocatable :: nodes(:,:),bc(:,:),id(:),work(:),fnode(:),telsd(:,:), &
     worku(:),workl(:),node_pos(:),node_neg(:),slip(:),perm(:),onlst(:,:),     &
     frc(:),slip_sum(:),slip0(:),idgp(:,:),idgp_loc(:,:),gpnlst(:,:),gpl2g(:), &
     nnd_fe2fd(:),fdact_loc(:),matFD(:,:)
  real(8),allocatable :: coords(:,:),mat(:,:),stress(:,:,:),vvec(:),cval(:,:), &
     fval(:,:),tval(:,:),vecf(:,:),fc(:),matf(:,:),st_init(:,:),xfnd(:,:),     &
     ocoord(:,:),oshape(:,:),fcd(:),dc(:),rsfb0(:),rsfV0(:),biot(:),           &
     rsfdtau0(:),rsfa(:),rsfb(:),rsfL(:),rsftheta(:),coh(:),dcoh(:),mu_hyb(:), &
     mu_cap(:),rsfv(:),ocoord_loc(:,:),xgp(:,:),gpshape(:,:)
  real(8),allocatable,target :: uu(:),tot_uu(:),uup(:),uu_dyn(:),tot_uu_dyn(:),&
     fl(:),ql(:),flc(:),fp(:),qu(:),ss(:),sh(:),f2s(:),dip(:),nrm(:),          &
     flt_slip(:),tot_flt_slip(:),res_flt_slip(:),qs_flt_slip(:)
  character(12) :: stype
  character(256) :: output_file
  logical :: poro,visco,fault,dyn,fail,write_dyn,crp,gf,kfv
  Vec :: Vec_F,Vec_U,Vec_Um,Vec_Up,Vec_lambda,Vec_I,Vec_lambda_tot,Vec_U_dyn,  &
     Vec_Um_dyn,Vec_U_dyn_tot,Vec_Up_dyn,Vec_I_dyn,Vec_fp,Vec_qu,Vec_Uu,Vec_Ul,&
     Vec_fl,Vec_flc,Vec_ql,Vec_SS,Vec_SH,Vec_f2s,Vec_dip,Vec_nrm,Vec_Cp,       &
     Vec_lambda_sta,Vec_lambda_sta0,Vec_lambda_bk,Vec_lm_pn,Vec_lm_pp,         &
     Vec_lm_f2s,Vec_Fm_dyn,Vec_F_dyn,Vec_Cp0,Vec_Up_fv,Vec_Up_hst,Vec_Cst 
  Vec,pointer :: Vec_W(:),Vec_Wlm(:)
  Mat :: Mat_K,Mat_M,Mat_Minv,Mat_Gt,Mat_G,Mat_GMinvGt,Mat_Kc,Mat_K_dyn,Mat_H, &
     Mat_Ht
  KSP :: Krylov
  PC :: PreCon
  ! Local element/side/node variables
  integer :: el,side,node
  real(8) :: E,nu,dns,visc,expn,H,B,phi,Kf
  integer,allocatable :: indx(:),indxp(:),enodes(:),indx_dyn(:)
  real(8),allocatable :: k(:,:),m(:,:),f(:),ecoords(:,:),kc(:,:),Hs(:),        &
     k_dyn(:,:),uu_obs(:,:),tot_uu_obs(:,:),uu_dyn_obs(:,:),                   &
     tot_uu_dyn_obs(:,:),uu_fd(:,:)
  ! Variables for parallel code
  integer :: nprcs,rank,ierr
  integer,allocatable :: epart(:),npart(:) ! Partitioning
  integer,allocatable :: nmap(:),emap(:),nl2g(:,:),indxmap(:,:),               &
     indxmap_u(:,:),FltMap(:,:),ol2g(:) ! L-G Mapping
  Vec :: Seq_U,Seq_U_dyn,Seq_fp,Seq_fl,Seq_flc,Seq_ql,Seq_qu,Seq_SS,Seq_SH,    &
     Seq_f2s,Seq_dip,Seq_nrm
  IS :: From,To,RI,From_u,To_u,RIu,From_p,To_p,RIl
  VecScatter :: Scatter,Scatter_dyn,Scatter_u,Scatter_q,Scatter_s2d,           &
     Scatter_pn2d,Scatter_pp2d
  real(8),pointer :: pntr(:)

contains

  ! Form local [K]
  subroutine FormLocalK(el,k,indx,strng)
    implicit none
    integer :: el,indx(:)
    real(8) :: k(:,:),estress(nip,cdmn)
    character(2) :: strng
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    if (visco) estress=stress(el,:,:)
    if (dyn) then
       E=mat(id(el),5+4*p+init+1); nu=mat(id(el),5+4*p+init+2)
    else
       E=mat(id(el),1); nu=mat(id(el),2)
    end if
    visc=mat(id(el),3); expn=mat(id(el),4)
    if ((.not. poro) .or. strng=="Ke") then
       call FormElK(ecoords,estress,E,nu,visc,expn,dt,k,strng)
    else
       H=mat(id(el),6)
       B=mat(id(el),7); phi=mat(id(el),8); Kf=mat(id(el),9)
       call FormElKp(ecoords,estress,E,nu,visc,expn,H,B,phi,Kf,1.0d0,scale,dt, &
          k,strng)
    end if
    call AddWinklerFdn(el,k)
    if (.not. dyn) call FixBCinLocalK(el,k)
    if (poro .and. kfv) call FixFVLocalK(el,k)
    if (dyn) then
       call FormLocalIndx_dyn(enodes,indx)
    else
       call FormLocalIndx(enodes,indx)
    end if
  end subroutine FormLocalK

  ! Custom permeability from FV
  subroutine FormLocalKPerm(el,k,indx,m_perm,strng)
    implicit none
    integer :: el,indx(:)
    real(8) :: k(:,:),estress(nip,cdmn),m_perm(dmn,dmn)
    character(2) :: strng
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    if (visco) estress=stress(el,:,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    B=mat(id(el),7); phi=mat(id(el),8); Kf=mat(id(el),9)
    call FormElKPerm(ecoords,estress,E,nu,visc,expn,m_perm,B,phi,Kf,1.0d0,     &
       scale,dt,k,strng)
    call AddWinklerFdn(el,k)
    call FixBCinLocalK(el,k)
    if (strng=="Kp") call FixFVLocalKBD(el,k)
    call FormLocalIndx(enodes,indx)
  end subroutine FormLocalKPerm

  ! Rescale local [Kv] for dt
  subroutine RscKv(el,k,indx,ddt)
    implicit none
    integer :: el,indx(:)
    real(8) :: k(:,:),estress(nip,cdmn),ddt
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    estress=stress(el,:,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    H=mat(id(el),6)
    call RscElKv(ecoords,estress,E,nu,visc,expn,dt,k,ddt)
    call FormLocalIndx(enodes,indx) 
  end subroutine RscKv

  ! Rescale local [Kp] for dt
  subroutine RscKp(el,k,indx,ddt)
    implicit none
    integer :: el,indx(:)
    real(8) :: k(:,:),estress(nip,cdmn),ddt
    character(2) :: strng
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    if (visco) estress=stress(el,:,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    H=mat(id(el),6)
    strng="Kp"
    if (visco) strng="Kv"
    call RscElKp(ecoords,estress,E,nu,visc,expn,H,1.0d0,scale,dt,k,ddt,strng)
    call FormLocalIndx(enodes,indx)
  end subroutine RscKp 

  ! Rescale [K] for new dt = fdt*dt (poro and/or linear visco) 
  subroutine Rscdt(fdt)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    integer :: i,ndof
    real(8) :: fdt,ddt
    ddt=(fdt-f1)*dt
    ndof=eldof+eldofp
    do i=1,nels
       if (visco .and. .not. poro) then
          call RscKv(i,k,indx,ddt)
          !call RscRHSv(i,f,indx,ddt)
       else 
          call RscKp(i,k,indx,ddt)
          !call RscRHSp(i,f,indx,ddt)
       end if
       indx=indxmap(indx,2)
       call MatSetValues(Mat_K,ndof,indx,ndof,indx,k,Add_Values,ierr)
       !call VecSetValues(Vec_F,eldof,indx,f,Add_Values,ierr)
    end do
    call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
    !call VecAssemblyBegin(Vec_F,ierr)
    !call VecAssemblyEnd(Vec_F,ierr)
    dt=fdt*dt
  end subroutine Rscdt

  ! Form local [M]
  subroutine FormLocalM(el,m,indx)
    implicit none
    integer :: el,indx(:)
    real(8) :: m(:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    dns=mat(id(el),5)
    call FormElM(ecoords,dns,m)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalM

  ! Account for constraint eqns
  subroutine ApplyConstraints
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: i,j,n,j1,j2,j3,j4
    open(15,file=trim(output_file(:index(output_file,"/",BACK=.TRUE.)))//      &
       "cnstrns.tmp",status='old')
    j4=0
    do i=1,nceqs
       if (poro .and. mod(i,dmn+1)>0) j4=j4+1
       read(15,*)n
       do j=1,n
          read(15,*)vvec,node
          if (stype/="explicit") then
             do j1=1,dmn+p
                j2=(dmn+p)*node-(dmn+p)+j1-1; j3=(dmn+p)*nnds+i-1
                call MatSetValue(Mat_K,j2,j3,wt*vvec(j1),Add_Values,ierr)
                call MatSetValue(Mat_K,j3,j2,wt*vvec(j1),Add_Values,ierr)
             end do
          end if
          if ((stype=="explicit" .and. .not. gf) .or. (fault .and.             &
             i<=nceqs_ncf .and. .not. kfv)) then ! Consider non-conforming grid
             do j1=1,dmn
                j2=dmn*node-dmn+j1-1; j3=i-1
                if (poro .and. mod(i,dmn+1)>0) then
                   call MatSetValue(Mat_Gt,j2,j4-1,vvec(j1),Add_Values,ierr)
                elseif (.not. poro) then
                   call MatSetValue(Mat_Gt,j2,j3,vvec(j1),Add_Values,ierr)
                end if
             end do
          end if
       end do
       if (stype/="explicit") then ! Constraint block diagonals have to be
          j1=(dmn+p)*nnds+i-1      ! explicitly set to zero (PETSc req)
          call MatSetValue(Mat_K,j1,j1,f0,Add_Values,ierr)
       end if
       read(15,*)cval(i,:)
    end do
    close(15)
  end subroutine ApplyConstraints

  ! Create full Mat_Gt for dynamic constraints
  subroutine GetMat_Gt
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,j4,j5
    if (rank==0) then
       do j=1,nfnd
          do j1=1,dmn
             j3=(j-1)*dmn+j1-1
             do j2=1,2
                if (j2==1) then 
                   node=node_pos(j)
                   vvec(:dmn)=vecf(j,(j1-1)*dmn+1:j1*dmn)
                else
                   node=node_neg(j)
                   vvec(:dmn)=-vecf(j,(j1-1)*dmn+1:j1*dmn) 
                end if
                do j4=1,dmn
                   j5=dmn*node-dmn+j4-1
                   call MatSetValue(Mat_Gt,j5,j3+nceqs_ncf*dmn/(dmn+p),vvec(j4)&
                      ,Add_Values,ierr)
                end do
             end do
          end do
       end do
    end if 
  end subroutine GetMat_Gt

  ! Scatter LMs to solution space
  subroutine GetVec_flambda
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j3,j4,u1,u2,workpos(dmn),workneg(dmn),row(1)
    real(8) :: lm(dmn),vecfl(dmn),lmq(1),q
    call VecGetOwnershipRange(Vec_Ul,u1,u2,ierr)
    do j=1,nfnd
       lm=f0; vecfl=f0; lmq=f0; q=f0
       do j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
          workneg(j1)=dmn*node_neg(j)-dmn+j1-1
          if (poro) then
             row=dmn*(j-1)+sum(perm(1:j-1))+j1-1
          else
             row=dmn*(j-1)+j1-1
          end if
          if (row(1)>=u1 .and. row(1)<u2) then
             call VecGetValues(Vec_Ul,1,row,lm(j1),ierr)
          end if
       end do
       call MPI_Reduce(lm,vecfl,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,  &
          ierr)
       if (poro .and. perm(j)==1) then
           row=row+1
           if (row(1)>=u1 .and. row(1)<u2) then
              call VecGetValues(Vec_Ul,1,row,lmq,ierr)
           end if
           call MPI_Reduce(lmq,q,1,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,   &
              ierr)
       end if
       if (rank==nprcs-1) then 
          call VecSetValues(Vec_fl,dmn,workpos,wt*vecfl,Insert_Values,ierr)
          call VecSetValues(Vec_fl,dmn,workneg,-wt*vecfl,Insert_Values,ierr)
          if (poro .and. perm(j)==1) then
             j3=node_pos(j)-1
             j4=node_neg(j)-1
             call VecSetValue(Vec_ql,j3,wt*q,Insert_Values,ierr)
             call VecSetValue(Vec_ql,j4,-wt*q,Insert_Values,ierr)
          end if
       end if
    end do
    call VecAssemblyBegin(Vec_fl,ierr)
    call VecAssemblyEnd(Vec_fl,ierr)
    if (poro) then
       call VecAssemblyBegin(Vec_ql,ierr)
       call VecAssemblyEnd(Vec_ql,ierr)
    end if
  end subroutine GetVec_flambda

  ! Extract the LM in fault's strike, dip and normal directions
  subroutine GetVec_fcoulomb
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,u1,u2,r1,r2,workpos(dmn),workneg(dmn),row(1),row_f2s(1)
    real(8) :: lm(dmn),vecfl(dmn),mattmp(dmn,dmn),vectmp(dmn,1),rvec(dmn),     &
       rf2s(dmn)
    integer,save :: k=0
    call VecGetOwnershipRange(Vec_Ul,u1,u2,ierr)
    if (k>0) call VecGetOwnershipRange(Vec_f2s,r1,r2,ierr)
    do j=1,nfnd
       lm=f0; vecfl=f0; rvec=f0; rf2s=f0
       do j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
          workneg(j1)=dmn*node_neg(j)-dmn+j1-1
          if (poro) then
             row=dmn*(j-1)+sum(perm(1:j-1))+j1-1
          else
             row=dmn*(j-1)+j1-1
          end if
          if (row(1)>=u1 .and. row(1)<u2) then
             call VecGetValues(Vec_Ul,1,row,lm(j1),ierr)
          end if
          if (k>0) then
             row_f2s=dmn*node_pos(j)-dmn+j1-1
             if (row_f2s(1)>=r1 .and. row_f2s(1)<r2) then
                call VecGetValues(Vec_f2s,1,row_f2s,rvec(j1),ierr)
             end if
          end if
       end do
       call MPI_Reduce(lm,vecfl,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,  &
          ierr)
       if (k>0) call MPI_Reduce(rvec,rf2s,dmn,MPI_Real8,MPI_Sum,nprcs-1,       &
          MPI_Comm_World,ierr)
       if (rank==nprcs-1) then
          vectmp=reshape(vecfl,(/dmn,1/))
          mattmp=transpose(reshape(vecf(j,:),(/dmn,dmn/)))
          vectmp=matmul(mattmp,vectmp)
          vecfl(:)=vectmp(:,1)
          if (k>0) vecfl=vecfl*rf2s
          call VecSetValues(Vec_flc,dmn,workpos,wt*vecfl,Insert_Values,ierr)
          call VecSetValues(Vec_flc,dmn,workneg,-wt*vecfl,Insert_Values,ierr)
       end if
    end do
    call VecAssemblyBegin(Vec_flc,ierr)
    call VecAssemblyEnd(Vec_flc,ierr)
    k=k+1
  end subroutine GetVec_fcoulomb

  ! Pass pseudo velocity to dynamic model
  subroutine Rsfv2dyn
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer:: j,j1,j2,j3,rw_loc(dmn)
    real(8) :: vec(dmn),flt_qs(dmn)
    real(8),target :: flt_ndf(n_lmnd*dmn)
    call VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_ndf=pntr
    call VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    do j=1,nfnd_loc
       j1=FltMap(j,1); j3=FltMap(j,2)
       rw_loc=(/((j1-1)*dmn+j2,j2=1,dmn)/)
       select case(dmn)
       case(2)
          vec(1)=rsfv(j3)
          vec(2)=f0 ! Zero normal velocity
       case(3)
          flt_qs=flt_ndf(rw_loc)+st_init(j3,:)
          vec(1)=rsfv(j3)*flt_qs(1)/sqrt(flt_qs(1)**2+flt_qs(2)**2)
          vec(2)=rsfv(j3)*flt_qs(2)/sqrt(flt_qs(1)**2+flt_qs(2)**2)
          vec(3)=f0 ! Zero normal velocity
       end select
       rw_loc=lmnd0*dmn+rw_loc-1
       call VecSetValues(Vec_I_dyn,dmn,rw_loc,vec*dt_dyn,Insert_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_I_dyn,ierr) 
    call VecAssemblyEnd(Vec_I_dyn,ierr)
  end subroutine Rsfv2Dyn

  ! Force to stress ratio and normal dip vector of the fault nodes
  subroutine GetVec_f2s
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,j4,j5,workpos(dmn),workneg(dmn),row_p(1)
    real(8) :: vecfl(dmn),vecss(dmn),vecsh(dmn),r,dip(dmn),nrm(dmn),pn(1),     &
       pp(1),ptmp,matrot(dmn,dmn),matst(dmn,dmn),st(dmn,dmn),vec(dmn)
    call VecGetOwnershipRange(Vec_flc,j2,j3,ierr)
    if (poro) call VecGetOwnershipRange(Vec_Um,j4,j5,ierr)
    do j=1,nfnd
       matrot=reshape(vecf(j,:),(/dmn,dmn/))
       vecfl=f0; vecss=f0; vecsh=f0
       do j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
          if (workpos(j1)>=j2 .and. workpos(j1)<j3) then
             call VecGetValues(Vec_flc,1,workpos(j1),vecfl(j1),ierr)
             call VecGetValues(Vec_SS,1,workpos(j1),vecss(j1),ierr)
             call VecGetValues(Vec_SH,1,workpos(j1),vecsh(j1),ierr)
          end if
       end do
       call MPI_Reduce(vecfl,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecfl=vec
       call MPI_Reduce(vecss,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecss=vec
       call MPI_Reduce(vecsh,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecsh=vec
       pn=f0; pp=f0
       if (poro) then
          row_p=(dmn+1)*node_neg(j)-1
          if (row_p(1)>=j4 .and. row_p(1)<j5) then
             call VecGetValues(Vec_Um,1,row_p,pn,ierr)
          end if
          row_p=(dmn+1)*node_pos(j)-1
          if (row_p(1)>=j4 .and. row_p(1)<j5) then
             call VecGetValues(Vec_Um,1,row_p,pp,ierr)
          end if
          call MPI_Reduce(pn,ptmp,1,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,  &
             ierr) 
          pn=ptmp
          call MPI_Reduce(pp,ptmp,1,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,  &
             ierr) 
          pp=ptmp
       end if
       if (rank==nprcs-1) then
          select case(dmn)
          case(2)
             st(1,1)=vecss(1); st(2,2)=vecss(2)
             st(1,2)=vecsh(1); st(2,1)=vecsh(2)
             dip(:)=matrot(:,1)
             nrm(:)=matrot(:,2)
          case(3)
             st(1,1)=vecss(1); st(2,2)=vecss(2); st(3,3)=vecss(3)
             st(1,2)=vecsh(1); st(2,3)=vecsh(2); st(1,3)=vecsh(3)
             st(2,1)=vecsh(1); st(3,2)=vecsh(2); st(3,1)=vecsh(3)
             dip(:)=matrot(:,2)
             nrm(:)=matrot(:,3)
          end select
          matst=matmul(matmul(transpose(matrot),st),matrot)
          vecss=matst(:,dmn)
          vecss(dmn)=vecss(dmn)-biot(j)*pp(1)*scale
          r=sqrt(sum(vecss*vecss)/sum(vecfl*vecfl))
          call VecSetValues(Vec_f2s,dmn,workpos,(/r,r,r/),Insert_Values,ierr)
          call VecSetValues(Vec_dip,dmn,workpos,dip,Insert_Values,ierr)
          call VecSetValues(Vec_nrm,dmn,workpos,nrm,Insert_Values,ierr)
       end if
       matrot=reshape(vecf(j,:),(/dmn,dmn/))
       vecfl=f0; vecss=f0; vecsh=f0
       do j1=1,dmn
          workneg(j1)=dmn*node_neg(j)-dmn+j1-1
          if (workneg(j1)>=j2 .and. workneg(j1)<j3) then
             call VecGetValues(Vec_flc,1,workneg(j1),vecfl(j1),ierr)
             call VecGetValues(Vec_SS,1,workneg(j1),vecss(j1),ierr)
             call VecGetValues(Vec_SH,1,workneg(j1),vecsh(j1),ierr)
          end if
       end do
       call MPI_Reduce(vecfl,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecfl=vec
       call MPI_Reduce(vecss,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecss=vec
       call MPI_Reduce(vecsh,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
          ierr)
       vecsh=vec
       if (rank==nprcs-1) then
          select case(dmn)
          case(2)
             st(1,1)=vecss(1); st(2,2)=vecss(2)
             st(1,2)=vecsh(1); st(2,1)=vecsh(2)
             dip(:)=matrot(:,1)
             nrm(:)=matrot(:,2)
          case(3)
             st(1,1)=vecss(1); st(2,2)=vecss(2); st(3,3)=vecss(3)
             st(1,2)=vecsh(1); st(2,3)=vecsh(2); st(1,3)=vecsh(3)
             st(2,1)=vecsh(1); st(3,2)=vecsh(2); st(3,1)=vecsh(3)
             dip(:)=matrot(:,2)
             nrm(:)=matrot(:,3)
          end select
          matst=matmul(matmul(transpose(matrot),st),matrot)
          vecss=matst(:,dmn)
          vecss(dmn)=vecss(dmn)-biot(j)*pn(1)*scale
          r=(r+sqrt(sum(vecss*vecss)/sum(vecfl*vecfl)))/f2
          ! Convert prestress to nodal force
          if (r>f0) then
             st_init(j,:)=st_init(j,:)/r
             coh(j)=coh(j)/r
             if (rsf==1) rsfdtau0(j)=rsfdtau0(j)/r
          end if
          call VecSetValues(Vec_f2s,dmn,workneg,-(/r,r,r/),Insert_Values,ierr)
          call VecSetValue(Vec_lm_f2s,nceqs_ncf/(dmn+1)+j-1,r,Insert_Values,   &
             ierr)
          call VecSetValues(Vec_dip,dmn,workneg,dip,Insert_Values,ierr)
          call VecSetValues(Vec_nrm,dmn,workneg,nrm,Insert_Values,ierr)
       end if
    end do
    call VecAssemblyBegin(Vec_f2s,ierr)
    call VecAssemblyEnd(Vec_f2s,ierr)
    call VecAssemblyBegin(Vec_lm_f2s,ierr)
    call VecAssemblyEnd(Vec_lm_f2s,ierr)
    call VecAssemblyBegin(Vec_dip,ierr)
    call VecAssemblyEnd(Vec_dip,ierr)
    call VecAssemblyBegin(Vec_nrm,ierr)
    call VecAssemblyEnd(Vec_nrm,ierr)
    call MPI_Bcast(st_init,nfnd*dmn,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
    call MPI_Bcast(coh,nfnd,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
    if (rsf==1) call MPI_Bcast(rsfdtau0,nfnd,MPI_Real8,nprcs-1,                &
       MPI_Comm_World,ierr)
  end subroutine GetVec_f2s

  ! Fault dip and normal vectors
  subroutine GetVec_ft
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,workpos(dmn)
    real(8):: dip(dmn),nrm(dmn),matrot(dmn,dmn)
    do j=1,nfnd
       do j1=1,dmn
          workpos(j1)=dmn*node_pos(j)-dmn+j1-1
       end do
       matrot=reshape(vecf(j,:),(/dmn,dmn/))
       select case(dmn)
       case(2)
          dip(:)=matrot(:,1)
          nrm(:)=matrot(:,2)
       case(3)
          dip(:)=matrot(:,2)
          nrm(:)=matrot(:,3)
       end select
       if (rank==nprcs-1) then
          call VecSetValues(Vec_dip,dmn,workpos,dip,Insert_Values,ierr)
          call VecSetValues(Vec_nrm,dmn,workpos,nrm,Insert_Values,ierr)
       end if
    end do
    call VecAssemblyBegin(Vec_dip,ierr)
    call VecAssemblyEnd(Vec_dip,ierr)
    call VecAssemblyBegin(Vec_nrm,ierr)
    call VecAssemblyEnd(Vec_nrm,ierr)
  end subroutine GetVec_ft

  ! RSF pseudo time update 
  subroutine RSF_QS_update(flt_ndf0,flt_ndf1,slip_loc)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,j4,rw_loc(dmn),nqs,ncut,slip_loc(nfnd),ntol,nslip
    real(8) :: theta,mu,a,b0,b,V0,L,dd,v_qs,dtpsd,flt_qs0(dmn),flt_qs1(dmn),   &
       flt_qs(dmn),vmax
    real(8),target :: flt_ndf0(n_lmnd*dmn),flt_ndf1(n_lmnd*dmn)
    real(8),allocatable :: rsfstate(:,:,:)
    integer,save :: k=0,rsflog=0
    dtpsd=f1
    nqs=int((dt)/dtpsd)
    allocate(rsfstate(nqs,nfnd_loc,3))
    rsfstate=f0; slip_loc=0
    ncut=nqs
    vtol=2.5D-4 ! Velocity threshold
    vmax=1.0D-3 ! Maximum velocity
    ntol=10 ! Slip node threshold 
    do j4=1,nqs
       do j=1,nfnd_loc
          j1=FltMap(j,1); j3=FltMap(j,2)
          ! RSF parameters
          a=rsfa(j3); b0=rsfb0(j3); b=rsfb(j3); V0=rsfV0(j3); L=rsfL(j3)
          if (j4==1) then
             theta=rsftheta(j3)
          else
             theta=rsfstate(j4-1,j,3)
          end if
          rw_loc=(/((j1-1)*dmn+j2,j2=1,dmn)/)
          flt_qs0=flt_ndf0(rw_loc)
          flt_qs1=flt_ndf1(rw_loc)
          ! Match shear with friction by updating v_qs
          flt_qs=flt_qs0+(flt_qs1-flt_qs0)*dble(j4)/dble(nqs)+st_init(j3,:)
          mu=sqrt(sum(flt_qs(:dmn-1)*flt_qs(:dmn-1)))/abs(flt_qs(dmn))
          ! Prevent pseudo velocity overflow
          v_qs=min(vmax,sinh(mu/a)*V0*f2/exp((b0+b*log(V0*theta/L))/a))
          dd=v_qs*dtpsd/L
          ! Normal stress dependent theta (Dieterich 2007, alpha = 0.5)
          dd=dd+0.5/b*(flt_qs1(dmn)-flt_qs0(dmn))/dble(nqs)/flt_qs(dmn)
          theta=dtpsd/(f1+0.5*dd)+theta*(f1-dd*0.5)/(f1+dd*0.5)
          rsfstate(j4,j,:)=(/v_qs,mu,theta/)
       end do
       call MPI_AllReduce(size(pack(rsfstate(j4,:,1),rsfstate(j4,:,1)>vtol)), &
          nslip,1,MPI_Integer,MPI_Sum,MPI_Comm_World,ierr)
       ! At least ntol fault nodes nucleate
       if (nslip>=ntol) then 
          ncut=j4
          exit
       end if
    end do
    do j=1,nfnd_loc
       j3=FltMap(j,2)
       rsftheta(j3)=rsfstate(ncut,j,3)
       mu_hyb(j3)=rsfstate(ncut,j,2)
       rsfv(j3)=rsfstate(ncut,j,1)
       if (rsfv(j3)>vtol .and. nslip>=ntol) slip_loc(j3)=1
    end do 
    if (ncut<nqs) then 
       call VecCopy(Vec_lambda_sta,Vec_lambda_bk,ierr)
       call VecScale(Vec_lambda_sta,dble(ncut)/dble(nqs),ierr)
       call VecAXPY(Vec_lambda_sta,f1-dble(ncut)/dble(nqs),Vec_lambda_sta0,ierr)
       trunc=dt*(f1-dble(ncut)/dble(nqs))
    else
       trunc=f0
    end if
    do j4=1,ncut
       if (mod(k,200)==0) then
          call WriteOutput_rsf(rsfstate(j4,:,:2))
          rsflog=rsflog+1
       end if
       k=k+1
    end do
    if (rank==0) call WriteOutput_log_rsf(rsflog,dtpsd*dble(200))
  end subroutine RSF_QS_update

  ! Write rate-state pseudo velocity and friction
  subroutine WriteOutput_rsf(state)
    implicit none
    integer :: j,j3
    integer,save :: k=0
    real(8) :: state(nfnd_loc,2)
    character(256) :: name,name0,name1
    character(64) :: fmt
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    if (nfnd_loc>0) then 
       write(name,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_rsf_",rank,".txt"
       fmt="(2(ES11.2E3,X))"
       if (k==0) then
          open(10,file=adjustl(name),status='replace')
          write(10,'(I0)')nfnd_loc
          do j=1,nfnd_loc
             j3=FltMap(j,2)
             select case(dmn)
                case(2); write(10,'(2(F0.6,X),I0)')xfnd(j3,:),j3
                case(3); write(10,'(3(F0.6,X),I0)')xfnd(j3,:),j3
             end select
          end do
       else
          open(10,file=adjustl(name),status='old',position='append',action=    &
             'write')
       end if
       do j=1,nfnd_loc
          write(10,fmt)state(j,:)
       end do
       close(10); k=k+1
    end if
  end subroutine WriteOutput_rsf

  subroutine WriteOutput_log_rsf(n,dt)
    implicit none
    integer :: n
    integer,save :: k=0
    real(8) :: dt
    character(256) :: name,name0,name1
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_rsf.log"
    if (k==0) then
       open(10,file=adjustl(name),status='replace')
       write(10,'(F0.3)')dt
    else 
       open(10,file=adjustl(name),status='old',position='append',action='write')
    end if
    write(10,'(I0)')n
    close(10); k=k+1
  end subroutine WriteOutput_log_rsf

  subroutine GetTracDat(dat_trac) 
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,rw_loc(dmn)
    real(8) :: dat_trac(:,:,:),flt_qs(dmn),flt_p
    real(8),target :: flt_ndf(n_lmnd*dmn),lm_pn(n_lmnd),lm_pp(n_lmnd),lm_f2s(n_lmnd)
    call VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_ndf=pntr
    call VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    lm_pn=f0; lm_pp=f0; lm_f2s=f0 
    if (poro) then
       call VecGetArrayF90(Vec_lm_pn,pntr,ierr)
       lm_pn=pntr
       call VecRestoreArrayF90(Vec_lm_pn,pntr,ierr)
       call VecGetArrayF90(Vec_lm_pp,pntr,ierr)
       lm_pp=pntr
       call VecRestoreArrayF90(Vec_lm_pp,pntr,ierr)
       call VecGetArrayF90(Vec_lm_f2s,pntr,ierr)
       lm_f2s=pntr
       call VecRestoreArrayF90(Vec_lm_f2s,pntr,ierr)
    end if
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       flt_qs=(flt_ndf(rw_loc)+st_init(j3,:dmn))*lm_f2s(j)
       if (poro) then
          flt_p=(lm_pp(j)+lm_pn(j))/f2
          ! Positive pressure affects normal stress
          flt_qs(dmn)=flt_qs(dmn)+max(f0,biot(j3)*flt_p)
          dat_trac(1,:,j1)=(/flt_qs,flt_p/)
       else
          dat_trac(1,:,j1)=flt_qs
       end if
    end do
  end subroutine GetTracDat 

  ! Update slip from the static model (from Vec_lambda_sta)
  subroutine GetSlip_sta 
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,slip_loc(nfnd),rw_loc(dmn)
    real(8) :: mu,theta102,flt_qs(dmn),fsh,fnrm,rsftau,d,fcoh, &
       a,b0,b,V0,L
    real(8),target :: flt_ndf(n_lmnd*dmn),flt_ndf0(n_lmnd*dmn),lm_pn(n_lmnd),  &
       lm_pp(n_lmnd),lm_f2s(n_lmnd)
    integer,save :: k=0
    call VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_ndf=pntr
    call VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    lm_pn=f0; lm_pp=f0; lm_f2s=f1 ! None-zero default f2s, denominator
    if (poro) then
       call VecGetArrayF90(Vec_lm_pn,pntr,ierr)
       lm_pn=pntr
       call VecRestoreArrayF90(Vec_lm_pn,pntr,ierr)
       call VecGetArrayF90(Vec_lm_pp,pntr,ierr)
       lm_pp=pntr
       call VecRestoreArrayF90(Vec_lm_pp,pntr,ierr)
       call VecGetArrayF90(Vec_lm_f2s,pntr,ierr)
       lm_f2s=pntr
       call VecRestoreArrayF90(Vec_lm_f2s,pntr,ierr)
    end if
    if (rsf==1 .and. k>0) then
       call VecGetArrayF90(Vec_lambda_sta0,pntr,ierr)
       flt_ndf0=pntr
       call VecRestoreArrayF90(Vec_lambda_sta0,pntr,ierr)
       call RSF_QS_update(flt_ndf0,flt_ndf,slip_loc)
       go to 250
    end if
    if (rsf==1 .and. k==0) rsfv=v_bg
    slip_loc=0
    rsftau=f0 
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       flt_qs=flt_ndf(rw_loc)+st_init(j3,:dmn)
       ! Positive pressure affects normal stress
       flt_qs(dmn)=flt_qs(dmn)+max(f0,biot(j3)*(lm_pp(j)+lm_pn(j))/lm_f2s(j)/f2)
       if (rsf==1) then  
          if (k==0) then
             ! RSF parameters
             a=rsfa(j3); b0=rsfb0(j3); b=rsfb(j3); V0=rsfV0(j3); L=rsfL(j3)
             mu_cap(j3)=a*asinh(v_bg/V0/f2*exp((b0+b*log(V0*rsftheta(j3)/L))/a))
             mu=min(mu_cap(j3),sqrt(sum(flt_qs(:dmn-1)*flt_qs(:dmn-1)))        &
                /abs(flt_qs(dmn)))
             mu_hyb(j3)=mu
             theta102=rsftheta(j3)
             rsftheta(j3)=L/V0*exp((a*log(f2*sinh(mu/a))-b0-a*log(v_bg         &
                /V0))/b)
             ! Only for SCEC102
             mu=mu_cap(j3)
             rsftheta(j3)=theta102
             call GetExSt(xfnd(j3,:),t_abs,rsfdtau0(j3),rsftau)
          else
             mu=mu_hyb(j3)
          end if
       else ! Slip weakening
          d=sqrt(sum(res_flt_slip(rw_loc(:dmn-1))*res_flt_slip(rw_loc(:dmn-1))))
          if (d<dc(j3)) then
             mu=(f1-d/dc(j3))*(fc(j3)-fcd(j3))+fcd(j3)
          else
             mu=fcd(j3) 
          end if
       end if
       flt_qs(1)=flt_qs(1)+rsftau
       select case(dmn)
       case(2)
          fsh=abs(flt_qs(1))
          fnrm=flt_qs(2)
       case(3)
          fsh=sqrt(flt_qs(1)**2+flt_qs(2)**2)
          fnrm=flt_qs(3)
       end select
       ! Cohesive stress if any
       if (coh(j3)>f0) then
          d=sqrt(sum(qs_flt_slip(rw_loc(:dmn-1))*qs_flt_slip(rw_loc(:dmn-1))))
          if (d<dcoh(j3)) then 
             fcoh=coh(j3)*(f1-d/dcoh(j3))
          else
             fcoh=f0
          end if
       else
           fcoh=f0
       end if
       if ((fnrm<f0 .and. abs((fsh-fcoh)/fnrm)>mu) .or. fnrm>f0) then
          slip_loc(j3)=1
       else
          slip_loc(j3)=0
       endif
    end do
250 call MPI_AllReduce(slip_loc,slip,nfnd,MPI_Integer,MPI_Sum,MPI_Comm_World,  &
       ierr)
    slip0=slip
    slip_sum=slip
    k=k+1
  end subroutine GetSlip_sta 

  ! Scatter from Vec_Ul to Vec_lambda_sta 
  subroutine LM_s2d
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,rw_dyn(dmn),rw_loc(dmn),rw_sta(dmn),                 &
       idxmp(nfnd_loc*dmn,2)
    integer,save :: k=0
    real(8) :: vec(dmn),mattmp(dmn,dmn),vectmp(dmn,1)
    real(8),target :: flt_sta(n_lmnd*dmn)
    if (k==0) then
       do j1=1,nfnd_loc
          j=FltMap(j1,1); j3=FltMap(j1,2)
          rw_loc=(/((j1-1)*dmn+j2,j2=1,dmn)/)
          rw_dyn=lmnd0*dmn+(/((j-1)*dmn+j2,j2=1,dmn)/)-1
          idxmp(rw_loc,1)=rw_dyn
          if (poro) then
             rw_sta=(/((j3-1)*dmn+sum(perm(1:j3-1))+j2-1,j2=1,dmn)/)
          else
             rw_sta=(/((j3-1)*dmn+j2-1,j2=1,dmn)/)
          end if
          idxmp(rw_loc,2)=rw_sta
       end do
       call ISCreateGeneral(Petsc_Comm_World,nfnd_loc*dmn,idxmp(:,2),          &
          Petsc_Copy_Values,From,ierr)
       call ISCreateGeneral(Petsc_Comm_World,nfnd_loc*dmn,idxmp(:,1),          &
          Petsc_Copy_Values,To,ierr)
       call VecScatterCreate(Vec_Ul,From,Vec_lambda_sta,To,Scatter_s2d,ierr)
    end if
    call VecScatterBegin(Scatter_s2d,Vec_Ul,Vec_lambda_sta,Insert_Values,      &
       Scatter_Forward,ierr)
    call VecScatterEnd(Scatter_s2d,Vec_Ul,Vec_lambda_sta,Insert_Values,        &
       Scatter_Forward,ierr)
    ! Rotate and scale
    call VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_sta=pntr
    call VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       rw_dyn=lmnd0*dmn+rw_loc-1
       vec=flt_sta(rw_loc)*wt
       vectmp=reshape(vec,(/dmn,1/))
       mattmp=transpose(reshape(vecf(j3,:),(/dmn,dmn/)))
       vectmp=matmul(mattmp,vectmp)
       vec=vectmp(:,1)
       call VecSetValues(Vec_lambda_sta,dmn,rw_dyn,vec,Insert_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_lambda_sta,ierr)
    call VecAssemblyEnd(Vec_lambda_sta,ierr)
    k=k+1
  end subroutine LM_s2d 

  ! Scatter two side fault pressure to dynamic LM space 
  subroutine Up_s2d
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,jn,jp,idxmp_n(nfnd_loc,2),idxmp_p(nfnd_loc,2)
    integer,save :: k=0
    if (k==0) then 
       do j1=1,nfnd_loc
          j=FltMap(j1,1)
          jn=node_neg(FltMap(j1,2))
          jp=node_pos(FltMap(j1,2))
          idxmp_n(j1,1)=lmnd0+j-1
          idxmp_p(j1,1)=lmnd0+j-1
          idxmp_n(j1,2)=jn-1
          idxmp_p(j1,2)=jp-1
       end do
       call ISCreateGeneral(Petsc_Comm_World,nfnd_loc,idxmp_n(:,2),            &
          Petsc_Copy_Values,From,ierr)
       call ISCreateGeneral(Petsc_Comm_World,nfnd_loc,idxmp_n(:,1),            &
          Petsc_Copy_Values,To,ierr)
       call VecScatterCreate(Vec_Up,From,Vec_lm_pn,To,Scatter_pn2d,ierr)
       call ISCreateGeneral(Petsc_Comm_World,nfnd_loc,idxmp_p(:,2),            &
          Petsc_Copy_Values,From,ierr)
       call ISCreateGeneral(Petsc_Comm_World,nfnd_loc,idxmp_p(:,1),            &
          Petsc_Copy_Values,To,ierr)
       call VecScatterCreate(Vec_Up,From,Vec_lm_pp,To,Scatter_pp2d,ierr)
    end if
    call VecScatterBegin(Scatter_pn2d,Vec_Up,Vec_lm_pn,Insert_Values,          &
       Scatter_Forward,ierr)
    call VecScatterEnd(Scatter_pn2d,Vec_Up,Vec_lm_pn,Insert_Values,            &
       Scatter_Forward,ierr)
    call VecScatterBegin(Scatter_pp2d,Vec_Up,Vec_lm_pp,Insert_Values,          &
       Scatter_Forward,ierr)
    call VecScatterEnd(Scatter_pp2d,Vec_Up,Vec_lm_pp,Insert_Values,            &
       Scatter_Forward,ierr)
    call VecScale(Vec_lm_pn,scale,ierr)
    call VecScale(Vec_lm_pp,scale,ierr)
    k=k+1
  end subroutine Up_s2d
  
  ! Cap dynamic LM by frictional laws 
  subroutine CapLM_dyn
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer ::j,j1,j2,j3,rw_loc(dmn)
    real(8) :: d,dd,fr,frs,frd,fsh,mu,fcoh,fnrm,rsftau,vec_init(dmn),vec(dmn), &
       vec_slip(dmn-1),lm_sta(dmn),lm_dyn(dmn),lm_dyn0(dmn)
    real(8),target :: flt_sta(n_lmnd*dmn),flt_dyn(n_lmnd*dmn),                 &
       flt_dyn0(n_lmnd*dmn),lm_pn(n_lmnd),lm_pp(n_lmnd),lm_f2s(n_lmnd)
    call VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_sta=pntr
    call VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    call VecGetArrayF90(Vec_lambda,pntr,ierr)
    flt_dyn=pntr
    call VecRestoreArrayF90(Vec_lambda,pntr,ierr)
    call VecGetArrayF90(Vec_lambda_tot,pntr,ierr)
    flt_dyn0=pntr
    call VecRestoreArrayF90(Vec_lambda_tot,pntr,ierr)
    lm_pn=f0; lm_pp=f0; lm_f2s=f1 ! None-zero default f2s, denominator
    if (poro) then
       call VecGetArrayF90(Vec_lm_pn,pntr,ierr)
       lm_pn=pntr
       call VecRestoreArrayF90(Vec_lm_pn,pntr,ierr)
       call VecGetArrayF90(Vec_lm_pp,pntr,ierr)
       lm_pp=pntr
       call VecRestoreArrayF90(Vec_lm_pp,pntr,ierr)
       call VecGetArrayF90(Vec_lm_f2s,pntr,ierr)
       lm_f2s=pntr
       call VecRestoreArrayF90(Vec_lm_f2s,pntr,ierr)
    end if
    rsftau=f0
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       vec_init=st_init(j3,:)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       lm_sta=flt_sta(rw_loc)
       lm_dyn0=flt_dyn0(rw_loc)
       lm_dyn=flt_dyn(rw_loc)
       vec_slip=tot_flt_slip(rw_loc(:dmn-1))+res_flt_slip(rw_loc(:dmn-1))
       d=sqrt(sum(vec_slip*vec_slip))
       vec=vec_init+lm_sta+lm_dyn0+lm_dyn 
       ! Positive pressure affects friction
       vec(dmn)=vec(dmn)+max(f0,biot(j3)*(lm_pp(j)+lm_pn(j))/lm_f2s(j)/f2)
       if (rsf==1) then ! Rate state friction
          ! RSF options
          !dd=max(v_bg*dt_dyn,sqrt(sum(flt_slip(rw_loc(:dmn-1))                 &
          !   *flt_slip(rw_loc(:dmn-1)))))
          !dd=v_bg*dt_dyn+sqrt(sum(flt_slip(rw_loc(:dmn-1))                     &
          !   *flt_slip(rw_loc(:dmn-1))))
          !dd=max(rsfv(j3)*dt_dyn,sqrt(sum(flt_slip(rw_loc(:dmn-1))             &
          !   *flt_slip(rw_loc(:dmn-1)))))
          !dd=rsfv(j3)*dt_dyn+sqrt(sum(flt_slip(rw_loc(:dmn-1))                 &
          !   *flt_slip(rw_loc(:dmn-1))))
          dd=max(v_bg*dt_dyn,flt_slip(rw_loc(1))) ! Only for SCEC102
          mu=rsfa(j3)*asinh(dd/dt_dyn/rsfV0(j3)/f2*exp((rsfb0(j3)              &
             +rsfb(j3)*log(rsfV0(j3)*rsftheta(j3)/rsfL(j3)))/rsfa(j3)))
          dd=dd/rsfL(j3)
          ! Normal stress dependent theta (Dieterich 2007, alpha = 0.5)
          !if (vec(dmn)<f0) then ! Only when under compression, not for SCEC102 
          !   dd=dd+0.5/rsfb(j3)*lm_dyn(dmn)/(vec(dmn)-lm_dyn(dmn)*0.5)
          !end if
          ! Prevent negative theta
          rsftheta(j3)=max(dt_dyn,dt_dyn/(f1+0.5*dd)+rsftheta(j3)*(f1-0.5*dd)/ &
             (f1+0.5*dd))
          ! Only for SCEC102
          call GetExSt(xfnd(j3,:),t_hyb,rsfdtau0(j3),rsftau)
       else ! Slip weakening
          if (d<dc(j3)) then
             mu=(f1-d/dc(j3))*(fc(j3)-fcd(j3))+fcd(j3)
          else
             mu=fcd(j3)
          end if
       end if
       mu_hyb(j3)=mu
       vec(1)=vec(1)+rsftau
       select case(dmn)
       case(2)
          fsh=abs(vec(1))
          fnrm=vec(2)
       case(3)
          fsh=sqrt(vec(1)**2+vec(2)**2)
          fnrm=vec(3)
       end select
       ! Slip weakening cohesion
       if (coh(j3)>f0) then
          ! Cumulative slip 
          d=d+sqrt(sum(qs_flt_slip(rw_loc(:dmn-1))*qs_flt_slip(rw_loc(:dmn-1))))
          if (d<dcoh(j3)) then
             fcoh=coh(j3)*(f1-d/dcoh(j3))
          else
             fcoh=f0
          end if
       else
          fcoh=f0
       end if
       ! If LM exceeds maximum friction (mu*fn)
       if (fnrm<f0 .and. abs(fsh)>mu*abs(fnrm)+fcoh) then
          fr=mu*abs(fnrm)+fcoh 
          frs=fr*abs(vec(1)/fsh)
          frs=sign(frs,vec(1))
          vec(1)=frs
          if (dmn==3) then
             frd=fr*abs(vec(2)/fsh)
             frd=sign(frd,vec(2))
             vec(2)=frd
          end if
       ! When fault faces detach under cohesion
       elseif (fnrm>=f0) then 
          vec=vec/sqrt(sum(vec*vec))*min(sqrt(sum(vec*vec)),fcoh)
       end if
       ! Subtract the static LM 
       vec(1)=vec(1)-rsftau
       lm_dyn=vec-vec_init-lm_sta-lm_dyn0
       lm_dyn(dmn)=lm_dyn(dmn)-max(f0,biot(j3)*(lm_pp(j)+lm_pn(j))/lm_f2s(j)/f2)
       rw_loc=lmnd0*dmn+rw_loc-1
       ! From the new dynamic LM
       call VecSetValues(Vec_lambda,dmn,rw_loc,lm_dyn,Insert_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_lambda,ierr)
    call VecAssemblyEnd(Vec_lambda,ierr)
  end subroutine CapLM_dyn

  ! Record latest moment fault stress from the hybrid run
  subroutine GetVec_lambda_hyb
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: i,j,j1,j2,j3,rw_loc(dmn),nqs
    real(8) :: vec(dmn),lm_sta(dmn),lm_bk(dmn),lm_dyn(dmn),dtpsd,a,b0,b,V0,L,  &
       mu,dd
    real(8),target :: flt_sta(n_lmnd*dmn),flt_bk(n_lmnd*dmn),flt_dyn(n_lmnd*dmn)
    call VecGetArrayF90(Vec_lambda_sta,pntr,ierr)
    flt_sta=pntr
    call VecRestoreArrayF90(Vec_lambda_sta,pntr,ierr)
    if (rsf==1) then
       call VecGetArrayF90(Vec_lambda_bk,pntr,ierr)
       flt_bk=pntr
       call VecRestoreArrayF90(Vec_lambda_bk,pntr,ierr)
    end if
    call VecGetArrayF90(Vec_lambda_tot,pntr,ierr)
    flt_dyn=pntr
    call VecRestoreArrayF90(Vec_lambda_tot,pntr,ierr)
    dtpsd=f1
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       lm_sta=flt_sta(rw_loc)
       lm_dyn=flt_dyn(rw_loc)
       if (rsf==1) then ! Pseudo time healing for RSF
          lm_bk=flt_bk(rw_loc)
          nqs=int(trunc/dtpsd)
          ! RSF parameters
          a=rsfa(j3); b0=rsfb0(j3); b=rsfb(j3); V0=rsfV0(j3); L=rsfL(j3)
          do i=1,nqs 
             vec=lm_sta+lm_dyn+(lm_bk-lm_sta)*dble(i)/dble(nqs)+st_init(j3,:)
             mu=sqrt(sum(vec(:dmn-1)*vec(:dmn-1)))/abs(vec(dmn))
             rsfv(j3)=min(vtol,sinh(mu/a)*V0*f2/exp((b0+b*log(V0*rsftheta(j3)/ &
                L))/a))
             dd=dtpsd*rsfv(j3)/L
             ! Normal stress dependent theta (Dieterich 2007, alpha = 0.5)
             dd=dd+0.5/b*(lm_bk(dmn)-lm_sta(dmn))/dble(nqs)/vec(dmn)
             rsftheta(j3)=dtpsd/(f1+0.5*dd)+rsftheta(j3)*(f1-0.5*dd)/(f1+0.5*dd)
          end do
          vec=vec-st_init(j3,:)
       else
          vec=lm_sta+lm_dyn
       end if
       rw_loc=lmnd0*dmn+rw_loc-1
       call VecSetValues(Vec_lambda_sta0,dmn,rw_loc,vec,Insert_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_lambda_sta0,ierr)
    call VecAssemblyEnd(Vec_lambda_sta0,ierr)
  end subroutine GetVec_lambda_hyb

  ! Time dependent initial stress (SCEC102)
  subroutine GetExSt(x,t,dst,st)
    implicit none
    real(8) :: t,x(dmn),rx2,dst,rsfF,rsfG,R2,rsfT,st
    rsfF=f0;rsfG=f1;R2=9.0;rsfT=f1;st=f0
    if (dmn==3) then
       rx2=x(2)**2+(x(3)+7.5)**2
    else
       rx2=x(1)**2
    end if
    if (rx2<R2) rsfF=exp(rx2/(rx2-R2))
    if (t>f0 .and. t<rsfT) rsfG=exp((t-rsfT)**2/t/(t-2*rsfT))
    st=dst*rsfF*rsfG
  end subroutine GetExSt

  ! Determine if the dynamic slip is stabilized
  subroutine GetSlip_dyn
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,i,nc,nr,rw_loc(dmn),slip_loc(nfnd),slip_sum_loc(nfnd)
    real(8) :: d0,d
    slip_loc=slip
    slip_sum_loc=slip_sum
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       d=sqrt(sum(flt_slip(rw_loc(:dmn-1))*flt_slip(rw_loc(:dmn-1))))
       d0=sqrt(sum(tot_flt_slip(rw_loc(:dmn-1))*tot_flt_slip(rw_loc(:dmn-1))))
       ! Stabilization tolerance 
       if (rsf==1) then ! Rate state friction
          if (d0>5.0D-2*rsfL(j3)) then
             if ((d/d0)<1.0D-4) then
                slip_loc(j3)=0
             else
                slip_loc(j3)=1
                slip_sum_loc(j3)=1
             end if
          end if
          if (d/dt_dyn>vtol) then
            slip_loc(j3)=1
            slip_sum_loc(j3)=1
          end if
       else ! Slip weakening
          if (d0>5.0D-2*dc(j3)) then 
             if ((d/d0)<1.0D-4) then
                slip_loc(j3)=0
             else
                slip_loc(j3)=1
                slip_sum_loc(j3)=1
             end if
          end if
       end if
    end do
    ! Zero off-rank slip
    if (nfnd_loc>0) then
       j1=FltMap(1,2)-1; j2=FltMap(nfnd_loc,2)+1
       if (j1>0) then
          slip_loc(:j1)=0
          slip_sum_loc(:j1)=0
       end if
       if (j2<=nfnd) then
          slip_loc(j2:)=0
          slip_sum_loc(j2:)=0
       end if
    else
       slip_loc=0
       slip_sum_loc=0
    end if
    call MPI_AllReduce(slip_loc,slip,nfnd,MPI_Integer,MPI_Sum,                 &
       MPI_Comm_World,ierr)
    call MPI_AllReduce(slip_sum_loc,slip_sum,nfnd,MPI_Integer,MPI_Sum,         &
       MPI_Comm_World,ierr)
    ! Identify aseismic slip, nc=10 for SCEC10/14 (slow weakening)
    nc=5; nr=15
    if (ih>nc+rsf*nr .and. sum((slip0-slip_sum)*(slip0-slip_sum))==0) then 
       crp=.true.
       slip=0
       do j1=1,nfnd_loc 
          j=FltMap(j1,1); j3=FltMap(j1,2)
          rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
          if (slip_loc(j3)>0) res_flt_slip(rw_loc)=res_flt_slip(rw_loc)+       &
             tot_flt_slip(rw_loc)
       end do
    elseif (ih>nc+rsf*nr) then
       i=0
       do j=1,nfnd
          if (slip0(j)==0 .and. slip(j)==1) i=i+1
       end do
       if (i==0) slip=0 ! No slip except nucleation patch
    elseif (ih<=nc+rsf*nr .and. sum(slip)==0) then 
       crp=.true. ! Zero slip within time nc+rsf*nr 
    end if
  end subroutine GetSlip_dyn 

  ! Get fault QS state 
!  subroutine GetVec_flt_qs
!    implicit none
!#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
!#include "petsc.h"
!#endif
!    integer :: l1,l2,j,j1,r1,r2,row_l(1),row_f2s(1),row_p(1) 
!    real(8) :: lm(dmn),vec(dmn),rvec(dmn),rf2s(dmn),mattmp(dmn,dmn),           &
!       vectmp(dmn,1),pn(1),pp(1),fltpn,fltpp
!    call VecGetOwnershipRange(Vec_Um,l1,l2,ierr)
!    call VecGetOwnershipRange(Vec_f2s,r1,r2,ierr)
!    do j=1,nfnd
!       lm=f0; vec=f0; rf2s=0; rvec=f0; pn=f0; pp=f0; fltpn=f0; fltpp=f0
!       do j1=1,dmn 
!          if (poro) then
!             row_l=(dmn+1)*nnds+nceqs_ncf+(j-1)*dmn+sum(perm(1:j-1))+j1-1
!          else
!             row_l=dmn*nnds+nceqs_ncf+(j-1)*dmn+j1-1
!          end if
!          if (row_l(1)>=l1 .and. row_l(1)<l2) then
!             call VecGetValues(Vec_Um,1,row_l,lm(j1),ierr)
!          end if
!          row_f2s=dmn*node_pos(j)-dmn+j1-1
!          if (row_f2s(1)>=r1 .and. row_f2s(1)<r2) then
!             call VecGetValues(Vec_f2s,1,row_f2s,rf2s(j1),ierr)
!          end if
!       end do
!       if (poro) then
!          row_p=(dmn+1)*node_neg(j)-1
!          if (row_p(1)>=l1 .and. row_p(1)<l2) then
!             call VecGetValues(Vec_Um,1,row_p,pn,ierr)
!          end if
!          row_p=(dmn+1)*node_pos(j)-1
!          if (row_p(1)>=l1 .and. row_p(1)<l2) then
!             call VecGetValues(Vec_Um,1,row_p,pp,ierr)
!          end if
!       end if
!       call MPI_Reduce(lm,vec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,ierr)
!       call MPI_Reduce(rf2s,rvec,dmn,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World, &
!          ierr)
!       call MPI_Reduce(pn,fltpn,1,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,ierr)
!       call MPI_Reduce(pp,fltpp,1,MPI_Real8,MPI_Sum,nprcs-1,MPI_Comm_World,ierr)
!       ! Rotate vec to fault coordinate
!       if (rank==nprcs-1) then
!          vectmp=reshape(vec,(/dmn,1/))
!          mattmp=transpose(reshape(vecf(j,:),(/dmn,dmn/)))
!          vectmp=matmul(mattmp,vectmp)
!          vec(:)=vectmp(:,1)
!          ! Nodal force to stress
!          flt_ss(j,:)=(vec*wt+st_init(j,:))*rvec(dmn)
!          if (poro) then 
!             flt_p(j)=0.5*(fltpn+fltpp)*scale 
!             flt_ss(j,dmn)=flt_ss(j,dmn)+max(f0,biot(j)*flt_p(j))
!          end if 
!       end if
!    end do
!    call MPI_Bcast(flt_ss,nfnd*dmn,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
!    if (poro) call MPI_Bcast(flt_p,nfnd,MPI_Real8,nprcs-1,MPI_Comm_World,ierr)
!  end subroutine GetVec_flt_qs

  ! Add fault slip from dynamic model to static model as constraint functions
  subroutine FaultSlip
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: j,j1,j2,j3,rw_loc(dmn),rw_sta(dmn)
    real(8) :: vec(dmn),vectmp(dmn,1),mattmp(dmn,dmn)
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2) 
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       if (poro) then
          rw_sta=(dmn+1)*nnds+nceqs_ncf+(/((j3-1)*dmn+sum(perm(1:j3-1))        &
             +j2-1,j2=1,dmn)/)
       else
          rw_sta=dmn*nnds+nceqs_ncf+(/((j3-1)*dmn+j2-1,j2=1,dmn)/)
       end if
       vec=tot_flt_slip(rw_loc)
       ! Zero non-slip components, and rotate to Cartesian 
       vec(dmn)=f0 
       vectmp=reshape(vec,(/dmn,1/))
       mattmp=reshape(vecf(j3,:),(/dmn,dmn/))
       vectmp=matmul(mattmp,vectmp)
       vec=vectmp(:,1)
       call VecSetValues(Vec_F,dmn,rw_sta,vec*wt,Add_Values,ierr)
     end do 
  end subroutine FaultSlip

  ! Extract solution at observation locations
  subroutine GetVec_obs
    implicit none
    integer :: ob,i,j,ind(npel)
    integer,allocatable :: row(:)
    real(8) :: vecshp(npel,1)
    real(8),allocatable :: vectmp(:,:),mattmp(:,:)
    if (dyn .and. .not. gf) then
       allocate(row(dmn*npel),vectmp(dmn,1),mattmp(dmn,npel))
    else
       allocate(row((dmn+p)*npel),vectmp(dmn+p,1),mattmp(dmn+p,npel))
    end if
    do ob=1,nobs_loc
       ind=onlst(ob,:)
       if (dyn .and. .not. gf) then
          do i=1,npel
             row((/((i-1)*dmn+j,j=1,dmn)/))=(/((ind(i)-1)*dmn+j,j=1,dmn)/)
          end do
          mattmp=reshape(uu_dyn(row),(/dmn,npel/))
          vecshp=reshape(oshape(ob,:),(/npel,1/))
          vectmp=matmul(mattmp,vecshp)
          uu_dyn_obs(ob,:)=vectmp(:,1)
       else
          do i=1,npel
             row((/((i-1)*(dmn+p)+j,j=1,dmn+p)/))=(/((ind(i)-1)*(dmn+p)+j,j=1, &
                dmn+p)/)
          end do
          mattmp=reshape(uu(row),(/dmn+p,npel/))
          vecshp=reshape(oshape(ob,:),(/npel,1/))
          vectmp=matmul(mattmp,vecshp) 
          uu_obs(ob,:)=vectmp(:,1)
          if (poro) uu_obs(ob,dmn+1)=vectmp(dmn+1,1)*scale
       end if
    end do
  end subroutine GetVec_obs

  ! Apply nodal force
  subroutine ApplyNodalForce(node,vvec)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: node,i,j
    real(8) :: vvec(:)
    do i=1,dmn+p
       j=(dmn+p)*node-(dmn+p)+i-1
       val=vvec(i); if (i==dmn+1) val=scale*val
       if(dyn) then
          call VecSetValue(Vec_F_dyn,j,val,Add_Values,ierr)
       else
          call VecSetValue(Vec_F,j,val,Add_Values,ierr)
       end if
    end do
  end subroutine ApplyNodalForce

  ! Apply traction (EbEAve)
  subroutine ApplyTraction(el,side,vvec)
    implicit none
    integer :: el,side,i,snodes(nps)
    real(8) :: vvec(:),area
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    call EdgeAreaNodes(enodes,ecoords,side,area,snodes)
    vvec=vvec*area/dble(nps)
    snodes=nl2g(snodes,2)
    do i=1,nps
       call ApplyNodalForce(snodes(i),vvec)
    end do
  end subroutine ApplyTraction

  ! Form RHS
  subroutine FormRHS
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer :: i,j
    real(8) :: t1,t2
    do i=1,nceqs
       t1=cval(i,2)/dt; t2=cval(i,3)/dt
       if (tstep>=nint(t1) .and. tstep<=nint(t2)) then
          j=(dmn+p)*nnds+i-1
          if (stype/="explicit" .and. rank==0) then
             val=wt*cval(i,1)
             call VecSetValue(Vec_F,j,val,Add_Values,ierr)
          end if
       end if
    end do
    do i=1,nfrcs
       t1=fval(i,dmn+p+1)/dt; t2=fval(i,dmn+p+2)/dt
       if (tstep>=nint(t1) .and. tstep<=nint(t2)) then
          node=fnode(i); vvec=fval(i,1:dmn+p)
          if (rank==0) call ApplyNodalForce(node,vvec)
       end if
    end do
    do i=1,ntrcs
       t1=tval(i,dmn+p+1)/dt; t2=tval(i,dmn+p+2)/dt
       if (tstep>=nint(t1) .and. tstep<=nint(t2)) then
          el=telsd(i,1); side=telsd(i,2); vvec=tval(i,1:dmn+p)
          if (el/=0) call ApplyTraction(el,side,vvec)
       end if
    end do
  end subroutine FormRHS

  ! Form local damping matrix for elements with viscous dampers
  subroutine FormLocalAbsC(el,side,dir,m,indx)
    implicit none
    integer :: el,side,dir,indx(:)
    real(8) :: m(:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    dns=mat(id(el),5)
    call FormElAbsC(enodes,ecoords,side,dir,E,nu,dns,m)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalAbsC

  ! Form local damping matrix for elements with viscous dampers
  subroutine FormLocalAbsC1(el,side,m,indx)
    implicit none
    integer :: el,side,indx(:)
    real(8) :: m(:,:),matabs(dmn,dmn),vec1(dmn),vec2(dmn),vec3(dmn)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    dns=mat(id(el),5)
    select case (eltype) 
    case("tri")
       select case(side) 
       case(1); vec1=ecoords(2,:)-ecoords(1,:)
       case(2); vec1=ecoords(3,:)-ecoords(2,:)
       case(3); vec1=ecoords(1,:)-ecoords(3,:)
       end select
       vec1=vec1/sqrt(sum(vec1*vec1))
       vec2(1)=vec1(2); vec2(2)=-vec1(1)
       matabs(:,1)=vec1; matabs(:,2)=vec2
    case("qua")
       select case(side) 
       case(1); vec1=ecoords(2,:)-ecoords(1,:)
       case(2); vec1=ecoords(3,:)-ecoords(2,:)
       case(3); vec1=ecoords(4,:)-ecoords(3,:)
       case(4); vec1=ecoords(1,:)-ecoords(4,:)
       end select
       vec1=vec1/sqrt(sum(vec1*vec1))
       vec2(1)=vec1(2); vec2(2)=-vec1(1)
       matabs(:,1)=vec1; matabs(:,2)=vec2
    case("tet")
       select case(side) 
       case(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          call Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
       case(2)
          vec1=ecoords(4,:)-ecoords(3,:)
          call Cross(vec1,ecoords(2,:)-ecoords(3,:),vec3)
       case(3)
          vec1=ecoords(4,:)-ecoords(1,:)
          call Cross(vec1,ecoords(3,:)-ecoords(1,:),vec3)
       case(4)
          vec1=ecoords(3,:)-ecoords(1,:)
          call Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
       end select
       vec1=vec1/sqrt(sum(vec1*vec1))
       vec3=vec3/sqrt(sum(vec3*vec3))
       call Cross(vec1,vec3,vec2)
       matabs(:,1)=vec1; matabs(:,2)=vec2; matabs(:,3)=vec3
    case("hex")
       select case(side) 
       case(1)
          vec1=ecoords(2,:)-ecoords(1,:)
          call Cross(vec1,ecoords(5,:)-ecoords(1,:),vec3)
       case(2)
          vec1=ecoords(3,:)-ecoords(2,:)
          call Cross(vec1,ecoords(6,:)-ecoords(2,:),vec3)
       case(3)
          vec1=ecoords(4,:)-ecoords(3,:)
          call Cross(vec1,ecoords(7,:)-ecoords(3,:),vec3)
       case(4)
          vec1=ecoords(5,:)-ecoords(1,:)
          call Cross(vec1,ecoords(4,:)-ecoords(1,:),vec3)
       case(5)
          vec1=ecoords(4,:)-ecoords(1,:)
          call Cross(vec1,ecoords(2,:)-ecoords(1,:),vec3)
       case(6)
          vec1=ecoords(6,:)-ecoords(5,:)
          call Cross(vec1,ecoords(8,:)-ecoords(5,:),vec3)
       end select
       vec1=vec1/sqrt(sum(vec1*vec1))
       vec3=vec3/sqrt(sum(vec3*vec3))
       call Cross(vec1,vec3,vec2)
       matabs(:,1)=vec1; matabs(:,2)=vec2; matabs(:,3)=vec3
    end select
    call FormElAbsC1(enodes,ecoords,side,matabs,E,nu,dns,m)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalAbsC1

  subroutine Cross(a,b,r)
    implicit none
    real(8) :: a(3),b(3),r(3)
    r(1)=a(2)*b(3)-a(3)*b(2)
    r(2)=a(3)*b(1)-a(1)*b(3)
    r(3)=a(1)*b(2)-a(2)*b(1)
  end subroutine Cross

  ! Form local index
  subroutine FormLocalIndx(enodes,indx)
    implicit none
    integer :: enodes(:),indx(:)
    if (.not. poro) then
       call FormElIndx(enodes,indx)
    else
       call FormElIndxp(enodes,indx)
    end if
  end subroutine FormLocalIndx

  subroutine FormLocalIndx_dyn(enodes,indx)
    implicit none
    integer :: enodes(:),indx(:)
    call FormElIndx(enodes,indx)
  end subroutine FormLocalIndx_dyn

  ! Recover stress
  subroutine RecoverStress(el,stress)
    implicit none
    integer :: el,indxl(eldof)
    real(8) :: stress(:,:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    call FormLocalIndx(enodes,indx)
    indxl=indx(1:eldof)
    call CalcElStress(ecoords,uu(indxl),E,nu,stress(el,:,:))
  end subroutine RecoverStress

  ! Scatter stress to vertices (normal: SS, shear: SH)
  subroutine GetVec_Stress
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    integer:: i,j,indx(eldof),row(dmn)
    real(8):: sigma(cdmn),stvec(cdmn),cst(dmn),ecoords(npel,dmn),detj
    integer,save :: k=0
    if (k==0) then
       call VecDuplicate(Vec_SS,Vec_Cst,ierr)
       call VecZeroEntries(Vec_Cst,ierr)
    end if
    do i=1,nels
       enodes=nodes(i,:)
       call FormElIndx(enodes,indx)
       indx=indxmap_u(indx,2)
       ecoords=coords(enodes,:)
       stvec=f0; cst=f0
       do j=1,nip
          sigma=stress(i,j,:)
          call FormdetJ(ipoint(j,:),ecoords,detj)
          stvec=stvec+sigma*weight(j)*detj
          cst=cst+weight(j)*detj
       end do
       do j=1,npel
          select case(dmn)
          case(2)
             row=(/indx((j-1)*dmn+1),indx((j-1)*dmn+2)/)
             call VecSetValues(Vec_SS,dmn,row,stvec(1:dmn),Add_Values,ierr)
             call VecSetValues(Vec_SH,dmn,row,(/stvec(dmn+1:cdmn),             &
                stvec(dmn+1:cdmn)/),Add_Values,ierr)
          case(3)
             row=(/indx((j-1)*dmn+1),indx((j-1)*dmn+2),indx((j-1)*dmn+3)/)
             call VecSetValues(Vec_SS,dmn,row,stvec(1:dmn),Add_Values,ierr)
             call VecSetValues(Vec_SH,dmn,row,stvec(dmn+1:cdmn),Add_Values,ierr)
          end select
          if (k==0) call VecSetValues(Vec_Cst,dmn,row,cst,Add_Values,ierr)
       end do
    end do
    call VecAssemblyBegin(Vec_SS,ierr)
    call VecAssemblyEnd(Vec_SS,ierr)
    call VecAssemblyBegin(Vec_SH,ierr)
    call VecAssemblyEnd(Vec_SH,ierr)
    if (k==0) then
       call VecAssemblyBegin(Vec_Cst,ierr)
       call VecAssemblyEnd(Vec_Cst,ierr)
    end if
    call VecPointwiseDivide(Vec_SS,Vec_SS,Vec_Cst,ierr)
    call VecPointwiseDivide(Vec_SH,Vec_SH,Vec_Cst,ierr)
    k=k+1
  end subroutine GetVec_Stress

  ! Form local Hs
  subroutine FormLocalHs(el,Hs,indxp)
    implicit none
    real(8) :: Hs(:)
    integer :: el,j,indxp(:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    call FormElIndxp(enodes,indx)
    indxp=indx(eldof+1:)
    E=mat(id(el),1); nu=mat(id(el),2)
    call FormElHs(ecoords,uu(indxp),E,nu,scale,Hs)
    do j=1,npel
       ! Fixed or synced to FV 
       if (bc(nodes(el,j),dmn+1)==0 .or. bc(nodes(el,j),dmn+1)==2) Hs(j)=f0
    end do
  end subroutine FormLocalHs

  ! Signed distance point to plane/line
  subroutine Mix3D(a,b,c,m)
    implicit none
    real(8) :: a(3),b(3),c(3),r(3),m 
    r(1)=a(2)*b(3)-a(3)*b(2)
    r(2)=a(3)*b(1)-a(1)*b(3)
    r(3)=a(1)*b(2)-a(2)*b(1)
    m=(r(1)*c(1)+r(2)*c(2)+r(3)*c(3))/(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
  end subroutine Mix3D 
  subroutine Mix2D(a,b,m)
    implicit none
    real(8) :: a(2),b(2),a3(3),b3(3),c3(3),m 
    a3=(/f0,f0,f1/); b3=(/a(:),f0/); c3=(/b(:),f0/)
    call Mix3D(a3,b3,c3,m) 
  end subroutine Mix2D

  ! Observation/FD nodal base 
  subroutine GetObsNd(strng)
    implicit none
    character(2) :: strng 
    integer :: neval,ob,el
    integer,allocatable :: nd_full(:,:),pick(:)
    real(8) :: xmin,xmax,ymin,ymax,zmin,zmax,xmind,xmaxd,ymind,ymaxd,zmind,    &
       zmaxd,dd,du,dl,dr,df,db,d,eta,nu,psi,xob(dmn),N(npel),c,vec12(dmn),     &
       vec13(dmn),vec14(dmn),vec23(dmn),vec24(dmn),vec34(dmn),vec1o(dmn),      &
       vec2o(dmn),vec3o(dmn),vec15(dmn),vec73(dmn),vec76(dmn),vec78(dmn),      &
       vec7o(dmn)
    real(8),allocatable :: N_full(:,:)
    logical :: p_in_dom, p_in_el
    c=0.125d0

    ! Type of the evaluation Obs or FD
    if (strng=="ob") neval=nobs
    if (strng=="fd") neval=ngp
    allocate(pick(neval))
    pick=0

    select case(eltype) 
    case("tri"); allocate(nd_full(neval,3),N_full(neval,3))
    case("qua"); allocate(nd_full(neval,4),N_full(neval,4))
    case("tet"); allocate(nd_full(neval,4),N_full(neval,4))
    case("hex"); allocate(nd_full(neval,8),N_full(neval,8))
    end select
    xmind=minval(coords(:,1)); xmaxd=maxval(coords(:,1))
    ymind=minval(coords(:,2)); ymaxd=maxval(coords(:,2))
    if (dmn>2) then
       zmind=minval(coords(:,3)); zmaxd=maxval(coords(:,3))
    end if

    do ob=1,neval ! Observation loop
       if (strng=="ob") then
          xob=ocoord(ob,:)*km2m
       elseif (strng=="fd") then
          xob=xgp(ob,:)
       end if
       p_in_dom=(xob(1)>=xmind .and. xob(1)<=xmaxd .and. xob(2)>=ymind .and.   &
          xob(2)<=ymaxd)
       if (dmn>2) p_in_dom=(p_in_dom .and. xob(3)>=zmind .and. xob(3)<=zmaxd)
       if (p_in_dom) then ! Point probably in domain
          do el=1,nels
             p_in_el=.false.
             enodes=nodes(el,:)
             ecoords=coords(enodes,:)
             select case(dmn)
             case(2)
                xmin=minval(ecoords(:,1)); xmax=maxval(ecoords(:,1))
                ymin=minval(ecoords(:,2)); ymax=maxval(ecoords(:,2))
                ! Point probably in 2D cell 
                if (xob(1)>=xmin .and. xob(1)<=xmax .and. xob(2)>=ymin .and.   &
                   xob(2)<=ymax) then
                   select case(eltype)
                   case("tri")
                      vec12=ecoords(3,:)-ecoords(1,:)
                      vec13=ecoords(3,:)-ecoords(1,:)
                      vec23=ecoords(3,:)-ecoords(2,:)
                      vec1o=xob-ecoords(1,:)
                      vec2o=xob-ecoords(2,:)
                      vec3o=xob-ecoords(3,:)
                      call Mix2D(vec12,vec1o,dd)
                      call Mix2D(-vec13,vec3o,dl)
                      call Mix2D(vec23,vec2o,d)
                      ! Point in tri
                      if (dd>=f0 .and. dl>=f0 .and. d>=f0) then
                         call Mix2D(-vec13,vec12,dr)
                         call Mix2D(vec12,vec13,du)
                         eta=dl/dr; nu=dd/du
                         N(1)=f1-eta-nu; N(2)=eta; N(3)=nu
                         p_in_el=.true.
                      end if
                   case("qua")
                     vec12=ecoords(2,:)-ecoords(1,:)
                     vec14=ecoords(4,:)-ecoords(1,:)
                     vec23=ecoords(3,:)-ecoords(2,:)
                     vec34=ecoords(4,:)-ecoords(3,:)
                     vec1o=xob-ecoords(1,:)
                     vec3o=xob-ecoords(3,:)
                     call Mix2D(vec12,vec1o,dd)
                     call Mix2D(-vec14,vec1o,dl)
                     call Mix2D(vec23,vec3o,dr)
                     call Mix2D(vec34,vec3o,du)
                     ! Point in quad
                     if (dd>=f0 .and. dl>=f0 .and. dr>=f0 .and. du>=f0) then
                        eta=(dl-dr)/(dl+dr); nu=(dd-du)/(du+dd)
                        N(1)=0.25d0*(f1-eta)*(f1-nu)
                        N(2)=0.25d0*(f1+eta)*(f1-nu)
                        N(3)=0.25d0*(f1+eta)*(f1+nu)
                        N(4)=0.25d0*(f1-eta)*(f1+nu)
                        p_in_el=.true.
                     end if
                   end select ! Tri/qua
                end if ! Point probably in 2D cell
             case(3) 
                xmin=minval(ecoords(:,1)); xmax=maxval(ecoords(:,1))
                ymin=minval(ecoords(:,2)); ymax=maxval(ecoords(:,2))
                zmin=minval(ecoords(:,3)); zmax=maxval(ecoords(:,3))
                ! Point probably in 3D cell
                if (xob(1)>=xmin .and. xob(1)<=xmax .and. xob(2)>=ymin .and.   &
                   xob(2)<=ymax .and. xob(3)>=zmin .and. xob(3)<=zmax) then
                   select case(eltype)
                   case("tet")
                      vec12=ecoords(2,:)-ecoords(1,:)
                      vec13=ecoords(3,:)-ecoords(1,:)
                      vec14=ecoords(4,:)-ecoords(1,:)
                      vec23=ecoords(3,:)-ecoords(2,:)
                      vec24=ecoords(4,:)-ecoords(2,:)
                      vec1o=xob-ecoords(1,:)
                      vec2o=xob-ecoords(2,:)
                      call Mix3D(vec12,vec13,vec1o,dd)
                      call Mix3D(vec13,vec14,vec1o,dl)
                      call Mix3D(vec14,vec12,vec1o,df)
                      call Mix3D(vec24,vec23,vec2o,d)
                      ! Point in tet
                      if (dd>=f0 .and. dl>=f0 .and. df>=f0 .and. d>=f0) then 
                         call Mix3D(vec12,vec13,vec14,du)
                         call Mix3D(vec13,vec14,vec12,dr)
                         call Mix3D(vec14,vec12,vec13,db)
                         eta=dl/dr; nu=df/db; psi=dd/du
                         N(1)=f1-eta-nu-psi; N(2)=eta; N(3)=nu; N(4)=psi
                         p_in_el=.true.
                      end if
                   case("hex")
                      vec12=ecoords(2,:)-ecoords(1,:)
                      vec14=ecoords(4,:)-ecoords(1,:)
                      vec15=ecoords(5,:)-ecoords(1,:)
                      vec73=ecoords(3,:)-ecoords(7,:)
                      vec76=ecoords(6,:)-ecoords(7,:)
                      vec78=ecoords(8,:)-ecoords(7,:)
                      vec1o=xob-ecoords(1,:)
                      vec7o=xob-ecoords(7,:)
                      call Mix3D(vec12,vec14,vec1o,dd)
                      call Mix3D(vec15,vec12,vec1o,df)
                      call Mix3D(vec14,vec15,vec1o,dl)
                      call Mix3D(vec76,vec78,vec7o,du)
                      call Mix3D(vec78,vec73,vec7o,db)
                      call Mix3D(vec73,vec76,vec7o,dr)
                      ! Point in hex
                      if (dd>=f0 .and. dl>=f0 .and. df>=f0 .and. du>=f0 .and.  &
                         dr>=f0 .and. db>=f0) then
                         eta=(dl-dr)/(dr+dl); nu=(df-db)/(df+db)
                         psi=(dd-du)/(dd+du)
                         N(1)=c*(f1-eta)*(f1-nu)*(f1-psi)
                         N(2)=c*(f1+eta)*(f1-nu)*(f1-psi)
                         N(3)=c*(f1+eta)*(f1+nu)*(f1-psi)
                         N(4)=c*(f1-eta)*(f1+nu)*(f1-psi)
                         N(5)=c*(f1-eta)*(f1-nu)*(f1+psi)
                         N(6)=c*(f1+eta)*(f1-nu)*(f1+psi)
                         N(7)=c*(f1+eta)*(f1+nu)*(f1+psi)
                         N(8)=c*(f1-eta)*(f1+nu)*(f1+psi)
                         p_in_el=.true.
                      end if
                   end select
                end if ! Point probably in 3D cell
             end select ! Dimension
             ! Update local pick list
             if (p_in_el) then
                pick(ob)=ob
                N_full(ob,:)=N
                nd_full(ob,:)=enodes
                exit
             end if
          end do ! Element loop
       end if ! Point probably in domain
    end do ! Observation loop
    if (strng=="ob") then
       nobs_loc=size(pack(pick,pick/=0))
       allocate(ol2g(nobs_loc),ocoord_loc(nobs_loc,dmn))
       select case(eltype) 
       case("tri"); allocate(onlst(nobs_loc,3),oshape(nobs_loc,3))
       case("qua"); allocate(onlst(nobs_loc,4),oshape(nobs_loc,4))
       case("tet"); allocate(onlst(nobs_loc,4),oshape(nobs_loc,4))
       case("hex"); allocate(onlst(nobs_loc,8),oshape(nobs_loc,8))
       end select
       ol2g=pack(pick,pick/=0)
       ocoord_loc=ocoord(ol2g,:)
       onlst=nd_full(ol2g,:)
       oshape=N_full(ol2g,:)
    elseif (strng=="fd") then
       ngp_loc=size(pack(pick,pick/=0))
       allocate(gpl2g(ngp_loc),idgp_loc(ngp_loc,dmn))
       select case(eltype) 
       case("tri"); allocate(gpnlst(ngp_loc,3),gpshape(ngp_loc,3))
       case("qua"); allocate(gpnlst(ngp_loc,4),gpshape(ngp_loc,4))
       case("tet"); allocate(gpnlst(ngp_loc,4),gpshape(ngp_loc,4))
       case("hex"); allocate(gpnlst(ngp_loc,8),gpshape(ngp_loc,8))
       end select
       gpl2g=pack(pick,pick/=0)
       gpnlst=nd_full(gpl2g,:)
       gpshape=N_full(gpl2g,:)
       idgp_loc=idgp(gpl2g,:)
    end if
  end subroutine GetObsNd

  ! Active fault node map loc->glb 
  subroutine GetFltMap
    implicit none
    integer :: j,j3,map(n_lmnd,2)
    integer, allocatable :: rw(:)
    map=0
    do j=1,n_lmnd
       if (lmnd0+j>nceqs_ncf/(dmn+p)) then ! Has on rank fault nodes
          j3=lmnd0+j-nceqs_ncf/(dmn+p)
          if (frc(j3)>0) then
             map(j,:)=(/j,j3/)
          end if
       end if
    end do
    nfnd_loc=size(pack(map(:,1),map(:,1)/=0))
    allocate(rw(nfnd_loc),FltMap(nfnd_loc,2))
    rw=pack(map(:,1),map(:,1)/=0)
    FltMap=map(rw,:)
  end subroutine GetFltMap

  ! Apply body force
  subroutine ApplyGravity
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    real(8) :: dns,gip,gsca,gvec(eldof),ecoords(npel,dmn),detj
    integer :: el,i,indx(eldof+eldofp),row(eldof)
    do el=1,nels
       enodes=nodes(el,:)
       ecoords=coords(enodes,:)
       call FormLocalIndx(enodes,indx)
       row=indxmap(indx(:eldof),2)
       dns=mat(id(el),5)
       gsca=dns*gravity/dble(npel)
       gip=f0; gvec=f0
       do i=1,nip
          call FormdetJ(ipoint(i,:),ecoords,detj)
          gip=gip+gsca*weight(i)*detj
       end do
       ! Assume last dim aligned with gravity
       do i=1,npel
          if (bc(nodes(el,i),dmn)/=0) gvec(i*dmn)=-gip
       end do
       call VecSetValues(Vec_F,dmn*npel,row,gvec,Add_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
  end subroutine ApplyGravity

  ! Apply fluid body source
  subroutine ApplySource
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    real(8) :: sdns,sip,ssca,svec(eldofp),ecoords(npel,dmn),detj
    integer :: el,i,indx(eldof+eldofp),row(eldofp)
    do el=1,nels
       enodes=nodes(el,:)
       ecoords=coords(enodes,:)
       call FormElIndxp(enodes,indx)
       row=indxmap(indx(eldof+1:),2)
       sdns=mat(id(el),10)
       ssca=sdns/dble(npel)
       sip=f0; svec=f0
       do i=1,nip
          call FormdetJ(ipoint(i,:),ecoords,detj)
          sip=sip+ssca*weight(i)*detj
       end do
       do i=1,npel
          if (bc(nodes(el,i),dmn+1)/=0) svec(i)=sip
       end do
       call VecSetValues(Vec_F,npel,row,svec*scale,Add_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
  end subroutine ApplySource

  ! Reform RHS
  subroutine ReformLocalRHS(el,f,indx)
    implicit none
    integer :: el,indx(:),j,j1
    real(8) :: f(:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    f=f0
    call ReformElRHS(ecoords,stress(el,:,:),E,nu,visc,expn,dt,f(1:eldof))
    call FormLocalIndx(enodes,indx)
    ! Fix BCs (i.e., zero out entries)
    do j=1,npel
       do j1=1,dmn
          if (bc(nodes(el,j),j1)==0) f(dmn*j-dmn+j1)=f0
       end do
    end do
  end subroutine ReformLocalRHS

  ! Recover Vstress during implicit time stepping
  subroutine RecoverVStress(el,stress)
    implicit none
    integer :: el,indxl(eldof)
    real(8) :: stress(:,:,:)
    enodes=nodes(el,:)
    ecoords=coords(enodes,:)
    E=mat(id(el),1); nu=mat(id(el),2)
    visc=mat(id(el),3); expn=mat(id(el),4)
    call FormLocalIndx(enodes,indx)
    indxl=indx(1:eldof)
    call CalcElVStress(ecoords,uu(indxl),stress(el,:,:),E,nu,visc,expn,dt)
  end subroutine RecoverVStress

  ! Account for Winkler foundation(s)
  subroutine AddWinklerFdn(el,k)
    implicit none
    integer :: el,j,j1,j2,elbc(npel,dmn),snodes(nps)
    real(8) :: k(:,:),area
    dns=mat(id(el),5)
    elbc=bc(enodes,1:dmn)
    do j1=1,dmn
       call GetWinklerEdge(elbc,j1,side)
       if (side==0) cycle
       call EdgeAreaNodes(enodes,ecoords,side,area,snodes)
       val=dns*gravity*area/dble(nps)
       do j=1,npel
          if (elbc(j,j1)==-1) then
             j2=dmn*j-dmn+j1; k(j2,j2)=k(j2,j2)+val
          end if
       end do
    end do
  end subroutine AddWinklerFdn

  ! Fix BCs (i.e., zero rows/columns) in local [K/Kp]
  subroutine FixBCinLocalK(el,k)
    implicit none
    integer :: el,j,j1,j2
    real(8) :: k(:,:)
    do j=1,npel
       do j1=1,dmn
          if (bc(nodes(el,j),j1)==0) then
             j2=dmn*j-dmn+j1
             val=k(j2,j2)
             k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
             k(j2,j2)=val
          end if
       end do
       if (poro) then
          if (bc(nodes(el,j),dmn+1)==0) then
             j2=dmn*npel+j
             val=k(j2,j2)
             k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
             k(j2,j2)=val
          end if
       end if
    end do
  end subroutine FixBCinLocalK

  ! Fixed pressure at local [k]
  subroutine FixFVLocalK(el,k)
    implicit none
    integer :: el,j,j2
    real(8) :: k(:,:),tmp
    do j=1,npel
       j2=dmn*npel+j
       tmp=k(j2,j2)
       k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
       k(j2,j2)=tmp
    end do
  end subroutine FixFVLocalK
 
  subroutine FixFVLocalKBD(el,k) ! Only for boundary condition
    implicit none
    integer :: el,j,j2
    real(8) :: k(:,:),tmp
    do j=1,npel
       if (bc(nodes(el,j),dmn+1)==2) then
          j2=dmn*npel+j
          tmp=k(j2,j2)
          k(j2,:)=f0; k(:,j2)=f0 ! Zero out rows and columns
          k(j2,j2)=tmp
       end if
    end do
  end subroutine FixFVLocalKBD

  ! Print message
  subroutine PrintMsg(msg)
    implicit none
    character(*) :: msg
    if (rank==0) print*,msg
  end subroutine PrintMsg

  ! Write results in ASCII VTK (legacy) format
  subroutine WriteOutput
    implicit none
    character(64) :: fmt
    character(256) :: name,name0,name1
    integer,save :: k=0
    integer :: i,j,j1,lnnds,lnels
    real(8),pointer :: field_val(:)
    if (dsp==0) field_val=>uu
    if (dsp==1) field_val=>tot_uu
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,I0,3A,I0.6,A)')trim(name0),rank,"_",trim(name1),"_",k,".vtk"
    if (write_dyn) then
       if (dsp_hyb==1) then
          field_val=>tot_uu_dyn
       else
          field_val=>uu_dyn
       end if
       write(name,'(A,I0,3A,I0.6,A)')trim(name0),rank,"_",trim(name1),"_dyn_", &
          k,".vtk"
    end if
    open(10,file=adjustl(name),status='replace')
    lnnds=size(coords,1)
    lnels=size(nodes,1)
    write(10,'(A)')"# vtk DataFile Version 2.0"
    write(10,'(A)')"File written by Defmod-dev"
    write(10,'(A)')"ASCII"
    write(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    write(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m),f0/)
       end do
    case(3)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m)/)
       end do
    end select
    write(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    select case(npel)
    case(3); fmt="(I0,3(1X,I0))"
    case(4); fmt="(I0,4(1X,I0))"
    case(8); fmt="(I0,8(1X,I0))"
    end select
    do i=1,lnels
       write(10,fmt)npel,nodes(i,:)-1
    end do
    write(10,'(A,I0)')"CELL_TYPES ",lnels
    do i=1,lnels
       write(10,'(I0)')vtkid
    end do
    write(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn+p
    if (write_dyn) j=dmn
    if (poro .and. .not. write_dyn) then
       write(10,'(A)')"SCALARS pressure double"
       write(10,'(A)')"LOOKUP_TABLE default"
       fmt="(1(F0.3,1X))"
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(scale*field_val(j1)) ! Pr
       end do
    end if
    write(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          j1=i*j-p
          if (write_dyn) j1=i*j
          write(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       end do
    case(3)
       do i=1,lnnds
          j1=i*j-p
          if (write_dyn) j1=i*j
          write(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       end do
    end select
    close(10); k=k+1
  end subroutine WriteOutput

  subroutine WriteOutput_x
    implicit none
    character(64) :: fmt
    character(256) :: name,name0,name1
    integer,save :: k=0
    integer :: i,j,j1,j2,j3,lnnds,lnels
    real(8),pointer :: field_val(:),field_val_fp(:),field_val_fl(:),           &
       field_val_qu(:),field_val_ql(:),field_val_flc(:),field_val_ss(:),       &
       field_val_sh(:)
    if (dsp==0) field_val=>uu
    if (dsp==1) field_val=>tot_uu
    if (nceqs>0) field_val_fl=>fl; field_val_ql=>ql; field_val_flc=>flc
    if (poro) field_val_fp=>fp; field_val_qu=>qu
    if (fault .and. lm_str==1) field_val_ss=>ss; field_val_sh=>sh
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,I0,3A,I0.6,A)')trim(name0),rank,"_",trim(name1),"_",k,".vtk"
    open(10,file=adjustl(name),status='replace')
    lnnds=size(coords,1)
    lnels=size(nodes,1)
    write(10,'(A)')"# vtk DataFile Version 2.0"
    write(10,'(A)')"File written by Defmod-dev"
    write(10,'(A)')"ASCII"
    write(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    write(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m),f0/)
       end do
    case(3)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m)/)
       end do
    end select
    write(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    select case(npel)
    case(3); fmt="(I0,3(1X,I0))"
    case(4); fmt="(I0,4(1X,I0))"
    case(8); fmt="(I0,8(1X,I0))"
    end select
    do i=1,lnels
       write(10,fmt)npel,nodes(i,:)-1
    end do
    write(10,'(A,I0)')"CELL_TYPES ",lnels
    do i=1,lnels
       write(10,'(I0)')vtkid
    end do
    write(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn+p; j2=dmn
    if (poro) then
       write(10,'(A)')"SCALARS pressure double"
       write(10,'(A)')"LOOKUP_TABLE default"
       fmt="(1(F0.3,1X))"
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(scale*field_val(j1)) ! Pr
       end do
       write(10,'(A)')"SCALARS qu double"
       write(10,'(A)')"LOOKUP_TABLE default"
       fmt="(1(F0.3,1X))"
       do i=1,lnnds
          write(10,fmt)(field_val_qu(i)) ! qu
       end do
       write(10,'(A)')"VECTORS fp double"
       fmt="(3(F0.6,1X))"
       select case(dmn)
       case(2)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_fp(j1-1),field_val_fp(j1),f0/) ! 2D U
          end do
       case(3)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_fp(j1-2),field_val_fp(j1-1),             &
                field_val_fp(j1)/) ! 3D U
          end do
       end select
       if (nceqs > 0) then
          write(10,'(A)')"SCALARS ql double"
          write(10,'(A)')"LOOKUP_TABLE default"
          fmt="(1(F0.3,1X))"
          do i=1,lnnds
             write(10,fmt)(field_val_ql(i)) ! ql
          end do
       end if
    end if
    if (fault .and. lm_str==1) then
       write(10,'(A)')"VECTORS ss double"
       fmt="(3(E10.4,1X))"
       select case(dmn)
       case(2)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_ss(j1-1),field_val_ss(j1),f0/) ! 2D U
          end do
       case(3)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_ss(j1-2),field_val_ss(j1-1),             &
                field_val_ss(j1)/) ! 3D U
          end do
       end select
       write(10,'(A)')"VECTORS sh double"
       fmt="(3(E10.4,1X))"
       select case(dmn)
       case(2)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_sh(j1-1),field_val_sh(j1),f0/) ! 2D U
          end do
       case(3)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_sh(j1-2),field_val_sh(j1-1),             &
                field_val_sh(j1)/) ! 3D U
          end do
       end select
    end if
    if (nceqs>0) then
       do j3=1,dmn
          if (j3==1) then
             write(10,'(A)')"SCALARS fls double"
             write(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             do i=1,lnnds
                j1=i*dmn-dmn+j3
                write(10,fmt)(field_val_flc(j1))
             end do
          elseif(j3==2) then
             if (dmn==2) then
                write(10,'(A)')"SCALARS fln double"
             else
                write(10,'(A)')"SCALARS fld double"
             end if
             write(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             do i=1,lnnds
                j1=i*dmn-dmn+j3
                write(10,fmt)(field_val_flc(j1))
             end do
          elseif(j3==3) then
             write(10,'(A)')"SCALARS fln double"
             write(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             do i=1,lnnds
                j1=i*dmn-dmn+j3
                write(10,fmt)(field_val_flc(j1))
             end do
          end if
       end do
       write(10,'(A)')"VECTORS fl double"
       fmt="(3(F0.6,1X))"
       select case(dmn)
       case(2)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_fl(j1-1),field_val_fl(j1),f0/) ! 2D U
          end do
       case(3)
          do i=1,lnnds
             j1=i*j2
             write(10,fmt)(/field_val_fl(j1-2),field_val_fl(j1-1),             &
                field_val_fl(j1)/) ! 3D U
          end do
       end select
    end if
    write(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          j1=i*j-p
          write(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       end do
    case(3)
       do i=1,lnnds
          j1=i*j-p
          write(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       end do
    end select
    close(10); k=k+1
  end subroutine WriteOutput_x

  ! Write out deformation and (absolute) Coulomb force due 
  !to pore pressure initialization
  subroutine WriteOutput_init
    implicit none
    character(64) :: fmt
    character(256) :: name,name0,name1
    integer :: i,j,j1,j2,lnnds,lnels
    real(8),pointer :: field_val(:),field_val_flc(:)
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,I0,A,A,A)')trim(name0),rank,"_",trim(name1),"_init.vtk"
    open(10,file=adjustl(name),status='replace')
    lnnds=size(coords,1)
    lnels=size(nodes,1)
    write(10,'(A)')"# vtk DataFile Version 2.0"
    write(10,'(A)')"File written by Defmod-dev"
    write(10,'(A)')"ASCII"
    write(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    write(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m),f0/)
       end do
    case(3)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m)/)
       end do
    end select
    write(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    select case(npel)
    case(3); fmt="(I0,3(1X,I0))"
    case(4); fmt="(I0,4(1X,I0))"
    case(8); fmt="(I0,8(1X,I0))"
    end select
    do i=1,lnels
       write(10,fmt)npel,nodes(i,:)-1
    end do
    write(10,'(A,I0)')"CELL_TYPES ",lnels
    do i=1,lnels
       write(10,'(I0)')vtkid
    end do
    write(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn+1
    field_val=>uu; field_val_flc=>flc
    write(10,'(A)')"SCALARS pressure double"
    write(10,'(A)')"LOOKUP_TABLE default"
    fmt="(1(F0.3,1X))"
    do i=1,lnnds
       j1=i*j
       write(10,fmt)(scale*field_val(j1)) ! Pr
    end do
    if (fault .and. nfnd>0) then    
       do j2=1,dmn
          if (j2==1) then
             write(10,'(A)')"SCALARS fls double"
             write(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             do i=1,lnnds
                j1=i*dmn-dmn+j2
                write(10,fmt)(field_val_flc(j1))
             end do
          elseif(j2==2) then
             if (dmn==2) then
                write(10,'(A)')"SCALARS fln double"
             else
                write(10,'(A)')"SCALARS fld double"
             end if
             write(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             do i=1,lnnds
                j1=i*dmn-dmn+j2
                write(10,fmt)(field_val_flc(j1))
             end do
          elseif(j2==3) then
             write(10,'(A)')"SCALARS fln double"
             write(10,'(A)')"LOOKUP_TABLE default"
             fmt="(1(F0.3,1X))"
             do i=1,lnnds
                j1=i*dmn-dmn+j2
                write(10,fmt)(field_val_flc(j1))
             end do
          end if
       end do
    end if
    write(10,'(A)')"VECTORS displacements double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          j1=i*j-p
          write(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
       end do
    case(3)
       do i=1,lnnds
          j1=i*j-p
          write(10,fmt)(/field_val(j1-2),field_val(j1-1),field_val(j1)/) ! 3D U
       end do
    end select
    close(10); k=k+1
  end subroutine WriteOutput_init

  ! Write out force to stress ratio on fault f2s
  subroutine WriteOutput_f2s
    implicit none
    character(64) :: fmt
    character(256) :: name,name0,name1
    integer :: i,j,j1,lnnds,lnels
    real(8),pointer :: field_val(:),field_val_dip(:),field_val_nrm(:)
    field_val=>f2s; field_val_dip=>dip; field_val_nrm=>nrm
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,I0,A,A,A)')trim(name0),rank,"_",trim(name1),"_f2s.vtk"
    open(10,file=adjustl(name),status='replace')
    lnnds=size(coords,1)
    lnels=size(nodes,1)
    write(10,'(A)')"# vtk DataFile Version 2.0"
    write(10,'(A)')"File written by Defmod-dev"
    write(10,'(A)')"ASCII"
    write(10,'(A)')"DATASET UNSTRUCTURED_GRID"
    write(10,'(A,I0,A)')"POINTS ",lnnds," double"
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m),f0/)
       end do
    case(3)
       do i=1,lnnds
          write(10,fmt)(/(coords(i,:)/km2m)/)
       end do
    end select
    write(10,'(A,I0,1X,I0)')"CELLS ",lnels,lnels*(npel+1)
    select case(npel)
    case(3); fmt="(I0,3(1X,I0))"
    case(4); fmt="(I0,4(1X,I0))"
    case(8); fmt="(I0,8(1X,I0))"
    end select
    do i=1,lnels
       write(10,fmt)npel,nodes(i,:)-1
    end do
    write(10,'(A,I0)')"CELL_TYPES ",lnels
    do i=1,lnels
       write(10,'(I0)')vtkid
    end do
    write(10,'(A,I0)')"POINT_DATA ",lnnds
    j=dmn
    ! Write force to stress ratio
    if (nceqs>0) then
       write(10,'(A)')"VECTORS f2s double"
       fmt="(3(F0.6,1X))"
       select case(dmn)
       case(2)
          do i=1,lnnds
             j1=i*j
             write(10,fmt)(/field_val(j1-1),field_val(j1),f0/) ! 2D U
          end do
       case(3)
          do i=1,lnnds
             j1=i*j
             write(10,fmt)(/field_val(j1-2),field_val(j1-1),                   &
                field_val(j1)/) ! 3D U
          end do
       end select
    end if
    ! Write fault's dip (strike for 2D) and normal vectors
    select case(dmn)
    case(2)
       write(10,'(A)')"VECTORS strk double"
    case(3)
       write(10,'(A)')"VECTORS dip double"
    end select
    fmt="(3(F0.6,1X))"
    select case(dmn)
    case(2)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val_dip(j1-1),field_val_dip(j1),f0/) ! 2D U
       end do
    case(3)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val_dip(j1-2),field_val_dip(j1-1),              &
             field_val_dip(j1)/) ! 3D U
       end do
    end select
    write(10,'(A)')"VECTORS nrm double"
    select case(dmn)
    case(2)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val_nrm(j1-1),field_val_nrm(j1),f0/) ! 2D U
       end do
    case(3)
       do i=1,lnnds
          j1=i*j
          write(10,fmt)(/field_val_nrm(j1-2),field_val_nrm(j1-1),              &
             field_val_nrm(j1)/) ! 3D U
       end do
    end select
    close(10)
  end subroutine WriteOutput_f2s

  ! Output observations
  subroutine WriteOutput_obs
    implicit none
    integer :: i,j
    character(64) :: fmt
    character(256) :: name,name0,name1
    integer,save :: k=0,k_dyn=0
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    if (dyn) then
       write(name,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_dyn_obs_",rank,   &
          ".txt"
       j=k_dyn
    else
       write(name,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_obs_",rank,".txt"
       j=k
    end if
    if (j==0) then
       open(10,file=adjustl(name),status='replace')
       select case(dmn)
       case(2); fmt="(2(F0.6,1X),I0)"
       case(3); fmt="(3(F0.6,1X),I0)"
       end select
       write(10,"(I0)")nobs_loc
       do i=1,nobs_loc
          write(10,fmt)ocoord_loc(i,:),ol2g(i)
       end do
    else
       open(10,file=adjustl(name),status='old',position='append',action='write')
    end if
    if (dyn) then
       select case(dmn)
          case(2); fmt="(2(F0.6,1X))"
          case(3); fmt="(3(F0.6,1X))"
       end select
    else
       select case(dmn)
       case(2)
          if (poro) then
             fmt="(3(F0.6,1X))"
          else
             fmt="(2(F0.6,1X))"
          end if
       case(3)
          if (poro) then
             fmt="(4(F0.6,1X))"
          else
             fmt="(3(F0.6,1X))"
          end if
       end select
    end if
    do i=1,nobs_loc
       if (dyn) then
          if (dsp_hyb==1) then
             write(10,fmt)tot_uu_dyn_obs(i,:)
          else
             write(10,fmt)uu_dyn_obs(i,:)
          end if
       else
          if (dsp==1) then
             write(10,fmt)tot_uu_obs(i,:)
          else
             write(10,fmt)uu_obs(i,:)
          end if
       end if
    end do
    close(10) 
    if (dyn) then
       k_dyn=k_dyn+1
    else
       k=k+1
    end if
  end subroutine WriteOutput_obs

  ! Write event log file for event info 
  subroutine WriteOutput_log
    implicit none
    character(256) :: name,name0,name1
    integer,save :: k=0
    integer :: seis
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),".log"
    if (k==0) then
       open(10,file=adjustl(name),status='replace')
       write(10,"(F0.6)")dt
    else
       open(10,file=adjustl(name),status='old',position='append',action='write')
    end if
    if (crp) then
       seis=0
    else
       seis=1
    end if
    write(10,"(2(I0,X))")n_log,seis
    close(10); k=k+1
  end subroutine WriteOutput_log

  ! Write event log file for seismic data
  subroutine WriteOutput_log_wave
    implicit none
    character(256) :: name,name0,name1
    integer,save :: k=0
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),"_dyn.log"
    if (k==0) then
       open(10,file=adjustl(name),status='replace')
       if (gf) then
          write(10,"(F0.6)")dt*frq_wave
       else
          write(10,"(F0.6)")dt_dyn*frq_wave
       end if
    else
       open(10,file=adjustl(name),status='old',position='append',action='write')
    end if
    write(10,"(I0)")n_log_wave
    close(10); k=k+1
  end subroutine WriteOutput_log_wave

  ! Write event log file for seismic data
  subroutine WriteOutput_log_slip
    implicit none
    character(256) :: name,name0,name1
    integer,save :: k=0
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),"_slip.log"
    if (k==0) then
       open(10,file=adjustl(name),status='replace')
       write(10,"(F0.6)")dt_dyn*frq_slip
    else
       open(10,file=adjustl(name),status='old',position='append',action='write')
    end if
    write(10,"(I0)")n_log_slip
    close(10); k=k+1
  end subroutine WriteOutput_log_slip

  ! Write temporal fault slip (and theta for rate state friction) 
  subroutine WriteOutput_slip
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    character(256) :: name,name0,name1
    character(64) :: fmt
    integer :: j,j1,j2,j3,rw_loc(dmn) 
    integer,save :: k=0
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    if (nfnd_loc>0) then ! Has on rank fault nodes
       write(name,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_slip_",rank,".txt"
       if (rsf<1) then
          select case(dmn)
             case(2); fmt="(2(ES11.2E3,X))"
             case(3); fmt="(3(ES11.2E3,X))"
          end select
       else
          select case(dmn)
             case(2); fmt="(2(ES11.2E3,X),E12.4)"
             case(3); fmt="(3(ES11.2E3,X),E12.4)"
          end select
       end if
       if (k==0) then
          open(10,file=adjustl(name),status='replace')
          write(10,'(I0)')nfnd_loc
          do j1=1,nfnd_loc
             j3=FltMap(j1,2) 
             select case(dmn)
                case(2); write(10,'(2(F0.6,X),I0)')xfnd(j3,:),j3
                case(3); write(10,'(3(F0.6,X),I0)')xfnd(j3,:),j3
             end select
          end do
       else
          open(10,file=adjustl(name),status='old',position='append',action=    &
             'write')
       end if
       if (rsf<1) then ! Slip weakening
          do j1=1,nfnd_loc
             j=FltMap(j1,1); j3=FltMap(j1,2)
             rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
             if (dsp_hyb==1) then
                write(10,fmt)tot_flt_slip(rw_loc(:dmn-1)),mu_hyb(j3)
             else
                write(10,fmt)flt_slip(rw_loc(:dmn-1))/dt_dyn,mu_hyb(j3)
             end if
          end do
       else ! Rate and state friction
          do j1=1,nfnd_loc
             j=FltMap(j1,1); j3=FltMap(j1,2)
             rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
             if (dsp_hyb==1) then 
                write(10,fmt)tot_flt_slip(rw_loc(:dmn-1)),mu_hyb(j3),          &
                   rsftheta(j3)
             else
                write(10,fmt)flt_slip(rw_loc(:dmn-1))/dt_dyn,mu_hyb(j3),       &
                   rsftheta(j3)
             end if
          end do
       end if ! Constitutive models
       close(10); k=k+1
    end if ! Has on rank fault nodes
  end subroutine WriteOutput_slip

  ! Write qs fault data 
!  subroutine WriteOutput_flt_qs
!    implicit none
!    character(256) :: name,name0,name1
!    character(64) :: fmt
!    integer :: i
!    integer,save :: k=0
!    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
!    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
!    write(name,'(A,A,A)')trim(name0),trim(name1),"_fqs.txt"
!    select case(dmn+p)
!       case(2); fmt="(2(F0.6,1X))"
!       case(3); fmt="(3(F0.6,1X))"
!       case(4); fmt="(4(F0.6,1X))"
!    end select
!    if (k==0) then
!       open(10,file=adjustl(name),status='replace')
!       write(10,'(I0)')sum(frc)
!       do i=1,nfnd
!          if (frc(i)>0) then
!             select case(dmn)
!                case(2); write(10,'(2(F0.6,1X))')xfnd(i,:)
!                case(3); write(10,'(3(F0.6,1X))')xfnd(i,:)
!             end select
!          end if
!       end do
!    else
!       open(10,file=adjustl(name),status='old',position='append',action='write')
!    end if
!    do i=1,nfnd
!       if (frc(i)>0) then
!          if (poro) then
!             write(10,fmt)flt_ss(i,:),flt_p(i)
!          else
!             write(10,fmt)flt_ss(i,:)
!          end if
!       end if
!    end do
!    close(10); k=k+1
!  end subroutine WriteOutput_flt_qs

end module global
