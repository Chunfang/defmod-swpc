! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module fvfe 

#include <petscversion.h>

  use global
  use HDF5
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
  implicit none
#else
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
#endif

  private
  save

  !! Public routines
  public :: MakeEl2g
  public :: FVReset
  public :: FVInit
  public :: FVInitUsg
  public :: FVReformKF
  public :: FVReformKFUsg
  public :: FVReformKPerm
  public :: FVReformKPermUsg
  public :: FVSyncBD
  public :: FVSyncBDUsg
 
  !! Local vars
  character(256) :: nameh5
  integer :: stp0
  integer(hsize_t) :: off_dom(3),dim_dom(3) ! H5 offset and dim
  real(8) :: xmin,ymin,zmin,dx,dy,dz,r_perm(3),p_top,dt_fv ! Local FV/FE vars 
  integer, allocatable :: el2g(:),bdnd(:) ! Boundary nodes
  !! Pressure and time
  real(8), allocatable :: p_fv0(:),p_fv_bd(:,:),t_fv(:),t_fe(:) 

contains

  ! Local to original element map
  subroutine MakeEl2g
    integer :: el 
    allocate(el2g(nels))
    do el=1,size(emap,1)  
       if(emap(el)/=0) el2g(emap(el))=el 
    end do
  end subroutine MakeEl2g

  subroutine FVInit ! Read <model>_fvfe.cfg and <model>.h5 files
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    character(256) :: name,name0,name1,namegrp,namedat
    integer :: i,row,iv,nbd,it,nt_fv,nx,ny,nz,idfile,err,idgrp,iddat,spc_dat,  &
       spc_dom,lnnds,ix,iy,iz
    real(8) :: xref,yref,zref,xmax,ymax,zmax,p_hst
    real(8),allocatable :: p_dom(:,:,:)
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),"_fvfe.cfg"
    open(153,file=adjustl(name),status='old')
    read(153,*)nt_fv,stp0,dt_fv,p_top
    read(153,*)nx,ny,nz
    read(153,*)dx,dy,dz
    dx=km2m*dx; dy=km2m*dy; dz=km2m*dz
    read(153,*)xref,yref,zref
    xref=km2m*xref; yref=km2m*yref; zref=km2m*zref
    read(153,*)r_perm ! Inverse viscosity with scales
    close(153)
    write(nameh5,'(A,A,A)')trim(name0),trim(name1),".h5"
    xmin=minval(coords(:,1)); xmax=maxval(coords(:,1))
    ymin=minval(coords(:,2)); ymax=maxval(coords(:,2))
    zmin=minval(coords(:,3)); zmax=maxval(coords(:,3))
    ! Reversed and zero based hdf5 offsets
    off_dom(3)=int((xmin-xref)/dx)
    off_dom(2)=int((ymin-yref)/dy)
    off_dom(1)=int((zmin-zref)/dz)
    ! One based
    dim_dom(3)=int((xmax-xmin)/dx)+1
    dim_dom(2)=int((ymax-ymin)/dy)+1
    dim_dom(1)=int((zmax-zmin)/dz)+1
    ! Ascertain domain containment
    if (off_dom(3)<0 .or. off_dom(2)<0 .or. off_dom(1)<0 .or.                  &
       off_dom(3)+dim_dom(3)>nx .or. off_dom(2)+dim_dom(2)>ny .or.             &
       off_dom(1)+dim_dom(1)>nz) then
       call PrintMsg("FV->FE dimension mismatch!")
    end if
    lnnds=size(coords,1)
    nbd=size(pack(bc(:,dmn+1),bc(:,dmn+1)==2))
    allocate(t_fv(nt_fv),t_fe(nt_fv),p_fv0(lnnds),bdnd(nbd),p_fv_bd(nt_fv,nbd))
    allocate(p_dom(dim_dom(1),dim_dom(2),dim_dom(3))) 
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    call h5screate_simple_f(dmn,dim_dom,spc_dom,err)
    ! Hydrostatic pressure from FV
    write(namegrp,'(A,ES12.5E2,A)')"/Time: ",f0," d" 
    write(namedat,'(A,A)')trim(namegrp),"/Liquid_Pressure [Pa]"
    call h5gopen_f(idfile,trim(namegrp),idgrp,err)
    call h5dopen_f(idfile,trim(namedat),iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off_dom,dim_dom,err)
    call h5dread_f(iddat,h5t_native_double,p_dom,dim_dom,err,spc_dom,spc_dat)
    call h5dclose_f(iddat,err) 
    call h5gclose_f(idgrp,err)
    ! Form hydrostatic pressure
    call VecGetSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecDuplicate(Vec_Up,Vec_Up_hst,ierr)
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecZeroEntries(Vec_Up_hst,ierr)
    do iv=1,lnnds
       ix=int((coords(iv,1)-xmin)/dx)+1 
       iy=int((coords(iv,2)-ymin)/dy)+1
       iz=int((coords(iv,3)-zmin)/dz)+1
       p_hst=p_dom(iz,iy,ix)
       row=nl2g(i,2)-1
       call VecSetValue(Vec_Up_hst,row,-(p_hst-p_top)/scale,Insert_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_Up_hst,ierr)
    call VecAssemblyEnd(Vec_Up_hst,ierr) 
    ! Sampling pressure from FV
    do it=1,nt_fv
       t_fv(it)=dble(stp0+it-1)*dt_fv
       write(namegrp,'(A,ES12.5E2,A)')"/Time: ",t_fv(it)," d" 
       write(namedat,'(A,A)')trim(namegrp),"/Liquid_Pressure [Pa]"
       call h5gopen_f(idfile,trim(namegrp),idgrp,err)
       call h5dopen_f(idfile,trim(namedat),iddat,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off_dom,dim_dom,err)
       call h5dread_f(iddat,h5t_native_double,p_dom,dim_dom,err,spc_dom,spc_dat)
       call h5dclose_f(iddat,err) 
       call h5gclose_f(idgrp,err)
       i=1
       do iv=1,lnnds
          ix=int((coords(iv,1)-xmin)/dx)+1 
          iy=int((coords(iv,2)-ymin)/dy)+1
          iz=int((coords(iv,3)-zmin)/dz)+1
          if (it==1) p_fv0(iv)=p_dom(iz,iy,ix)
          if (bc(iv,dmn+1)==2) then 
             p_fv_bd(it,i)=p_dom(iz,iy,ix)
             bdnd(i)=iv
             i=i+1
          end if
       end do 
    end do
    call h5fclose_f(idfile,err)
    call h5close_f(err)
    t_fe=t_fv*24.d0*3600.d0 ! Day to sec
    ! Zero based incremental
    t_fe=t_fe-t_fe(1)
    p_fv_bd(2:nt_fv,:)=p_fv_bd(2:nt_fv,:)-p_fv_bd(1:nt_fv-1,:)
    p_fv_bd(1,:)=f0
  end subroutine FVInit

  subroutine FVInitUsg
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    character(256) :: name,name0,name1,namegrp,namedat
    integer :: i,j,row,nt_fv,idfile,idgrp,iddat,err,spc_dat,spc_cell,nbd
    integer(hsize_t) :: off(1),dim_cell(1)
    real(8) :: p_cell,cp,detj
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),"_fvfe.cfg"
    open(153,file=adjustl(name),status='old')
    read(153,*)nt_fv,stp0,dt_fv,p_top
    read(153,*)r_perm ! Inverse viscosity with scales
    close(153)
    write(nameh5,'(A,A,A)')trim(name0),trim(name1),"_usg.h5"
    allocate(t_fv(nt_fv),t_fe(nt_fv))
    do i=1,nt_fv
       t_fv(i)=dble(stp0+i-1)*dt_fv
    end do
    t_fe=t_fv*24.d0*3600.d0 ! Day to sec
    ! Zero based incremental
    t_fe=t_fe-t_fe(1)
    ! Form hydrostatic pressure
    call VecGetSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecDuplicate(Vec_Up,Vec_Up_hst,ierr)
    call VecDuplicate(Vec_Up,Vec_Cp0,ierr)
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecZeroEntries(Vec_Up_hst,ierr)
    call VecZeroEntries(Vec_Cp0,ierr)
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    write(namegrp,'(A,I4,A,ES12.5E2,A)')"/",0," Time ",f0," d"
    write(namedat,'(A,A)')trim(namegrp),"/Liquid Pressure [Pa]"
    call h5gopen_f(idfile,trim(namegrp),idgrp,err)
    call h5dopen_f(idfile,trim(namedat),iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    dim_cell=1
    call h5screate_simple_f(1,dim_cell,spc_cell,err)
    do i=1,nels
       off=el2g(i)-1
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off,dim_cell,err)
       call h5dread_f(iddat,h5t_native_double,p_cell,dim_cell,err,spc_cell,    &
          spc_dat)
       ecoords=coords(nodes(i,:),:)
       cp=f0
       do j=1,nip
          call FormdetJ(ipoint(j,:),ecoords,detj)
          cp=cp+weight(j)*detj
       end do
       do j=1,npel
          row=nl2g(nodes(i,j),2)-1
          call VecSetValue(Vec_Cp0,row,cp,Add_Values,ierr)
          call VecSetValue(Vec_Up_hst,row,-cp*(p_cell-p_top)/scale,Add_Values, &
             ierr)
       end do
    end do
    call h5dclose_f(iddat,err)
    call h5gclose_f(idgrp,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err) 
    call VecAssemblyBegin(Vec_Cp0,ierr)
    call VecAssemblyEnd(Vec_Cp0,ierr)
    call VecAssemblyBegin(Vec_Up_hst,ierr)
    call VecAssemblyEnd(Vec_Up_hst,ierr)
    ! Take cell-wise average
    call VecPointwiseDivide(Vec_Up_hst,Vec_Up_hst,Vec_Cp0,ierr)
    nbd=size(pack(bc(:,dmn+1),bc(:,dmn+1)==2))
    allocate(bdnd(nbd))
    j=1
    do i=1,size(coords,1)
       if (bc(i,dmn+1)==2) then
          bdnd(j)=i 
          j=j+1
       end if
    end do
  end subroutine FVInitUsg

  ! Reform [Kp] and [F] to impose FE domain pressure
  subroutine FVReformKF(ef_eldof) 
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    integer :: i,j,j2,row,ef_eldof
    call MatZeroEntries(Mat_K,ierr)
    call VecGetSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecDuplicate(Vec_Up,Vec_Cp,ierr)
    call VecDuplicate(Vec_Up,Vec_Up_fv,ierr)
    call VecZeroEntries(Vec_Cp,ierr)
    call VecZeroEntries(Vec_Up_fv,ierr)
    do i=1,nels 
       call FormLocalK(i,k,indx,"Kp") 
       indx=indxmap(indx,2)
       call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)
       ! Record pressure coefficient and pressure
       do j=1,npel
          j2=dmn*npel+j      
          row=nl2g(nodes(i,j),2)-1
          call VecSetValue(Vec_Cp,row,k(j2,j2),Add_Values,ierr)
          call VecSetValue(Vec_Up,row,p_fv0(nodes(i,j))/scale,Insert_Values,   &
             ierr)
          call VecSetValue(Vec_Up_fv,row,p_fv0(nodes(i,j))/scale,Insert_Values,&
             ierr)
       end do
    end do
    ! Account for constraint eqn's
    if (rank==0 .and. nceqs>0) call ApplyConstraints
    call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
    call VecAssemblyBegin(Vec_Cp,ierr)
    call VecAssemblyEnd(Vec_Cp,ierr)
    call VecAssemblyBegin(Vec_Up_fv,ierr)
    call VecAssemblyEnd(Vec_Up_fv,ierr)
    call VecAssemblyBegin(Vec_Up,ierr)
    call VecAssemblyEnd(Vec_Up,ierr)
    ! Rescale RHS by pressure coefficient
    call VecPointwiseMult(Vec_up,Vec_Up,Vec_Cp,ierr)
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
    ! Add FV pressure contribution to nodal force
    call MatMult(Mat_H,Vec_Up_fv,Vec_fp,ierr)
    call VecGetSubVector(Vec_F,RIu,Vec_Uu,ierr)
    call VecAXPY(Vec_Uu,-f1,Vec_fp,ierr)
    call VecRestoreSubVector(Vec_F,RIu,Vec_Uu,ierr)
    deallocate(p_fv0)
  end subroutine FVReformKF

  subroutine FVReformKFUsg(ef_eldof) 
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    character(256) :: namegrp,namedat
    integer :: i,j,j2,row,ef_eldof,idfile,idgrp,iddat,spc_dat,spc_cell,err
    integer(hsize_t) :: off(1),dim_cell(1)
    real(8) :: p_cell,cp,detj
    call MatZeroEntries(Mat_K,ierr)
    call VecGetSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecDuplicate(Vec_Up,Vec_Cp,ierr)
    call VecDuplicate(Vec_Up,Vec_Up_fv,ierr)
    call VecZeroEntries(Vec_Cp,ierr)
    call VecZeroEntries(Vec_Up_fv,ierr)
    ! Initial pressure
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    write(namegrp,'(A,I4,A,ES12.5E2,A)')"/",stp0," Time ",t_fv(1)," d"
    write(namedat,'(A,A)')trim(namegrp),"/Liquid Pressure [Pa]"
    call h5gopen_f(idfile,trim(namegrp),idgrp,err)
    call h5dopen_f(idfile,trim(namedat),iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    dim_cell=1
    call h5screate_simple_f(1,dim_cell,spc_cell,err)
    do i=1,nels 
       off=el2g(i)-1
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off,dim_cell,err)
       call h5dread_f(iddat,h5t_native_double,p_cell,dim_cell,err,spc_cell,    &
          spc_dat)
       call FormLocalK(i,k,indx,"Kp") 
       indx=indxmap(indx,2)
       call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)
       ecoords=coords(nodes(i,:),:)
       cp=f0
       do j=1,nip
          call FormdetJ(ipoint(j,:),ecoords,detj)
          cp=cp+weight(j)*detj
       end do
       do j=1,npel
          j2=dmn*npel+j      
          row=nl2g(nodes(i,j),2)-1
          call VecSetValue(Vec_Cp,row,k(j2,j2),Add_Values,ierr)
          call VecSetValue(Vec_Up,row,cp*p_cell/scale,Add_Values,ierr)
          call VecSetValue(Vec_Up_fv,row,cp*p_cell/scale,Add_Values,ierr)
       end do
    end do
    call h5dclose_f(iddat,err) 
    call h5gclose_f(idgrp,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
    ! Account for constraint eqn's
    if (rank==0 .and. nceqs>0) call ApplyConstraints
    call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
    call VecAssemblyBegin(Vec_Cp,ierr)
    call VecAssemblyEnd(Vec_Cp,ierr)
    call VecAssemblyBegin(Vec_Up,ierr)
    call VecAssemblyEnd(Vec_Up,ierr)
    call VecAssemblyBegin(Vec_Up_fv,ierr)
    call VecAssemblyEnd(Vec_Up_fv,ierr)
    ! Average by weight
    call VecPointwiseDivide(Vec_Up,Vec_Up,Vec_Cp0,ierr)
    call VecPointwiseDivide(Vec_Up_fv,Vec_Up_fv,Vec_Cp0,ierr)
    ! Rescale RHS by pressure coefficient
    call VecPointwiseMult(Vec_up,Vec_Up,Vec_Cp,ierr)
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
    ! Add FV pressure contribution to nodal force
    call MatMult(Mat_H,Vec_Up_fv,Vec_fp,ierr)
    call VecGetSubVector(Vec_F,RIu,Vec_Uu,ierr)
    call VecAXPY(Vec_Uu,-f1,Vec_fp,ierr)
    call VecRestoreSubVector(Vec_F,RIu,Vec_Uu,ierr)
  end subroutine FVReformKFUsg 

  subroutine FVReformKPerm(t_sync,ef_eldof)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    character(256) :: namegrp,namedat
    integer :: it,i,j,j2,ix,iy,iz,row,ef_eldof,idfile,err,idgrp,iddat,spc_dat, &
       spc_dom
    real(8) :: t_sync,perm1(dim_dom(1),dim_dom(2),dim_dom(3)),m_perm(dmn,dmn), &
       perm2(dim_dom(1),dim_dom(2),dim_dom(3)),wfv
    it=size(pack(t_fe,t_fe<=t_sync))
    wfv=(t_sync-t_fe(it))/(t_fe(it+1)-t_fe(it))
    ! Read permeability
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    call h5screate_simple_f(dmn,dim_dom,spc_dom,err)
    ! Time 1
    write(namegrp,'(A,ES12.5E2,A)')"/Time: ",t_fv(it)," d"
    call h5gopen_f(idfile,trim(namegrp),idgrp,err)
    write(namedat,'(A,A)')trim(namegrp),"/Permeability_X [m^2]"
    call h5dopen_f(idfile,trim(namedat),iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off_dom,dim_dom,err)
    call h5dread_f(iddat,h5t_native_double,perm1,dim_dom,err,spc_dom,spc_dat)
    call h5dclose_f(iddat,err) 
    call h5gclose_f(idgrp,err)
    ! Time 2
    write(namegrp,'(A,ES12.5E2,A)')"/Time: ",t_fv(it+1)," d"
    call h5gopen_f(idfile,trim(namegrp),idgrp,err)
    write(namedat,'(A,A)')trim(namegrp),"/Permeability_X [m^2]"
    call h5dopen_f(idfile,trim(namedat),iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off_dom,dim_dom,err)
    call h5dread_f(iddat,h5t_native_double,perm2,dim_dom,err,spc_dom,spc_dat)
    call h5dclose_f(iddat,err) 
    call h5gclose_f(idgrp,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
    ! Interpolate in time and space
    perm1=(f1-wfv)*perm1+wfv*perm2
    call MatZeroEntries(Mat_K,ierr)
    call MatZeroEntries(Mat_Kc,ierr)
    call VecZeroEntries(Vec_Cp,ierr)
    do i=1,nels
       enodes=nodes(el,:)
       ecoords=coords(enodes,:)
       ix=int((sum(ecoords(:,1))/dble(size(ecoords,1))-xmin)/dx)+1
       iy=int((sum(ecoords(:,2))/dble(size(ecoords,1))-ymin)/dy)+1
       iz=int((sum(ecoords(:,3))/dble(size(ecoords,1))-zmin)/dz)+1
       m_perm=f0
       do j=1,dmn ! Tensor valued permeability
          m_perm(j,j)=perm1(iz,iy,ix)*r_perm(j) ! HDF5 reverse order
       end do
       call FormLocalKPerm(i,k,indx,m_perm,"Kp") 
       indx=indxmap(indx,2)
       call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)
       ! Record pressure coefficient
       do j=1,npel
          j2=dmn*npel+j      
          row=nl2g(nodes(i,j),2)-1
          if (bc(nodes(i,j),dmn+1)==2) call VecSetValue(Vec_Cp,row,k(j2,j2),   &
             Add_Values,ierr)
       end do
       ! Update [Kc]
       call FormLocalKPerm(i,k,indx,m_perm,"Kc")
       kc=k(eldof+1:,eldof+1:)
       indxp=(indxmap(indx(eldof+1:),2)+1)/(dmn+1)-1
       call MatSetValues(Mat_Kc,eldofp,indxp,eldofp,indxp,kc,Add_Values,ierr)
    end do
    ! Account for constraint eqn's
    if (rank==0 .and. nceqs>0) call ApplyConstraints
    call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyBegin(Mat_Kc,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_Kc,Mat_Final_Assembly,ierr)
    call VecAssemblyBegin(Vec_Cp,ierr)
    call VecAssemblyEnd(Vec_Cp,ierr)
    do i=1,size(coords,1)
       row=nl2g(i,2)-1
       if (bc(i,dmn+1)/=2) call VecSetValue(Vec_Cp,row,f1,Insert_Values,ierr)  
    end do
    call VecAssemblyBegin(Vec_Cp,ierr)
    call VecAssemblyEnd(Vec_Cp,ierr)
  end subroutine FVReformKPerm

  subroutine FVReformKPermUsg(t_sync,ef_eldof)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    character(256) :: namegrp,namedat
    integer :: it,i,j,j2,row,ef_eldof,idfile,err,idgrp1,idgrp2,iddat1,iddat2,  &
        spc_dat1,spc_dat2,spc_cell
    integer(hsize_t) :: off(1),dim_cell(1)
    real(8) :: t_sync,perm1,perm2,m_perm(dmn,dmn),wfv
    it=size(pack(t_fe,t_fe<=t_sync))
    wfv=(t_sync-t_fe(it))/(t_fe(it+1)-t_fe(it))
    call MatZeroEntries(Mat_K,ierr)
    call MatZeroEntries(Mat_Kc,ierr)
    call VecZeroEntries(Vec_Cp,ierr)
    ! Read permeability
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    ! Time 1
    write(namegrp,'(A,I4,A,ES12.5E2,A)')"/",it-1," Time ",t_fv(it)," d"
    write(namedat,'(A,A)')trim(namegrp),"/Permeability X [m^2]"
    call h5gopen_f(idfile,trim(namegrp),idgrp1,err)
    call h5dopen_f(idfile,trim(namedat),iddat1,err)
    call h5dget_space_f(iddat1,spc_dat1,err)
    ! Time 2
    write(namegrp,'(A,I4,A,ES12.5E2,A)')"/",it," Time ",t_fv(it+1)," d"
    write(namedat,'(A,A)')trim(namegrp),"/Permeability X [m^2]"
    call h5gopen_f(idfile,trim(namegrp),idgrp2,err)
    call h5dopen_f(idfile,trim(namedat),iddat2,err)
    call h5dget_space_f(iddat2,spc_dat2,err)
    dim_cell=1
    call h5screate_simple_f(1,dim_cell,spc_cell,err)
    do i=1,nels
       off=el2g(i)-1
       ! Time 1 
       call h5sselect_hyperslab_f(spc_dat1,H5S_SELECT_SET_F,off,dim_cell,err)
       call h5dread_f(iddat1,h5t_native_double,perm1,dim_cell,err,spc_cell,    &
          spc_dat1)
       ! Time 2
       call h5sselect_hyperslab_f(spc_dat2,H5S_SELECT_SET_F,off,dim_cell,err)
       call h5dread_f(iddat2,h5t_native_double,perm2,dim_cell,err,spc_cell,    &
          spc_dat2)
       ! Interpolate
       perm1=(f1-wfv)*perm1+wfv*perm2
       m_perm=f0
       do j=1,dmn ! Tensor valued permeability
          m_perm(j,j)=perm1*r_perm(j) 
       end do
       call FormLocalKPerm(i,k,indx,m_perm,"Kp") 
       indx=indxmap(indx,2)
       call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)
       ! Record pressure coefficient
       do j=1,npel
          j2=dmn*npel+j      
          row=nl2g(nodes(i,j),2)-1
          if (bc(nodes(i,j),dmn+1)==2) call VecSetValue(Vec_Cp,row,k(j2,j2),   &
             Add_Values,ierr)
       end do
       ! Update [Kc]
       call FormLocalKPerm(i,k,indx,m_perm,"Kc")
       kc=k(eldof+1:,eldof+1:)
       indxp=(indxmap(indx(eldof+1:),2)+1)/(dmn+1)-1
       call MatSetValues(Mat_Kc,eldofp,indxp,eldofp,indxp,kc,Add_Values,ierr)
    end do
    call h5dclose_f(iddat1,err) 
    call h5gclose_f(idgrp1,err)
    call h5dclose_f(iddat2,err) 
    call h5gclose_f(idgrp2,err) 
    call h5fclose_f(idfile,err)
    call h5close_f(err)
    ! Account for constraint eqn's
    if (rank==0 .and. nceqs>0) call ApplyConstraints
    call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyBegin(Mat_Kc,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_Kc,Mat_Final_Assembly,ierr)
    call VecAssemblyBegin(Vec_Cp,ierr)
    call VecAssemblyEnd(Vec_Cp,ierr)
    do i=1,size(coords,1)
       row=nl2g(i,2)-1
       if (bc(i,dmn+1)/=2) call VecSetValue(Vec_Cp,row,f1,Insert_Values,ierr)  
    end do
    call VecAssemblyBegin(Vec_Cp,ierr)
    call VecAssemblyEnd(Vec_Cp,ierr)
  end subroutine FVReformKPermUsg

  ! Remove hydrostatic pressure gradient
  subroutine FVReset
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
    call VecAXPY(Vec_Up,f1,Vec_Up_hst,ierr)
    call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
    call VecDestroy(Vec_Up_hst,ierr)
  end subroutine FVReset

  ! Impose boundary pore pressure via RHS
  subroutine FVSyncBD(t_sync)  
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    integer :: i,it,row
    real(8) :: t_sync,wfv,pfv
    call VecZeroEntries(Vec_Up_fv,ierr)
    do it=1,size(t_fe)-1
       if (t_sync>=t_fe(it) .and. t_sync<t_fe(it+1)) then 
          wfv=(t_sync-t_fe(it))/(t_fe(it+1)-t_fe(it)) 
          do i=1,size(bdnd)
             pfv=(p_fv_bd(it,i)*(f1-wfv)+p_fv_bd(it+1,i)*wfv)/scale
             row=indxmap(bdnd(i)*(dmn+1),2)
             call VecSetValue(Vec_F,row,pfv,Insert_Values,ierr)
             row=nl2g(bdnd(i),2)-1
             call VecSetValue(Vec_Up_fv,row,pfv,Insert_Values,ierr)
          end do
       end if
    end do
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
    call VecAssemblyBegin(Vec_Up_fv,ierr)
    call VecAssemblyEnd(Vec_Up_fv,ierr)
    call VecGetSubVector(Vec_F,RI,Vec_Up,ierr) 
    call VecPointwiseMult(Vec_Up,Vec_Up,Vec_Cp,ierr) 
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
    ! Add FV pressure contribution to nodal force
    call MatMult(Mat_H,Vec_Up_fv,Vec_fp,ierr)
    call VecGetSubVector(Vec_F,RIu,Vec_Uu,ierr)
    call VecAXPY(Vec_Uu,-f1,Vec_fp,ierr)
    call VecRestoreSubVector(Vec_F,RIu,Vec_Uu,ierr)  
  end subroutine FVSyncBD 

  subroutine FVSyncBDUsg(t_sync)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    character(256) :: namegrp,namedat
    integer :: i,j,it,row,idfile,idgrp1,iddat1,idgrp2,iddat2,err,spc_dat1,     &
       spc_dat2,spc_cell
    integer(hsize_t) :: off(1),dim_cell(1)
    real(8) :: t_sync,pfv1,pfv2,cp,detj
    it=size(pack(t_fe,t_fe<t_sync))
    call VecGetSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecSet(Vec_Cp0,f1,ierr)
    do i=1,size(bdnd)
       row=nl2g(bdnd(i),2)-1 
       call VecSetValue(Vec_Up,row,f0,Insert_Values,ierr)
       call VecSetValue(Vec_Cp0,row,f0,Insert_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_Up,ierr)
    call VecAssemblyEnd(Vec_Up,ierr)
    call VecAssemblyBegin(Vec_Cp0,ierr)
    call VecAssemblyEnd(Vec_Cp0,ierr)
    call VecZeroEntries(Vec_Up_fv,ierr)
    ! Read pressure
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    ! Time 1
    write(namegrp,'(A,I4,A,ES12.5E2,A)')"/",it-1," Time ",t_fv(it)," d"
    write(namedat,'(A,A)')trim(namegrp),"/Liquid Pressure [Pa]"
    call h5gopen_f(idfile,trim(namegrp),idgrp1,err)
    call h5dopen_f(idfile,trim(namedat),iddat1,err)
    call h5dget_space_f(iddat1,spc_dat1,err)
    ! Time 2
    write(namegrp,'(A,I4,A,ES12.5E2,A)')"/",it," Time ",t_fv(it+1)," d"
    write(namedat,'(A,A)')trim(namegrp),"/Liquid Pressure [Pa]"
    call h5gopen_f(idfile,trim(namegrp),idgrp2,err)
    call h5dopen_f(idfile,trim(namedat),iddat2,err)
    call h5dget_space_f(iddat2,spc_dat2,err)
    dim_cell=1
    call h5screate_simple_f(1,dim_cell,spc_cell,err)
    do i=1,nels
       off=el2g(i)-1 
       ! Time 1  
       call h5sselect_hyperslab_f(spc_dat1,H5S_SELECT_SET_F,off,dim_cell,err)
       call h5dread_f(iddat1,h5t_native_double,pfv1,dim_cell,err,              &
          spc_cell,spc_dat1)
       ! Time 2
       call h5sselect_hyperslab_f(spc_dat2,H5S_SELECT_SET_F,off,dim_cell,err)
       call h5dread_f(iddat2,h5t_native_double,pfv2,dim_cell,err,spc_cell,     &
          spc_dat2)
       pfv1=(pfv2-pfv1)*dt/dt_fv/3600/24 ! Incremental pressure
       do j=1,nip
          call FormdetJ(ipoint(j,:),ecoords,detj)
          cp=cp+weight(j)*detj
       end do
       do j=1,npel
          if (bc(nodes(i,j),dmn+1)==2) then 
             row=nl2g(nodes(i,j),2)-1 
             call VecSetValue(Vec_Up,row,cp*pfv1/scale,Add_Values,ierr) 
             call VecSetValue(Vec_Up_fv,row,cp*pfv1/scale,Add_Values,ierr)
             call VecSetValue(Vec_Cp0,row,cp,Add_Values,ierr)
          end if
       end do 
    end do
    call h5dclose_f(iddat1,err) 
    call h5gclose_f(idgrp1,err)
    call h5dclose_f(iddat2,err) 
    call h5gclose_f(idgrp2,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
    call VecAssemblyBegin(Vec_Up,ierr)
    call VecAssemblyEnd(Vec_Up,ierr)
    call VecAssemblyBegin(Vec_Up_fv,ierr)
    call VecAssemblyEnd(Vec_Up_fv,ierr)
    call VecAssemblyBegin(Vec_Cp0,ierr)
    call VecAssemblyEnd(Vec_Cp0,ierr)
    call VecPointwiseMult(Vec_Up,Vec_Up,Vec_Cp,ierr) 
    call VecPointwiseDivide(Vec_Up,Vec_Up,Vec_Cp0,ierr)
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
    call VecPointwiseDivide(Vec_Up_fv,Vec_Up_fv,Vec_Cp0,ierr)
    ! Add FV pressure contribution to nodal force
    call MatMult(Mat_H,Vec_Up_fv,Vec_fp,ierr)
    call VecGetSubVector(Vec_F,RIu,Vec_Uu,ierr)
    call VecAXPY(Vec_Uu,-f1,Vec_fp,ierr)
    call VecRestoreSubVector(Vec_F,RIu,Vec_Uu,ierr)
  end subroutine FVSyncBDUsg 

end module fvfe 
