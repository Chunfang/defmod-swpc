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
  public :: FVInit
  public :: FVInitUsg
  public :: FVReformKF
  public :: FVReformKPerm
  public :: FVReset
  public :: FVSyncBD
 
  !! Local vars
  character(256) :: nameh5
  integer(hsize_t) :: off_dom(3),dim_dom(3) ! H5 offset and dim
  real(8) :: xmin,ymin,zmin,dx,dy,dz,r_perm(3),p_top ! Local FV/FE domain vars 
  integer, allocatable :: bdnd(:) ! Boundary nodes
  !! Pressure and time
  real(8), allocatable :: p_fv0(:),p_fv_sta(:),p_fv_bd(:,:),t_fv(:),t_fe(:) 

contains

  subroutine FVInitUsg
    implicit none
    character(256) :: name,name0,name1
    integer :: it,nt_fv
    real(8) :: t0_fv,dt_fv
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),"_fvfe.cfg"
    open(153,file=adjustl(name),status='old')
    read(153,*)nt_fv,t0_fv,dt_fv,p_top
    read(153,*)r_perm ! Inverse viscosity with scales
    close(153)
    do it=1,nt_fv
       t_fv(it)=(t0_fv+dble(it-1)*dt_fv)
    end do
  end subroutine FVInitUsg

  subroutine FVInit ! Read <model>_fvfe.cfg and <model>.h5 files
    implicit none
    character(256) :: name,name0,name1,namegrp,namedat
    integer :: i,iv,nbd,it,nt_fv,nx,ny,nz,idfile,err,idgrp,iddat,spc_dat,      &
       spc_dom,lnnds
    integer(hsize_t) :: ix,iy,iz
    real(8) :: t0_fv,dt_fv,xref,yref,zref,xmax,ymax,zmax
    real(8),allocatable :: p_dom(:,:,:)
    name0=output_file(:index(output_file,"/",BACK=.TRUE.))
    name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
    write(name,'(A,A,A)')trim(name0),trim(name1),"_fvfe.cfg"
    open(153,file=adjustl(name),status='old')
    read(153,*)nt_fv,t0_fv,dt_fv,p_top
    read(153,*)nx,ny,nz
    read(153,*)dx,dy,dz
    read(153,*)xref,yref,zref
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
    allocate(t_fv(nt_fv),t_fe(nt_fv),p_fv0(lnnds),bdnd(nbd),p_fv_bd(nt_fv,nbd),&
       p_fv_sta(lnnds))
    allocate(p_dom(dim_dom(1),dim_dom(2),dim_dom(3))) 
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    call h5screate_simple_f(dmn,dim_dom,spc_dom,err)
    ! Hydrostatic pressure from FV
    write(namegrp,'(AES12.5E2A)')"/Time: ",f0," d" 
    write(namedat,'(A,A)')trim(namegrp),"/Liquid_Pressure [Pa]"
    call h5gopen_f(idfile,trim(namegrp),idgrp,err)
    call h5dopen_f(idfile,trim(namedat),iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off_dom,dim_dom,err)
    call h5dread_f(iddat,h5t_native_double,p_dom,dim_dom,err,spc_dom,spc_dat)
    call h5dclose_f(iddat,err) 
    call h5gclose_f(idgrp,err)
    do iv=1,lnnds
       ix=int((coords(iv,1)-xmin)/dx)+1 
       iy=int((coords(iv,2)-ymin)/dy)+1
       iz=int((coords(iv,3)-zmin)/dz)+1
       p_fv_sta(iv)=p_dom(iz,iy,ix)
    end do
    ! Sampling pressure from FV
    do it=1,nt_fv
       t_fv(it)=(t0_fv+dble(it-1)*dt_fv)
       write(namegrp,'(AES12.5E2A)')"/Time: ",t_fv(it)," d" 
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
    call VecZeroEntries(Vec_Cp,ierr)
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
       end do
    end do
    ! Account for constraint eqn's
    if (rank==0 .and. nceqs>0) call ApplyConstraints
    call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
    call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
    call VecAssemblyBegin(Vec_Cp,ierr)
    call VecAssemblyEnd(Vec_Cp,ierr)
    call VecAssemblyBegin(Vec_Up,ierr)
    call VecAssemblyEnd(Vec_Up,ierr)
    ! Rescale RHS by pressure coefficient
    call VecPointwiseMult(Vec_up,Vec_Up,Vec_Cp,ierr)
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
  end subroutine FVReformKF
  
  subroutine FVReformKPerm(t_sync,ef_eldof)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    implicit none
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
    write(namegrp,'(AES12.5E2A)')"/Time: ",t_fv(it)," d"
    call h5gopen_f(idfile,trim(namegrp),idgrp,err)
    write(namedat,'(A,A)')trim(namegrp),"/Permeability_X [m^2]"
    call h5dopen_f(idfile,trim(namedat),iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,off_dom,dim_dom,err)
    call h5dread_f(iddat,h5t_native_double,perm1,dim_dom,err,spc_dom,spc_dat)
    call h5dclose_f(iddat,err) 
    call h5gclose_f(idgrp,err)
    ! Time 2
    write(namegrp,'(AES12.5E2A)')"/Time: ",t_fv(it+1)," d"
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
          if (bc(nodes(i,j),dmn+1)==2) call VecSetValue(Vec_Cp,row,k(j2,j2),  &
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

  subroutine FVReset
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    implicit none
    integer :: i,row
    do i=1,size(coords,1)
       row=i*(dmn+1)
       row=indxmap(row,2)
       call VecSetValue(Vec_Um,row,-(p_fv_sta(i)-p_top)/scale,Add_Values,ierr)
    end do
    call VecAssemblyBegin(Vec_Um,ierr)
    call VecAssemblyEnd(Vec_Um,ierr)
    deallocate(p_fv_sta,p_fv0)
  end subroutine FVReset

  subroutine FVSyncBD(t_sync) ! Impose boundary pore pressure via RHS
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h" 
#endif
    integer :: i,it,row
    real(8) :: t_sync,wfv,pfv
    do it=1,size(t_fe)-1
       if (t_sync>=t_fe(it) .and. t_sync<t_fe(it+1)) then 
          wfv=(t_sync-t_fe(it))/(t_fe(it+1)-t_fe(it)) 
          do i=1,size(bdnd)
             pfv=(p_fv_bd(it,i)*(f1-wfv)+p_fv_bd(it+1,i)*wfv)/scale
             row=indxmap(bdnd(i)*(dmn+1),2)
             call VecSetValue(Vec_F,row,pfv,Insert_Values,ierr)
          end do
       end if
    end do
    call VecAssemblyBegin(Vec_F,ierr)
    call VecAssemblyEnd(Vec_F,ierr)
    call VecGetSubVector(Vec_F,RI,Vec_Up,ierr) 
    call VecPointwiseMult(Vec_Up,Vec_Up,Vec_Cp,ierr) 
    call VecRestoreSubVector(Vec_F,RI,Vec_Up,ierr)
  end subroutine FVSyncBD 

end module fvfe 
