!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!! Set fault rupture source; output obs defined by Defmod.
!!
!! @copyright
!!   Copyright 2013-2019 Takuto Maeda, 2016-2018 Chunfang Meng. All rights reserved. This project is released under the MIT license.
!<
!! ----
#include "m_debug.h"
module m_rup

  !! dependency
  use m_std
  use m_global
  use m_debug
  use HDF5
  implicit none

  private
  save

  !! public routines
  public :: rup__Init
  public :: rup__getTime
  public :: rup__getFE2FD
  public :: rup__setSrc
  public :: rup__getObsFD
  public :: rup__evalObsFD
  public :: rup__writeObsFD

  !! local vars
  integer :: i0_rup ! start frame, number rupture samples
  integer :: nobsFE ! number of obs passed by FE
  integer :: nobsFD ! number on-rank obs for FD
  integer :: nviz ! obs sample interval
  integer :: nproc_fe ! number of FE MPI
  integer :: np_tot ! total point sources held by local FD rank
  integer, allocatable :: xobsFD(:,:)  ! Local obs grid indices
  integer, allocatable :: idx_src(:,:) ! Local source grid indices
  integer, allocatable :: idobsFE(:) ! global obs id alignment
  integer, allocatable :: idobsFD(:) ! obs id alignment
  real(MP) :: dt_rup ! rupture sample interval
  real(SP) :: xref,yref,zref ! FD corner at FE domain
  real(sp) :: km2m
  real(SP), allocatable :: xobsFE(:,:) ! Global obs coordinates, by FE
  real(SP), allocatable :: xobsFE_loc(:,:) ! Global obs coordinates, by FE
  real(MP), allocatable :: u_src(:,:,:) ! Local rupture source
  real(MP), allocatable :: obsFD(:,:)   ! Local obs passed from FE

contains

  subroutine rup__Init
    implicit none
    integer :: i
    character(256) :: name0
    write(name0,'(A,A)')trim(name_fe),"_fefd.cfg"
    km2m=1.0d3
    open(248,file=adjustl(name0),status='old')
    do i=1,3
      read(248,*)
    end do
    read(248,*)xref,yref,zref
    xref=km2m*xref; yref=km2m*yref; zref=km2m*zref
    read(248,*) nviz,nobsFE
    allocate(xobsFE(nobsFE,3),idobsFE(nobsFE))
    do i=1,nobsFE
      read(248,*)xobsFE(i,:),idobsFE(i)
    end do
    close(248)
  end subroutine rup__Init

  ! Read time series from FE output -> dt_rup, nt_rup
  subroutine rup__getTime
    implicit none
    integer :: i,i1_rup
    character(256) :: name0
    write(name0,'(A,A)')trim(name_fe),"_dyn.log"
    open(249,file=adjustl(name0),status='old')
    read(249,*)dt_rup
    if (eid>1) then
      ! skip unwanted events
      do i=1,eid-2
        read(249,*)
      end do
      read(249,*)i0_rup
      read(249,*)i1_rup
      nt_rup=i1_rup-i0_rup
    else
      i1_rup=0
      read(249,*)nt_rup ! only read first event
    end if
    close(249)
  end subroutine rup__getTime

  ! Read rupture source from FE outputs -> idx_src(:,:), u_src(:,:,:)
  subroutine rup__getFE2FD
    implicit none
    integer :: j,j2,nfile,nproc_fe2fd,ierr
    integer, allocatable :: idfile_pt(:),np_pt(:),idx_src_loc(:,:),rw(:),      &
      fe2fd(:),idfile(:),mattmp(:,:),fdact_loc(:),fdact(:),idx_tmp(:,:)
    real(MP), allocatable :: u_src_loc(:,:,:),u_tmp(:,:,:)
    character(256) :: name0,name1
    write(name0,'(A,A)')trim(name_fe),"_fe2fd.txt"
    open(250,file=adjustl(name0),status='old')
    read(250,*)nproc_fe,nproc_fe2fd
    call assert( nproc == nproc_fe2fd )
    allocate(fe2fd(nproc_fe),idfile(nproc_fe),mattmp(nproc_fe,nproc))
    do j=1,nproc_fe
      read(250,*)mattmp(j,:)
    end do
    close(250)
    fe2fd=0
    do j=1,nproc
      if (myid==j-1) fe2fd(:)=mattmp(:,j)
    end do
    do j=1,nproc_fe
      if(fe2fd(j)/=0) then
        idfile(j)=j-1
      else
        idfile(j)=-1
      end if
    end do
    ! Loop over files to read sources
    nfile=size(pack(fe2fd,fe2fd/=0))
    np_tot=sum(fe2fd)
    call MPI_Allreduce(np_tot,np_rup,1,MPI_INTEGER,MPI_Sum,MPI_Comm_World,ierr)
    allocate(idfile_pt(nfile),np_pt(nfile),idx_tmp(np_tot,3),fdact(np_tot),    &
      u_tmp(np_tot,nt_rup,3))
    idfile_pt=pack(idfile,idfile>-1)
    np_pt=pack(fe2fd,fe2fd/=0)
    do j=1,nfile
      allocate(rw(np_pt(j)),idx_src_loc(np_pt(j),3),u_src_loc(np_pt(j),        &
        nt_rup,3),fdact_loc(np_pt(j)))
      rw=(/(sum(np_pt(:j))-np_pt(j)+j2,j2=1,np_pt(j))/)
      write(name1,'(A,A,I0.6,A)')trim(name_fe),"_fe2fd_",idfile_pt(j),".h5"
      call RupSrc(name1,idx_src_loc,u_src_loc,fdact_loc,size(rw),nt_rup)
      idx_tmp(rw,:)=idx_src_loc
      u_tmp(rw,:,:)=u_src_loc
      fdact(rw)=fdact_loc
      deallocate(rw,idx_src_loc,u_src_loc,fdact_loc)
    end do
    np_tot=size(pack(fdact,fdact>0))
    allocate(idx_src(np_tot,3),u_src(np_tot,nt_rup,3))
    j2=0
    do j=1,size(fdact)
       if (fdact(j)>0) then
          j2=j2+1
          idx_src(j2,:)=idx_tmp(j,:)
          u_src(j2,:,:)=u_tmp(j,:,:)
       end if
    end do
  end subroutine rup__getFE2FD

  ! Read rupture source from FE patch
  subroutine RupSrc(nameh5,idx_src,u_src,fdact_loc,nrw,nt_rup)
    implicit none
    integer(hid_t) :: idfile,iddat,spc_dat,spc_src
    integer(hsize_t) :: dim_dat(3),dim_tmp(3),offset(3)
    character(256) :: nameh5
    integer :: rankxy,isrc,hit,i,j,k,n_glb,nrw,nt_rup,idx_src(nrw,3),          &
      fdact_loc(nrw),act,err
    integer,allocatable :: ID_glb(:,:)
    real(MP) :: vx,vy,vz,u_src(nrw,nt_rup,3)
    real(MP),allocatable :: src_glb(:,:,:)
    call h5open_f(err)
    call h5fopen_f(trim(nameh5),H5F_ACC_RDWR_F,idfile,err)
    call h5dopen_f(idfile,"ID",iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    call h5sget_simple_extent_dims_f(spc_dat,dim_dat(2:3),dim_tmp(2:3),err)
    n_glb=dim_dat(3)
    allocate(ID_glb(dim_dat(2),n_glb),src_glb(nt_rup,3,n_glb))
    call h5dread_f(iddat,h5t_native_integer,ID_glb,dim_dat(2:3),err)
    call h5dclose_f(iddat,err)
    call h5dopen_f(idfile,"disp",iddat,err)
    call h5dget_space_f(iddat,spc_dat,err)
    ! Skip unwanted frames
    offset=(/i0_rup,0,0/)
    dim_tmp(1:2)=(/nt_rup,3/)
    call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_tmp,err)
    call h5screate_simple_f(3,dim_tmp,spc_src,err)
    call h5dread_f(iddat,h5t_native_double,src_glb,dim_tmp,err,spc_src,spc_dat)
    call h5dclose_f(iddat,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
    hit=0
    do isrc=1,n_glb
      i=ID_glb(1,isrc); j=ID_glb(2,isrc); k=ID_glb(3,isrc)
      rankxy=((j-1)/nyp)*nproc_x+(i-1)/nxp
      if (myid==rankxy .and. hit<nrw) then
        hit=hit+1
        idx_src(hit,:)=(/i,j,k/)
        fdact_loc(hit)=ID_glb(3+max(1,eid),isrc)
        u_src(hit,:,:)=src_glb(:,:,isrc)
      end if
    end do
  end subroutine RupSrc

  ! Add rupture source
  subroutine rup__setSrc(it)
    implicit none
    integer :: it,i,i_rup,ii,jj,kk
    real(MP) :: t_rup,t,vxs,vys,vzs
    t = dble(it-1)*dt
    i_rup=int(t/dt_rup)
    if (i_rup>0 .and. i_rup<nt_rup) then
      t_rup=i_rup*dt_rup
      do i=1,np_tot
        ii=idx_src(i,1)
        jj=idx_src(i,2)
        kk=idx_src(i,3)
        ! Linearly interpolate the source velocity
        vx(kk,ii,jj)=(u_src(i,i_rup+1,1)-u_src(i,i_rup,1))/dt_rup
        vy(kk,ii,jj)=(u_src(i,i_rup+1,2)-u_src(i,i_rup,2))/dt_rup
        vz(kk,ii,jj)=(u_src(i,i_rup+1,3)-u_src(i,i_rup,3))/dt_rup
      end do
    elseif (i_rup>=nt_rup) then
      do i=1,np_tot
        ii=idx_src(i,1)
        jj=idx_src(i,2)
        kk=idx_src(i,3)
        ! Zero velocity (hard source) to prevent back slip
        vx(kk,ii,jj)=0.d0
        vy(kk,ii,jj)=0.d0
        vz(kk,ii,jj)=0.d0
      end do
    endif
  end subroutine rup__setSrc

  ! FE observation to FD->nobsFD, xobsFD, obsFD, xobsFE_loc
  subroutine rup__getObsFD
    implicit none
    integer :: i,xid,yid,zid,rankxy,xobsFD_full(nobsFE,3),hit(nobsFE)
    real(MP) :: xfe,yfe,zfe,x0,y0,z0,x1,y1,z1
    do i=1,nobsFE
      xfe= xobsFE(i,1)*km2m-xref
      yfe= xobsFE(i,2)*km2m-yref
      zfe=-xobsFE(i,3)*km2m+zref ! deep negative by FE
      xid=int(xfe/dx/km2m)+1
      yid=int(yfe/dy/km2m)+1
      zid=int(zfe/dz/km2m)+1
      x0=(xid-1)*dx*km2m; x1=xid*dx*km2m
      y0=(yid-1)*dy*km2m; y1=yid*dy*km2m
      z0=(zid-1)*dz*km2m; z1=zid*dz*km2m
      if (xfe-x0>=x1-xfe) xid=xid+1
      if (yfe-y0>=y1-yfe) yid=yid+1
      if (zfe-z0>=z1-zfe) zid=zid+1
      rankxy=((yid-1)/nyp)*nproc_x+(xid-1)/nxp
      if (myid==rankxy) then
        xobsFD_full(i,:)=(/xid,yid,zid/)
        hit(i)=i
      else
        hit(i)=0
      end if
    end do
    nobsFD=size(pack(hit,hit/=0))
    allocate(xobsFD(nobsFD,3),obsFD(nobsFD,3),xobsFE_loc(nobsFD,3),idobsFD(nobsFD))
    xobsFD=xobsFD_full(pack(hit,hit/=0),:)
    xobsFE_loc=xobsFE(pack(hit,hit/=0),:)
    idobsFD=idobsFE(pack(hit,hit/=0))
  end subroutine rup__getObsFD

  ! Evaluate FD observation -> obsFD
  subroutine rup__evalObsFD(it)
    implicit none
    integer :: i,it
    if (mod(it,nviz)==0 .and. nobsFD>0) then
      do i=1,nobsFD
        obsFD(i,:)=(/vx(xobsFD(i,3),xobsFD(i,1),xobsFD(i,2)),&
                     vy(xobsFD(i,3),xobsFD(i,1),xobsFD(i,2)),&
                    -vz(xobsFD(i,3),xobsFD(i,1),xobsFD(i,2))/) ! deep negative by FE
      end do
    end if
  end subroutine rup__evalObsFD

  ! Output FD observations
  subroutine rup__writeObsFD(it)
    implicit none
    character(256) :: name_fd
    character(64) :: fmt
    integer :: i,it
    integer,save :: k=0
    if (mod(it,nviz)==0 .and. nobsFD>0) then
      write(name_fd,'(A,A,I0,A)')trim(name_fe),"_",myid,"_ofd.txt"
      if (k==0) then
        open(252,file=adjustl(name_fd),status='replace')
        write(252,'(I0,1X,I0,1X,F0.6)')nobsFD,nt/nviz,dt*dble(nviz)
        fmt='(3(F0.3,1X),I0)'
        do i=1,nobsFD
           write(252,fmt)xobsFE_loc(i,:),idobsFD(i)
        end do
      else
        open(252,file=adjustl(name_fd),status='old',position='append',action='write')
      endif
      fmt='(3(F0.6,1X))'
      do i=1,nobsFD
        write(252,fmt)obsFD(i,:)
      end do
      close(252)
      k=k+1
    end if
  end subroutine rup__writeObsFD

end module m_rup
