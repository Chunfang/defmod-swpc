! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module h5io

  use global
  use HDF5

  private
  save

  !! Public routines
  public :: Write_fault
  public :: Write_obs
  public :: Write_fd

  !! Local vars
  character(256) :: file_flt,file_obs,file_fd
  integer(hsize_t) :: dim_flt(3),dim_st(3),dim_fd(3)
contains

  ! Output fault slip qs_flt_slip (sta) or (tot_)flt_slip (dyn)
  subroutine Write_fault(strng)
    implicit none
    integer(hid_t) :: idfile,iddat,spc_dat,spc_flt,spc_trc,config
    integer(hsize_t) :: dim_flts(3),dim_trc(3),dim_trcs(3),offset(3),          &
       dim_fndx(2),limit(3)
    integer :: j1,j3,err
    integer,save :: k=0,kd=0
    character(3) :: strng
    character(256) :: name0,name1,namedat,nametrc
    real(8) :: dat_slip(1,dmn,nfnd_loc),dat_trac(1,2*dmn,nfnd_loc),            &
       dat_fndx(dmn*3,nfnd_loc) !H5 reverse order
    call h5open_f(err)
    if (k==0 .and. kd==0) then
       name0=output_file(:index(output_file,"/",BACK=.TRUE.))
       name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
       write(file_flt,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_flt_",rank,   &
          ".h5"
       call h5fcreate_f(trim(file_flt),H5F_ACC_TRUNC_F,idfile,err)
       ! Create initial date space and close
       do j1=1,nfnd_loc
          j3=FltMap(j1,2)
          dat_fndx(:dmn,j1)=xfnd(j3,:)
          dat_fndx(dmn+1:,j1)=vecf(j3,:dmn*2)
       end do
       ! Write coordinates
       dim_fndx=(/dmn*3,nfnd_loc/)
       call h5screate_simple_f(2,dim_fndx,spc_dat,err)
       call h5dcreate_f(idfile,"fault_x",h5t_native_double,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_double,dat_fndx,dim_fndx,err)
       call h5dclose_f(iddat,err)
       call h5fclose_f(idfile,err)
       dim_flt=(/1,dmn,nfnd_loc/)
    end if
    if (strng=="sta") then
       namedat="slip_sta"
       nametrc="trac_sta"
       call GetSlipDat(dat_slip,qs_flt_slip)
       dim_trc=(/1,dmn+p,nfnd_loc/)
       k=k+1; dim_flts=(/k,dmn,nfnd_loc/); dim_trcs=(/k,dmn+p,nfnd_loc/)
    elseif (strng=="dyn") then
       namedat="slip_dyn"
       nametrc="trac_dyn"
       if (dsp_hyb==0) then
          call GetSlipDat(dat_slip,flt_slip/dt_dyn)
       elseif (dsp_hyb==1) then
          call GetSlipDat(dat_slip,tot_flt_slip)
       end if
       dim_trc=(/1,2*dmn,nfnd_loc/)
       kd=kd+1; dim_flts=(/kd,dmn,nfnd_loc/); dim_trcs=(/kd,2*dmn,nfnd_loc/)
    end if
    call GetTracDat(dat_trac,strng)
    call h5fopen_f(trim(file_flt),H5F_ACC_RDWR_F,idfile,err)
    call h5screate_simple_f(3,dim_flt,spc_flt,err)
    call h5screate_simple_f(3,dim_trc,spc_trc,err)
    if ((k==1 .and. strng=="sta") .or. (kd==1 .and. strng=="dyn")) then !Create
       limit(1)=H5S_UNLIMITED_F; limit(2:3)=(/dmn,nfnd_loc/)
       call h5screate_simple_f(3,dim_flt,spc_dat,err,limit)
       call h5pcreate_f(H5P_DATASET_CREATE_F,config,err)
       call h5pset_chunk_f(config,3,dim_flt,err)
       call h5dcreate_f(idfile,namedat,h5t_native_double,spc_dat,iddat,err,    &
          config)
       call h5dwrite_f(iddat,h5t_native_double,dat_slip,dim_flt,err)
       if (strng=="sta") then ! Static traction and pressure
          limit(2:3)=(/dmn+p,nfnd_loc/)
          call h5screate_simple_f(3,dim_trc,spc_dat,err,limit)
          call h5pcreate_f(H5P_DATASET_CREATE_F,config,err)
          call h5pset_chunk_f(config,3,dim_trc,err)
          call h5dcreate_f(idfile,nametrc,h5t_native_double,spc_dat,iddat,err, &
             config)
          call h5dwrite_f(iddat,h5t_native_double,dat_trac(:,:dmn+p,:),        &
             dim_trc,err)
       else ! Dynamic traction
          limit(2:3)=(/2*dmn,nfnd_loc/)
          call h5screate_simple_f(3,dim_trc,spc_dat,err,limit)
          call h5pcreate_f(H5P_DATASET_CREATE_F,config,err)
          call h5pset_chunk_f(config,3,dim_trc,err)
          call h5dcreate_f(idfile,nametrc,h5t_native_double,spc_dat,iddat,err, &
             config)
          call h5dwrite_f(iddat,h5t_native_double,dat_trac(:,:2*dmn,:),dim_trc,&
             err)
       end if
    else ! Extend data
       ! Slip
       offset(1)=dim_flts(1)-1; offset(2:3)=0
       call h5dopen_f(idfile,namedat,iddat,err)
       call h5dset_extent_f(iddat,dim_flts,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_flt,err)
       call H5dwrite_f(iddat,h5t_native_double,dat_slip,dim_flt,err,spc_flt,   &
          spc_dat)
       ! Traction
       offset(1)=dim_trcs(1)-1; offset(2:3)=0
       call h5dopen_f(idfile,nametrc,iddat,err)
       call h5dset_extent_f(iddat,dim_trcs,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_trc,err)
       if (strng=="sta") then ! Static traction and pressure
          call H5dwrite_f(iddat,h5t_native_double,dat_trac(:,:dmn+p,:),        &
             dim_trc,err,spc_trc,spc_dat)
       else ! Dynamic traction
          call H5dwrite_f(iddat,h5t_native_double,dat_trac(:,:2*dmn,:),dim_trc,&
             err,spc_trc,spc_dat)
       end if
    end if
    call h5dclose_f(iddat,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
  end subroutine Write_fault

  ! Write observation data (tot_)uu(_dyn)_obs, (tot_)st_obs
  subroutine Write_obs(strng)
    implicit none
    integer(hid_t) :: idfile,iddat,spc_dat,spc_uu,spc_st,config
    integer(hsize_t) :: dim_uu(3),dim_uus(3),dim_sts(3),offset(3),dim_obsx(2), &
       dim_idx(1), limit(3)
    integer :: j1,j3,err,dat_idx(nobs_loc)
    integer,save :: k=0,kd=0
    character(3) :: strng
    character(256) :: name0,name1,nameuu,namest
    real(8) :: dat_uu(1,dmn+p,nobs_loc),dat_st(1,cdmn,nobs_loc),               &
       dat_obsx(dmn,nobs_loc) !H5 reverse order
    call h5open_f(err)
    if (k==0 .and. kd==0) then
       name0=output_file(:index(output_file,"/",BACK=.TRUE.))
       name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
       write(file_obs,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_ofe_",rank,   &
          ".h5"
       call h5fcreate_f(trim(file_obs),H5F_ACC_TRUNC_F,idfile,err)
       ! Create initial date space and close
       do j1=1,nobs_loc
          j3=ol2g(j1)
          dat_idx(j1)=j3
          dat_obsx(:,j1)=ocoord(j3,:)
       end do
       ! Write coordinates
       dim_obsx=(/dmn,nobs_loc/)
       call h5screate_simple_f(2,dim_obsx,spc_dat,err)
       call h5dcreate_f(idfile,"obs_x",h5t_native_double,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_double,dat_obsx,dim_obsx,err)
       dim_idx=(/nobs_loc/)
       call h5screate_simple_f(1,dim_idx,spc_dat,err)
       call h5dcreate_f(idfile,"idx",h5t_native_integer,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_integer,dat_idx,dim_idx,err)
       call h5dclose_f(iddat,err)
       call h5fclose_f(idfile,err)
    end if
    limit(1)=H5S_UNLIMITED_F
    if (strng=="sta") then
       dim_uu=(/1,dmn+p,nobs_loc/)
       dim_st=(/1,cdmn,nobs_loc/)
       nameuu="uu_sta"
       namest="st_sta"
       do j1=1,nobs_loc
          if (dsp==0) then
             dat_uu(1,:,j1)=uu_obs(j1,:)
             dat_st(1,:,j1)=st_obs(j1,:)
          else
             dat_uu(1,:,j1)=tot_uu_obs(j1,:)
             dat_st(1,:,j1)=tot_st_obs(j1,:)
          end if
       end do
       k=k+1; dim_uus=(/k,dmn+p,nobs_loc/); dim_sts=(/k,cdmn,nobs_loc/)
       limit(2:3)=(/dmn+p,nobs_loc/)
    else
       dim_uu=(/1,dmn,nobs_loc/)
       nameuu="uu_dyn"
       do j1=1,nobs_loc
          if (dsp_hyb==0) then
             dat_uu(1,:dmn,j1)=uu_dyn_obs(j1,:)
          else
             dat_uu(1,:dmn,j1)=tot_uu_dyn_obs(j1,:)
          end if
       end do
       kd=kd+1; dim_uus=(/kd,dmn,nobs_loc/)
       limit(2:3)=(/dmn,nobs_loc/)
    end if
    call h5fopen_f(trim(file_obs),H5F_ACC_RDWR_F,idfile,err)
    call h5screate_simple_f(3,dim_uu,spc_uu,err)
    call h5screate_simple_f(3,dim_st,spc_st,err)
    if ((k==1 .and. strng=="sta") .or. (kd==1 .and. strng=="dyn")) then !Create
       call h5screate_simple_f(3,dim_uu,spc_dat,err,limit)
       call h5pcreate_f(H5P_DATASET_CREATE_F,config,err)
       call h5pset_chunk_f(config,3,dim_uu,err)
       call h5dcreate_f(idfile,nameuu,h5t_native_double,spc_dat,iddat,err,     &
          config)
       if (strng=="sta") then ! Static stress
          call h5dwrite_f(iddat,h5t_native_double,dat_uu,dim_uu,err)
          limit(2:3)=(/cdmn,nobs_loc/)
          call h5screate_simple_f(3,dim_st,spc_dat,err,limit)
          call h5pcreate_f(H5P_DATASET_CREATE_F,config,err)
          call h5pset_chunk_f(config,3,dim_st,err)
          call h5dcreate_f(idfile,namest,h5t_native_double,spc_dat,iddat,err,  &
             config)
          call h5dwrite_f(iddat,h5t_native_double,dat_st,dim_st,err)
       else
          call h5dwrite_f(iddat,h5t_native_double,dat_uu(1,:dmn,:),dim_uu,err)
       end if
    else ! Extend data
       offset(1)=dim_uus(1)-1; offset(2:3)=0
       call h5dopen_f(idfile,nameuu,iddat,err)
       call h5dset_extent_f(iddat,dim_uus,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_uu,err)
       if (strng=="sta") then ! Static traction and pressure
          call H5dwrite_f(iddat,h5t_native_double,dat_uu,dim_uu,err,spc_uu,    &
             spc_dat)
          offset(1)=dim_sts(1)-1; offset(2:3)=0
          call h5dopen_f(idfile,namest,iddat,err)
          call h5dset_extent_f(iddat,dim_sts,err)
          call h5dget_space_f(iddat,spc_dat,err)
          call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_st,err)
          call H5dwrite_f(iddat,h5t_native_double,dat_st,dim_st,err,spc_st,    &
             spc_dat)
       else
          call H5dwrite_f(iddat,h5t_native_double,dat_uu(:,:dmn,:),dim_uu,err, &
             spc_uu,spc_dat)
       end if
    end if
    call h5dclose_f(iddat,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
  end subroutine Write_obs

  subroutine Write_fd(strng)
    implicit none
    integer(hid_t) :: idfile,iddat,spc_dat,spc_fd,config
    integer(hsize_t) :: dim_dat(3),dim_tmp(3),dim_ID(2),offset(3),limit(3)
    integer :: dat_fdact(1,ngp_loc),err
    integer,save :: k=0
    character(256) :: name0,name1
    character(3) :: strng
    real(8) :: dat_fd(1,dmn,ngp_loc) ! H5 reverse order
    call h5open_f(err)
    dim_fd=(/1,dmn,ngp_loc/)
    dat_fd(1,:,:)=transpose(uu_fd)
    dat_fd(1,3,:)=-dat_fd(1,3,:) ! deep positive in swpc
    if (strng=='mat') then
       name0=output_file(:index(output_file,"/",BACK=.TRUE.))
       name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
       write(file_fd,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_fe2fd_",rank,  &
          ".h5"
       call h5fcreate_f(trim(file_fd),H5F_ACC_TRUNC_F,idfile,err)
       dim_dat(1:2)=(/size(matFD,2),size(matFD,1)/)
       call h5screate_simple_f(2,dim_dat(1:2),spc_fd,err)
       call h5dcreate_f(idfile,"mat",h5t_native_integer,spc_fd,iddat,err)
       call h5dwrite_f(iddat,h5t_native_integer,transpose(matFD),dim_dat(1:2)  &
          ,err)
    elseif (k==1) then
       dim_ID=(/dmn,ngp_loc/)
       limit(1)=H5S_UNLIMITED_F; limit(2)=ngp_loc
       call h5screate_simple_f(2,dim_ID,spc_fd,err,limit(1:2))
       call h5pcreate_f(H5P_DATASET_CREATE_F,config,err)
       call h5pset_chunk_f(config,2,dim_ID,err)
       call h5fopen_f(trim(file_fd),H5F_ACC_RDWR_F,idfile,err)
       call h5dcreate_f(idfile,"ID",h5t_native_integer,spc_fd,iddat,err,config)
       call h5dwrite_f(iddat,h5t_native_integer,transpose(idgp_loc),dim_ID,err)
       call h5dclose_f(iddat,err)
       ! Create data
       limit(1)=H5S_UNLIMITED_F; limit(2:3)=(/dmn,ngp_loc/)
       call h5screate_simple_f(3,dim_fd,spc_dat,err,limit)
       call h5pcreate_f(H5P_DATASET_CREATE_F,config,err)
       call h5pset_chunk_f(config,3,dim_fd,err)
       call h5dcreate_f(idfile,"disp",h5t_native_double,spc_dat,iddat,err,     &
          config)
       call h5dwrite_f(iddat,h5t_native_double,dat_fd,dim_fd,err)
    elseif (strng=="dyn") then ! Extend data
       call h5fopen_f(trim(file_fd),H5F_ACC_RDWR_F,idfile,err)
       call h5screate_simple_f(3,dim_fd,spc_fd,err)
       call h5dopen_f(idfile,"disp",iddat,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sget_simple_extent_dims_f(spc_dat,dim_dat,dim_tmp,err)
       offset(1)=dim_dat(1); offset(2:3)=0
       dim_dat(1)=dim_dat(1)+1
       call h5dset_extent_f(iddat,dim_dat,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_fd,err)
       call H5dwrite_f(iddat,h5t_native_double,dat_fd,dim_fd,err,spc_fd,spc_dat)
    elseif (strng=="act") then ! Extend node activity
       call h5fopen_f(trim(file_fd),H5F_ACC_RDWR_F,idfile,err)
       call h5dopen_f(idfile,"ID",iddat,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sget_simple_extent_dims_f(spc_dat,dim_dat(1:2),dim_tmp(1:2),err)
       offset(1)=dim_dat(1); offset(2)=0
       dim_dat(1)=dim_dat(1)+1
       dim_tmp(1)=1
       call h5dset_extent_f(iddat,dim_dat(1:2),err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset(1:2),        &
          dim_tmp(1:2),err)
       call h5screate_simple_f(2,dim_tmp(1:2),spc_fd,err)
       dat_fdact(1,:)=fdact_loc ! H5 Reverse
       call H5dwrite_f(iddat,h5t_native_integer,dat_fdact,dim_tmp(1:2),err,    &
          spc_fd,spc_dat)
    end if
    call h5dclose_f(iddat,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
    k=k+1
  end subroutine Write_fd

  ! Create fault slip data to output
  subroutine GetSlipDat(dat_slip,dat_src)
    implicit none
    integer :: j,j1,j2,rw_loc(dmn)
    real(8) :: dat_slip(:,:,:),dat_src(:)
    do j1=1,nfnd_loc
       j=FltMap(j1,1)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       dat_slip(1,:,j1)=dat_src(rw_loc)
    end do
  end subroutine GetSlipDat

  ! Create fault traction/pressure data to output
  subroutine GetTracDat(dat_trac,strng)
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    character(3) :: strng
    integer :: j,j1,j2,j3,rw_loc(dmn)
    real(8) :: dat_trac(:,:,:),flt_qs(dmn),flt_p,flt_dyn(dmn)
    real(8),target :: flt_ndf(n_lmnd*dmn),lm_pn(n_lmnd),lm_pp(n_lmnd),         &
       lm_f2s(n_lmnd),flt_ndf_dyn(n_lmnd*dmn)
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
    end if
    call VecGetArrayF90(Vec_lm_f2s,pntr,ierr)
    lm_f2s=pntr
    call VecRestoreArrayF90(Vec_lm_f2s,pntr,ierr)
    if (strng=="dyn") then
       call VecGetArrayF90(Vec_lambda_tot,pntr,ierr)
       flt_ndf_dyn=pntr
       call VecRestoreArrayF90(Vec_lambda_tot,pntr,ierr)
    end if
    do j1=1,nfnd_loc
       j=FltMap(j1,1); j3=FltMap(j1,2)
       rw_loc=(/((j-1)*dmn+j2,j2=1,dmn)/)
       ! Rotate to fault coordinate
       flt_qs=flt_ndf(rw_loc)
       call Cart2Flt(vecf(j3,:),flt_qs,1)
       flt_qs=(flt_qs+st_init(j3,:dmn))*lm_f2s(j)
       if (poro) then
          flt_p=(lm_pp(j)+lm_pn(j))/f2
          ! Positive pressure affects normal stress
          flt_qs(dmn)=flt_qs(dmn)+max(f0,biot(j3)*flt_p)
       end if
       if (poro .and. strng=="sta") then
          dat_trac(1,:dmn+p,j1)=(/flt_qs,flt_p/)
       else if (strng=="sta") then
          dat_trac(1,:dmn,j1)=flt_qs
       else ! Dynamic traction
          dat_trac(1,:dmn,j1)=trac_dyn(j,:dmn)*lm_f2s(j)
          flt_dyn=flt_ndf_dyn(rw_loc)*lm_f2s(j)
          call Cart2Flt(vecf(j3,:),flt_dyn,1)
          dat_trac(1,dmn+1:,j1)=flt_qs+flt_dyn
       end if
    end do
  end subroutine GetTracDat

end module h5io
