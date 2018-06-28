! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module h5io

  use global
  use HDF5

  private
  save

  !! Public routines
  public :: Write_fe
  public :: Write_fd

  !! Local vars
  character(256) :: fileh5, file_fd
  integer(hsize_t) :: dim_flt(3),dim_trc(3),dim_fd(3)
contains
  
  ! Output fault slip qs_flt_slip (sta) or (tot_)flt_slip (dyn)
  subroutine Write_fe(strng)
    implicit none
    integer(hid_t) :: idfile,iddat,spc_dat,spc_flt,spc_trc,config
    integer(hsize_t) :: dim_slip(3),dim_strs(3),offset(3),dim_fndx(2),limit(3)
    integer :: j1,j3,err
    integer,save :: k=0,kd=0
    character(3) :: strng
    character(256) :: name0,name1,namedat,nametrc
    real(8) :: dat_slip(1,dmn,nfnd_loc),dat_trac(1,dmn+p,nfnd_loc),            &
       dat_fndx(dmn,nfnd_loc) !H5 reverse order
    call h5open_f(err)
    if (k==0 .and. kd==0) then
       name0=output_file(:index(output_file,"/",BACK=.TRUE.))
       name1=output_file(index(output_file,"/",BACK=.TRUE.)+1:)
       write(fileh5,'(A,A,A,I0.6,A)')trim(name0),trim(name1),"_fe_",rank,".h5"
       call h5fcreate_f(trim(fileh5),H5F_ACC_TRUNC_F,idfile,err)
       ! Create initial date space and close 
       do j1=1,nfnd_loc
          j3=FltMap(j1,2)
          dat_fndx(:,j1)=xfnd(j3,:)
       end do
       ! Write coordinates
       dim_fndx=(/dmn,nfnd_loc/)
       call h5screate_simple_f(2,dim_fndx,spc_dat,err)
       call h5dcreate_f(idfile,"fault_x",h5t_native_double,spc_dat,iddat,err)
       call h5dwrite_f(iddat,h5t_native_double,dat_fndx,dim_fndx,err)
       call h5dclose_f(iddat,err)
       call h5fclose_f(idfile,err)
    end if
    dim_flt=(/1,dmn,nfnd_loc/)
    dim_trc=(/1,dmn+p,nfnd_loc/)
    if (strng=="sta") then
       namedat="slip_sta" 
       nametrc="trac_sta"
       call GetSlipDat(dat_slip,qs_flt_slip)
       call GetTracDat(dat_trac)
       k=k+1; dim_slip=(/k,dmn,nfnd_loc/); dim_strs=(/k,dmn+p,nfnd_loc/)
    elseif (strng=="dyn") then
       namedat="slip_dyn" 
       if (dsp_hyb==0) then
          call GetSlipDat(dat_slip,flt_slip/dt_dyn)
       elseif (dsp_hyb==1) then 
          call GetSlipDat(dat_slip,tot_flt_slip)
       end if
       kd=kd+1; dim_slip=(/kd,dmn,nfnd_loc/)
    end if
    call h5fopen_f(trim(fileh5),H5F_ACC_RDWR_F,idfile,err)
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
          call h5dwrite_f(iddat,h5t_native_double,dat_trac,dim_trc,err)
       end if
    else ! Extend data data
       offset(1)=dim_slip(1)-1; offset(2:3)=0
       call h5dopen_f(idfile,namedat,iddat,err)
       call h5dset_extent_f(iddat,dim_slip,err)
       call h5dget_space_f(iddat,spc_dat,err)
       call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_flt,err)
       call H5dwrite_f(iddat,h5t_native_double,dat_slip,dim_flt,err,spc_flt,   &
          spc_dat)
       if (strng=="sta") then ! Static traction and pressure
          offset(1)=dim_strs(1)-1; offset(2:3)=0
          call h5dopen_f(idfile,nametrc,iddat,err)
          call h5dset_extent_f(iddat,dim_strs,err)
          call h5dget_space_f(iddat,spc_dat,err)
          call h5sselect_hyperslab_f(spc_dat,H5S_SELECT_SET_F,offset,dim_trc,  &
             err)
          call H5dwrite_f(iddat,h5t_native_double,dat_trac,dim_trc,err,spc_trc,&
             spc_dat)
       end if
    end if
    call h5dclose_f(iddat,err)
    call h5fclose_f(idfile,err)
    call h5close_f(err)
  end subroutine Write_fe

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

end module h5io 
