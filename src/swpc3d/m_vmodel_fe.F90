!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!!  Overwrite the FD velocity with FE velocity
!!
!! @copyright
!!   Copyright 2013-2018 Takuto Maeda, 2016-2018 Chunfang Meng. All rights reserved. This project is released under the MIT license.
!<
!! ----
module m_vmodel_fe

  !! dependency
  use m_std
  use m_global
  use m_debug
  use HDF5
  implicit none

  private
  save

  !! public routines
  public :: vmodel__fe
  !! local vars
contains

  subroutine vmodel__fe
    implicit none
    character(256) :: name0,name1
    integer :: nproc_fe,nproc_fe2fd,nmat,cmat,imat,nels,irank
    integer :: i0,i1,j0,j1,k0,k1,i,j,k,el,err ! i,j,k regions
    integer,allocatable :: mattmp(:,:),dat_mat(:,:)
    real(SP),allocatable :: mat(:,:)
    real(SP) :: E,nu,lam0,mu0,rho0
    integer(hid_t) :: idfile,iddat,spc_dat
    integer(hsize_t) :: dim_dat(2),dim_tmp(2)
    write(name0,'(A,A)')trim(name_fe),"_mfe2fd.txt"
    ! Read material data and rank match matrix
    open(247,file=adjustl(name0),status='old')
    read(247,*)nmat,cmat,nproc_fe,nproc_fe2fd
    call assert( nproc == nproc_fe2fd )
    allocate(mat(nmat,cmat),mattmp(nproc_fe,nproc_fe2fd))
    do i=1,nmat
      read(247,*)mat(i,:)
    end do
    do i=1,nproc_fe
      read(247,*)mattmp(i,:)
    end do
    close(247)
    do irank=0,nproc_fe-1
      if (mattmp(irank+1,myid+1)==1) then ! has on-rank pixel
        write(name1,'(A,A,I0.6,A)')trim(name_fe),"_fe2fd_",irank,".h5"
        call h5open_f(err)
        call h5fopen_f(trim(name1),H5F_ACC_RDWR_F,idfile,err)
        call h5dopen_f(idfile,"mat",iddat,err)
        call h5dget_space_f(iddat,spc_dat,err)
        call h5sget_simple_extent_dims_f(spc_dat,dim_dat,dim_tmp,err)
        nels=dim_dat(2)
        allocate(dat_mat(7,nels))
        call h5dread_f(iddat,h5t_native_integer,dat_mat,dim_dat,err)
        call h5dclose_f(iddat,err)
        call h5fclose_f(idfile,err)
        call h5close_f(err)
        do el=1,nels
          i0=max(dat_mat(1,el),ibeg);i1=min(dat_mat(2,el),iend)
          j0=max(dat_mat(3,el),jbeg);j1=min(dat_mat(4,el),jend)
          k0=max(dat_mat(5,el),kbeg);k1=min(dat_mat(6,el),kend)
          if (i0<=i1 .and. j0<=j1 .and. k0<=k1) then ! positive volume
            imat=dat_mat(7,el)
            rho0=mat(imat,5)
            select case(cmat)
            case(7)
              E=mat(imat,6); nu=mat(imat,7)
            case(11)
              E=mat(imat,10); nu=mat(imat,11)
            case(12)
              E=mat(imat,11); nu=mat(imat,12)
            end select
            lam0=E*nu/(1+nu)/(1-2*nu)
            mu0=.5d0*E/(1.d0+nu)
            do j = j0, j1
              do i = i0, i1
                do k = k0, k1
                  rho(k,i,j) = rho0*1.d-3 ! g/cm^3
                  mu (k,i,j) = mu0*1.d-9  ! g/cm^3)*(km/s)^2
                  lam(k,i,j) = lam0*1.d-9 ! g/cm^3)*(km/s)^2
                end do !k
              end do !i
            end do !j
          end if ! positive volume
        end do
        deallocate(dat_mat)
      end if ! has on-rank pixel
    end do
  end subroutine vmodel__fe

end module m_vmodel_fe
