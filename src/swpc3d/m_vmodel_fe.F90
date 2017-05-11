! Overwrite the FD velocity with FE velocity 
module m_vmodel_fe
  
  !! dependency
  use m_std
  use m_global
  use m_debug
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
    integer :: i0,i1,j0,j1,k0,k1,i,j,k,el ! i,j,k regions
    integer,allocatable :: mattmp(:,:)
    real(SP),allocatable :: mat(:,:)
    real(SP) :: E,nu,lam0,mu0,rho0
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
          write(name1,'(A,A,I0.6,A)')trim(name_fe),"_",irank,"_ife2fd.txt"
          open(248,file=adjustl(name1),status='old')
          read(248,*)nels
          do el=1,nels
             read(248,*)i0,i1,j0,j1,k0,k1,imat
             i0=max(i0,ibeg);i1=min(i1,iend)
             j0=max(j0,jbeg);j1=min(j1,jend)
             k0=max(k0,kbeg);k1=min(k1,kend)
             if (i0<=i1 .and. j0<=j1 .and. k0<=k1) then ! positive volume
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
          close(248)
       end if ! has on-rank pixel
    end do
  end subroutine vmodel__fe
  
end module m_vmodel_fe
