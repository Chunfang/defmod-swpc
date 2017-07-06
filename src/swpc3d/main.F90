!! ----------------------------------------------------------------------------------------------------------------------------- !!
!>
!!  SWPC: Seismic Wave Propagation Code
!!
!! @detail
!!   This software simulate seismic wave propagation, by solving equations of motion with constitutive equations of
!!   elastic/visco-elastic medium by finite difference method (FDM).
!!
!! @copyright
!!   Copyright 2013-2017 Takuto Maeda. All rights reserved. This project is released under the MIT license.
!!
!! @see
!!   - Noguchi et al.     (2016) GJI     doi:10.1093/gji/ggw074
!!   - Maeda et al.       (2014) EPSL    doi:10.1016/j.epsl.2014.04.037
!!   - Maeda et al.       (2013) BSSA    doi:10.1785/0120120118
!!   - Maeda and Furumura (2013) PAGEOPH doi:10.1007/s00024-011-0430-z
!!   - Noguchi et al.     (2013) PAGEOPH doi:10.1007/s00024-011-0412-1
!!   - Furumura et al.    (2008) PAGEOPH doi:10.1007/s00024-008-0318-8
!!   - Furumura ad Chen   (2005) PARCO   doi:10.1016/j.parco.2005.02.003
!<
!! ----------------------------------------------------------------------------------------------------------------------------- !!
#include "m_debug.h"
program SWPC_3D

  !! -- Dependency
  use m_std
  use m_debug
  use m_global
  use m_kernel
  use m_getopt
  use m_source
  use m_vmodel_fe
  use m_rup
  use m_medium
  use m_report
  use m_pwatch
  use m_output
  use m_absorb
  use m_ckprst
  use m_green
  use m_readini
  use mpi

  !! -- Declarations
  implicit none
  !! --

  character(256) :: fn_prm, str_eid
  integer :: it
  integer :: ierr
  logical :: is_opt
  logical :: stopwatch_mode
  integer :: io_prm, io_watch
  integer :: it0
  !! ----

  !!
  !! option processing
  !!
  call getopt('i', is_opt, fn_prm, './in/input.inf' )
  call std__getio( io_prm )
  open( io_prm, file=trim(fn_prm), action='read', status='old' )
  call readini( io_prm, 'stopwatch_mode', stopwatch_mode, .true.  )
    
  !!
  !! Point/rupture source
  !!
  call getopt('r', is_opt, name_fe, 'RS' )
  have_rs=.false.; have_ps=.false.
  if (name_fe/='RS' .and. is_opt) have_rs=.true.
  if (.not. have_rs) have_ps=.true.
  
  !!
  !! Select event, see XXX_dyn.log for event numbers  
  !!
  if (have_rs) then
    call getopt('e', is_opt,str_eid, 'ID' )
    if (str_eid/='ID' .and. is_opt) then
      str_eid=trim(str_eid)
      read(str_eid,'(I3.0)')eid ! max 999 events
    else
      eid=1
    end if
  end if

  !!
  !! Launch MPI process
  !!
  call mpi_init( ierr )

  !!
  !! Read control parameters
  !!
  call global__setup( io_prm )

  !!
  !! check if there are restart files
  !!
  call ckprst__setup( io_prm )
  call ckprst__restart( it0 )

  !!
  !! setup for new simulation
  !!
  if( it0 == 1 ) then

    !!
    !! stopwatch start
    !!
    call pwatch__setup( stopwatch_mode )

    !!
    !! set-up each module
    !!
    call global__setup2( )
    call medium__setup( io_prm )

     ! FE vmodel overwrite
     call vmodel__fe
     ! wait until deallocation by medium__setup
     call mpi_barrier( mpi_comm_world, ierr ) 

     call kernel__setup( )
     if (have_ps) then 
       call source__setup( io_prm )
     else if (have_rs) then ! rupture source no rescale
       UC=1D0; M0=1D0
     end if 

     ! Rupture source setup
     if (have_rs) then
       call rup__Init
       call rup__getTime
       call rup__getFE2FD 
       call rup__getObsFD
       call mpi_barrier( mpi_comm_world, ierr )
     end if

     call absorb__setup( io_prm )
     call output__setup( io_prm )
     call green__setup( io_prm )
     call report__setup( io_prm )
  end if

  close( io_prm )

#ifdef _FX
  call fipp_start()  !! performance measurement
#endif
  !! mainloop
  do it = it0, nt

    call report__progress(it)

    call green__store( it )
    call output__store_wav ( it )
    call output__write_snap( it )

    call kernel__update_stress()
    call absorb__update_stress()

    if (have_ps) call source__stressdrop(it)
    call global__comm_stress()

    call kernel__update_vel()
    call absorb__update_vel()

    if (have_ps) call source__bodyforce(it)

     ! Set rupture source
     if (have_rs) then
       call rup__setSrc(it) 
       call rup__evalObsFD(it)
       call rup__writeObsFD(it)
     end if 

    call green__source( it )

    call global__comm_vel()

    call ckprst__checkpoint( it )


  end do
#ifdef _FX
  call fipp_stop()  !! performance measurement
#endif

  call green__export()
  call output__export_wav()
  call output__closefiles()

  !! ending message
  call report__terminate()

  !!
  !! Checkpoint file delete
  !!
  call ckprst__filedelete()

  !!
  !! stopwatch report from 0-th node
  !!
  if( stopwatch_mode ) then
    if( myid == 0 ) then
      call std__getio( io_watch )
      open( io_watch, file=trim( odir ) // '/' // trim( title ) //  '.tim', action='write', status='unknown' )
    end if
    call pwatch__report( io_watch, 0 )
  end if

  !!
  !! Program termination
  !!
  call mpi_finalize( ierr )

end program SWPC_3D

!! ----------------------------------------------------------------------------------------------------------------------------- !!
