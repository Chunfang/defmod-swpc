! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module galpha 

#include <petscversion.h>

  use global
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
  implicit none
#include "petscdef.h"
#else
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
#endif
  ! Generalized alpha parameters
  real(8) :: f1m,f2m,f3m,f1d,f2d,f3d,alpha_m,alpha_f,bt,gm 
  Mat :: Mat_Ka,Mat_D
  Vec :: Vec_V,Vec_A,Vec_Fa

contains

  ! Output initial Mat_K_dyn, Vec_Fa for generalized-alpha method
  subroutine AlphaInit
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    Mat :: Mat_tmp
    alpha_m=0.2 !0.2
    alpha_f=0.4 !0.4
    bt     =0.25*(f1-alpha_m+alpha_f)**2 !0.36
    gm     =0.5-alpha_m+alpha_f !0.7
    f1m=(f1-alpha_m)/bt/dt/dt
    f2m=(f1-alpha_m)/bt/dt
    f3m=(f1-alpha_m-f2*bt)/f2/bt
    f1d=(f1-alpha_f)*gm/bt/dt
    f2d=((f1-alpha_f)*gm-bt)/bt
    f3d=(f1-alpha_f)*(gm-f2*bt)/f2/bt*dt

    ! Create Mat/Vec spaces
    call VecDuplicate(Vec_U_dyn,Vec_V,ierr)
    call VecDuplicate(Vec_U_dyn,Vec_A,ierr)
    call VecZeroEntries(Vec_V,ierr)
    call VecZeroEntries(Vec_A,ierr)
    call VecDuplicate(Vec_F_dyn,Vec_Fa,ierr) 

    ! [K]-alpha
    call MatScale(Mat_M,f1m+alpha*f1d,ierr)
    call MatDuplicate(Mat_K_dyn,Mat_Copy_Values,Mat_tmp,ierr)
    call MatAYPX(Mat_tmp,f1-alpha_f+beta*f1d,Mat_M,Different_Nonzero_Pattern   &
       ,ierr) 
    call MatAXPY(Mat_tmp,f1d,Mat_D,Different_Nonzero_Pattern,ierr)
    call MatDuplicate(Mat_tmp,Mat_Copy_Values,Mat_Ka,ierr)

    ! Initial RHS
    !call MatCopy(Mat_K_dyn,Mat_tmp,Same_Nonzero_Pattern,ierr)
    !call MatAYPX(Mat_tmp,beta*f1d-alpha_f,Mat_M,Same_Nonzero_Pattern,ierr)
    !call MatAXPY(Mat_tmp,f1d,Mat_D,Different_Nonzero_Pattern,ierr)
    !call MatMult(Mat_tmp,Vec_U_dyn,Vec_Fa,ierr) 

    ! Add initial force (dyn=true)
    !call FormRHS
    !call VecAssemblyBegin(Vec_F_dyn,ierr)
    !call VecAssemblyEnd(Vec_F_dyn,ierr)
    !call VecCopy(Vec_F_dyn,Vec_Fa,ierr)
    
    ! Restore
    call MatScale(Mat_M,f1/(f1m+alpha*f1d),ierr)
    call MatDestroy(Mat_tmp,ierr) 
  end subroutine AlphaInit

  ! RHS of implicit dynamics 
  subroutine AlphaRHS
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    Mat :: Mat_tmp
    Vec,pointer :: Vec_tmp(:)
    real(8) :: fct(3)
    call MatDuplicate(Mat_K_dyn,Mat_Copy_Values,Mat_tmp,ierr)
    call VecDuplicateVecsF90(Vec_U_dyn,3,Vec_tmp,ierr)
    
    ! New force term (dyn=true)
    call VecCopy(Vec_F_dyn,Vec_Fm_dyn,ierr)
    call VecZeroEntries(Vec_F_dyn,ierr)
    call FormRHS
    call VecAssemblyBegin(Vec_F_dyn,ierr)
    call VecAssemblyEnd(Vec_F_dyn,ierr)

    ! Displacement term
    call MatScale(Mat_M,f1m+alpha*f1d,ierr)
    call MatAYPX(Mat_tmp,beta*f1d-alpha_f,Mat_M,Different_Nonzero_Pattern,ierr)
    call MatMult(Mat_tmp,Vec_U_dyn,Vec_tmp(1),ierr)

    ! Velocity term
    call MatCopy(Mat_K_dyn,Mat_tmp,Same_Nonzero_Pattern,ierr)
    call MatScale(Mat_M,(f2m+alpha*f2d)/(f1m+alpha*f1d),ierr)
    call MatAYPX(Mat_tmp,beta*f2d,Mat_M,Different_Nonzero_Pattern,ierr)
    call MatMult(Mat_tmp,Vec_V,Vec_tmp(2),ierr)
    call VecAXPY(Vec_tmp(1),f1,Vec_tmp(2),ierr)

    ! Acceleration term
    call MatCopy(Mat_K_dyn,Mat_tmp,Same_Nonzero_Pattern,ierr)
    call MatScale(Mat_M,(f3m+alpha*f3d)/(f2m+alpha*f2d),ierr)
    call MatAYPX(Mat_tmp,beta*f3d,Mat_M,Different_Nonzero_Pattern,ierr)
    call MatMult(Mat_tmp,Vec_A,Vec_tmp(2),ierr)
    call VecWAXPY(Vec_Fa,f1,Vec_tmp(1),Vec_tmp(2),ierr)

    ! Force terms at time n+1-alpha_f 
    call VecCopy(Vec_Fm_dyn,Vec_tmp(1),ierr)
    call VecCopy(vec_F_dyn,Vec_tmp(2),ierr)
    fct(:2)=(/alpha_f,f1-alpha_f/)
    call VecMAXPY(Vec_Fa,2,fct(:2),Vec_tmp(:2),ierr)

    ! Absorbing damping terms
    call MatMult(Mat_D,Vec_U_dyn,Vec_tmp(1),ierr)
    call MatMult(Mat_D,Vec_V,    Vec_tmp(2),ierr)
    call MatMult(Mat_D,Vec_A,    Vec_tmp(3),ierr)
    fct=(/f1d,f2d,f3d/)
    call VecMAXPY(Vec_Fa,3,fct,Vec_tmp,ierr)
    
    ! Restore 
    call MatScale(Mat_M,f1/(f3m+alpha*f3d),ierr)
    call MatDestroy(Mat_tmp,ierr)
    call VecDestroyVecsF90(3,Vec_tmp,ierr)
  end subroutine AlphaRHS

  ! Update velocity and acceleration
  subroutine AlphaUpdate
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    Vec,pointer :: Vec_tmp(:)
    real(8) :: fct(2)
    !real(8),target :: flt_dyn(n_lmnd*dmn)

    ! Solve and constrain
    call VecCopy(Vec_U_dyn,Vec_Um_dyn,ierr)
    call KSPSolve(Krylov,Vec_Fa,Vec_U_dyn,ierr)
    call AlphaCnstr

    ! Extract and output temporal fault slip
    call MatMult(Mat_G,Vec_U_dyn,Vec_Wlm(1),ierr)
    call VecGetArrayF90(Vec_Wlm(1),pntr,ierr)
    flt_slip=pntr
    tot_flt_slip=tot_flt_slip+flt_slip
    call VecRestoreArrayF90(Vec_Wlm(1),pntr,ierr)
    if (mod(tstep,frq)==0) then
       call GetSlip_dyn
       rslip=real(sum(slip))/real(size(slip))
       if (rank==0) print'(A11,I0,X,F0.2,A)'," Time Step ",tstep,rslip*100.0,  &
          "% fault slipping."
    end if

    ! Velocity update
    call VecDuplicateVecsF90(Vec_U_dyn,3,Vec_tmp,ierr)
    call VecWAXPY(Vec_tmp(1),-f1,Vec_Um_dyn,Vec_U_dyn,ierr)
    call VecCopy(Vec_V,Vec_tmp(2),ierr)
    call VecCopy(Vec_A,Vec_tmp(3),ierr) 
    call VecScale(Vec_tmp(3),dt*(f1-gm/f2/bt),ierr)
    fct=(/gm/bt/dt,f1-gm/bt/)
    call VecMAXPY(Vec_tmp(3),2,fct,Vec_tmp(:2),ierr) 

    ! Acceleration update
    call VecScale(Vec_A,f1-f1/f2/bt,ierr)
    fct=(/f1/bt/dt/dt,-f1/bt/dt/)
    call VecMAXPY(Vec_A,2,fct,Vec_tmp(:2),ierr)

    ! Cleanup
    call VecCopy(Vec_tmp(3),Vec_V,ierr)
    call VecDestroyVecsF90(3,Vec_tmp,ierr)
  end subroutine AlphaUpdate

  ! Constrain solution
  subroutine AlphaCnstr
    implicit none
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
#include "petsc.h"
#endif
    ! Form lambda=(GUp-Flm)/(dt^2GMinvGt)
    call MatMult(Mat_G,Vec_U_dyn,Vec_Wlm(1),ierr)
    call VecWAXPY(Vec_Wlm(2),-f1,Vec_I_dyn,Vec_Wlm(1),ierr)
    call MatMult(Mat_GMinvGt,Vec_Wlm(2),Vec_lambda,ierr)
    call CapLM_dyn
    call VecAXPY(Vec_lambda_tot,f1,Vec_lambda,ierr)
    ! Form Up=Up-dt^2(Minv(Gtlambda))
    call MatMult(Mat_Gt,Vec_lambda,Vec_W(1),ierr)
    call MatMult(Mat_Minv,Vec_W(1),Vec_W(2),ierr)
    ! FIXME Scale with 1/f1m, dt**2*beta/(1-alpha_m), instead of dt**2 seems to
    ! replicate explicit solver, why?
    !call VecAXPY(Vec_U_dyn,-dt**2,Vec_W(2),ierr)
    call VecAXPY(Vec_U_dyn,-f1/f1m,Vec_W(2),ierr)
  end subroutine AlphaCnstr 

end module galpha
