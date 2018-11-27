! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

program main

#include <petscversion.h>

  use global
  use galpha
  use fefd 
  use fvfe
  use h5io
  use esh3d
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7 && PETSC_VERSION_SUBMINOR<5)
  implicit none
#include "petsc.h"
#else
#include <petsc/finclude/petscksp.h>
  use petscksp
  implicit none
#endif
  character(256) :: input_file,viz,fd,fv
  logical :: l,v,w,pbc
  integer,pointer :: null_i=>null()
  real(8),pointer :: null_r=>null()
  real(8) :: fdt
  integer :: i,j,j1,j2,j3,j4,j5,j6,n,n_dyn,nodal_bw,ef_eldof,ng,vout,fdout
  
  call PetscInitialize(Petsc_Null_Character,ierr)

  call MPI_Comm_Rank(MPI_Comm_World,rank,ierr)
  call MPI_Comm_Size(MPI_Comm_World,nprcs,ierr)

#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=6)
  call PetscOptionsGetString(Petsc_Null_Character,'-f',input_file,l,ierr)
  call PetscOptionsGetString(Petsc_Null_Character,'-ss',viz,v,ierr)
  call PetscOptionsGetString(Petsc_Null_Character,'-fd',fd,w,ierr)
  call PetscOptionsGetString(Petsc_Null_Character,'-fv',fv,pbc,ierr)
#elif (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR==7 && PETSC_VERSION_SUBMINOR<5)
  call PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-f',      &
     input_file,l,ierr)
  call PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-ss',     &
     viz,v,ierr)
  call PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-fd',     &
     fd,w,ierr)
  call PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-fv',     &
     fv,pbc,ierr)
#else
  call PetscOptionsGetString(Petsc_Null_Options,Petsc_Null_Character,'-f',     &
     input_file,l,ierr)
  call PetscOptionsGetString(Petsc_Null_Options,Petsc_Null_Character,'-ss',    &
     viz,v,ierr)
  call PetscOptionsGetString(Petsc_Null_Options,Petsc_Null_Character,'-fd',    &
     fd,w,ierr)
  call PetscOptionsGetString(Petsc_Null_Options,Petsc_Null_Character,'-fv',    &
     fv,pbc,ierr)
#endif
  if (.not. l) then
     call PrintMsg("Usage: [mpiexec -n <np>] defmod -f <input_filename>")
     go to 9
  end if
  vout=0
  if (v) then
     read (viz,'(I1.0)')vout 
     if (vout==1) then 
        call PrintMsg("Snapshot output will slow the run!")
     else
        call PrintMsg("Shosen NOT to output snapshot.")
     end if
  else
     call PrintMsg("Use -ss 1 to turn on the snapshot output.")
  end if
  fdout=0
  if (w) then
     read (fd,'(I1.0)')fdout
     if (fdout==1) call PrintMsg("Running in FE-FD mixed mode.")
  else
     call PrintMsg("Use -fd 1 to turn on FE-FD mixd mode.")
  end if

  ! Read input file parameters
  open(10,file=input_file,status='old')

  ! Output file name is the same as the input
  if (index(input_file,".inp")>0) then
     output_file=input_file(1:index(input_file,".inp")-1)
  else
     output_file=input_file
  end if

  call PrintMsg("Reading input ...")
  call ReadParameters
  ! Bounding box pressure (from pflotran) option
  if (poro) then
     fvin=0
     if (pbc) then
        read (fv,'(I1.0)')fvin
        if (fvin/=0) call PrintMsg("FV model expected")
     else
        call PrintMsg("Use -fv 1 (or 2) to load pflotran model.")
     end if
  end if
  ! Set element specific constants
  call InitializeElement
  p=0; ef_eldof=eldof
  if (poro) then
     p=1; ef_eldof=eldof+eldofp
     if (eltype=="tri") nip=3; if (eltype=="tet") nip=4
  elseif (fault .or. galf) then
     if (eltype=="tri") nip=3; if (eltype=="tet") nip=4
  end if

  ! Partition mesh using METIS, create mappings, and read on-rank mesh data
  allocate(npart(nnds),epart(nels)); epart=0; npart=0
  if (nprcs>1) then
     call PrintMsg("Partitioning mesh ...")
     if (rank==0) then
        allocate(nodes(1,npel*nels),work(nels+1)); work(1)=0
        do i=1,nels
           j=npel*(i-1)+1; n=npel*i; read(10,*)nodes(1,j:n); work(i+1)=n
        end do
        nodes=nodes-1
        call METIS_PartMeshNodal(nels,nnds,work,nodes,null_i,null_i,nprcs,     &
           null_r,null_i,n,epart,npart)
        deallocate(nodes,work)
        rewind(10); call ReadParameters
     end if
     call MPI_Bcast(npart,nnds,MPI_Integer,0,MPI_Comm_World,ierr)
     call MPI_Bcast(epart,nels,MPI_Integer,0,MPI_Comm_World,ierr)
  end if
  call PrintMsg("Reading mesh data ...")
  allocate(emap(nels),nmap(nnds)); emap=0; nmap=0
  ! Create original to local element mappings and read on-rank element data
  j=1
  do i=1,nels
     if (epart(i)==rank) then
        epart(i)=1; emap(i)=j; j=j+1
     else
        epart(i)=0
     end if
  end do
  n=sum(epart); allocate(nodes(n,npel),id(n)) ! id -> mtrl flag
  j=1
  do i=1,nels
     if (epart(i)==1) then
        read(10,*)nodes(j,:),id(j); j=j+1
     else
        read(10,*)val
     end if
  end do
  nels=n
  ! Create original to global nodal mappings and read on-rank + ghost nodes
  allocate(work(0:nprcs-1))
  j=0
  do i=1,nnds
     if (npart(i)==rank) j=j+1
  end do
  call MPI_AllGather(j,1,MPI_Integer,work,1,MPI_Integer,MPI_Comm_World,ierr)
  if (rank==0) then
     n=1
  else
     n=sum(work(0:rank-1))+1
  end if
  do i=1,nnds
     if (npart(i)==rank) then
        nmap(i)=n; n=n+1
     end if
  end do
  deallocate(work)
  allocate(work(nnds))
  call MPI_AllReduce(nmap,work,nnds,MPI_Integer,MPI_Sum,MPI_Comm_World,ierr)
  nmap=work
  npart=0
  do i=1,nels
     do j=1,npel
        npart(nodes(i,j))=1
     end do
  end do
  j=1; work=0
  do i=1,nnds
     if (npart(i)==1) then
        work(i)=j; j=j+1
     end if
  end do
  n=sum(npart); allocate(coords(n,dmn),bc(n,dmn+p))
  j=1
  do i=1,nnds
     if (npart(i)==1) then
        read(10,*)coords(j,:),bc(j,:)
        j=j+1
     else
        read(10,*)val
     end if
  end do
  coords=km2m*coords
  ! Re-number on-rank nodes and create local to global node and dof mappings
  do i=1,nels
     do j=1,npel
        nodes(i,j)=work(nodes(i,j))
     end do
  end do
  n=sum(npart); allocate(nl2g(n,2),indxmap((dmn+p)*n,2))
  if (fault .or. galf) allocate(indxmap_u(dmn*n,2))
  j=1
  do i=1,nnds
     if (work(i)==j) then
        nl2g(j,1)=j; nl2g(j,2)=nmap(i); j=j+1
     end if
  end do
  do i=1,n
     do j=1,dmn+p
        indxmap((dmn+p)*i-j+1,:)=(dmn+p)*nl2g(i,:)-j ! 0 based index
     end do
  end do
  if (fault .or. galf) then
     do i=1,n
        do j=1,dmn
           indxmap_u(dmn*i-j+1,:)=dmn*nl2g(i,:)-j
        end do
     end do
  end if
  deallocate(work)
  ! Read material data assuming nels >> nmts, i.e., all ranks store all data
  if (fault .or. galf) then
     allocate(mat(nmts,5+4*p+init+2+2*rve))
  else
     allocate(mat(nmts,5+4*p+2*rve))
  end if
  do i=1,nmts
     read(10,*)mat(i,:)
  end do
  if (rve>0 .and. eltype=="hex") then ! Only for 8-node hex
     allocate(incls(nincl,11)) 
     do i=1,nincl
        read(10,*)incls(i,:) 
     end do
     allocate(matDstr(nmts,8,6,6))
     call RVEHex8Dstr(mat,incls,matDstr)
     ! Inspect Dstr(l,k,:,:) for l-th material at k-th Gauss point
     !do i=1,6
     !   print'(6(F0.6X))',matDstr(1,5,i,:) 
     !end do
  end if
  deallocate(epart,npart)

  ! Initialize local element variables and global U
  allocate(ipoint(nip,dmn),weight(nip),k(ef_eldof,ef_eldof),m(eldof,eldof),    &
     f(ef_eldof),indx(ef_eldof),enodes(npel),ecoords(npel,dmn),vvec(dmn+p))
  if (fault .or. galf) allocate(k_dyn(eldof,eldof),indx_dyn(eldof))
  call SamPts(ipoint,weight)
  n=(dmn+p)*nnds; if (.not. explicit) n=n+nceqs
  call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n,Vec_U,ierr)
  ! Global U for dynamic run
  if (fault .or. galf) then
     n_dyn=dmn*nnds
     call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n_dyn,Vec_U_dyn,ierr)
  end if
  if (visco .or. ((fault .or. galf) .and. lm_str))                             &
     allocate(stress(nels,nip,cdmn))

  ! Set scaling constants
  wt=f2*exp((log(maxval(mat(:,1)))+log(minval(mat(:,1))))/f2)
  if (dmn==3) wt=wt*km2m
  if (poro) then
     scale=sqrt(exp((log(maxval(mat(:,1)))+log(minval(mat(:,1))))/f2)/         &
                exp((log(maxval(mat(:,6)))+log(minval(mat(:,6))))/f2)/         &
                dt/km2m)
  else
     scale=f1 
  end if

  ! Form stiffness matrix
  call PrintMsg("Forming [K] ...")
  nodal_bw=(dmn+p)*(nodal_bw+1)
  call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,nodal_bw,   &
     Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K,ierr)
  call MatSetOption(Mat_K,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)
  if (fault .or. galf) then ! Dynamic K
     nodal_bw=nodal_bw/(dmn+p)*dmn
     call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n_dyn,n_dyn, &
        nodal_bw,Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K_dyn,ierr)
     call MatSetOption(Mat_K_dyn,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
        ierr)
  end if

  kfv=.false. ! Full matrix without FV Bc
  do i=1,nels
     if (rve>0 .and. eltype=="hex") then  
        call FormLocalK(i,k,indx,"Ke",matDstr)
     else 
        if (poro) then
           call FormLocalK(i,k,indx,"Kp")
        else
           call FormLocalK(i,k,indx,"Ke")
        end if
     end if
     indx=indxmap(indx,2)
     call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)
     if (fault .or. galf) then
        dyn=.true.
        if (rve>0 .and. eltype=="hex") then
           call FormLocalK(i,k_dyn,indx_dyn,"Ke",matDstr)
        else
           call FormLocalK(i,k_dyn,indx_dyn,"Ke")
        end if
        indx_dyn=indxmap_u(indx_dyn,2)
        call MatSetValues(Mat_K_dyn,eldof,indx_dyn,eldof,indx_dyn,k_dyn,       &
           Add_Values,ierr)
        dyn=.false.
     end if
  end do

  ! Initialize and form mass matrix and its inverse
  if (explicit .or. galf .or. fault) then
     call PrintMsg("Forming [M] & [M]^-1 ...")
     if (fault .or. galf) then
        call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n_dyn,    &
           n_dyn,3,Petsc_Null_Integer,3,Petsc_Null_Integer,Mat_M,ierr)
     elseif (gf) then
        call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,3,    &
           Petsc_Null_Integer,3,Petsc_Null_Integer,Mat_M,ierr)
     else
        call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,1,    &
           Petsc_Null_Integer,0,Petsc_Null_Integer,Mat_M,ierr)
     end if
     call MatSetOption(Mat_M,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)
     do i=1,nels
        if (fault .or. galf) then
           call FormLocalM(i,m,indx_dyn)
           indx_dyn=indxmap_u(indx_dyn,2)
           do j=1,eldof
              val=m(j,j)
              call MatSetValue(Mat_M,indx_dyn(j),indx_dyn(j),val,Add_Values,   &
                 ierr)
           end do
        else
           call FormLocalM(i,m,indx)
           indx=indxmap(indx,2)
           do j=1,eldof
              val=m(j,j)
              call MatSetValue(Mat_M,indx(j),indx(j),val,Add_Values,ierr)
           end do
        end if
     end do
     call MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
     call MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)
     call MatDuplicate(Mat_M,Mat_Do_Not_Copy_Values,Mat_Minv,ierr)
     if (fault .or. galf) then
        call MatGetDiagonal(Mat_M,Vec_U_dyn,ierr) ! Vec_U_dyn -> work vector
        call VecReciprocal(Vec_U_dyn,ierr)
        call MatDiagonalSet(Mat_Minv,Vec_U_dyn,Insert_Values,ierr)
        call VecZeroEntries(Vec_U_dyn,ierr)
     else
        call MatGetDiagonal(Mat_M,Vec_U,ierr) ! Vec_U is used as a work vector
        call VecReciprocal(Vec_U,ierr)
        call MatDiagonalSet(Mat_Minv,Vec_U,Insert_Values,ierr)
        call VecZeroEntries(Vec_U,ierr)
     end if
     if (.not. galf) call MatScale(Mat_M,alpha,ierr) ! M -> alpha x M
  end if

  ! Allocate arrays to store loading history
  allocate(cval(nceqs,3),fnode(nfrcs),fval(nfrcs,dmn+p+2),telsd(ntrcs,2),      &
     tval(ntrcs,dmn+p+2)); cval=f0; fval=f0; tval=f0

  ! Account for constraint eqn's
  if (nceqs>0) then
     call PrintMsg("Applying constraints ...")
     if (explicit) then
        call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,        &
           nceqs,3,Petsc_Null_Integer,3,Petsc_Null_Integer,Mat_Gt,ierr)
        call MatSetOption(Mat_Gt,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
           ierr)
     elseif ((fault .or. galf) .and. nceqs-nceqs_ncf>0 .and. hyb>0) then
        call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,nceqs_ncf/(dmn+p)+nfnd,&
           Vec_lm_pn,ierr)
        call VecGetLocalSize(Vec_lm_pn,n_lmnd,ierr)
        call VecGetOwnershipRange(Vec_lm_pn,lmnd0,j,ierr)
        call VecDuplicate(Vec_lm_pn,Vec_lm_f2s,ierr)
        if (poro) then
           call VecDuplicate(Vec_lm_pn,Vec_lm_pp,ierr)
        else
           call VecDestroy(Vec_lm_pn,ierr)
        end if
        ! Dofs of one fault node are not split by different ranks
        call VecCreateMPI(Petsc_Comm_World,n_lmnd*dmn,(nceqs_ncf/(dmn+p)+nfnd)*&
           dmn,Vec_lambda_sta,ierr)
        call VecDuplicate(Vec_lambda_sta,Vec_lambda_sta0,ierr)
        call VecZeroEntries(Vec_lambda_sta,ierr)
        call VecZeroEntries(Vec_lambda_sta0,ierr)
        if (rsf==1) call VecDuplicate(Vec_lambda_sta,Vec_lambda_bk,ierr)
        call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,n_lmnd*dmn,dmn*nnds,   &
           (nceqs_ncf/(dmn+p)+nfnd)*dmn,5,Petsc_Null_Integer,5,                &
           Petsc_Null_Integer,Mat_Gt,ierr)
        call MatSetOption(Mat_Gt,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
           ierr)
        call MatZeroEntries(Mat_Gt,ierr)
        allocate(flt_slip(n_lmnd*dmn),tot_flt_slip(n_lmnd*dmn),                &
           res_flt_slip(n_lmnd*dmn),qs_flt_slip(n_lmnd*dmn))
        qs_flt_slip=f0; tot_flt_slip=f0; res_flt_slip=f0; flt_slip=f0
     end if
     if (rank==0) open(15,file=trim(output_file(:index(output_file,"/",        &
        BACK=.TRUE.)))//"cnstrns.tmp",status="replace")
     do i=1,nceqs
        read(10,*)n
        if (rank==0) write(15,*)n
        do j=1,n
           read(10,*)vvec,node; node=nmap(node)
           if (rank==0) write(15,*)real(vvec),node
        end do
        read(10,*)cval(i,:)
        if (poro) then
           if (vvec(dmn+1)/=f0) cval(i,1)=cval(i,1)/scale
        end if
        if (rank==0) write(15,*)real(cval(i,:))
     end do
     if (rank==0) close(15)

     ! Read fault orientation vectors and friction parameters
     if (fault .or. gf .or. galf) then
        allocate(node_pos(nfnd),node_neg(nfnd),vecf(nfnd,dmn*dmn),fc(nfnd),    &
           fcd(nfnd),dc(nfnd),perm(nfnd),st_init(nfnd,dmn),xfnd(nfnd,dmn),     &
           frc(nfnd),coh(nfnd),dcoh(nfnd),biot(nfnd))    
        ! Have frictional fault intersection 
        ! xpair global fault node id paired with intersection
        ! vecxf alternative fault (strike,dip,normal) vectors at intersection 
        ! XFltMap map global fault id to intersection id, 0 off-intersection
        if (Xflt>1) allocate(XFltMap(nfnd),xpair(nxflt),vecxf(nxflt,dmn*dmn),  &
           stx_init(nxflt,dmn),xlock(nxflt))
        biot=f1
        if (rsf==1) allocate(rsfb0(nfnd),rsfV0(nfnd),rsfdtau0(nfnd),rsfa(nfnd),&
            rsfb(nfnd),rsfL(nfnd),rsftheta(nfnd))
        do i=1,nfnd
           if (poro .and. init==1) then
              if (rsf==1) then
                 read(10,*) node_pos(i),node_neg(i),vecf(i,:),rsfb0(i),        &
                    rsfV0(i),rsfdtau0(i),rsfa(i),rsfb(i),rsfL(i),rsftheta(i),  &
                    perm(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i)
              else
                 read(10,*) node_pos(i),node_neg(i),vecf(i,:),fc(i),fcd(i),    &
                    dc(i),perm(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i)
              end if
           else if (poro) then
              if (rsf==1) then
                 read(10,*) node_pos(i),node_neg(i),vecf(i,:),rsfb0(i),        &
                    rsfV0(i),rsfdtau0(i),rsfa(i),rsfb(i),rsfL(i),rsftheta(i),  &
                    perm(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i),biot(i)
              else
                 read(10,*) node_pos(i),node_neg(i),vecf(i,:),fc(i),fcd(i),    &
                    dc(i),perm(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i),&
                    biot(i)
              end if
           else
              if (rsf==1) then
                 read(10,*) node_pos(i),node_neg(i),vecf(i,:),rsfb0(i),        &
                    rsfV0(i),rsfdtau0(i),rsfa(i),rsfb(i),rsfL(i),rsftheta(i),  &
                    st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i)
              else
                read(10,*) node_pos(i),node_neg(i),vecf(i,:),fc(i),fcd(i),     &
                   dc(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),dcoh(i)
              end if
           end if
           node_pos(i)=nmap(node_pos(i)); node_neg(i)=nmap(node_neg(i))
        end do
        if (Xflt>1) then
           do i=1,nxflt
              read(10,*) xpair(i),vecxf(i,:),stx_init(i,:)
           end do
        end if
     end if ! Has frictional fault
     if (rank==0) call ApplyConstraints
     if ((explicit .or. fault .or. galf) .and. nceqs>0) then
        call MatAssemblyBegin(Mat_Gt,Mat_Final_Assembly,ierr)
        call MatAssemblyEnd(Mat_Gt,Mat_Final_Assembly,ierr)
     end if
  end if ! Has constraint equations

  call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
  call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
  if ((fault .and. lm_str) .or. galf) then
     call MatAssemblyBegin(Mat_K_dyn,Mat_Final_Assembly,ierr)
     call MatAssemblyEnd(Mat_K_dyn,Mat_Final_Assembly,ierr)
  end if

  ! Form RHS
  call PrintMsg("Forming RHS ...")
  call VecDuplicate(Vec_U,Vec_F,ierr)
  tstep=0
  do i=1,nfrcs
     read(10,*)fnode(i),fval(i,:)
  end do
  do i=1,ntrcs
     read(10,*)telsd(i,:),tval(i,:)
  end do
  fnode=nmap(fnode); telsd(:,1)=emap(telsd(:,1)) ! Remap nodes/els
  call FormRHS
  call VecAssemblyBegin(Vec_F,ierr)
  call VecAssemblyEnd(Vec_F,ierr)

  ! Observation solution space
  if (nobs>0) then
     allocate(ocoord(nobs,dmn))
     do i=1,nobs
        read(10,*) ocoord(i,:)
     end do
     call GetObsNd("ob")
     if (nobs_loc>0) then
        allocate(uu_obs(nobs_loc,dmn+p),tot_uu_obs(nobs_loc,dmn+p),            &
           uu_dyn_obs(nobs_loc,dmn),tot_uu_dyn_obs(nobs_loc,dmn))
        tot_uu_obs=f0; tot_uu_dyn_obs=f0
     end if
     n_log_dyn=0
  end if
  
  ! FD domain grid bocks containing fault nodes, xgp, idgp(_loc), gpnlst, 
  ! gpshape
  if (nceqs-nceqs_ncf>0 .and. fdout==1) then
     call PrintMsg("Locating FD grid points ...")
     call FDInit
     call GetFDFnd 
     call GetObsNd("fd")
     !deallocate(xgp,idgp)  
     deallocate(xgp)
     if (ngp_loc>0) allocate(uu_fd(ngp_loc,dmn))
     call NndFE2FD
     allocate(matFD(nels,dmn*2+1))
     call MatFE2FD
     call Write_fd("mat")
     deallocate(matFD)
  end if

  ! Account for absorbing boundaries (FormLocalAbsC1 for arbitrary boundary)
  if (explicit .or. ((fault .and. lm_str) .or. galf)) then
     call PrintMsg("Absorbing boundary ...")
     if (galf) then 
        call MatDuplicate(Mat_M,Mat_Do_Not_Copy_Values,Mat_D,ierr)
        call MatZeroEntries(Mat_D,ierr)
     end if
     do i=1,nabcs
        ! For axis aligned absorbing boundaries
        !read(10,*)el,side,j; el=emap(el)
        ! j (dir ID) has no use for arbitrarily facing AbsC 
        read(10,*)el,side; el=emap(el)
        if (el/=0) then
           if (fault) then
              !call FormLocalAbsC(el,side,abs(j),m,indx_dyn)
              call FormLocalAbsC1(el,side,m,indx_dyn)
              indx_dyn=indxmap_u(indx_dyn,2)
              !do j=1,eldof
              !   val=m(j,j)
              !   call MatSetValue(Mat_M,indx_dyn(j),indx_dyn(j),val,           &
              !      Add_Values,ierr)
              !end do
              ! [C] being sparse, bw=3, instead of 
              !call MatSetValues(Mat_M,eldof,indx_dyn,eldof,indx_dyn,m,         &
              !   Add_Values,ierr)
              ! , follow 
              do j1=1,eldof
                 do j2=1,eldof
                    val=m(j1,j2)
                    if (abs(val)>f0) call MatSetValue(Mat_M,indx_dyn(j1),      &
                       indx_dyn(j2),val,Add_Values,ierr)
                 end do
              end do
           elseif (galf) then
              call FormLocalAbsC1(el,side,m,indx_dyn)
              indx_dyn=indxmap_u(indx_dyn,2)
              ! [D] being sparse ...
              !call MatSetValues(Mat_D,eldof,indx_dyn,eldof,indx_dyn,m,         &
              !   Add_Values,ierr)
              do j1=1,eldof
                 do j2=1,eldof
                    val=m(j1,j2)
                    if (abs(val)>f0) call MatSetValue(Mat_D,indx_dyn(j1),      &
                       indx_dyn(j2),val,Add_Values,ierr)
                 end do
              end do
           else
              !call FormLocalAbsC(el,side,abs(j),m,indx)
              call FormLocalAbsC1(el,side,m,indx)
              indx=indxmap(indx,2)
              !do j=1,eldof
              !   val=m(j,j)
              !   call MatSetValue(Mat_M,indx(j),indx(j),val,Add_Values,ierr)
              !end do
              ! [C] being sparse ...
              !call MatSetValues(Mat_M,eldof,indx,eldof,indx,m,Add_Values,ierr)
              do j1=1,eldof
                 do j2=1,eldof
                    val=m(j1,j2)
                    if (abs(val)>f0) call MatSetValue(Mat_M,indx(j1),indx(j2), &
                       val,Add_Values,ierr)
                 end do
              end do
           end if
        end if
     end do
     if (galf) then
        call MatAssemblyBegin(Mat_D,Mat_Final_Assembly,ierr)
        call MatAssemblyEnd(Mat_D,Mat_Final_Assembly,ierr)
     else
        call MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
        call MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)
     end if
  end if
  close(10)
  if (fvin>2) call MakeEl2g
  deallocate(nmap,emap) ! End of input reading

  ! Initialize arrays to communicate ghost node values
  call PrintMsg("Setting up solver ...")
  j=size(indxmap,1)
  call VecCreateSeq(Petsc_Comm_Self,j,Seq_U,ierr)
  call ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,2),Petsc_Copy_Values,From,  &
     ierr)
  call ISCreateGeneral(Petsc_Comm_Self,j,indxmap(:,1),Petsc_Copy_Values,To,    &
     ierr)
  call VecScatterCreate(Vec_U,From,Seq_U,To,Scatter,ierr)
  allocate(uu(j),tot_uu(j)); tot_uu=f0
  if (fault .or. galf) then
     if (lm_str) then
        j=size(indxmap_u,1)
        call VecCreateSeq(Petsc_Comm_Self,j,Seq_SS,ierr)
        call VecCreateSeq(Petsc_Comm_Self,j,Seq_SH,ierr)
        if (nceqs>0) call VecCreateSeq(Petsc_Comm_Self,j,Seq_f2s,ierr)
        call VecCreateSeq(Petsc_Comm_Self,j,Seq_dip,ierr)
        call VecCreateSeq(Petsc_Comm_Self,j,Seq_nrm,ierr)
     end if
     j=size(indxmap_u,1)
     call VecCreateSeq(Petsc_Comm_Self,j,Seq_U_dyn,ierr)
     call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),Petsc_Copy_Values,  &
        From,ierr)
     call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),Petsc_Copy_Values,  &
        To,ierr)
     call VecScatterCreate(Vec_U_dyn,From,Seq_U_dyn,To,Scatter_dyn,ierr)
     allocate(uu_dyn(j),tot_uu_dyn(j)) 
     if (lm_str) then
        allocate(ss(j),sh(j),f2s(j),dip(j),nrm(j))
     end if
  end if

  ! Implicit Solver
  if (.not. (explicit .or. fault .or. galf)) then
     call KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     call KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     call SetupKSPSolver
     call PrintMsg("Solving ...")
     call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
     call GetVec_U; tot_uu=tot_uu+uu
     if (nobs_loc>0) then
        call GetVec_obs
        tot_uu_obs=tot_uu_obs+uu_obs
     end if
     if (visco) then
        ! Recover stress
        call PrintMsg("Recovering stress ...")
        do i=1,nels
           call RecoverStress(i,stress)
        end do
     end if
     ! Write output
     call WriteOutput
     ! Prepare for time stepping
     if (t>f0 .and. dt>f0 .and. t>=dt) then
        if (poro) then
           ! Form Kc and Up
           call VecDuplicate(Vec_U,Vec_Um,ierr) ! U->du & Um->u
           call VecCopy(Vec_U,Vec_Um,ierr)
           call VecGetLocalSize(Vec_U,j,ierr)
           call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
           j2=0
           do i=1,j
              if (mod(j1+i,dmn+1)==0 .and. j1+i-1<(dmn+1)*nnds) then
                 j2=j2+1
              end if
           end do
           allocate(work(j2))
           j2=0
           do i=1,j
              if (mod(j1+i,dmn+1)==0 .and. j1+i-1<(dmn+1)*nnds) then
                 j2=j2+1
                 work(j2)=j1+i-1
              end if
           end do
           j=size(work)
           call ISCreateGeneral(Petsc_Comm_World,j,work,Petsc_Copy_Values,     &
              RI,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
           call MatGetSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
#else
           call MatCreateSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
#endif
           call MatZeroEntries(Mat_Kc,ierr)
           allocate(kc(eldofp,eldofp),indxp(eldofp),Hs(eldofp))
           do i=1,nels
              call FormLocalK(i,k,indx,"Kc")
              kc=k(eldof+1:,eldof+1:)
              indxp=indx(eldof+1:); indxp=indxmap(indxp,2)
              indxp=(indxp+1)/(dmn+1)-1
              call MatSetValues(Mat_Kc,eldofp,indxp,eldofp,indxp,kc,           &
                 Add_Values,ierr)
           end do
           call MatAssemblyBegin(Mat_Kc,Mat_Final_Assembly,ierr)
           call MatAssemblyEnd(Mat_Kc,Mat_Final_Assembly,ierr)
           call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
           call VecDuplicate(Vec_Up,Vec_I,ierr) ! I->KcUp
           call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
           allocate(uup(size(work)))
        end if
        steps=int(ceiling(t/dt))
        ! Start time stepping
        do tstep=1,steps
           if (rank==0) print'(A11,I0)'," Time Step ",tstep
           ! Reform stiffness matrix, if needed
           if (visco .and. (tstep==1 .or. maxval(mat(:,4))>f1)) then
              call PrintMsg(" Reforming [K] ...")
              call MatZeroEntries(Mat_K,ierr)
              do i=1,nels
                 call FormLocalK(i,k,indx,"Kv")
                 indx=indxmap(indx,2)
                 call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,        &
                    Add_Values,ierr)
              end do
              ! Account for constraint eqn's
              if (rank==0 .and. nceqs>0) call ApplyConstraints
              call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
              call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
           end if
           ! Reform RHS
           call VecZeroEntries(Vec_F,ierr)
           call PrintMsg(" Reforming RHS ...")
           if (poro) then
              call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
              call MatMult(Mat_Kc,Vec_Up,Vec_I,ierr)
              call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
              call VecScale(Vec_I,-dt,ierr)
              call VecGetArrayF90(Vec_I,pntr,ierr)
              uup=pntr
              call VecRestoreArrayF90(Vec_I,pntr,ierr)
              j=size(uup)
              call VecSetValues(Vec_F,j,work,uup,Add_Values,ierr)
              ! Stabilize RHS
              do i=1,nels
                 call FormLocalHs(i,Hs,indxp)
                 indxp=indxmap(indxp,2)
                 call VecSetValues(Vec_F,npel,indxp,Hs,Add_Values,ierr)
              end do
           end if
           if (visco) then
              do i=1,nels
                 call ReformLocalRHS(i,f,indx)
                 indx=indxmap(indx,2)
                 call VecSetValues(Vec_F,eldof,indx,f,Add_Values,ierr)
              end do
           end if
           call FormRHS
           call VecAssemblyBegin(Vec_F,ierr)
           call VecAssemblyEnd(Vec_F,ierr)
           ! Solve
           call PrintMsg(" Solving ...")
           call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
           call GetVec_U; tot_uu=tot_uu+uu
           if (poro) call VecAXPY(Vec_Um,f1,Vec_U,ierr)
           if (visco) then
              ! Recover stress
              call PrintMsg(" Recovering stress ...")
              do i=1,nels
                 call RecoverVStress(i,stress)
              end do
           end if
           ! Write output
           if (mod(tstep,frq)==0) call WriteOutput
           if (nobs_loc>0) then
              call GetVec_obs
              tot_uu_obs=tot_uu_obs+uu_obs
              call WriteOutput_obs
           end if
           !end if
        end do
        if (poro) deallocate(uup,kc,indxp,Hs,work)
     end if
     if (poro) then
        call VecDestroy(Vec_I,ierr)
        call VecDestroy(Vec_Up,ierr)
        call MatDestroy(Mat_Kc,ierr)
        call ISDestroy(RI,ierr)
        call VecDestroy(Vec_Um,ierr)
     end if
     call KSPDestroy(Krylov,ierr)
  end if

  ! Generalized-alpha implicit Solver
  if (galf) then
     ! Local to global fault node map
     if (nceqs-nceqs_ncf>0) then 
        call GetFltMap
        if (Xflt>1) call GetXfltMap ! Intersection
     end if
     if (bod_frc==1) then
        call PrintMsg("Applying gravity ...")
        call ApplyGravity
     end if
     call KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     call KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     call SetupKSPSolver
     call PrintMsg("Static solving ...")
     call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     if (rank==nprcs-1) print'(I0,A,I0,A)',j2," dofs on ",nprcs," processors."
     call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
     call GetVec_U; tot_uu=tot_uu+uu
     ! Get observation
     if (nobs_loc>0) then
        call GetVec_obs
        tot_uu_obs=tot_uu_obs+uu_obs
        dyn=.false.; call WriteOutput_obs 
     end if
     if (lm_str) then
        ! Recover stress
        call PrintMsg("Recovering stress ...")
        do i=1,nels
           if (rve>0 .and. eltype=="hex") then 
              call RecoverStress(i,stress,matDstr)
           else    
              call RecoverStress(i,stress)
           end if
        end do
     end if
     ! Initialize spaces, and stress-to-force.
     call VecGetLocalSize(Vec_U,j,ierr)
     call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     j3=0; j4=0
     do i=1,j
        if (j1+i-1>=dmn*nnds+nceqs_ncf) then
           j3=j3+1
        elseif (j1+i-1<dmn*nnds) then
           j4=j4+1
        end if
     end do
     allocate(workl(j3)); allocate(worku(j4))
     j3=0; j4=0
     do i=1,j
        if (j1+i-1>=dmn*nnds+nceqs_ncf) then
           j3=j3+1
           workl(j3)=j1+i-1
        elseif (j1+i-1<dmn*nnds) then
           j4=j4+1
           worku(j4)=j1+i-1
        end if
     end do
     j=size(worku)
     call ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,RIu,ierr)
     j=size(workl)
     call ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values,RIl,ierr)
     call VecGetSubVector(Vec_U,RIu,Vec_Uu,ierr)
     call VecDuplicate(Vec_Uu,Vec_fl,ierr) ! Ifl->-Gtuul
     if (nceqs>0) call VecDuplicate(Vec_Uu,Vec_flc,ierr)
     if (lm_str) then
        call VecDuplicate(Vec_Uu,Vec_SS,ierr)
        call VecDuplicate(Vec_Uu,Vec_SH,ierr)
        call VecZeroEntries(Vec_SS,ierr)
        call VecZeroEntries(Vec_SH,ierr)
        if (nceqs>0) then
           call VecDuplicate(Vec_Uu,Vec_f2s,ierr)
           call VecZeroEntries(Vec_f2s,ierr)
        end if
        call VecDuplicate(Vec_Uu,Vec_dip,ierr)
        call VecZeroEntries(Vec_dip,ierr)
        call VecDuplicate(Vec_Uu,Vec_nrm,ierr)
        call VecZeroEntries(Vec_nrm,ierr)
     end if
     call VecRestoreSubVector(Vec_U,RIu,Vec_Uu,ierr)
     j=size(indxmap_u,1)
     if (nceqs>0) then
        call VecCreateSeq(Petsc_Comm_Self,j,Seq_fl,ierr)
        call VecCreateSeq(Petsc_Comm_Self,j,Seq_flc,ierr)
        call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),                 &
           Petsc_Copy_Values,From_u,ierr)
        call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),                 &
           Petsc_Copy_Values,To_u,ierr)
        call VecScatterCreate(Vec_fl,From_u,Seq_fl,To_u,Scatter_u,ierr)
     elseif (lm_str) then
        call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),                 &
           Petsc_Copy_Values,From_u,ierr)
        call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),                 &
           Petsc_Copy_Values,To_u,ierr)
        call VecScatterCreate(Vec_U,From_u,Seq_U,To_u,Scatter_u,ierr)
     end if
     if (nceqs>0) then
        allocate(fl(size(indxmap_u,1)))
        allocate(flc(size(indxmap_u,1)))
        ! Vector to communicate with dynamic LM
        call VecGetSubVector(Vec_U,RIl,Vec_Ul,ierr)
        ! Extract Coulomb force (needed by force-stress translation).
        ! TODO call LM_s2d, and have GetVec_fcoulomb GetVec_f2s in parallel
        call GetVec_fcoulomb
        call VecRestoreSubVector(Vec_U,RIl,Vec_Ul,ierr)
     end if
     if (lm_str) then ! Scatter stress to nodes
        call GetVec_Stress
        call GetVec_S
        if (nceqs>0) then
           call GetVec_f2s
           call GetVec_f2s_seq
        else
           call GetVec_ft
        end if
        call GetVec_dip_nrm
        if (vout==1) call WriteOutput_f2s
        ! Cleanup
        call VecDestroy(Vec_dip,ierr)
        call VecDestroy(Seq_dip,ierr)
        call VecDestroy(Vec_nrm,ierr)
        call VecDestroy(Seq_nrm,ierr)
        call VecDestroy(Seq_f2s,ierr)
        deallocate(f2s,dip,nrm)
     end if

     ! Initialize slip indicator and slip history
     if (nceqs>0 .and. hyb>0) then
        allocate(slip(nfnd),slip0(nfnd),slip_sum(nfnd)) 
        slip(:)=0;slip_sum=0
        n_log_wave=0
        n_log_slip=0
        if (rsf==1) then 
           allocate(mu_cap(nfnd),rsfv(nfnd))
           trunc=f0
        end if
        allocate(mu_hyb(nfnd))
     end if
     ! Prepare dynamic spaces
     call VecZeroEntries(Vec_U_dyn,ierr)
     call VecDuplicate(Vec_U_dyn,Vec_Um_dyn,ierr)
     call VecDuplicate(Vec_U_dyn,Vec_U_dyn_tot,ierr)
     call VecZeroEntries(Vec_Um_dyn,ierr)
     call VecZeroEntries(Vec_U_dyn_tot,ierr)
     call MatCreateTranspose(Mat_Gt,Mat_G,ierr)

     ! Cleanup one-time static space
     call MatDestroy(Mat_K,ierr)
     call VecDestroy(Vec_Ul,ierr)

     ! Prepare implicit dynamic run 
     dyn=.true.
     write_dyn=.true.
     call VecDuplicate(Vec_U_dyn,Vec_F_dyn,ierr)
     call VecDuplicate(Vec_U_dyn,Vec_Fm_dyn,ierr)
     call VecZeroEntries(Vec_F_dyn,ierr)
     call VecZeroEntries(Vec_Fm_dyn,ierr)
     call AlphaInit

     if (nceqs>0) then
        ! Scatter static LM to dynamic Vec_lambda_sta
        call VecGetSubVector(Vec_U,RIl,Vec_Ul,ierr)
        call LM_s2d  
        call VecRestoreSubVector(Vec_U,RIl,Vec_Ul,ierr)
        call VecDuplicate(Vec_lambda_sta,Vec_lambda,ierr)
        ! GMinvGt is replaced by its inverse, assuming GMinvGt is diagonal
        call MatPtAP(Mat_Minv,Mat_Gt,Mat_Initial_Matrix,f1,Mat_GMinvGt,ierr)
        call VecDuplicateVecsF90(Vec_lambda,2,Vec_Wlm,ierr)
        call MatGetDiagonal(Mat_GMinvGt,Vec_Wlm(1),ierr)
        call MatSetOption(Mat_GMinvGt,Mat_New_Nonzero_Allocation_Err,          &
           PETSC_FALSE,ierr)
        call VecReciprocal(Vec_Wlm(1),ierr)
        call MatDiagonalSet(Mat_GMinvGt,Vec_Wlm(1),Insert_Values,ierr)
        call VecZeroEntries(Vec_Wlm(1),ierr)
        ! Form 1/dt^2GMinvGt assuming it doesnt change with time
        ! FIXME Scale with f1m, (1-alpha_m)/beta/dt**2, instead of 1/dt**2 
        ! seems to replicate explicit solver, why? 
        !call MatScale(Mat_GMinvGt,f1/dt**2,ierr)
        call MatScale(Mat_GMinvGt,f1m,ierr)
        call VecDuplicateVecsF90(Vec_U_dyn,2,Vec_W,ierr)
        ! Dynamic constraint I=0 
        call VecDuplicate(Vec_lambda,Vec_I_dyn,ierr)
        call VecZeroEntries(Vec_I_dyn,ierr)
        call VecDuplicate(Vec_lambda,Vec_lambda_tot,ierr)
        call VecZeroEntries(Vec_lambda_tot,ierr)
        ! Record static slip and traction  
        if (nfnd_loc>0) call Write_fe("sta")
     end if

     call KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(Krylov,Mat_Ka,Mat_Ka,Different_Nonzero_Pattern,ierr)
#else
     call KSPSetOperators(Krylov,Mat_Ka,Mat_Ka,ierr)
#endif
     call SetupKSPSolver    
     call PrintMsg("Alpha-solving ...")
     call VecGetOwnershipRange(Vec_U_dyn,j1,j2,ierr)
     if (rank==nprcs-1) print'(I0,A,I0,A)',j2," dofs on ",nprcs," processors."
     ! Dummy static event log
     n_log=1; crp=.false.
     if (rank==0) call WriteOutput_log 
     ! Start time stepping
     steps=int(ceiling(t/dt)); t_hyb=f0
     do tstep=1,steps
        t_hyb=t_hyb+dt
        ! Reform RHS
        call AlphaRHS
        ! Solve
        call AlphaUpdate
        call GetVec_U; tot_uu_dyn=tot_uu_dyn+uu_dyn
        ! Output observation
        if (nobs_loc>0) then
           call GetVec_obs
           tot_uu_dyn_obs=tot_uu_dyn_obs+uu_dyn_obs
           if (mod(tstep,frq_wave)==0) then
              uu_dyn_obs=uu_dyn_obs/dt
              call WriteOutput_obs
           end if
        end if
        if (mod(tstep,frq_wave)==0) n_log_wave=n_log_wave+1
        ! Output fault slip
        if (mod(tstep,frq_slip)==0 .and. nfnd>0) then
           if (nfnd_loc>0) call Write_fe("dyn") ! H5 file
           n_log_slip=n_log_slip+1
        end if
        ! Output snapshots
        if (vout==1 .and. mod(tstep,frq)==0) then 
           uu_dyn=uu_dyn/dt
           call WriteOutput
        end if
     end do
     if (rank==0) then
        call WriteOutput_log_wave
        if (nfnd>0) call WriteOutput_log_slip
        ! Dummy static event log
        crp=.true.; call WriteOutput_log
     end if
     call KSPDestroy(Krylov,ierr)
  end if

  ! Fault/hybrid solver
  if (fault) then
     ! Local to global fault node map
     if (nceqs-nceqs_ncf>0) then 
        call GetFltMap
        if (Xflt>1) call GetXfltMap ! Intersection
     end if
     if (bod_frc==1) then
        call PrintMsg("Applying gravity ...")
        call ApplyGravity
     end if
     ! Create cumulative static solution space Vec_Um 
     call VecDuplicate(Vec_U,Vec_Um,ierr) ! U->du & Um->u
     call VecZeroEntries(Vec_Um,ierr)
     call VecGetLocalSize(Vec_U,j,ierr)
     call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     if (rank==nprcs-1) print'(I0,A,I0,A)',j2," dofs on ", nprcs," processors."
     ! Create pressure, force and traction IDs (RI RIu RIl) 
     if (poro) then
        j2=0; j3=0; j4=0; j5=0; j6=0
        do i=1,j
           if (mod(j1+i,dmn+1)==0 .and. j1+i-1<(dmn+1)*nnds) then
              j2=j2+1
           end if
           if (nceqs>0) then
              if (j1+i-1>=(dmn+1)*nnds+nceqs_ncf) then
                 j3=j3+1
              end if
           end if
           if (mod(j1+i,dmn+1)>0 .and. j1+i-1<(dmn+1)*nnds) then
              j4=j4+1
           end if
           if (j1+i-1<(dmn+1)*nnds) then
              j5=j5+1
           end if
        end do
        allocate(work(j2),workl(j3),worku(j4))
        j2=0; j3=0; j4=0; j5=0
        do i=1,j
           if (mod(j1+i,dmn+1)==0 .and. j1+i-1<(dmn+1)*nnds) then
              j2=j2+1
              work(j2)=j1+i-1
           end if
           if (mod(j1+i,dmn+1)>0 .and. j1+i-1<(dmn+1)*nnds) then
              j4=j4+1
              worku(j4)=j1+i-1
           end if
           if (nceqs>0) then
              if (j1+i-1>=(dmn+1)*nnds+nceqs_ncf) then
                 j3=j3+1
                 workl(j3)=j1+i-1
              end if
           end if
        end do
        j=size(work)
        allocate(uup(j))
        call ISCreateGeneral(Petsc_Comm_World,j,work,Petsc_Copy_Values,RI,ierr)
        j=size(worku)
        call ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,RIu,   &
           ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
        call MatGetSubMatrix(Mat_K,RIu,RI,Mat_Initial_Matrix,Mat_H,ierr)
#else
        call MatCreateSubMatrix(Mat_K,RIu,RI,Mat_Initial_Matrix,Mat_H,ierr)
#endif
        if (nceqs > 0) then
           j=size(workl)
           call ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values,    &
              RIl,ierr)
        end if
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=7)
        call MatGetSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
#else
        call MatCreateSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
#endif
        call MatZeroEntries(Mat_Kc,ierr)
        allocate(kc(eldofp,eldofp),indxp(eldofp),Hs(eldofp))
        do i=1,nels
           call FormLocalK(i,k,indx,"Kc")
           kc=k(eldof+1:,eldof+1:)
           indxp=indx(eldof+1:); indxp=indxmap(indxp,2)
           indxp=((indxp+1)/(dmn+1))-1
           call MatSetValues(Mat_Kc,eldofp,indxp,eldofp,indxp,kc,Add_Values,   &
              ierr)
        end do
        call MatAssemblyBegin(Mat_Kc,Mat_Final_Assembly,ierr)
        call MatAssemblyEnd(Mat_Kc,Mat_Final_Assembly,ierr)
        ! Initialize space for lambda, p related nodal force
        call VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
        call VecDuplicate(Vec_Uu,Vec_fp,ierr) ! fp->Hp (required by FV)
        call VecDuplicate(Vec_Uu,Vec_fl,ierr) ! Ifl->-Gtuul
        call VecCopy(Vec_Uu,Vec_fl,ierr) ! Hold Uu
        call VecDuplicate(Vec_Uu,Vec_flc,ierr)
        if (lm_str) then
           call VecDuplicate(Vec_Uu,Vec_SS,ierr)
           call VecDuplicate(Vec_Uu,Vec_SH,ierr)
           call VecZeroEntries(Vec_SS,ierr)
           call VecZeroEntries(Vec_SH,ierr)
           if (nceqs>0) then
              call VecDuplicate(Vec_Uu,Vec_f2s,ierr)
              call VecZeroEntries(Vec_f2s,ierr)
           end if
           call VecDuplicate(Vec_Uu,Vec_dip,ierr)
           call VecZeroEntries(Vec_dip,ierr)
           call VecDuplicate(Vec_Uu,Vec_nrm,ierr)
           call VecZeroEntries(Vec_nrm,ierr)
        end if
        call VecRestoreSubVector(Vec_Um,RIu,Vec_Uu,ierr)
        ! Initialize pore pressure from FV model
        if (fvin>0) then
           call PrintMsg("Initializing FV pressure ...")
           call VecZeroEntries(Vec_Um,ierr) ! Zero absolute U
           kfv=.true. ! FV bc 
           if (fvin<3) then
              call FVInit
              call FVReformKF(ef_eldof)
           else
              call FVInitUsg
              call FVReformKFUsg(ef_eldof)
           end if
        end if
     else ! Not poroelastic
        j3=0; j4=0
        do i=1,j
           if (j1+i-1>=dmn*nnds+nceqs_ncf) then
              j3=j3+1
           elseif (j1+i-1<dmn*nnds) then
              j4=j4+1
           end if
        end do
        allocate(workl(j3)); allocate(worku(j4))
        j3=0; j4=0
        do i=1,j
           if (j1+i-1>=dmn*nnds+nceqs_ncf) then
              j3=j3+1
              workl(j3)=j1+i-1
           elseif (j1+i-1<dmn*nnds) then
              j4=j4+1
              worku(j4)=j1+i-1
           end if
        end do
        j=size(worku)
        call ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,       &
           RIu,ierr)
        j=size(workl)
        call ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values,       &
           RIl,ierr)
     end if ! Is poro? 

     call KSPCreate(Petsc_Comm_World,Krylov,ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
     call KSPSetOperators(Krylov,Mat_K,Mat_K,Different_Nonzero_Pattern,ierr)
#else
     call KSPSetOperators(Krylov,Mat_K,Mat_K,ierr)
#endif
     call SetupKSPSolver
     call PrintMsg("Static solving (step 0) ...")
     call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
     call GetVec_U; tot_uu=tot_uu+uu
     call VecAXPY(Vec_Um,f1,Vec_U,ierr)
     ! Get observation
     if (nobs_loc>0) then
        call GetVec_obs
        tot_uu_obs=tot_uu_obs+uu_obs
     end if
     if (visco .or. lm_str) then
        ! Recover stress
        call PrintMsg("Recovering stress ...")
        do i=1,nels
           if (rve>0 .and. eltype=="hex") then  
              call RecoverStress(i,stress,matDstr)
           else
              call RecoverStress(i,stress)
           end if
        end do
     end if
     
     ! Have more then one static step
     if (t>f0 .and. dt>f0 .and. t>=dt) then
        if (poro) then
           ! Initialize space for u, lambda related nodal flux 
           call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
           call VecDuplicate(Vec_Up,Vec_I,ierr) ! I->KcUp
           call VecCopy(Vec_Up,Vec_I,ierr) ! Hold Up
           call VecDuplicate(Vec_Up,Vec_qu,ierr) ! qu->Htu
           call VecDuplicate(Vec_Up,Vec_ql,ierr)
           call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
           ! Sequential vectors for output 
           j=size(indxmap_u,1)
           if (nceqs>0) then
              allocate(fl(size(indxmap_u,1)))
              allocate(ql(size(nl2g,1)))
              allocate(flc(size(indxmap_u,1)))
              call VecCreateSeq(Petsc_Comm_Self,j,Seq_fl,ierr)
              call VecCreateSeq(Petsc_Comm_Self,j,Seq_flc,ierr)
           end if
           call VecCreateSeq(Petsc_Comm_Self,j,Seq_fp,ierr)
           call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),              &
              Petsc_Copy_Values,From_u,ierr)
           call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),              &
              Petsc_Copy_Values,To_u,ierr)
           call VecScatterCreate(Vec_fp,From_u,Seq_fp,To_u,Scatter_u,ierr)
           j=size(nl2g,1)
           if (nceqs>0) call VecCreateSeq(Petsc_Comm_Self,j,Seq_ql,ierr)
           call VecCreateSeq(Petsc_Comm_Self,j,Seq_qu,ierr)
           call ISCreateGeneral(Petsc_Comm_Self,j,nl2g(:,2)-1,                 &
              Petsc_Copy_Values,From_p,ierr)
           call ISCreateGeneral(Petsc_Comm_Self,j,nl2g(:,1)-1,                 &
              Petsc_Copy_Values,To_p,ierr)
           call VecScatterCreate(Vec_qu,From_p,Seq_qu,To_p,Scatter_q,ierr)
           allocate(fp(size(indxmap_u,1)))
           allocate(qu(size(nl2g,1)))
           ! Extract nodal force by p, and fluid source by u
           if (vout==1) then
              call MatMult(Mat_H,Vec_I,Vec_fp,ierr)
              call VecScale(Vec_fp,-f1,ierr)
              call GetVec_fp
              call MatCreateTranspose(Mat_H,Mat_Ht,ierr)
              call MatMult(Mat_Ht,Vec_fl,Vec_qu,ierr)
              call GetVec_qu
           end if
           ! Extract lambda
           if (nceqs>0) then
              ! Vector to communicate with dynamic LM
              call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              ! Extract lambda induced nodal force
              ! TODO call LM_s2d, and have GetVec_fcoulomb in parallel 
              if (vout==1) then
                 call VecZeroEntries(Vec_fl,ierr)
                 call VecZeroEntries(Vec_flc,ierr)
                 call VecZeroEntries(Vec_ql,ierr)
                 call GetVec_ql
                 call GetVec_flambda
                 call GetVec_fl
                 call GetVec_fcoulomb
                 call GetVec_flc
              else
                 ! Extract Coulomb force (needed by force-stress translation).  
                 call GetVec_fcoulomb
              end if
              call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
           end if
           if (lm_str) then ! Scatter stress to nodes
              call GetVec_Stress
              call GetVec_S
              if (nceqs>0) then
                 call GetVec_f2s
                 call GetVec_f2s_seq
              else
                 call GetVec_ft
              end if
              call GetVec_dip_nrm
              if (vout==1) call WriteOutput_f2s
              call VecDestroy(Vec_dip,ierr)
              call VecDestroy(Seq_dip,ierr)
              call VecDestroy(Vec_nrm,ierr)
              call VecDestroy(Seq_nrm,ierr)
              call VecDestroy(Seq_f2s,ierr)
              deallocate(f2s,dip,nrm)
           end if
           ! Write output
           if (vout==1) call WriteOutput_x
           if (nobs_loc>0) call WriteOutput_obs 
           if (init==1) then
              call PrintMsg("Applying one step (24 hr) fluid source ...")
              ! Zero initial pressure 
              call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
              call VecZeroEntries(Vec_Up,ierr)
              call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
              call VecZeroEntries(Vec_F,ierr)
              call ApplySource
              call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
              call GetVec_U; tot_uu=uu
              call VecAXPY(Vec_Um,f1,Vec_U,ierr)
              call Rscdt(fdt) ! Scale [K] with dt = fdt*24hr (pv)
              if (vout==1) then
                 if (nceqs>0) then
                    call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                    call VecZeroEntries(Vec_flc,ierr)
                    call GetVec_fcoulomb
                    call GetVec_flc
                    call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                 end if
                 ! Should be analyzed to see if any initial slip 
                 call WriteOutput_init
              end if
           end if 
           if (fvin/=0) then
              if (vout==1) then
                 call WriteOutput_init
              end if
              ! Single phase FV model 
              if (fvin==1) call FVReformKPerm(f0,ef_eldof) 
              if (fvin==3) call FVReformKPermUsg(f0,ef_eldof)
           end if
        else ! Not poro
           call VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           if (nceqs>0) then 
              call VecDuplicate(Vec_Uu,Vec_fl,ierr) ! Ifl->-Gtuul
              call VecDuplicate(Vec_Uu,Vec_flc,ierr)
           end if
           if (lm_str) then
              call VecDuplicate(Vec_Uu,Vec_SS,ierr)
              call VecDuplicate(Vec_Uu,Vec_SH,ierr)
              call VecZeroEntries(Vec_SS,ierr)
              call VecZeroEntries(Vec_SH,ierr)
              if (nceqs>0) then
                 call VecDuplicate(Vec_Uu,Vec_f2s,ierr)
                 call VecZeroEntries(Vec_f2s,ierr)
              end if
              call VecDuplicate(Vec_Uu,Vec_dip,ierr)
              call VecZeroEntries(Vec_dip,ierr)
              call VecDuplicate(Vec_Uu,Vec_nrm,ierr)
              call VecZeroEntries(Vec_nrm,ierr)
           end if
           call VecRestoreSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           ! Sequential vectors for output
           j=size(indxmap_u,1)
           if (nceqs>0) then
              allocate(fl(size(indxmap_u,1)))
              allocate(flc(size(indxmap_u,1)))
              call VecCreateSeq(Petsc_Comm_Self,j,Seq_fl,ierr)
              call VecCreateSeq(Petsc_Comm_Self,j,Seq_flc,ierr)
              call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),           &
                 Petsc_Copy_Values,From_u,ierr)
              call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),           &
                 Petsc_Copy_Values,To_u,ierr)
              call VecScatterCreate(Vec_fl,From_u,Seq_fl,To_u,Scatter_u,ierr)
              ! Vector to communicate with dynamic LM
              call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              ! Extract lambda induced nodal force
              if (vout==1) then
                 call VecZeroEntries(Vec_fl,ierr)
                 call VecZeroEntries(Vec_flc,ierr)
                 call GetVec_flambda
                 call GetVec_fl
                 call GetVec_fcoulomb
                 call GetVec_flc
              else
                 ! Extract Coulomb force (needed by force-stress translation).
                 call GetVec_fcoulomb
              end if
              call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
           end if
           if (lm_str) then ! Scatter stress to nodes
              if (nceqs==0) then ! Sequential vectors for stress
                 call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),        &
                     Petsc_Copy_Values,From_u,ierr)
                 call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),        &
                    Petsc_Copy_Values,To_u,ierr)
                 call VecScatterCreate(Vec_U,From_u,Seq_U,To_u,Scatter_u,ierr)
              end if
              call GetVec_Stress
              call GetVec_S
              if (nceqs>0) then
                 call GetVec_f2s
                 call GetVec_f2s_seq
              else
                 call GetVec_ft
              end if
              call GetVec_dip_nrm
              if (vout==1) call WriteOutput_f2s
              ! Cleanup
              call VecDestroy(Vec_dip,ierr)
              call VecDestroy(Seq_dip,ierr)
              call VecDestroy(Vec_nrm,ierr)
              call VecDestroy(Seq_nrm,ierr)
              call VecDestroy(Seq_f2s,ierr)
              deallocate(f2s,dip,nrm)
           end if
           ! Write output
           if (vout==1) call WriteOutput_x
           if (nobs_loc>0) call WriteOutput_obs
        end if ! Poro or not

        ! Solution space is allocated differently for static and dynamic runs, 
        ! so we keep mat_K and Mat_K_dyn separate instead of
        !call MatGet(Create)SubMatrix(Mat_K,RIu,RIu,Mat_Initial_Matrix,         &
        !Mat_K_dyn,ierr)

        ! Initialize slip indicator and slip history
        if (nceqs>0 .and. hyb>0) then
           allocate(slip(nfnd),slip0(nfnd),slip_sum(nfnd)) 
           slip(:)=0;slip_sum=0
           n_log_wave=0
           n_log_slip=0
           if (rsf==1) then
              allocate(mu_cap(nfnd),rsfv(nfnd))
              trunc=f0
           end if
           allocate(mu_hyb(nfnd))
        end if

        ! Start implicit time step
        steps=int(ceiling(t/dt)); t_abs=f0
        dyn=.false.; fail=.false.; n_log=0
        do tstep=1,steps
           t_abs=t_abs+dt
           if (rank==0) print'(A11,I0)'," Time Step ",tstep
           ! Reform stiffness matrix, if needed
           if (visco .and. (tstep==1 .or. maxval(mat(:,4))>f1) .and. .not.     &
              (poro .and. (fvin==2 .or. fvin==4))) then
              call PrintMsg(" Reforming [K] ...")
              call MatZeroEntries(Mat_K,ierr)
              do i=1,nels
                 call FormLocalK(i,k,indx,"Kv")
                 indx=indxmap(indx,2)
                 call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,        &
                    Add_Values,ierr)
              end do
              ! Account for constraint eqn's
              if (rank==0 .and. nceqs>0) call ApplyConstraints
              call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
              call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
           end if
           if (poro .and. (fvin==2 .or. fvin==4)) then ! Multiphase FV model
              call PrintMsg(" Reforming [Kp] ...")
              if (fvin==2) call FVReformKPerm(t_abs,ef_eldof)
              if (fvin==4) call FVReformKPermUsg(t_abs,ef_eldof)
           end if
           ! Reform RHS
           call VecZeroEntries(Vec_F,ierr)
           call PrintMsg(" Reforming RHS ...")
           if (poro) then
              call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
              call MatMult(Mat_Kc,Vec_Up,Vec_I,ierr)
              call VecScale(Vec_I,-dt,ierr)
              call VecGetArrayF90(Vec_I,pntr,ierr)
              uup=pntr
              call VecRestoreArrayF90(Vec_I,pntr,ierr)
              j=size(uup)
              call VecSetValues(Vec_F,j,work,uup,Add_Values,ierr)
              !call Up_s2d
              call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
              if (fvin>0) then! Remove hydrostatic pressure gradient from RHS
                 call MatMult(Mat_Kc,Vec_Up_hst,Vec_I,ierr)
                 call VecScale(Vec_I,dt,ierr)
                 call VecGetArrayF90(Vec_I,pntr,ierr)
                 uup=pntr
                 call VecRestoreArrayF90(Vec_I,pntr,ierr)
                 j=size(uup)
                 call VecSetValues(Vec_F,j,work,uup,Add_Values,ierr)
              end if
              ! Stabilize RHS
              do i=1,nels
                 call FormLocalHs(i,Hs,indxp)
                 indxp=indxmap(indxp,2)
                 call VecSetValues(Vec_F,npel,indxp,Hs,Add_Values,ierr)
              end do
           end if
           if (visco) then
              do i=1,nels
                 call ReformLocalRHS(i,f,indx)
                 indx=indxmap(indx,2)
                 call VecSetValues(Vec_F,eldof,indx,f,Add_Values,ierr)
              end do
           end if
           call FormRHS
           if (fail) then 
              call FaultSlip_sta
              ! Backup QS slip 
              qs_flt_slip=qs_flt_slip+tot_flt_slip
              fail=.false.
              if (.not. crp) res_flt_slip=f0
           end if
           call VecAssemblyBegin(Vec_F,ierr)
           call VecAssemblyEnd(Vec_F,ierr)
           ! Impose FV pressure
           if (fvin==1 .or. fvin==2) call FVSyncBD(t_abs)
           if (fvin==3 .or. fvin==4) call FVSyncBDUsg(t_abs)
           ! Solve
           call PrintMsg(" Solving ...")
           call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
           ! Reset dynamic (slip) solutions 
           if (nceqs>0 .and. hyb>0) then
              tot_uu_dyn=f0
              if (nobs_loc>0) tot_uu_dyn_obs=f0
              flt_slip=f0; tot_flt_slip=f0
           end if
           n_log=n_log+1
           call GetVec_U; tot_uu=tot_uu+uu
           if (nobs_loc>0) then 
              call GetVec_obs
              tot_uu_obs=tot_uu_obs+uu_obs
           end if
           ! Record absolute solution
           if (poro .or. nceqs>0) call VecAXPY(Vec_Um,f1,Vec_U,ierr)
           if (visco .or. vout==1) then
              ! Recover stress
              call PrintMsg(" Recovering stress ...")
              do i=1,nels
                 if (visco) then   
                    call RecoverVStress(i,stress)
                 else if (rve>0 .and. eltype=="hex") then 
                    call RecoverStress(i,stress,matDstr)
                 else 
                    call RecoverStress(i,stress)
                 end if
              end do
              call GetVec_Stress
              call GetVec_S
           end if
           ! Extract nodal force by p
           if (poro) then
              ! Reset tot_uu at step 1
              if (tstep==1 .and. fvin>0) then
                 tot_uu=f0 
                 if (nobs_loc>0) tot_uu_obs=f0
              end if
              if (vout==1) then
                 ! Extract force by p
                 call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
                 call MatMult(Mat_H,Vec_Up,Vec_fp,ierr)
                 call VecScale(Vec_fp,-f1,ierr)
                 call GetVec_fp
                 call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
                 ! Extract fluid source by u
                 call VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
                 call MatMult(Mat_Ht,Vec_Uu,Vec_qu,ierr)
                 call GetVec_qu
                 call VecRestoreSubVector(Vec_Um,RIu,Vec_Uu,ierr)
                 ! Extract lambda and induced nodal f and q
                 if (nceqs>0) then
                    call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                    call VecZeroEntries(Vec_fl,ierr)
                    call VecZeroEntries(Vec_flc,ierr)
                    call VecZeroEntries(Vec_ql,ierr)
                    call GetVec_flambda
                    call GetVec_fl
                    call GetVec_ql
                    call GetVec_fcoulomb
                    call GetVec_flc
                    call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                 end if
                 ! Write output
                 if (mod(tstep,frq)==0) call WriteOutput_x
              end if
              if (nobs_loc>0) call WriteOutput_obs
           else
              if (vout==1) then
                 ! Extract lambda and induced nodal f
                 if (nceqs>0) then
                    call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                    call VecZeroEntries(Vec_fl,ierr)
                    call VecZeroEntries(Vec_flc,ierr)
                    call GetVec_flambda
                    call GetVec_fl
                    call GetVec_fcoulomb
                    call GetVec_flc
                    call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                 end if
                 ! Write output
                 if (mod(tstep,frq)==0) call WriteOutput_x
              end if
              if (nobs_loc>0) call WriteOutput_obs
           end if

           ! Determine if the fault shall fail by LM 
           if (nceqs-nceqs_ncf>0 .and. hyb>0) then
              call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              call LM_s2d  
              call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              if (poro) then 
                 call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
                 call Up_s2d
                 call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
              end if  
              ! Record static slip and traction  
              if (mod(tstep,frq)==0 .and. nfnd_loc>0) call Write_fe("sta")
              ! Determine if the fault shall fail
              call GetSlip_sta
              rslip=real(sum(slip))/real(size(slip))
              if (rank==0) print'(F0.2,A,I0,A)',rslip*100.0,"% (",sum(slip),   &
                 ") fault nodes critical."
              if (rslip>f0 .and. sum(slip)>=1) then ! Failure threshold 
                 dyn=.true. 
              end if
           end if

           ! Prepare variables for hybrid run
           if (dyn) then
              ! Prepare the working space for the dynamic model
              call VecZeroEntries(Vec_U_dyn,ierr)
              call VecDuplicate(Vec_U_dyn,Vec_Um_dyn,ierr)
              call VecDuplicate(Vec_U_dyn,Vec_Up_dyn,ierr)
              call VecDuplicate(Vec_U_dyn,Vec_U_dyn_tot,ierr)
              call VecZeroEntries(Vec_Um_dyn,ierr)
              call VecZeroEntries(Vec_Up_dyn,ierr)
              call VecZeroEntries(Vec_U_dyn_tot,ierr)
              call VecDuplicateVecsF90(Vec_U_dyn,6,Vec_W,ierr)
              call MatCreateTranspose(Mat_Gt,Mat_G,ierr)
              call VecDuplicate(Vec_lambda_sta,Vec_lambda,ierr)
              ! GMinvGt is replaced by its inverse, assuming GMinvGt is diagonal
              call MatPtAP(Mat_Minv,Mat_Gt,Mat_Initial_Matrix,f1,              &
                 Mat_GMinvGt,ierr)
              call VecDuplicateVecsF90(Vec_lambda,2,Vec_Wlm,ierr)
              call MatGetDiagonal(Mat_GMinvGt,Vec_Wlm(1),ierr)
              call MatSetOption(Mat_GMinvGt,Mat_New_Nonzero_Allocation_Err,    &
                 PETSC_FALSE,ierr)
              call VecReciprocal(Vec_Wlm(1),ierr)
              call MatDiagonalSet(Mat_GMinvGt,Vec_Wlm(1),Insert_Values,ierr)
              call VecZeroEntries(Vec_Wlm(1),ierr)
              ! Dynamic constraint I=0 
              call VecDuplicate(Vec_lambda,Vec_I_dyn,ierr)
              call VecZeroEntries(Vec_I_dyn,ierr)
              ! Pass pseudo velocity to U_dyn
              if (rsf==1 .and. tstep>1) call Rsfv2Dyn 
              call VecDuplicate(Vec_lambda,Vec_lambda_tot,ierr)
              call VecZeroEntries(Vec_lambda_tot,ierr)
              ! Form 1/dt^2GMinvGt assuming it doesnt change with time
              call MatScale(Mat_GMinvGt,f1/dt_dyn**2,ierr)
              call PrintMsg("Hybrid Solving ...")
              steps_dyn=int(ceiling(t_dyn/dt_dyn)); t_hyb=t_abs
              ih=0; t_sta=f0; fail=.true.;crp=.false.
           end if
           
           ! Explicit/implicit hybrid step for rupture propagation
           do while (dyn .or. (t_sta>f0 .and. (t_sta<t_lim) .and. .not. crp))
              ! Explicit time step
              do tstep_dyn=0,steps_dyn
                 t_hyb=t_hyb+dt_dyn
                 ! Up=Minv(dt^2(F-KU)-dt(C(U-Um)))+2U-Um
                 call MatMult(Mat_K_dyn,Vec_U_dyn,Vec_W(1),ierr)
                 call VecScale(Vec_W(1),-f1*dt_dyn**2,ierr)
                 call VecWAXPY(Vec_W(2),-f1,Vec_Um_dyn,Vec_U_dyn,ierr)
                 call MatMult(Mat_K_dyn,Vec_W(2),Vec_W(3),ierr)
                 call MatMult(Mat_M,Vec_W(2),Vec_W(4),ierr)
                 call VecAXPY(Vec_W(4),beta,Vec_W(3),ierr)
                 call VecScale(Vec_W(4),dt_dyn,ierr)
                 call VecWAXPY(Vec_W(5),-f1,Vec_W(4),Vec_W(1),ierr)
                 call MatMult(Mat_Minv,Vec_W(5),Vec_W(6),ierr)
                 call VecWAXPY(Vec_Up_dyn,f2,Vec_U_dyn,Vec_W(6),ierr)
                 call VecAXPY(Vec_Up_dyn,-f1,Vec_Um_dyn,ierr)
                 ! Form lambda=(GUp-Flm)/(dt^2GMinvGt)
                 call MatMult(Mat_G,Vec_Up_dyn,Vec_Wlm(1),ierr)
                 call VecWAXPY(Vec_Wlm(2),-f1,Vec_I_dyn,Vec_Wlm(1),ierr)
                 call MatMult(Mat_GMinvGt,Vec_Wlm(2),Vec_lambda,ierr)
                 if (rsf==1 .and. tstep>1 .and. tstep_dyn+ih==0) then 
                     ! Skip friction law, and reset constraint 
                     call VecZeroEntries(Vec_I_dyn,ierr)
                 else 
                    ! Cap the nodal LM not to exceed max friction
                    call CapLM_dyn
                 end if
                 call VecAXPY(Vec_lambda_tot,f1,Vec_lambda,ierr)
                 ! Form Up=Up-dt^2(Minv(Gtlambda))
                 call MatMult(Mat_Gt,Vec_lambda,Vec_W(1),ierr)
                 call MatMult(Mat_Minv,Vec_W(1),Vec_W(2),ierr)
                 call VecAXPY(Vec_Up_dyn,-dt_dyn**2,Vec_W(2),ierr)
                 ! Update solution
                 call VecCopy(Vec_U_dyn,Vec_Um_dyn,ierr)
                 call VecCopy(Vec_Up_dyn,Vec_U_dyn,ierr)
                 call VecAXPY(Vec_U_dyn_tot,f1,Vec_U_dyn,ierr)
                 ! Extract dynamic solution
                 if (.not. dyn) then
                    dyn=.true.
                    call GetVec_U
                    tot_uu_dyn=tot_uu_dyn+uu_dyn
                    if (nobs_loc>0) then
                       call GetVec_obs
                       tot_uu_dyn_obs=tot_uu_dyn_obs+uu_dyn_obs
                       uu_dyn_obs=uu_dyn_obs/dt_dyn
                       if (mod(n_log_dyn,frq_wave)==0) call WriteOutput_obs 
                    end if
                    if (mod(n_log_dyn,frq_wave)==0) n_log_wave=n_log_wave+1
                    dyn=.false.
                 else
                    call GetVec_U
                    tot_uu_dyn=tot_uu_dyn+uu_dyn
                    if (nobs_loc>0) then
                       call GetVec_obs
                       tot_uu_dyn_obs=tot_uu_dyn_obs+uu_dyn_obs
                       uu_dyn_obs=uu_dyn_obs/dt_dyn
                       if (mod(n_log_dyn,frq_wave)==0) call WriteOutput_obs
                    end if
                    if (mod(n_log_dyn,frq_wave)==0) n_log_wave=n_log_wave+1
                 end if
                 ! Evaluated FD grid movement
                 if (ngp_loc>0 .and. mod(n_log_dyn,frq_wave)==0 .and.          &
                    fdout==1) then  
                    call GetVec_fd
                    call Write_fd("dyn")
                 end if
                 ! Extract and output temporal fault slip: flt_slip
                 call FaultSlip_dyn
                 tot_flt_slip=tot_flt_slip+flt_slip
                 if (mod(n_log_dyn,frq_slip)==0) then 
                    if (nfnd_loc>0) call Write_fe("dyn") ! H5 file
                    n_log_slip=n_log_slip+1
                 end if
                 ! Export dynamic snapshot
                 write_dyn=.true.
                 if (vout==1 .and. mod(n_log_dyn,frq_dyn)==0) then
                    uu_dyn=uu_dyn/dt_dyn
                    call WriteOutput
                 end if
                 write_dyn=.false.
                 n_log_dyn=n_log_dyn+1
              end do ! Explicit loop 
              ! Assess the fault status 
              if (dyn) call GetSlip_dyn
              rslip=real(sum(slip))/real(size(slip))
              if ((.not. sum(slip)>0) .or. (.not. dyn)) then
                 dyn=.false. ! Failure threshold
                 t_sta=t_sta+t_dyn
              end if
              if (rank==0 .and. dyn) print'(F0.2,A)',rslip*100.0,              &
                 "% fault slipping."
              if (rank==0 .and. crp) print'(A)',"Aseismic"
              if (rank==0 .and. .not. (dyn .or. crp)) print'(A,F0.4,A,F0.4)',  &
                 "Wave traveling ",t_sta,"/",t_lim
              ! Turn off "dyn" for debug if the fault cannot stabilize
              if (dble(ih+1)*t_dyn>=t_lim .and. dyn) then
                 dyn=.false.
                 slip=0
                 t_sta=t_sta+t_dyn
              end if
              ! Hybrid iteration count
              ih=ih+1
           end do ! Hybrid run
           if (fail) then
              if (rank==0 .and. nobs>0) then 
                 call WriteOutput_log
                 call WriteOutput_log_wave
                 if (nceqs>0) call WriteOutput_log_slip
              end if
              if (fdout==1 .and. ngp_loc>0) then 
                 call GetFDAct
                 call Write_fd("act")
              end if
              ! Latest fault stress (Vec_lambda_sta0)
              call GetVec_lambda_hyb
              ! Cleanup dynamics
              call VecDestroy(Vec_Um_dyn,ierr)
              call VecDestroy(Vec_Up_dyn,ierr)
              call VecDestroy(Vec_U_dyn_tot,ierr)
              call VecDestroy(Vec_I_dyn,ierr)
              call VecDestroyVecsF90(6,Vec_W,ierr)
              call VecDestroy(Vec_lambda,ierr)
              call VecDestroy(Vec_lambda_tot,ierr)
              call VecDestroyVecsF90(2,Vec_Wlm,ierr)
              call MatDestroy(Mat_G,ierr)
              call MatDestroy(Mat_GMinvGt,ierr)
           else
              ! Backup fault stress
              call VecCopy(Vec_lambda_sta,Vec_lambda_sta0,ierr)
           end if
        end do ! Implicit time steps
     elseif (.not. poro) then ! One time static
        call VecDuplicate(Vec_U,Vec_SS,ierr)
        call VecDuplicate(Vec_U,Vec_SH,ierr)
        call VecZeroEntries(Vec_SS,ierr)
        call VecZeroEntries(Vec_SH,ierr)
        call GetVec_Stress
        Scatter_u=Scatter
        call GetVec_S
        if (vout==1) call WriteOutput_x
        if (nobs_loc>0) call WriteOutput_obs
     end if ! Assert implicit time
     crp=.true.
     if (rank==0) call WriteOutput_log 
     ! Cleanup hybrid
     if (poro) deallocate(uup,kc,indxp,Hs,work,fp,qu)
     if (nceqs>0) then 
        deallocate(fl,flc,workl,worku)
        if (lm_str) deallocate(ss,sh)
        if (poro) deallocate(ql)
     end if
     call MatDestroy(Mat_Kc,ierr)
     call MatDestroy(Mat_H,ierr)
     call MatDestroy(Mat_Ht,ierr)
     call MatDestroy(Mat_K_dyn,ierr)
     call MatDestroy(Mat_Gt,ierr)
     call MatDestroy(Mat_M,ierr)
     call MatDestroy(Mat_Minv,ierr)
     call VecDestroy(Vec_I,ierr)
     call VecDestroy(Vec_fp,ierr)
     call VecDestroy(Vec_Up,ierr)
     call VecDestroy(Vec_Uu,ierr)
     call VecDestroy(Vec_Um,ierr)
     call VecDestroy(Vec_qu,ierr)
     call VecDestroy(Vec_Ul,ierr)
     call VecDestroy(Vec_lambda_sta,ierr)
     call VecDestroy(Vec_lambda_sta0,ierr)
     if (rsf==1) call VecDestroy(Vec_lambda_bk,ierr)
     call VecDestroy(Vec_fl,ierr)
     call VecDestroy(Vec_flc,ierr)
     call VecDestroy(Vec_f2s,ierr)
     call VecDestroy(Vec_ql,ierr)
     call VecDestroy(Vec_U_dyn,ierr)
     call VecDestroy(Vec_SS,ierr)
     call VecDestroy(Vec_SH,ierr)
     call VecDestroy(Vec_Cst,ierr)
     call VecDestroy(Seq_fp,ierr)
     call VecDestroy(Seq_qu,ierr)
     call VecDestroy(Seq_fl,ierr)
     call VecDestroy(Seq_flc,ierr)
     call VecDestroy(Seq_ql,ierr)
     call VecDestroy(Seq_U_dyn,ierr)
     call VecDestroy(Seq_SS,ierr)
     call VecDestroy(Seq_SH,ierr)
     call ISDestroy(RI,ierr)
     call ISDestroy(RIu,ierr)
     call ISDestroy(RIl,ierr)
     call ISDestroy(From_u,ierr)
     call ISDestroy(To_u,ierr)   
     call ISDestroy(From_p,ierr)
     call ISDestroy(To_p,ierr)
     call VecScatterDestroy(Scatter_u,ierr)
     call VecScatterDestroy(Scatter_q,ierr)
     call VecScatterDestroy(Scatter_s2d,ierr)
     call VecScatterDestroy(Scatter_dyn,ierr)
     call KSPDestroy(Krylov,ierr)
  end if ! Fault solver

  ! Explicit (Green's function) Solver
  if (explicit .and. t>f0 .and. dt>f0 .and. t>=dt) then
     steps=int(ceiling(t/dt))
     ! Initialize work vectors
     call VecDuplicate(Vec_U,Vec_Um,ierr)
     call VecDuplicate(Vec_U,Vec_Up,ierr)
     call VecDuplicateVecsF90(Vec_U,6,Vec_W,ierr)
     if (nceqs>0) then
        call MatCreateTranspose(Mat_Gt,Mat_G,ierr)
        call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,nceqs,Vec_lambda,ierr)
        call VecDuplicate(Vec_lambda,Vec_I,ierr)
        call VecDuplicateVecsF90(Vec_lambda,2,Vec_Wlm,ierr)
     end if
     ! Start time stepping
     call PrintMsg("Solving ...") ! Up=Minv(dt^2(F-KU)-dt(C(U-Um)))+2U-Um
     call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     if (rank==nprcs-1) print'(I0,A,I0,A)',j2+nceqs," dofs on ", nprcs,        &
        " processors."
     ng=1
     if (nobs_loc>0) tot_uu_dyn_obs=f0
     do igf=1,nceqs
        if (ng>10) exit ! Max mumber of GFs
        j=(igf-1)/dmn+1
        do tstep=0,steps
           if (gf .and. frc(j)<1) then
              exit
           elseif (gf .and. tstep==0) then
              if (mod(igf,dmn)==0) then
                 exit
              else
                ng=ng+1
                if (rank==0) print'(A,I0,A,I0)',"Source ",j," component ",     &
                   mod(igf,dmn)
              end if
           end if
           if (rank==0 .and. .not. gf) print'(A12,I0)'," Time Step ",tstep
           ! Solve
           call MatMult(Mat_K,Vec_U,Vec_W(1),ierr)
           call VecAYPX(Vec_W(1),-f1,Vec_F,ierr)
           call VecScale(Vec_W(1),dt**2,ierr)
           call VecWAXPY(Vec_W(2),-f1,Vec_Um,Vec_U,ierr)
           call MatMult(Mat_K,Vec_W(2),Vec_W(3),ierr)
           call MatMult(Mat_M,Vec_W(2),Vec_W(4),ierr)
           call VecAXPY(Vec_W(4),beta,Vec_W(3),ierr)
           call VecScale(Vec_W(4),dt,ierr)
           call VecWAXPY(Vec_W(5),-f1,Vec_W(4),Vec_W(1),ierr)
           call MatMult(Mat_Minv,Vec_W(5),Vec_W(6),ierr)
           call VecWAXPY(Vec_Up,f2,Vec_U,Vec_W(6),ierr)
           call VecAXPY(Vec_Up,-f1,Vec_Um,ierr)
           if (nceqs>0) then
              if (tstep==0) then
                 call MatPtAP(Mat_Minv,Mat_Gt,Mat_Initial_Matrix,f1,           &
                    Mat_GMinvGt,ierr)
                 ! GMinvGt is replaced by its inverse, assuming GMinvGt 
                 ! is diagonal
                 call MatGetDiagonal(Mat_GMinvGt,Vec_Wlm(1),ierr)
                 call VecReciprocal(Vec_Wlm(1),ierr)
                 call MatDiagonalSet(Mat_GMinvGt,Vec_Wlm(1),Insert_Values,ierr)
                 call VecZeroEntries(Vec_Wlm(1),ierr)
                 ! Form 1/dt^2GMinvGt assuming it doesnt change with time
                 call MatScale(Mat_GMinvGt,f1/dt**2,ierr)
              end if
              ! Apply constraint values using Flm
              if (gf) then
                 if (tstep==0) then
                    call VecZeroEntries(Vec_I,ierr)
                    if (rank==0) call VecSetValue(Vec_I,igf-1,f1,Insert_Values,&
                       ierr) 
                    call VecAssemblyBegin(Vec_I,ierr)
                    call VecAssemblyEnd(Vec_I,ierr)
                 elseif (tstep==5) then
                    call VecZeroEntries(Vec_I,ierr)  
                 end if
              else
                call VecZeroEntries(Vec_I,ierr)
                if (rank==0) then
                   do i=1,nceqs
                      j1=nint(cval(i,2)/dt); j2=nint(cval(i,3)/dt)
                      if (tstep>=j1 .and. tstep<=j2) then
                         val=cval(i,1)
!                        ! Use a decaying function instead of boxcar, e.g., ...
!                        if (j2>j1) val=3.5d0*val*(dble(j2-tstep)/dble(j2-j1)) &
!                           **2.5d0
                         call VecSetValue(Vec_I,i-1,val,Insert_Values,ierr)
                      end if
                   end do
                end if
                call VecAssemblyBegin(Vec_I,ierr)
                call VecAssemblyEnd(Vec_I,ierr)
              end if
              ! Form lambda=(GUp-Flm)/(dt^2GMinvGt)
              call MatMult(Mat_G,Vec_Up,Vec_Wlm(1),ierr)
              call VecWAXPY(Vec_Wlm(2),-f1,Vec_I,Vec_Wlm(1),ierr)
              call MatMult(Mat_GMinvGt,Vec_Wlm(2),Vec_lambda,ierr)
              ! Form Up=Up-dt^2(Minv(Gtlambda))
              call MatMult(Mat_Gt,Vec_lambda,Vec_W(1),ierr)
              call MatMult(Mat_Minv,Vec_W(1),Vec_W(2),ierr)
              call VecAXPY(Vec_Up,-dt**2,Vec_W(2),ierr)
           end if
           call VecCopy(Vec_U,Vec_Um,ierr)
           call VecCopy(Vec_Up,Vec_U,ierr)
           if (gf) then
              if (nobs_loc>0) then
                 call GetVec_obs
                 tot_uu_dyn_obs=tot_uu_dyn_obs+uu_obs
                 uu_dyn_obs=uu_obs/dt
                 if (mod(n_log_dyn,frq_wave)==0) then
                    dyn=.true.
                    call WriteOutput_obs;
                    dyn=.false.
                 end if
              end if
              if (mod(n_log_dyn,frq_wave)==0) n_log_wave=n_log_wave+1
              n_log_dyn=n_log_dyn+1
              ! Write output
              call GetVec_U; tot_uu=tot_uu+uu
              if (mod(tstep,frq)==0) call WriteOutput
           else
              call VecZeroEntries(Vec_F,ierr)
              call FormRHS
              call VecAssemblyBegin(Vec_F,ierr)
              call VecAssemblyEnd(Vec_F,ierr)
              ! Write output
              call GetVec_U; tot_uu=tot_uu+uu
              if (mod(tstep,frq)==0) call WriteOutput
              if (nceqs>0) call VecZeroEntries(Vec_I,ierr)
           end if
        end do ! Explicit run
        if (.not. gf) exit 
        if (rank==0 .and. frc(j)==1 .and. mod(igf,dmn)/=0)                     &
           call WriteOutput_log_wave
        ! Reset solutions
        call VecZeroEntries(Vec_U,ierr) 
        call VecZeroEntries(Vec_Um,ierr)
        call VecZeroEntries(Vec_Up,ierr)
        if (nobs_loc>0) tot_uu_dyn_obs=f0
     end do ! Unit rupture loop
     call PrintMsg("Done")
     if (nceqs>0) then
        call MatDestroy(Mat_Gt,ierr)
        call MatDestroy(Mat_G,ierr)
        call MatDestroy(Mat_GMinvGt,ierr)
        call VecDestroy(Vec_lambda,ierr)
        call VecDestroy(Vec_I,ierr)
        call VecDestroyVecsF90(2,Vec_Wlm,ierr)
     end if
     call MatDestroy(Mat_M,ierr)
     call MatDestroy(Mat_Minv,ierr)
     call VecDestroy(Vec_Um,ierr)
     call VecDestroy(Vec_Up,ierr)
     call VecDestroyVecsF90(6,Vec_W,ierr)
  end if ! Explicit solver

  ! Cleanup
  call PrintMsg("Cleaning up ...")
  call VecScatterDestroy(Scatter,ierr)
  call ISDestroy(To,ierr)
  call ISDestroy(From,ierr)
  call VecDestroy(Seq_U,ierr)
  call VecDestroy(Vec_F,ierr)
  call VecDestroy(Vec_U,ierr)
  if (.not. galf) call MatDestroy(Mat_K,ierr) ! [K]-alpha destroyed early
  if (visco .or. ((fault .or. galf) .and. lm_str)) deallocate(stress)
  deallocate(coords,nodes,bc,mat,id,k,m,f,indx,ipoint,weight,enodes,ecoords,   &
     vvec,indxmap,tot_uu,uu,cval,fnode,fval,telsd,tval,nl2g)
  ! Delete cnstr.tmp 
  if (rank==0 .and. nceqs>0) then
    open(15, file=trim(output_file(:index(output_file,"/",BACK=.TRUE.)))//     &
       "cnstrns.tmp",status='old')
    close(15, status='delete')
  end if
  call PrintMsg("Finished")

9 call PetscFinalize(ierr)

contains

  ! Read simulation parameters
  subroutine ReadParameters
    implicit none
    read(10,*)stype,eltype,nodal_bw
    fault=.false.
    if (stype=="fault" .or. stype=="fault-p" .or. stype=="fault-pv" .or.       &
       stype=="fault-v" .or. stype=="fault-rve") fault=.true.
    explicit=.false.; gf=.false. 
    if (stype=="explicit" .or. stype=="explicit-gf" .or. stype=="explicit-rve")&
       then
       explicit=.true.
       if (stype=="explicit-gf") gf=.true.
    end if
    rve=0
    if (stype=="implcit-rve" .or. stype=="explcit-rve" .or. stype=="fault-rve" &
       .or. stype=="alpha-rve") rve=1
    galf=.false.   
    if (stype=="alpha" .or. stype=="alpha-rve") galf=.true.
    if (fault .or. gf .or. galf) then
       if (rve>0) then 
          read(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nfnd,nobs,nceqs_ncf,&
             nincl
       else
          read(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nfnd,nobs,nceqs_ncf
       end if
    else
       if (rve>0) then 
          read(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nobs,nincl
       else
          read(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nobs
       end if
    end if
    read(10,*)t,dt,frq,dsp
    poro=.false.; visco=.false.
    if (stype=="implicit-p" .or. stype=="implicit-pv") poro=.true.
    if (stype=="implicit-v" .or. stype=="implicit-pv") visco=.true.
    if (stype=="fault-p" .or. stype=="fault-pv") poro=.true.
    if (stype=="fault-v" .or. stype=="fault-pv") visco=.true.
    ! Dynamic run time before jumping back to static model
    if (fault .or. galf) then 
       lm_str=.false. ! Calculate LM to stress ratio: lm*r_f2s = str
       read(10,*)t_dyn,dt_dyn,frq_dyn,t_lim,dsp_hyb,Xflt,bod_frc,hyb,rsf,init
       if (nceqs-nceqs_ncf>0) lm_str=.true.
       if (galf) dt_dyn=dt
       ! One time (dt = 24hr) fluid source for initial pressure condition 
       if (init==1 .and. poro) then
          fdt=dt/dble(3600*24)
          dt=dble(3600*24)
       else
          init=0 ! Only for poro
       end if
    end if
    if (hyb==1 .and. rsf==0) read(10,*)frq_wave,frq_slip 
    if (hyb==1 .and. rsf==1) read(10,*)frq_wave,frq_slip,v_bg
    if (dt==f0) dt=f1
    if (explicit .or. fault .or. galf) then
       if (Xflt>1) then
          read(10,*)alpha,beta,nxflt
       else
          read(10,*)alpha,beta
      end if
    end if   
  end subroutine ReadParameters

  ! Setup implicit solver
  subroutine SetupKSPSolver
    implicit none
    call KSPSetType(Krylov,"gmres",ierr)
    call KSPGetPC(Krylov,PreCon,ierr)
    call PCSetType(PreCon,"asm",ierr)
#if (PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR<=4)
    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Double_Precision,        &
       Petsc_Default_Double_Precision,Petsc_Default_Integer,ierr)
#else
    call KSPSetTolerances(Krylov,1.0D-9,Petsc_Default_Real,Petsc_Default_Real, &
       Petsc_Default_Integer,ierr)
#endif
    call KSPSetFromOptions(Krylov,ierr)
  end subroutine SetupKSPSolver

  ! Scatter U and get all local values
  subroutine GetVec_U
    implicit none
    if (dyn) then
       call GetVec(Scatter_dyn,Vec_U_dyn,Seq_U_dyn,uu_dyn)
    else
       call GetVec(Scatter,Vec_U,Seq_U,uu)
    end if
  end subroutine GetVec_U

  subroutine GetVec_S
    implicit none
    call GetVec(Scatter_u,Vec_SS,Seq_SS,ss)
    call GetVec(Scatter_u,Vec_SH,Seq_SH,sh)
  end subroutine GetVec_S

  subroutine GetVec_dip_nrm
    implicit none
    call GetVec(Scatter_u,Vec_dip,Seq_dip,dip)
    call GetVec(Scatter_u,Vec_nrm,Seq_nrm,nrm)
  end subroutine GetVec_dip_nrm

  subroutine GetVec_f2s_seq
    implicit none
    call GetVec(Scatter_u,Vec_f2s,Seq_f2s,f2s)
  end subroutine GetVec_f2s_seq

  subroutine GetVec_fp
    implicit none
    call GetVec(Scatter_u,Vec_fp,Seq_fp,fp)
  end subroutine GetVec_fp

  subroutine GetVec_fl
    implicit none
    call GetVec(Scatter_u,Vec_fl,Seq_fl,fl)
  end subroutine GetVec_fl

  subroutine GetVec_flc
    implicit none
    call GetVec(Scatter_u,Vec_flc,Seq_flc,flc)
  end subroutine GetVec_flc

  subroutine GetVec_qu
    implicit none
    call GetVec(Scatter_q,Vec_qu,Seq_qu,qu)
  end subroutine GetVec_qu

  subroutine GetVec_ql
    implicit none
    call GetVec(Scatter_q,Vec_ql,Seq_ql,ql)
  end subroutine GetVec_ql

  subroutine GetVec(Scatter,Vec_X,Vec_Y,array)
    implicit none
    Vec :: Vec_X, Vec_Y
    VecScatter :: Scatter
    real(8) :: array(:)
    real(8),pointer :: pntr(:)
    call VecScatterBegin(Scatter,Vec_X,Vec_Y,Insert_Values,Scatter_Forward,ierr)
    call VecScatterEnd(Scatter,Vec_X,Vec_Y,Insert_Values,Scatter_Forward,ierr)
    call VecGetArrayF90(Vec_Y,pntr,ierr)
    array=pntr
    call VecRestoreArrayF90(Vec_Y,pntr,ierr)
  end subroutine GetVec

end program main
