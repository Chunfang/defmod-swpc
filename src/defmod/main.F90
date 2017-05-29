! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

program main

#include <petscversion.h>

  use global
  use FD2FE
  implicit none
#include "petsc.h"
  character(256) :: input_file,viz,fd
  logical :: l,v,w
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
#else
  call PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-f',      &
     input_file,l,ierr)
  call PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-ss',     &
     viz,v,ierr)
  call PetscOptionsGetString(Petsc_Null_Object,Petsc_Null_Character,'-fd',     &
     fd,w,ierr)
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
     if (fdout==1) call PrintMsg("Runing in FE-FD mixed mode.")
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

  ! Set element specific constants
  call InitializeElement
  p=0; ef_eldof=eldof
  if (poro) then
     p=1; ef_eldof=eldof+eldofp
     if (eltype=="tri") nip=3; if (eltype=="tet") nip=4
  elseif (fault) then
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
  if (rank==0) n=1
  if (rank/=0) n=sum(work(0:rank-1))+1
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
        read(10,*)coords(j,:),bc(j,:); j=j+1
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
  if (fault) allocate(indxmap_u(dmn*n,2))
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
  if (fault) then
     do i=1,n
        do j=1,dmn
           indxmap_u(dmn*i-j+1,:)=dmn*nl2g(i,:)-j
        end do
     end do
  end if
  deallocate(work)
  ! Read material data assuming nels >> nmts, i.e., all ranks store all data
  if (fault) then
     allocate(mat(nmts,5+4*p+init+2))
  else
     allocate(mat(nmts,5+4*p))
  end if
  do i=1,nmts
     read(10,*)mat(i,:)
  end do
  deallocate(epart,npart)

  ! Initialize local element variables and global U
  allocate(ipoint(nip,dmn),weight(nip),k(ef_eldof,ef_eldof),m(eldof,eldof),    &
     f(ef_eldof),indx(ef_eldof),enodes(npel),ecoords(npel,dmn),vvec(dmn+p))
  if (fault) allocate(k_dyn(eldof,eldof),indx_dyn(eldof))
  call SamPts(ipoint,weight)
  n=(dmn+p)*nnds; if (stype/="explicit") n=n+nceqs
  call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n,Vec_U,ierr)
  ! Global U for dynamic run
  if (fault) then
     n_dyn=dmn*nnds
     call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,n_dyn,Vec_U_dyn,ierr)
  end if
  if (visco .or. (fault .and. lm_str==1)) allocate(stress(nels,nip,cdmn))

  ! Set scaling constants
  wt=f2*exp((log(maxval(mat(:,1)))+log(minval(mat(:,1))))/f2)
  if (dmn==3) wt=wt*km2m
  if (poro) then
     scale=sqrt(exp((log(maxval(mat(:,1)))+log(minval(mat(:,1))))/f2)/         &
                exp((log(maxval(mat(:,6)))+log(minval(mat(:,6))))/f2)/         &
                dt/km2m)
  end if

  ! Form stiffness matrix
  call PrintMsg("Forming [K] ...")
  nodal_bw=(dmn+p)*(nodal_bw+1)
  call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,n,nodal_bw,   &
     Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K,ierr)
  call MatSetOption(Mat_K,Mat_New_Nonzero_Allocation_Err,Petsc_False,ierr)
  if (fault) then ! Dynamic K
     nodal_bw=nodal_bw/(dmn+p)*dmn
     call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n_dyn,n_dyn, &
        nodal_bw,Petsc_Null_Integer,nodal_bw,Petsc_Null_Integer,Mat_K_dyn,ierr)
     call MatSetOption(Mat_K_dyn,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
        ierr)
  end if
  do i=1,nels
     if (poro) then
        call FormLocalK(i,k,indx,"Kp")
     else
        call FormLocalK(i,k,indx,"Ke")
     end if
     indx=indxmap(indx,2)
     call MatSetValues(Mat_K,ef_eldof,indx,ef_eldof,indx,k,Add_Values,ierr)
     if (fault) then
        dyn=.true.
        call FormLocalK(i,k_dyn,indx_dyn,"Ke")
        indx_dyn=indxmap_u(indx_dyn,2)
        call MatSetValues(Mat_K_dyn,eldof,indx_dyn,eldof,indx_dyn,k_dyn,       &
           Add_Values,ierr)
        dyn=.false.
     end if
  end do

  ! Initialize and form mass matrix and its inverse
  if (stype=="explicit" .or. fault) then
     call PrintMsg("Forming [M] & [M]^-1 ...")
     if (fault) then
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
        if (fault) then
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
     if (fault) then
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
     call MatScale(Mat_M,alpha,ierr)
  end if

  ! Allocate arrays to store loading history
  allocate(cval(nceqs,3),fnode(nfrcs),fval(nfrcs,dmn+p+2),telsd(ntrcs,2),      &
     tval(ntrcs,dmn+p+2)); cval=f0; fval=f0; tval=f0

  ! Account for constraint eqn's
  if (nceqs>0) then
     call PrintMsg("Applying constraints ...")
     if (stype=="explicit") then
        call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,Petsc_Decide,n,        &
           nceqs,3,Petsc_Null_Integer,3,Petsc_Null_Integer,Mat_Gt,ierr)
        call MatSetOption(Mat_Gt,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
           ierr)
     elseif(fault .and. nceqs-nceqs_ncf>0 .and. hyb>0) then
        call VecCreateMPI(Petsc_Comm_World,Petsc_Decide,nceqs_ncf/(dmn+p)+nfnd,&
           Vec_lmnd,ierr)
        call VecGetLocalSize(Vec_lmnd,n_lmnd,ierr)
        call VecGetOwnershipRange(Vec_lmnd,lmnd0,j,ierr)
        Call VecDestroy(Vec_lmnd,ierr)
        ! Dofs of one fault node are not split by different ranks
        call VecCreateMPI(Petsc_Comm_World,n_lmnd*dmn,(nceqs_ncf/(dmn+p)+nfnd)*&
           dmn,Vec_lambda_sta,ierr)
        call VecDuplicate(Vec_lambda_sta,Vec_lambda_sta0,ierr)
        call VecZeroEntries(Vec_lambda_sta,ierr)
        call VecZeroEntries(Vec_lambda_sta0,ierr)
        call MatCreateAIJ(Petsc_Comm_World,Petsc_Decide,n_lmnd*dmn,dmn*nnds,   &
           (nceqs_ncf/(dmn+p)+nfnd)*dmn,5,Petsc_Null_Integer,5,                &
           Petsc_Null_Integer,Mat_Gt,ierr)
        call MatSetOption(Mat_Gt,Mat_New_Nonzero_Allocation_Err,Petsc_False,   &
           ierr)
        call MatZeroEntries(Mat_Gt,ierr)
        allocate(flt_slip(n_lmnd*dmn),tot_flt_slip(n_lmnd*dmn),                &
           qs_flt_slip(n_lmnd*dmn))
        qs_flt_slip=f0
     end if
     if (rank==0) open(15,file="cnstrns.tmp",status="replace")
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

     ! Read fault orientation vectors
     if (fault .or. gf) then
        allocate(node_pos(nfnd),node_neg(nfnd),vecf(nfnd,dmn*dmn),fc(nfnd),    &
           fcd(nfnd),dc(nfnd),perm(nfnd),vvec_all(2*nfnd*dmn,dmn),             &
           node_all(2*nfnd*dmn),st_init(nfnd,dmn),xfnd(nfnd,dmn),frc(nfnd),    &
           coh(nfnd),dcoh(nfnd))    
        if (rsf==1) allocate(rsfb0(nfnd),rsfV0(nfnd),rsfdtau0(nfnd),rsfa(nfnd),&
            rsfb(nfnd),rsfL(nfnd),rsftheta(nfnd))
        do i=1,nfnd
           if (poro) then
              if (rsf==1) then
                 read(10,*) node_pos(i),node_neg(i),vecf(i,:),rsfb0(i),        &
                    rsfV0(i),rsfdtau0(i),rsfa(i),rsfb(i),rsfL(i),rsftheta(i),  &
                    perm(i),st_init(i,:),xfnd(i,:),frc(i),coh(i),              &
                    dcoh(i)
              else
                 read(10,*) node_pos(i),node_neg(i),vecf(i,:),fc(i),fcd(i),    &
                    dc(i),perm(i),st_init(i,:),xfnd(i,:),frc(i),               &
                    coh(i),dcoh(i)
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
     end if

     ! Create variables for rotated constrain matrix
     if (fault .or. gf) then
        do i=1,nfnd*dmn
           do j=1,2
              read(10,*)vvec(1:dmn),node; node=nmap(node)
              vvec_all(2*(i-1)+j,:)=vvec(1:dmn); node_all(2*(i-1)+j)=node
           end do
        end do
     end if
     if (rank==0) call ApplyConstraints
     if (stype=="explicit" .and. .not. gf) then
        call MatAssemblyBegin(Mat_Gt,Mat_Final_Assembly,ierr)
        call MatAssemblyEnd(Mat_Gt,Mat_Final_Assembly,ierr)
     end if

  ! For dummy fault without split nodes
  elseif (fault) then
     call PrintMsg("Dummy fault placed ...")
     allocate(node_pos(nfnd),vecf(nfnd,dmn*dmn),fc(nfnd),perm(nfnd))
     do i=1,nfnd
        if (poro) then
           read(10,*) node_pos(i),vecf(i,:),fc(i),perm(i)
        else
           read(10,*) node_pos(i),vecf(i,:),fc(i)
        end if
        node_pos(i)=nmap(node_pos(i))
     end do
  end if

  call MatAssemblyBegin(Mat_K,Mat_Final_Assembly,ierr)
  call MatAssemblyEnd(Mat_K,Mat_Final_Assembly,ierr)
  if (fault .and. nceqs>0) then
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
        tot_uu_obs=f0
     end if
     n_log_dyn=0
  end if
  
  ! FD domain grid bocks containing fault nodes, xgp, idgp(_loc), gpl2g, gpnlst, gpshape
  if (nceqs-nceqs_ncf>0 .and. fdout==1) then
     call PrintMsg("Locating FD grid points ...")
     call FDInit
     call GetFDFnd 
     call GetObsNd("fd")
     deallocate(xgp,idgp)  
     if (ngp_loc>0) allocate(uu_fd(ngp_loc,dmn))
     call NndFE2FD
     call MatFE2FD
  end if

  ! Account for absorbing boundaries
  call PrintMsg("Absorbing bundary ...")
  if (stype=="explicit" .or. (fault .and. nceqs>0)) then
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
              !call MatSetValues(Mat_M,eldof,indx_dyn,eldof,indx_dyn,m,         &
              !   Add_Values,ierr)
              do j1=1,eldof
                 do j2=1,eldof
                    val=m(j1,j2)
                    if (abs(val)>f0) call MatSetValue(Mat_M,indx_dyn(j1),      &
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
     call MatAssemblyBegin(Mat_M,Mat_Final_Assembly,ierr)
     call MatAssemblyEnd(Mat_M,Mat_Final_Assembly,ierr)
  end if
  close(10)
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
  if (fault) then
     if (lm_str==1) then
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
     if (lm_str==1) then
        allocate(ss(j),sh(j),f2s(j),dip(j),nrm(j))
     end if
  end if

  ! Implicit Solver
  if (stype/="explicit" .and. (.not. fault)) then
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
           call MatGetSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
           call MatZeroEntries(Mat_Kc,ierr)
           allocate(kc(eldofp,eldofp),indxp(eldofp),Hs(eldofp))
           do i=1,nels
              call FormLocalK(i,k,indx,"Kc")
              kc=k(eldof+1:,eldof+1:)
              indxp=indx(eldof+1:); indxp=indxmap(indxp,2)
              indxp=((indxp+1)/(dmn+1))-1
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

  ! Fault/hybrid solver
  if (fault) then
     ! Local to global fault node map
     if (nceqs-nceqs_ncf>0) call GetFltMap
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
     call PrintMsg("Solving ...")
     call VecGetOwnershipRange(Vec_U,j1,j2,ierr)
     !print'(A,I0,A,I0,A,I0)',"  Rank ",rank," has dofs ",j1+1," to ",j2
     if (rank==nprcs-1) print'(I0,A,I0,A)',j2+nceqs," dofs on ", nprcs,        &
        " processors."
     call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
     call GetVec_U; tot_uu=tot_uu+uu
     ! Get observation
     if (nobs_loc>0) then
        call GetVec_obs
        tot_uu_obs=tot_uu_obs+uu_obs
     end if
     if (visco .or. lm_str==1) then
        ! Recover stress
        call PrintMsg("Recovering stress ...")
        do i=1,nels
           call RecoverStress(i,stress)
        end do
     end if
     if (t>f0 .and. dt>f0 .and. t>=dt) then
        call VecDuplicate(Vec_U,Vec_Um,ierr) ! U->du & Um->u
        call VecCopy(Vec_U,Vec_Um,ierr)
        call VecGetLocalSize(Vec_U,j,ierr)
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
           call ISCreateGeneral(Petsc_Comm_World,j,work,Petsc_Copy_Values,RI,  &
              ierr)
           j=size(worku)
           call ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,    &
              RIu,ierr)
           call MatGetSubMatrix(Mat_K,RIu,RI,Mat_Initial_Matrix,Mat_H,ierr)
           if (nceqs > 0) then
              j=size(workl)
              call ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values, &
                 RIl,ierr)
           end if
           call MatGetSubMatrix(Mat_K,RI,RI,Mat_Initial_Matrix,Mat_Kc,ierr)
           call MatZeroEntries(Mat_Kc,ierr)
           allocate(kc(eldofp,eldofp),indxp(eldofp),Hs(eldofp))
           do i=1,nels
              call FormLocalK(i,k,indx,"Kc")
              kc=k(eldof+1:,eldof+1:)
              indxp=indx(eldof+1:); indxp=indxmap(indxp,2)
              indxp=((indxp+1)/(dmn+1))-1
              call MatSetValues(Mat_Kc,eldofp,indxp,eldofp,indxp,kc,           &
                 Add_Values,ierr)
           end do
           call MatAssemblyBegin(Mat_Kc,Mat_Final_Assembly,ierr)
           call MatAssemblyEnd(Mat_Kc,Mat_Final_Assembly,ierr)
           call VecGetSubVector(Vec_Um,RI,Vec_Up,ierr)
           if (init==1) call VecDuplicate(Vec_Up,Vec_Up0,ierr) ! Initial p
           call VecDuplicate(Vec_Up,Vec_I,ierr) ! I->KcUp
           call VecCopy(Vec_Up,Vec_I,ierr) ! Hold Up
           call VecDuplicate(Vec_Up,Vec_qu,ierr) ! qu->Htu
           call VecDuplicate(Vec_Up,Vec_ql,ierr)
           call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
           allocate(uup(size(work)))
           ! Initialize space for lambda, p related nodal force
           call VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           call VecDuplicate(Vec_Uu,Vec_fp,ierr) ! fp->Hp
           call VecDuplicate(Vec_Uu,Vec_fl,ierr) ! Ifl->-Gtuul
           call VecCopy(Vec_Uu,Vec_fl,ierr) ! Hold Uu
           call VecDuplicate(Vec_Uu,Vec_flc,ierr)
           if (lm_str==1) then
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
           j=size(indxmap_u,1)
           if (nceqs>0) then
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
           if (nceqs>0) then
              allocate(fl(size(indxmap_u,1)))
              allocate(ql(size(nl2g,1)))
              allocate(flc(size(indxmap_u,1)))
           end if
           if (vout==1) then ! Extract nodal force by p, and fluid source by u
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
              !call LM_s2d 
              ! Extract lambda induced nodal force
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
           if (lm_str==1) then ! Scatter stress to nodes
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
              call PrintMsg("Pore fluid initialization...")
              call VecGetSubVector(Vec_Um,RI,Vec_Up0,ierr)
              ! Zero initial pressure 
              call VecZeroEntries(Vec_Up0,ierr)
              call VecRestoreSubVector(Vec_Um,RI,Vec_Up0,ierr)
              call VecDestroy(Vec_Up0,ierr)
              tot_uu=f0
              call VecZeroEntries(Vec_F,ierr)
              call ApplySource
              call KSPSolve(Krylov,Vec_F,Vec_U,ierr)
              call GetVec_U; tot_uu=tot_uu+uu
              call VecAXPY(Vec_Um,f1,Vec_U,ierr)
              if (vout==1) then
                 if (nceqs>0) then
                    call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                    call VecZeroEntries(Vec_flc,ierr)
                    call GetVec_fcoulomb
                    call GetVec_flc
                    call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
                 end if
                 ! Initial state should be analyzed to see if any initial slip 
                 call WriteOutput_init
              end if
           end if 
           if (nceqs-nceqs_ncf>0) then
              allocate(flt_ss(nfnd,dmn),flt_p(nfnd))
              call GetVec_flt_qs 
              if (rank==0) call WriteOutput_flt_qs
           end if
           call Rscdt(fdt)
        else ! Not poro
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
           call ISCreateGeneral(Petsc_Comm_World,j,worku,Petsc_Copy_Values,    &
              RIu,ierr)
           j=size(workl)
           call ISCreateGeneral(Petsc_Comm_World,j,workl,Petsc_Copy_Values,    &
              RIl,ierr)
           call VecGetSubVector(Vec_Um,RIu,Vec_Uu,ierr)
           call VecDuplicate(Vec_Uu,Vec_fl,ierr) ! Ifl->-Gtuul
           if (nceqs>0) call VecDuplicate(Vec_Uu,Vec_flc,ierr)
           if (lm_str==1) then
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
           j=size(indxmap_u,1)
           if (nceqs>0) then
              call VecCreateSeq(Petsc_Comm_Self,j,Seq_fl,ierr)
              call VecCreateSeq(Petsc_Comm_Self,j,Seq_flc,ierr)
              call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),           &
                 Petsc_Copy_Values,From_u,ierr)
              call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),           &
                 Petsc_Copy_Values,To_u,ierr)
              call VecScatterCreate(Vec_fl,From_u,Seq_fl,To_u,Scatter_u,ierr)
           elseif (lm_str==1) then
              call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,2),           &
                 Petsc_Copy_Values,From_u,ierr)
              call ISCreateGeneral(Petsc_Comm_Self,j,indxmap_u(:,1),           &
                 Petsc_Copy_Values,To_u,ierr)
              call VecScatterCreate(Vec_U,From_u,Seq_U,To_u,Scatter_u,ierr)
           end if
           if (nceqs>0) then
              allocate(fl(size(indxmap_u,1)))
              allocate(flc(size(indxmap_u,1)))
              ! Vector to communicate with dynamic LM
              call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              !call LM_s2d
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
           if (lm_str==1) then ! Scatter stress to nodes
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
              if (nceqs-nceqs_ncf>0) then
                 allocate(flt_ss(nfnd,dmn))
                 call GetVec_flt_qs
                 if (rank==0) call WriteOutput_flt_qs
              end if
           end if
           ! Write output
           if (vout==1) call WriteOutput_x
           !if (rank==0 .and. nobs>0) call WriteOutput_obs
           if (nobs_loc>0) call WriteOutput_obs
        end if ! Poro or not
        ! Solution space is allocated differently for static and dynamic runs, 
        ! so we keep mat_K and Mat_K_dyn separate instead of
        ! call MatGetSubMatrix(Mat_K,RIu,RIu,Mat_Initial_Matrix,Mat_K_dyn,ierr)

        ! Initialize slip indicator and slip history
        if (nceqs>0 .and. hyb>0) then
           allocate(slip(nfnd),slip0(nfnd),slip_sum(nfnd)) 
           slip(:)=0;slip_sum=0
           n_log_wave=0
           n_log_slip=0
           if (rsf==1) allocate(mu_cap(nfnd),rsfv(nfnd),rsf_sta(nfnd_loc))
           trunc=f0
           allocate(mu_hyb(nfnd))
        end if

        ! Start implicit time step
        steps=int(ceiling(t/dt)); t_abs=f0
        dyn=.false.; fail=.false.; n_log=0
        do tstep=1,steps
           t_abs=t_abs+dt
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
              call VecScale(Vec_I,-dt,ierr)
              call VecGetArrayF90(Vec_I,pntr,ierr)
              uup=pntr
              call VecRestoreArrayF90(Vec_I,pntr,ierr)
              j=size(uup)
              call VecSetValues(Vec_F,j,work,uup,Add_Values,ierr)
              call VecRestoreSubVector(Vec_Um,RI,Vec_Up,ierr)
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
              call FaultSlip
              ! Backup QS slip 
              qs_flt_slip=qs_flt_slip+tot_flt_slip
              fail=.false.
           end if
           ! Solve
           call VecAssemblyBegin(Vec_F,ierr)
           call VecAssemblyEnd(Vec_F,ierr)
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
           if (visco .or. lm_str==1) then
              ! Recover stress
              call PrintMsg(" Recovering stress ...")
              do i=1,nels
                 call RecoverVStress(i,stress)
              end do
              call GetVec_Stress
              call GetVec_S
              if (fault) then
                 call GetVec_flt_qs
                 if (rank==0) call WriteOutput_flt_qs
              end if 
           end if
           ! Extract nodal force by p
           if (poro) then
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
              !if (rank==0 .and. nobs>0) call WriteOutput_obs
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
              !if (rank==0 .and. nobs>0) call WriteOutput_obs
              if (nobs_loc>0) call WriteOutput_obs
           end if

           ! Determine if the fault shall fail by LM 
           if (nceqs-nceqs_ncf>0 .and. hyb>0) then
              call VecGetSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              call LM_s2d  
              call VecRestoreSubVector(Vec_Um,RIl,Vec_Ul,ierr)
              ! Determine if the fault shall fail
              call GetSlip_sta
              rslip=real(sum(slip))/real(size(slip))
              if (rank==0) print'(F0.2,A)',rslip*100.0,"% fault critical."
              if (rslip>f0) then ! Failure threshold
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
              if (n_log_dyn==0) then! Create full Gt
                 call GetMat_Gt
                 call MatAssemblyBegin(Mat_Gt,Mat_Final_Assembly,ierr)
                 call MatAssemblyEnd(Mat_Gt,Mat_Final_Assembly,ierr)
              end if
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
              if (rsf==1 .and. tstep>1) then
                 call Rsfv2Dyn 
                 rsf_sta=0
              end if
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
                 if (ngp_loc>0 .and. mod(n_log_dyn,frq_wave)==0 .and. fdout==1) then  
                    call GetVec_fd
                    call WriteOutput_fd
                 end if
                 ! Extract and output temporal fault slip
                 call MatMult(Mat_G,Vec_U_dyn,Vec_Wlm(1),ierr)
                 call VecGetArrayF90(Vec_Wlm(1),pntr,ierr)
                 flt_slip=pntr
                 tot_flt_slip=tot_flt_slip+flt_slip
                 call VecRestoreArrayF90(Vec_Wlm(1),pntr,ierr)
                 if (mod(n_log_dyn,frq_slip)==0) then 
                    call WriteOutput_slip
                    n_log_slip=n_log_slip+1
                 end if
                 ! Export dynamic snapshot
                 dsp_dyn=.true.
                 if (vout==1) then
                    uu_dyn=uu_dyn/dt_dyn
                    if (mod(n_log_dyn,frq_dyn)==0) call WriteOutput
                 end if
                 dsp_dyn=.false.
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
              ! Latest fault stress (Vec_lambda_sta0)
              call GetVec_lambda_hyb(trunc) ! Last time truncation
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
        if (lm_str==1) deallocate(ss,sh)
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
     call VecDestroy(Vec_fl,ierr)
     call VecDestroy(Vec_flc,ierr)
     call VecDestroy(Vec_f2s,ierr)
     call VecDestroy(Vec_ql,ierr)
     call VecDestroy(Vec_U_dyn,ierr)
     call VecDestroy(Vec_SS,ierr)
     call VecDestroy(Vec_SH,ierr)
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
  if (stype=="explicit" .and. t>f0 .and. dt>f0 .and. t>=dt) then
     steps=int(ceiling(t/dt))
     ! Initialize work vectors
     call VecDuplicate(Vec_U,Vec_Um,ierr)
     call VecDuplicate(Vec_U,Vec_Up,ierr)
     call VecDuplicateVecsF90(Vec_U,6,Vec_W,ierr)
     if (nceqs>0) then
        if (gf) call GetMat_Gt
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
                         ! Use a decaying function instead of boxcar ...
                         !if (j2>j1) val=val*(dble(j2-tstep)/dble(j2-j1))**2.5
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
  call MatDestroy(Mat_K,ierr)
  if (visco) deallocate(stress)
  deallocate(coords,nodes,bc,mat,id,k,m,f,indx,ipoint,weight,enodes,ecoords,   &
     vvec,indxmap,tot_uu,uu,cval,fnode,fval,telsd,tval,nl2g)
  ! Delete cnstr.tmp 
  if (rank==0 .and. nceqs>0) then
    open(15, file="cnstrns.tmp",status='old')
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
       stype=="fault-v") fault=.true.
    gf=.false.
    if (stype=="explicit-gf") then
       stype="explicit"; gf=.true.
    end if
    if (fault .or. gf) then
       read(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nfnd,nobs,nceqs_ncf
    else
       read(10,*)nels,nnds,nmts,nceqs,nfrcs,ntrcs,nabcs,nobs
    end if
    read(10,*)t,dt,frq,dsp
    ! Dynamic run time before jumping back to static model
    if (fault) then 
       read(10,*)t_dyn,dt_dyn,frq_dyn,t_lim,dsp_hyb,lm_str,bod_frc,            &
       hyb,rsf,init
       if (init==1) then
          fdt=dt/dble(3600*24)
          dt=dble(3600*24)
       end if
    end if
    if (hyb==1 .and. rsf==0) read(10,*)frq_wave,frq_slip 
    if (hyb==1 .and. rsf==1) read(10,*)frq_wave,frq_slip,v_bg
    poro=.false.; visco=.false.
    if (stype=="implicit-p" .or. stype=="implicit-pv") poro=.true.
    if (stype=="implicit-v" .or. stype=="implicit-pv") visco=.true.
    if (dt==f0) dt=f1
    if (stype=="explicit") read(10,*)alpha,beta
    if (stype=="fault-p" .or. stype=="fault-pv") poro=.true.
    if (stype=="fault-v" .or. stype=="fault-pv") visco=.true.
    if (fault) read(10,*)alpha,beta
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
