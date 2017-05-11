! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module elems

  use utils
  implicit none
  character(3) :: eltype
  integer :: dmn,npel,nps,nip,vtkid,eldof,eldofp,cdmn

contains

  ! Sets element specific constants
  subroutine InitializeElement
    implicit none
    select case(eltype)
    case("tri"); dmn=2; npel=3; nps=2; nip=1; vtkid=5
    case("qua"); dmn=2; npel=4; nps=2; nip=4; vtkid=9
    case("tet"); dmn=3; npel=4; nps=3; nip=1; vtkid=10
    case("hex"); dmn=3; npel=8; nps=4; nip=8; vtkid=12
    end select
    eldof=dmn*npel; eldofp=npel
    if (dmn==2) cdmn=3; if (dmn==3) cdmn=6
  end subroutine InitializeElement

  ! Returns quadrature ipoints and weights
  subroutine SamPts(ipoint,weight)
    implicit none
    real(8) :: ipoint(:,:),weight(:)
    if (eltype=="tri") call SamPtsTri(ipoint,weight)
    if (eltype=="qua") call SamPtsQad(ipoint,weight)
    if (eltype=="tet") call SamPtsTet(ipoint,weight)
    if (eltype=="hex") call SamPtsHex(ipoint,weight)
  end subroutine SamPts

  ! Computes 'N', the shape functions
  subroutine ShapeFunc(N,coord)
    implicit none
    real(8) :: N(:),coord(:)
    if (eltype=="tri") call ShapeFuncTri(N,coord)
    if (eltype=="qua") call ShapeFuncQad(N,coord)
    if (eltype=="tet") call ShapeFuncTet(N,coord)
    if (eltype=="hex") call ShapeFuncHex(N,coord)
  end subroutine ShapeFunc

  ! Computes 'dN', derivative of shape functions
  subroutine ShapeFuncd(dN,coord)
    implicit none
    real(8) :: dN(:,:),coord(:)
    if (eltype=="tri") call ShapeFuncdTri(dN,coord)
    if (eltype=="qua") call ShapeFuncdQad(dN,coord)
    if (eltype=="tet") call ShapeFuncdTet(dN,coord)
    if (eltype=="hex") call ShapeFuncdHex(dN,coord)
  end subroutine ShapeFuncd

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodes(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(:),side,snodes(:)
    real(8) :: ecoords(:,:),area
    if (eltype=="tri") call EdgeAreaNodesTri(enodes,ecoords,side,area,snodes)
    if (eltype=="qua") call EdgeAreaNodesQad(enodes,ecoords,side,area,snodes)
    if (eltype=="tet") call EdgeAreaNodesTet(enodes,ecoords,side,area,snodes)
    if (eltype=="hex") call EdgeAreaNodesHex(enodes,ecoords,side,area,snodes)
  end subroutine EdgeAreaNodes

  ! Get 'edge' with Winkler foundation
  subroutine GetWinklerEdge(elbc,dir,side)
    implicit none
    integer :: elbc(:,:),dir,side
    if (eltype=="tri") call GetWinklerEdgeTri(elbc,dir,side)
    if (eltype=="qua") call GetWinklerEdgeQad(elbc,dir,side)
    if (eltype=="tet") call GetWinklerEdgeTet(elbc,dir,side)
    if (eltype=="hex") call GetWinklerEdgeHex(elbc,dir,side)
  end subroutine GetWinklerEdge

! Linear Tri ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsTri(ipoint,weight)
    implicit none
    real(8) :: ipoint(:,:),weight(:)
    if (nip==1) then
       ipoint(1,:)=(/f1/f3,f1/f3/)
       weight=0.5d0 ! Weight = Wi = wi*wj
    end if
    if (nip==3) then
       ipoint(1,:)=(/f1/f6,f1/f6/)
       ipoint(2,:)=(/f4/f6,f1/f6/)
       ipoint(3,:)=(/f1/f6,f4/f6/)
       weight=f1/f6 ! Weight = Wi = wi*wj
    end if
  end subroutine SamPtsTri

  ! Computes 'N', the shape functions
  subroutine ShapeFuncTri(N,coord)
    implicit none
    real(8) :: N(3),coord(2),eta,nu
    eta=coord(1); nu=coord(2)
    N(1)=f1-eta-nu; N(2)=eta; N(3)=nu
  end subroutine ShapeFuncTri

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta' and 'nu'
  subroutine ShapeFuncdTri(dN,coord)
    implicit none
    real(8) :: dN(2,3),coord(2),e,n
    e=coord(1); n=coord(2)
    dN(1,:)=(/-f1,f1,f0/) ! dN/de
    dN(2,:)=(/-f1,f0,f1/) ! dN/dn
  end subroutine ShapeFuncdTri

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesTri(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(3),side,snodes(2),j(2)
    real(8) :: ecoords(3,2),area
    if (side==1) then; j(1)=1; j(2)=2
    else if (side==2) then; j(1)=2; j(2)=3
    else if (side==3) then; j(1)=3; j(2)=1
    end if
    snodes=enodes(j)
    area=sqrt((ecoords(j(1),1)-ecoords(j(2),1))**2+                            &
              (ecoords(j(1),2)-ecoords(j(2),2))**2)
  end subroutine EdgeAreaNodesTri

  ! Get 'edge' with Winkler foundation
  subroutine GetWinklerEdgeTri(elbc,dir,side)
    implicit none
    integer :: elbc(3,2),dir,side
    side=0
    if (elbc(1,dir)==-1 .and. elbc(2,dir)==-1) side=1
    if (elbc(2,dir)==-1 .and. elbc(3,dir)==-1) side=2
    if (elbc(1,dir)==-1 .and. elbc(3,dir)==-1) side=3
  end subroutine GetWinklerEdgeTri

! Linear Quad ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsQad(ipoint,weight)
    implicit none
    real(8) :: ipoint(4,2),weight(4)
    ipoint(1,:)=(/-sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoint(2,:)=(/ sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoint(3,:)=(/ sqrt(f1/f3), sqrt(f1/f3)/)
    ipoint(4,:)=(/-sqrt(f1/f3), sqrt(f1/f3)/)
    weight=f1 ! Weight = Wi = wi*wj
  end subroutine SamPtsQad

  ! Computes 'N', the shape functions
  subroutine ShapeFuncQad(N,coord)
    implicit none
    real(8) :: N(4),coord(2),eta,nu
    eta=coord(1); nu=coord(2)
    N(1)=0.25d0*(f1-eta)*(f1-nu); N(2)=0.25d0*(f1+eta)*(f1-nu)
    N(3)=0.25d0*(f1+eta)*(f1+nu); N(4)=0.25d0*(f1-eta)*(f1+nu)
  end subroutine ShapeFuncQad

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta' and 'nu'
  subroutine ShapeFuncdQad(dN,coord)
    implicit none
    real(8) :: dN(2,4),coord(2),e,n
    e=coord(1); n=coord(2)
    dN(1,:)=0.25d0*(/-(f1-n), (f1-n),(f1+n),-(f1+n)/) ! dN/de
    dN(2,:)=0.25d0*(/-(f1-e),-(f1+e),(f1+e), (f1-e)/) ! dN/dn
  end subroutine ShapeFuncdQad

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesQad(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(4),side,snodes(2),j(2)
    real(8) :: ecoords(4,2),area
    if (side==1) then; j(1)=1; j(2)=2
    else if (side==2) then; j(1)=2; j(2)=3;
    else if (side==3) then; j(1)=3; j(2)=4;
    else if (side==4) then; j(1)=4; j(2)=1;
    end if
    snodes=enodes(j)
    area=sqrt((ecoords(j(1),1)-ecoords(j(2),1))**2+                            &
              (ecoords(j(1),2)-ecoords(j(2),2))**2)
  end subroutine EdgeAreaNodesQad

  ! Get 'edge' with Winkler foundation
  subroutine GetWinklerEdgeQad(elbc,dir,side)
    implicit none
    integer :: elbc(4,2),dir,side
    side=0
    if (elbc(1,dir)==-1 .and. elbc(2,dir)==-1) side=1
    if (elbc(2,dir)==-1 .and. elbc(3,dir)==-1) side=2
    if (elbc(3,dir)==-1 .and. elbc(4,dir)==-1) side=3
    if (elbc(4,dir)==-1 .and. elbc(1,dir)==-1) side=4
  end subroutine GetWinklerEdgeQad

! Linear Tet ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsTet(ipoint,weight)
    implicit none
    real(8) :: ipoint(:,:),weight(:)
    real(8),parameter :: f24=24.0d0,a=0.13819660d0,b=0.58541020d0
    if (nip==1) then
       ipoint(1,:)=(/f1/f4,f1/f4,f1/f4/)
       weight(1)=f1/f6 ! Weight = Wi = wi*wj*wk
    end if
    if (nip==4) then
       ipoint(1,:)=(/a,a,a/) ! a=(5.0-sqrt(5.0))/20.0
       ipoint(2,:)=(/b,a,a/) ! b=(5.0+3.0*sqrt(5.0))/20.0
       ipoint(3,:)=(/a,b,a/)
       ipoint(4,:)=(/a,a,b/)
       weight=f1/f24 ! Weight = Wi = wi*wj*wk
    end if
  end subroutine SamPtsTet

  ! Computes 'N', the shape functions
  subroutine ShapeFuncTet(N,coord)
    implicit none
    real(8) :: N(4),coord(3),eta,nu,psi
    eta=coord(1); nu=coord(2); psi=coord(3)
    N(1)=f1-eta-nu-psi; N(2)=eta; N(3)=nu; N(4)=psi
  end subroutine ShapeFuncTet

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta','nu' and 'psi'
  subroutine ShapeFuncdTet(dN,coord)
    implicit none
    real(8) :: dN(3,4),coord(3),e,n,s
    e=coord(1); n=coord(2); s=coord(3)
    dN(1,:)=(/-f1,f1,f0,f0/) ! dN/de
    dN(2,:)=(/-f1,f0,f1,f0/) ! dN/dn
    dN(3,:)=(/-f1,f0,f0,f1/) ! dN/ds
  end subroutine ShapeFuncdTet

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesTet(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(4),side,snodes(3),j(3)
    real(8) :: ecoords(4,3),area
    if (side==1) then; j(1)=1; j(2)=2; j(3)=4
    else if (side==2) then; j(1)=2; j(2)=3; j(3)=4
    else if (side==3) then; j(1)=1; j(2)=3; j(3)=4
    else if (side==4) then; j(1)=1; j(2)=2; j(3)=3
    end if
    snodes=enodes(j)
    call TriArea(ecoords(j(1),1),ecoords(j(1),2),ecoords(j(1),3),              &
                 ecoords(j(2),1),ecoords(j(2),2),ecoords(j(2),3),              &
                 ecoords(j(3),1),ecoords(j(3),2),ecoords(j(3),3),area)
  end subroutine EdgeAreaNodesTet

  ! Get 'edge' with Winkler foundation
  subroutine GetWinklerEdgeTet(elbc,dir,side)
    implicit none
    integer :: elbc(4,3),dir,side
    side=0
    if (elbc(1,dir)==-1 .and. elbc(2,dir)==-1 .and. elbc(4,dir)==-1) side=1
    if (elbc(2,dir)==-1 .and. elbc(3,dir)==-1 .and. elbc(4,dir)==-1) side=2
    if (elbc(1,dir)==-1 .and. elbc(3,dir)==-1 .and. elbc(4,dir)==-1) side=3
    if (elbc(1,dir)==-1 .and. elbc(2,dir)==-1 .and. elbc(3,dir)==-1) side=4
  end subroutine GetWinklerEdgeTet

! Linear Hex ...

  ! Returns quadrature ipoints and weights
  subroutine SamPtsHex(ipoint,weight)
    implicit none
    real(8) :: ipoint(8,3),weight(8)
    ipoint(1,:)=(/-sqrt(f1/f3),-sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoint(2,:)=(/ sqrt(f1/f3),-sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoint(3,:)=(/ sqrt(f1/f3), sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoint(4,:)=(/-sqrt(f1/f3), sqrt(f1/f3),-sqrt(f1/f3)/)
    ipoint(5,:)=(/-sqrt(f1/f3),-sqrt(f1/f3), sqrt(f1/f3)/)
    ipoint(6,:)=(/ sqrt(f1/f3),-sqrt(f1/f3), sqrt(f1/f3)/)
    ipoint(7,:)=(/ sqrt(f1/f3), sqrt(f1/f3), sqrt(f1/f3)/)
    ipoint(8,:)=(/-sqrt(f1/f3), sqrt(f1/f3), sqrt(f1/f3)/)
    weight=f1 ! Weight = Wi = wi*wj*wk
  end subroutine SamPtsHex

  ! Computes 'N', the shape functions
  subroutine ShapeFuncHex(N,coord)
    implicit none
    real(8) :: N(8),coord(3),eta,nu,psi
    real(8),parameter :: c=0.125d0
    eta=coord(1); nu=coord(2); psi=coord(3)
    N(1)=c*(f1-eta)*(f1-nu)*(f1-psi); N(2)=c*(f1+eta)*(f1-nu)*(f1-psi)
    N(3)=c*(f1+eta)*(f1+nu)*(f1-psi); N(4)=c*(f1-eta)*(f1+nu)*(f1-psi)
    N(5)=c*(f1-eta)*(f1-nu)*(f1+psi); N(6)=c*(f1+eta)*(f1-nu)*(f1+psi)
    N(7)=c*(f1+eta)*(f1+nu)*(f1+psi); N(8)=c*(f1-eta)*(f1+nu)*(f1+psi)
  end subroutine ShapeFuncHex

  ! Computes 'dN', derivative of shape functions w.r.t. 'eta','nu' and 'psi'
  subroutine ShapeFuncdHex(dN,coord)
    implicit none
    real(8) :: dN(3,8),coord(3),e,n,s
    real(8),parameter :: c=0.125d0
    e=coord(1); n=coord(2); s=coord(3)
    dN(1,:)=c*(/-(f1-n)*(f1-s), (f1-n)*(f1-s), (f1+n)*(f1-s),-(f1+n)*(f1-s),   &
       -(f1-n)*(f1+s), (f1-n)*(f1+s), (f1+n)*(f1+s),-(f1+n)*(f1+s)/) ! dN/de
    dN(2,:)=c*(/-(f1-e)*(f1-s),-(f1+e)*(f1-s), (f1+e)*(f1-s), (f1-e)*(f1-s),   &
       -(f1-e)*(f1+s),-(f1+e)*(f1+s), (f1+e)*(f1+s), (f1-e)*(f1+s)/) ! dN/dn
    dN(3,:)=c*(/-(f1-e)*(f1-n),-(f1+e)*(f1-n),-(f1+e)*(f1+n),-(f1-e)*(f1+n),   &
        (f1-e)*(f1-n), (f1+e)*(f1-n), (f1+e)*(f1+n), (f1-e)*(f1+n)/) ! dN/ds
  end subroutine ShapeFuncdHex

  ! Computes edge 'area' and list of nodes 'snodes' in the edge
  subroutine EdgeAreaNodesHex(enodes,ecoords,side,area,snodes)
    implicit none
    integer :: enodes(8),side,snodes(4),j(4)
    real(8) :: ecoords(8,3),area
    if (side==1) then; j(1)=1; j(2)=2; j(3)=6; j(4)=5
    else if (side==2) then; j(1)=2; j(2)=3; j(3)=7; j(4)=6
    else if (side==3) then; j(1)=3; j(2)=4; j(3)=8; j(4)=7
    else if (side==4) then; j(1)=4; j(2)=1; j(3)=5; j(4)=8
    else if (side==5) then; j(1)=1; j(2)=2; j(3)=3; j(4)=4
    else if (side==6) then; j(1)=5; j(2)=6; j(3)=7; j(4)=8
    end if
    snodes=enodes(j)
    call QuadArea(ecoords(j(1),1),ecoords(j(1),2),ecoords(j(1),3),             &
                  ecoords(j(2),1),ecoords(j(2),2),ecoords(j(2),3),             &
                  ecoords(j(3),1),ecoords(j(3),2),ecoords(j(3),3),             &
                  ecoords(j(4),1),ecoords(j(4),2),ecoords(j(4),3),area)
  end subroutine EdgeAreaNodesHex

  ! Get 'edge' with Winkler foundation
  subroutine GetWinklerEdgeHex(elbc,dir,side)
    implicit none
    integer :: elbc(8,3),dir,side
    side=0
    if (elbc(1,dir)==-1 .and. elbc(2,dir)==-1 .and. elbc(6,dir)==-1 .and.      &
        elbc(5,dir)==-1) side=1
    if (elbc(2,dir)==-1 .and. elbc(3,dir)==-1 .and. elbc(7,dir)==-1 .and.      &
        elbc(6,dir)==-1) side=2
    if (elbc(3,dir)==-1 .and. elbc(4,dir)==-1 .and. elbc(8,dir)==-1 .and.      &
        elbc(7,dir)==-1) side=3
    if (elbc(4,dir)==-1 .and. elbc(1,dir)==-1 .and. elbc(5,dir)==-1 .and.      &
        elbc(8,dir)==-1) side=4
    if (elbc(1,dir)==-1 .and. elbc(2,dir)==-1 .and. elbc(3,dir)==-1 .and.      &
        elbc(4,dir)==-1) side=5
    if (elbc(5,dir)==-1 .and. elbc(6,dir)==-1 .and. elbc(7,dir)==-1 .and.      &
        elbc(8,dir)==-1) side=6
  end subroutine GetWinklerEdgeHex

end module elems
