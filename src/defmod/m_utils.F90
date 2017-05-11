! Copyright (C) 2010-2015 ../AUTHORS. All rights reserved.
! This file is part of Defmod. See ../COPYING for license information.

module utils

  implicit none
  real(8),parameter :: f0=0.0d0,f1=1.0d0,f2=2.0d0,f3=3.0d0,f4=4.0d0,f6=6.0d0,  &
     km2m=1000.0d0,gravity=9.80665d0

contains

  ! Computes area of a triangle in 3D space
  subroutine TriArea(x1,y1,z1,x2,y2,z2,x3,y3,z3,area)
    implicit none
    real(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,area,m1(3,3),m2(3,3),m3(3,3),d1,d2,d3
    m1=reshape((/x1,x2,x3,y1,y2,y3,f1,f1,f1/),shape=(/3,3/))
    m2=reshape((/y1,y2,y3,z1,z2,z3,f1,f1,f1/),shape=(/3,3/))
    m3=reshape((/z1,z2,z3,x1,x2,x3,f1,f1,f1/),shape=(/3,3/))
    call MatDet(m1,d1); call MatDet(m2,d2); call MatDet(m3,d3)
    area=0.5d0*sqrt((abs(d1))**2+(abs(d2))**2+(abs(d3))**2)
  end subroutine TriArea

  ! Computes area of a quad in 3D space
  subroutine QuadArea(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,area)
    implicit none
    real(8) :: x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,area,ta1,ta2
    call TriArea(x1,y1,z1,x2,y2,z2,x3,y3,z3,ta1) ! 1st triangle
    call TriArea(x1,y1,z1,x3,y3,z3,x4,y4,z4,ta2) ! 2nd triangle
    area=abs(ta1)+abs(ta2)
  end subroutine QuadArea

  ! Computes inverse of a matrix
  subroutine MatInv(A,Ainv)
    implicit none
    real(8) :: A(:,:),Ainv(:,:)
    real(8),allocatable :: work(:)
    integer :: n,info
    integer,allocatable :: ipiv(:)
    Ainv=A
    n=size(Ainv,1)
    allocate(ipiv(n),work(n))
    call dgetrf(n,n,Ainv,n,ipiv,info)
    call dgetri(n,Ainv,n,ipiv,work,n,info)
    deallocate(ipiv,work)
  end subroutine MatInv

  ! Computes determinant of a matrix of size 2/3
  subroutine MatDet(A,det)
    implicit none
    real(8) :: A(:,:),det
    if (size(A,1)==2) det=A(1,1)*A(2,2)-A(1,2)*A(2,1)
    if (size(A,1)==3) det=A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2))-                &
       A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1))+A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))
  end subroutine MatDet

end module utils
