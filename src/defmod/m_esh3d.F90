module esh3d

  use utils  
  implicit none

contains 

  ! Eshelby's tensor S(6,6)
  subroutine EshS2(vm,a,S2,PIvec)  
    implicit none  
    real(8) :: vm,a(3),Ifir(3),Isec(3,3),rate,theta,m,B,D,F,E,denom,           &
       S4(3,3,3,3),S2(6,6),PIvec(3)
    rate=a(1)/a(3)
    ! Check geometry category, assuming a(1)>a(2)>a(3)
    if (a(1)-a(2)<1.d-6*a(1) .and. a(2)-a(3)<1.d-6*a(1)) then ! Spherical case 
       Ifir=(f4/f3)*pi
       Isec=0.8*pi*a(1)**2
    elseif (a(1)-a(2)>1.d-6*a(1) .and. a(2)-a(3)<1.d-6*a(1)) then ! Prolate case 
       Ifir(2)=f2*pi*a(1)*a(3)**2/((a(1)**2-a(3)**2)**1.5)                     &
          *(rate*sqrt(rate**2-f1)-acosh(rate))
       Ifir(3)=Ifir(2) 
       Ifir(1)=f4*pi-2*Ifir(2) 
       Isec=f0
       Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)**2-a(2)**2)
       Isec(1,3)=Isec(1,2) 
       Isec(2,1)=Isec(1,2)
       Isec(3,1)=Isec(1,3)
       Isec(1,1)=(f4*pi/a(1)**2-f2*Isec(1,2))/f3
       Isec(2,3)=pi/(a(2)**2)-(Ifir(2)-Ifir(1))/(f4*(a(1)**2-a(2)**2))
       Isec(3,2)=Isec(2,3)
       Isec(2,2)=Isec(2,3)
       Isec(3,3)=Isec(2,3)
    elseif (a(1)-a(2)<1.d-6*a(1) .and. a(2)-a(3)>1.d-6*a(2)) then ! Oblate case
       rate=a(3)/a(1)
       Ifir(1)=f2*pi*a(1)**2*a(3)/((a(1)**2-a(3)**2)**1.5)*(acos(rate)         &
               -rate*sqrt(1-rate**2))
       Ifir(2)=Ifir(1)
       Ifir(3)=f4*pi-2*Ifir(1)
       Isec=f0
       Isec(1,3)=(Ifir(1)-Ifir(3))/(a(3)**2-a(1)**2)
       Isec(3,1)=Isec(1,3)
       Isec(2,3)=Isec(1,3)
       Isec(3,2)=Isec(2,3)
       Isec(1,2)=pi/a(1)**2-Isec(1,3)*0.25
       Isec(2,1)=Isec(1,2)
       Isec(1,1)=Isec(1,2)
       Isec(2,2)=Isec(1,2)
       Isec(3,3)=(f4*pi/a(3)**2-f2*Isec(1,3))/f3
    else ! Ellipsoid case
       theta=asin(sqrt(1-(a(3)/a(1))**2))
       m=(a(1)**2-a(2)**2)/(a(1)**2-a(3)**2)
       call elbd(theta,0.5*pi-theta,f1-m,B,D)
       F=B+D; E=B+(f1-m)*D
       Ifir(1)=(f4*pi*product(a)/((a(1)**2-a(2)**2)*sqrt(a(1)**2-a(3)**2)))    &
          *(F-E)
       Ifir(3)=(f4*pi*product(a)/((a(2)**2-a(3)**2)*sqrt((a(1)**2-a(3)**2))))  &
          *(a(2)*sqrt((a(1)**2-a(3)**2))/(a(1)*a(3))-E)
       Ifir(2)=f4*pi-Ifir(1)-Ifir(3)
       Isec=f0
       Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)**2-a(2)**2)
       Isec(2,3)=(Ifir(3)-Ifir(2))/(a(2)**2-a(3)**2)
       Isec(3,1)=(Ifir(1)-Ifir(3))/(a(3)**2-a(1)**2)
       Isec(2,1)=Isec(1,2)
       Isec(3,2)=Isec(2,3)
       Isec(1,3)=Isec(3,1)
       Isec(1,1)=(f4*pi/a(1)**2-Isec(1,2)-Isec(1,3))/f3
       Isec(2,2)=(f4*pi/a(2)**2-Isec(2,3)-Isec(2,1))/f3
       Isec(3,3)=(f4*pi/a(3)**2-Isec(3,1)-Isec(3,2))/f3
    end if
    denom=f8*pi*(f1-vm);
    S4=f0 ! Sparse tensor
    S4(1,1,1,1)=(f3*a(1)**2*Isec(1,1)+(f1-2*vm)*Ifir(1))/denom
    S4(2,2,2,2)=(f3*a(2)**2*Isec(2,2)+(f1-2*vm)*Ifir(2))/denom
    S4(3,3,3,3)=(f3*a(3)**2*Isec(3,3)+(f1-2*vm)*Ifir(3))/denom
    S4(1,1,2,2)=(a(2)**2*Isec(1,2)-(f1-2*vm)*Ifir(1))/denom
    S4(2,2,3,3)=(a(3)**2*Isec(2,3)-(f1-2*vm)*Ifir(2))/denom
    S4(3,3,1,1)=(a(1)**2*Isec(3,1)-(f1-2*vm)*Ifir(3))/denom
    S4(1,1,3,3)=(a(3)**2*Isec(1,3)-(f1-2*vm)*Ifir(1))/denom
    S4(2,2,1,1)=(a(1)**2*Isec(2,1)-(f1-2*vm)*Ifir(2))/denom
    S4(3,3,2,2)=(a(2)**2*Isec(3,2)-(f1-2*vm)*Ifir(3))/denom
    S4(1,2,1,2)=((a(1)**2+a(2)**2)*Isec(1,2)+(f1-f2*vm)*(Ifir(1)+Ifir(2)))   &
       /(2*denom)
    S4(2,3,2,3)=((a(2)**2+a(3)**2)*Isec(2,3)+(f1-f2*vm)*(Ifir(2)+Ifir(3)))   &
       /(2*denom)
    S4(3,1,3,1)=((a(3)**2+a(1)**2)*Isec(3,1)+(f1-f2*vm)*(Ifir(3)+Ifir(1)))   &
       /(2*denom)
    S4(1,3,1,3)=S4(3,1,3,1)
    ! Original order 
    !S2(1,:)=(/S4(1,1,1,1),f0,f0,S4(1,1,2,2),f0,S4(1,1,3,3)/) 
    !S2(2,:)=(/f0,2*S4(1,2,1,2),f0,f0,f0,f0/)
    !S2(3,:)=(/f0,f0,2*S4(1,3,1,3),f0,f0,f0/)
    !S2(4,:)=(/S4(2,2,1,1),f0,f0,S4(2,2,2,2),f0,S4(2,2,3,3)/)
    !S2(5,:)=(/f0,f0,f0,f0,2*S4(2,3,2,3),f0/)
    !S2(6,:)=(/S4(3,3,1,1),f0,f0,S4(3,3,2,2),f0,S4(3,3,3,3)/)
    ! FE order
    S2(1,:)=(/S4(1,1,1,1),S4(1,1,2,2),S4(1,1,3,3),f0,f0,f0/)
    S2(2,:)=(/S4(2,2,1,1),S4(2,2,2,2),S4(2,2,3,3),f0,f0,f0/)
    S2(3,:)=(/S4(3,3,1,1),S4(3,3,2,2),S4(3,3,3,3),f0,f0,f0/)
    S2(4,:)=(/f0,f0,f0,f2*S4(1,2,1,2),f0,f0/)
    S2(5,:)=(/f0,f0,f0,f0,f2*S4(2,3,2,3),f0/)
    S2(6,:)=(/f0,f0,f0,f0,f0,f2*S4(1,3,1,3)/)
    PIvec(2)=-f2*(Ifir(1)-Ifir(3))/(f8*pi)
    PIvec(3)=-f2*(Ifir(2)-Ifir(1))/(f8*pi)
    PIvec(1)=-f2*(Ifir(3)-Ifir(2))/(f8*pi)
  end subroutine EshS2    
  
  ! Eshelby's tensor D(3,3,3,3)
  subroutine EshD4(vm,a,x,D4,fderphi,tderpsi)
    implicit none
    integer :: i,j,k,l,q,p,r
    real(8) :: vm,a(3),x(3),D4(3,3,3,3),fderphi(3),tderpsi(3,3,3) ! In(out)put
    real(8) :: a_21(3,3),a_22(3,3),coef(1:4,0:4),root(1:4),root_im2,root_im3,  &
               lambda,theta,m,B,D,F,F_21(3,3),F_22(3,3),E,Ifir(3),Isec(3,3),   &
               del,bbar,dbar,ultadelfir(3),ultadelfir_21(3,3),                 &
               ultadelfir_22(3,3),ultadelsec(3,3),fderlambda(3),               &
               fderlambda_21(3,3),c1,c2,c3,fderlambda_22(3,3),diagvals(3,3),   &
               nondiagvals(3,3),fderc1(3),fderc1_21(3,3),fderc1_22(3,3),       &
               fderF(3,3),sderc1(3,3),sderlambda(3,3),fderIfir(3,3),           &
               sderF(3,3,3),zeefir(3),zeesec(3,3),sderIfir(3,3,3),             &
               fderIsec(3,3,3),sderIsec(3,3,3,3),tderlambda(3,3,3),            &
               sderVfir(3,3,3),tderVfir(3,3,3,3),sderphi(3,3),tderphi(3,3,3),  &
               foderpsi(3,3,3,3),premult1,delta1,delta2,delta3,delta4,delta5,  &
               Fvec(3),fderc2(3),fderc2_21(3,3),fderc2_22(3,3)
    coef=f0
    coef(3,3)=f1 ! coefficient of lambds**3 term
    coef(3,2)=a(1)**2+a(2)**2+a(3)**2-(x(1)**2+x(2)**2+x(3)**2) 
    coef(3,1)=a(1)**2*a(2)**2+a(1)**2*a(3)**2+a(2)**2*a(3)**2-((a(2)**2+       &
              a(3)**2)*x(1)**2+(a(1)**2+a(3)**2)*x(2)**2+(a(1)**2+a(2)**2)     &
              *x(3)**2) 
    coef(3,0)=a(1)**2*a(2)**2*a(3)**2-(a(2)**2*a(3)**2*x(1)**2+a(1)**2*a(3)**2*&
              x(2)**2+a(1)**2*a(2)**2*x(3)**2) 
    call Root_4(3,coef,root,root_im2,root_im3)
    lambda=f0
    if (x(1)**2/a(1)**2+x(2)**2/a(2)**2+x(3)**2/a(3)**2>f1) then
       if (root_im2==f0 .and. root_im3==f0) then 
          lambda=max(f0,maxval(root))
       else 
          lambda=max(f0,root(3))
       end if 
    end if
    theta=asin(sqrt((a(1)**2-a(3)**2)/(a(1)**2+lambda))) ! the amplitude
    m=(a(1)**2-a(2)**2)/(a(1)**2-a(3)**2) ! m=k**2 is the parameter
    call elbd(theta,0.5*pi-theta,f1-m,B,D)
    F=B+D; E=B+(f1-m)*D
    ! Calculation of Is 
    if (a(1)==a(2) .and. a(1)==a(3)) then
       Ifir=(4/3)*pi*a(1)**3/(a(1)**2+lambda)**(3/2)
       Isec=(4/5)*pi*a(1)**3/(a(1)**2+lambda)**(1/2)
    elseif (a(1)>a(2) .and. a(3)==a(2)) then
       del=sqrt((a(1)**2+lambda)*(a(2)**2+lambda)*(a(3)**2+lambda))
       bbar=sqrt(a(1)**2+lambda)/sqrt(a(3)**2+lambda)
       dbar=sqrt(a(1)**2-a(3)**2)/sqrt(a(3)**2+lambda)
       Ifir(1)=4*pi*a(1)*a(2)**2*(acosh(bbar)-dbar/bbar)/(a(1)**2-a(2)**2)**1.5
       Ifir(2)=2*pi*a(1)*a(2)**2*(-acosh(bbar)+dbar*bbar)/(a(1)**2-a(2)**2)**1.5
       Ifir(3)=Ifir(2)
       Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)**2-a(2)**2)
       Isec(1,3)=Isec(1,2)
       Isec(2,1)=Isec(1,2)
       Isec(3,1)=Isec(1,3)
       Isec(2,3)=pi*product(a)/((a(3)**2+lambda)*del)-Isec(1,3)*0.25
       Isec(3,2)=Isec(2,3)
       Isec(1,1)=((f4*pi*product(a))/((a(1)**2+lambda)*del)-Isec(1,2)-         &
                 Isec(1,3))/f3
       Isec(2,2)=Isec(2,3)
       Isec(3,3)=Isec(2,3)
    elseif (a(1)==a(2) .and. a(2)>a(3)) then
       del=sqrt((a(1)**2+lambda)*(a(2)**2+lambda)*(a(3)**2+lambda))
       bbar=sqrt(a(3)**2+lambda)/sqrt(a(1)**2+lambda)
       dbar=sqrt(a(1)**2-a(3)**2)/sqrt(a(1)**2+lambda)
       Ifir(1)=2*pi*a(1)**2*a(3)*(acos(bbar)-dbar*bbar)/(a(1)**2-a(3)**2)**1.5
       Ifir(2)=Ifir(1)
       Ifir(3)=4*pi*product(a)/del-2*Ifir(1)
       Isec(1,3)=(Ifir(3)-Ifir(1))/(a(1)**2-a(3)**2)
       Isec(3,1)=Isec(1,3)
       Isec(2,3)=Isec(1,3)
       Isec(3,2)=Isec(2,3)
       Isec(1,1)=pi*product(a)/((a(1)**2+lambda)*del)-Isec(1,3)*0.25
       Isec(1,2)=Isec(1,1)
       Isec(2,1)=Isec(1,2)
       Isec(2,2)=Isec(1,1)
       Isec(3,3)=((f4*pi*product(a))/((a(3)**2+lambda)*del)-Isec(1,3)-         &
                 Isec(2,3))/f3
    else
       del=sqrt((a(1)**2+lambda)*(a(2)**2+lambda)*(a(3)**2+lambda))
       Ifir(1)=f4*pi*product(a)*F/sqrt(a(1)**2-a(3)**2)*(f1-E/F)/(a(1)**2-     &
               a(2)**2)
       Ifir(2)=f4*pi*product(a)*(E*sqrt(a(1)**2-a(3)**2)/((a(1)**2-a(2)**2)*   &
               (a(2)**2-a(3)**2))-F/((a(1)**2-a(2)**2)*sqrt(a(1)**2-a(3)**2))- &
               (f1/(a(2)**2-a(3)**2))*sqrt((a(3)**2+lambda)/((a(1)**2+lambda)* &
               (a(2)**2+lambda))))
       Ifir(3)=f4*pi*product(a)/del-Ifir(1)-Ifir(2)
       Isec(1,2)=(Ifir(2)-Ifir(1))/(a(1)**2-a(2)**2)
       Isec(2,1)=Isec(1,2)
       Isec(1,3)=(Ifir(3)-Ifir(1))/(a(1)**2-a(3)**2)
       Isec(3,1)=Isec(1,3)
       Isec(2,3)=(Ifir(3)-Ifir(2))/(a(2)**2-a(3)**2)
       Isec(3,2)=Isec(2,3)
       Isec(1,1)=((4*pi*product(a))/((a(1)**2+lambda)*del)-Isec(1,2)-          &
                 Isec(1,3))/f3
       Isec(2,2)=((4*pi*product(a))/((a(2)**2+lambda)*del)-Isec(1,2)-          &
                 Isec(2,3))/f3
       Isec(3,3)=((4*pi*product(a))/((a(3)**2+lambda)*del)-Isec(1,3)-          &
                 Isec(2,3))/f3
    end if

    ! I derivatives
    call buildtensors(a,a_21,a_22) 
    ultadelfir=-f2*pi*product(a)/((a**2+lambda)*del)
    call buildtensors(ultadelfir,ultadelfir_21,ultadelfir_22) 
    ultadelsec=-f2*pi*product(a)/((a_21**2+lambda)*(a_22**2+lambda)*del)
    ! derivatives of lambda 
    c1=sum((x**2)/((a**2+lambda)**2))
    c2=sum((x**2)/((a**2+lambda)**3))
    c3=sum((x**2)/((a**2+lambda)**4))
    Fvec=f2*x/(a**2+lambda)
    call buildtensors(Fvec,F_21,F_22)
    if (lambda==f0) then 
       fderlambda=f0 
    else
       fderlambda=Fvec/c1
    end if 
    call buildtensors(fderlambda,fderlambda_21,fderlambda_22)
    diagvals=f0;nondiagvals=f0
    do i=1,3
       do j=1,3
          if (i==j) then
             diagvals(i,j)=f1
          else
             nondiagvals(i,j)=f1
          end if
       end do
    end do
    fderF=nondiagvals*(f1/(a_21**2+lambda))*(-F_21*fderlambda_22)+diagvals*    &
          (f1/(a_21**2+lambda))*(f2-F_21*fderlambda_22)
    fderc1=Fvec/(a**2+lambda)-f2*c2*fderlambda
    call buildtensors(fderc1,fderc1_21,fderc1_22) 
    fderc2=Fvec/(a**2+lambda)**2-f3*c3*fderlambda
    call buildtensors(fderc2,fderc2_21,fderc2_22) 
    if (lambda==f0) then 
       sderlambda=f0 
    else
       sderlambda=(fderF-fderlambda_21*fderc1_22)/c1
    end if 
    sderc1=(f1/(a_21**2+lambda))*(fderF-fderlambda_22*F_21/(a_21**2+lambda))-  &
       f2*(fderc2_22*fderlambda_21+c2*sderlambda)
    fderIfir=ultadelfir_21*fderlambda_22
    do q=1,3 
       do p=1,3 
          do r=1,3 
             sderF(q,p,r)=-(fderF(q,p)*fderlambda(r)+fderF(q,r)*fderlambda(p)+ &
                          Fvec(q)*sderlambda(p,r))/(a(q)**2+lambda)
          end do 
       end do 
    end do
    zeefir=f1/(a**2+lambda)+0.5*sum(f1/(a**2+lambda))
    zeesec=f1/(a_21**2+lambda)+f1/(a_22**2+lambda)+0.5*sum(f1/(a**2+lambda))
    do i=1,3 
       do j=1,3 
          do k=1,3 
             sderIfir(i,j,k)=ultadelfir(i)*(sderlambda(j,k)-fderlambda(j)*     &
                             fderlambda(k)*zeefir(i))
          end do 
       end do 
    end do 
    do i=1,3 
       do j=1,3 
          do k=1,3 
             fderIsec(i,j,k)=ultadelsec(i,j)*fderlambda(k)
          end do 
       end do 
    end do 
    do i=1,3 
       do j=1,3 
          do k=1,3 
             do l=1,3 
                sderIsec(i,j,k,l)=ultadelsec(i,j)*(sderlambda(k,l)-            &
                                  fderlambda(k)*fderlambda(l)*zeesec(i,j))
             end do 
          end do 
       end do 
    end do 
    do q=1,3 
       do p=1,3 
          do r=1,3 
             if (lambda==f0) then 
                tderlambda(q,p,r)=f0
             else
                tderlambda(q,p,r)=(-f1/c1)*(sderlambda(q,p)*fderc1(r)-         &
                                  sderF(q,p,r)+sderlambda(q,r)*fderc1(p)+      &
                                  fderlambda(q)*sderc1(p,r))
             end if 
          end do 
       end do 
    end do 
    
    !Calculation of V-potentials 
    do i=1,3 
       do p=1,3 
          do q=1,3 
             call kdelta(p,q,delta1) 
             sderVfir(i,p,q)=-(delta1*Isec(p,i)+x(p)*fderIsec(p,i,q))
          end do 
       end do 
    end do 
    do i=1,3 
       do p=1,3 
          do q=1,3 
             do r=1,3 
                call kdelta(p,q,delta1); call kdelta(p,r,delta2)
                tderVfir(i,p,q,r)=-(delta1*fderIsec(p,i,r)+delta2*             &
                                  fderIsec(p,i,q)+x(p)*sderIsec(p,i,q,r))
             end do 
          end do 
       end do 
    end do 

    !calculation of phi derivatives 
    do p=1,3 
       do q=1,3 
          call kdelta(p,q,delta1) 
          sderphi(p,q)=-(delta1*Ifir(p)+x(p)*fderIfir(p,q))
       end do 
    end do 
    do p=1,3 
       do q=1,3 
          do r=1,3 
             call kdelta(p,q,delta1); call kdelta(p,r,delta2)
             tderphi(p,q,r)=-(delta1*fderIfir(p,r)+delta2*fderIfir(p,q)+x(p)*  &
                            sderIfir(p,q,r))
          end do 
       end do 
    end do 

    !psi's 
    do i=1,3 
       do j=1,3 
          do k=1,3 
             do l=1,3 
                call kdelta(i,j,delta1); call kdelta(i,k,delta2)
                call kdelta(i,l,delta3)
                foderpsi(i,j,k,l)=delta1*(sderphi(k,l)-a(i)**2*                &
                                  sderVfir(i,k,l))+delta2*(sderphi(j,l)-       &
                                  a(i)**2*sderVfir(i,j,l))+delta3*             &
                                  (sderphi(j,k)-a(i)**2*sderVfir(i,j,k))+      &
                                  x(i)*(tderphi(j,k,l)-a(i)**2*                &
                                  tderVfir(i,j,k,l))
             end do 
          end do 
       end do 
    end do 

    !calculation of D4 
    premult1=f1/(f8*pi*(f1-vm))
    do i=1,3 
       do j=1,3 
          do k=1,3 
             do l=1,3 
                call kdelta(k,l,delta1); call kdelta(i,l,delta2)
                call kdelta(j,l,delta3); call kdelta(i,k,delta4)
                call kdelta(j,k,delta5)
                D4(i,j,k,l)=premult1*(foderpsi(k,l,i,j)-f2*vm*delta1*          &
                            sderphi(i,j)-(f1-vm)*(sderphi(k,j)*delta2+         &
                            sderphi(k,i)*delta3+sderphi(l,j)*delta4+           &
                            sderphi(l,i)*delta5))
             end do 
          end do 
       end do 
    end do

    ! fderphi tderpsi
    fderphi=-x*Ifir
    do i=1,3
       do j=1,3
          do l=1,3  
             call kdelta(i,j,delta1); call kdelta(i,l,delta2)
             call kdelta(j,l,delta3)
             tderpsi(i,j,l)=-delta1*x(l)*(Ifir(l)-a(i)**2*Isec(i,l))-          &  
                            x(i)*x(j)*(fderIfir(j,l)-a(i)**2*fderIsec(i,j,l))- &
                            (delta2*x(j)+delta3*x(i))*(Ifir(j)-                &
                            a(i)**2*Isec(i,j))
          end do
       end do
    end do
  end subroutine EshD4   

  ! Form RVE Dstr(6,6) at hex Gauss point, 8 nodes element only
  subroutine Hex8Dstr(ellip,Em,vm,Dstr)
    ! ellip(nellip,11), nellip inclusions in one element, center, semi-axes, 
    ! rotation angels, Eh, vh;  
    ! Dstr(nip,6,6), Strain=>elastic strain at Gauss points;
    ! Dstr(i,:,:)=[I+sum_e(Dei(Cm-dCeSe)^-1dCe)]^-1, i inp, e ellip 
    ! S2(6,6), Rank 2 first Eshelby tensors of each inclusion;
    ! D2(6,6), Rank 2 second Eshelby tensors of each inclusion and gauss point;
    ! D2e(6,6), elastic stran=>eigen strain (Cm-dCeSe)^-1dCe of each inclusion
    ! Cm(6,6), Rank 2 host stiffness matrix;
    ! Ch(6,6), Rank 2 stiffness matrix of each inclusion;
    ! dC(6,6), Cm-Ch;
    implicit none  
    integer :: i,j,k,l,nellip 
    real(8) :: ellip(:,:),Dstr(8,6,6),D4(3,3,3,3),a(3),Em,vm,Eh,vh,ang(3),     &
       R_init(3,3),Rb_init(3,3),R(3,3),Rb(3,3),xint(8,3),xobs(3,1),PIvec(3),   &
       fderphi(3),tderpsi(3,3,3),S2(6,6),D2(6,6),R2(3,3),R2b(3,3),Cm(6,6),     &
       Ch(6,6),D2e(6,6),D2R(6,6),dC(6,6),DInv(8,6,6),exh(3,3),tmp
    nellip=size(ellip,1)
    DInv=f0
    do i=1,6
       DInv(:,i,i)=f1
    end do
    ! Gauss points of 8-node element
    xint(1,:)=(/-sqrt(f1/f3),-sqrt(f1/f3),-sqrt(f1/f3)/)
    xint(2,:)=(/ sqrt(f1/f3),-sqrt(f1/f3),-sqrt(f1/f3)/)
    xint(3,:)=(/ sqrt(f1/f3), sqrt(f1/f3),-sqrt(f1/f3)/)
    xint(4,:)=(/-sqrt(f1/f3), sqrt(f1/f3),-sqrt(f1/f3)/)
    xint(5,:)=(/-sqrt(f1/f3),-sqrt(f1/f3), sqrt(f1/f3)/)
    xint(6,:)=(/ sqrt(f1/f3),-sqrt(f1/f3), sqrt(f1/f3)/)
    xint(7,:)=(/ sqrt(f1/f3), sqrt(f1/f3), sqrt(f1/f3)/)
    xint(8,:)=(/-sqrt(f1/f3), sqrt(f1/f3), sqrt(f1/f3)/)
    call CMat(Em,vm,Cm); 
    do i=1,nellip
       a=ellip(i,4:6)
       ! Stage a1>=a2>=a3
       exh=f0
       do k=1,2
          do l=2,3
             if (a(k)<a(l)) then
                exh(k,l)=f1
                tmp=a(k)
                a(k)=a(l)
                a(l)=tmp
             end if 
          end do
       end do
       ! Initial rotation matrices due to axis exchange
       ang=pi/f2*(/exh(2,3),exh(1,3),exh(1,2)/)
       call Ang2Mat(ang,R_init,f1)
       call Ang2Mat(ang,Rb_init,-f1)
       ! Rotation matrices w.r.t the ellipsoid
       ang=ellip(i,7:9)
       call Ang2Mat(ang,R,f1)
       call Ang2Mat(ang,Rb,-f1) 
       R2=matmul(R_init,Rb)  ! Gauss=>ellipsoid 
       R2b=matmul(R,Rb_init) ! Ellipsoid=>Gauss
       call EshS2(vm,a,S2,PIvec)
       Eh=ellip(i,10); vh=ellip(i,11)
       call CMat(Eh,vh,Ch); dC=Cm-Ch
       call MatInv(Cm-matmul(dC,S2),D2e)
       D2e=matmul(D2e,dC)
       do j=1,8
          xobs(:,1)=xint(j,:)-ellip(i,:3) ! Relative coordinate
          xobs=matmul(R2,xobs) ! Element=>ellipsoid
          if (xobs(1,1)**2/a(1)**2+xobs(2,1)**2/a(2)**2+xobs(3,1)**2/a(3)**2<= &
             f1) then
             print*, "Warining, Gauss point inside inclusion!"
             D2=S2
             go to 250
          end if  
          call EshD4(vm,a,xobs(:,1),D4,fderphi,tderpsi)
          call T4T2(D4,D2)
250       D2=matmul(D2,D2e)
          call T2Rot(D2,R2b,D2R) ! Rotate to Cartesian before add
          DInv(j,:,:)=DInv(j,:,:)+D2R 
       end do
    end do
    do j=1,8
       call MatInv(DInv(j,:,:),Dstr(j,:,:))
    end do
  end subroutine Hex8Dstr

  ! Wrapper of Hex8Dstr for all RVEs
  subroutine RVEHex8Dstr(mat,ellip,matDstr)
    ! ellip(nellip,11), total ellipsoid shapes  
    ! i1, i2, start and end of the sampling for i-th material
    ! matDstr(i,j,6,6), Dstr for i-th material at j-th gauss point
    implicit none
    integer :: nmt,ncol,i1,i2,i
    real(8) :: mat(:,:),ellip(:,:),matDstr(:,:,:,:),Dstr(8,6,6),I8(8,6,6) 
    nmt=size(mat,1)
    ncol=size(mat,2)
    I8=f0
    do i=1,6
       I8(:,i,i)=f1 
    end do
    do i=1,nmt
       i1=int(mat(i,ncol-1)); i2=int(mat(i,ncol))
       if (i1>0 .and. i2>=i1) then
          call Hex8Dstr(ellip(i1:i2,:),mat(i,1),mat(i,2),Dstr)
          matDstr(i,:,:,:)=Dstr
       else
          matDstr(i,:,:,:)=I8
       end if
    end do
  end subroutine RVEHex8Dstr 

  subroutine EshDisp(vm,eigen,fderphi,tderpsi,u)
    implicit none 
    integer :: i,j,k
    real(8) :: vm,u(3),eigen(6),fderphi(3),tderpsi(3,3,3),MatEigen(3,3),       &
        ut(3,1),SumDiag,fderphit(3,1),premult
    call Vec2Mat(eigen,MatEigen)
    SumDiag=f0; ut=f0
    do i=1,3
       do j=1,3
          do k=1,3
             ut(i,1)=ut(i,1)+tderpsi(i,j,k)*MatEigen(j,k);
          end do
       end do
       fderphit(i,1)=fderphi(i)
       SumDiag=SumDiag+MatEigen(i,i)
    end do
    premult=f1/(f8*pi*(f1-vm))
    ut=premult*(ut-f2*vm*SumDiag*fderphit-f4*(f1-vm)*matmul(MatEigen,fderphit))
    u=ut(:,1)
  end subroutine EshDisp

  subroutine EshSol(Em,vm,stress,ellip,ocoord,sol)
    implicit none
    integer :: i,j,k,l,m,n,nobs,nellip
    real(8) :: Em,vm,Eh,vh,ocoord(:,:),ellip(:,:),stress(6),sol(:,:) !In(out)put
    ! ellip(nellip,17): 1-3 ellipsoid centroid coordinate, 4-6 semi-axises, 7-9
    ! rotation angles around x,y and z axises, 10,11 inclusion Young's modulus 
    ! and Poisson's ratio, 12-17 eigen strain 
    ! stress(6): remote stress
    ! sol(nobs,9): 1-3 displacement, 4-9 stress 
    real(8) :: ang(3),a(3),tmp,exh(3,3),R_init(3,3),Rb_init(3,3),R(3,3),       &
       Rb(3,3),PIvec(3),Tstress(3,3),Cm(6,6),straint(6,1),stresst(6,1),        &
       eigent(6,1),Ch(6,6),dC(6,6),vert(3,1),D4(3,3,3,3),fderphi(3),           &
       tderpsi(3,3,3),disp(3),dispt(3,1),Ttmp(3,3),Vtmp(6,1),Teigen(3,3),      &
       S2(6,6)
    nobs=size(ocoord,1); nellip=size(ellip,1)
    !sol=f0 ! Initial solution space
    do i=1,nellip
       a=ellip(i,4:6)
       ! Stage a1>=a2>=a3
       exh=f0
       do k=1,2
          do l=2,3
             if (a(k)<a(l)) then
                exh(k,l)=f1
                tmp=a(k)
                a(k)=a(l)
                a(l)=tmp
             end if 
          end do
       end do
       ! Initial rotation matrices due to axis exchange
       ang=pi/f2*(/exh(2,3),exh(1,3),exh(1,2)/)
       call Ang2Mat(ang,R_init,f1)
       call Ang2Mat(ang,Rb_init,-f1)
       ! Rotation matrices w.r.t the ellipsoid
       ang=ellip(i,7:9)
       call Ang2Mat(ang,R,f1)
       call Ang2Mat(ang,Rb,-f1) 
       ! Eshelby's tensor
       call EshS2(vm,a,S2,PIvec) 
       ! Rotate stress and initial eigenstrain against oblique ellipsoid
       call Vec2Mat(stress,Tstress) 
       Tstress=matmul(matmul(matmul(R_init,Rb),Tstress),                       &
               transpose(matmul(R_init,Rb)))
       stresst(:,1)=(/Tstress(1,1),Tstress(2,2),Tstress(3,3),Tstress(1,2),     &
                    Tstress(2,3),Tstress(1,3)/)    
       call Vec2Mat(ellip(i,12:17),Teigen)    
       !Teigen=matmul(matmul(matmul(R_init,Rb),Teigen),                         &
       !       transpose(matmul(R_init,Rb)))
       !print*, matmul(S2,eigent)
       Teigen=matmul(matmul(R_init,Teigen),transpose(R_init))
       eigent(:,1)=(/Teigen(1,1),Teigen(2,2),Teigen(3,3),Teigen(1,2),          &
                   Teigen(2,3),Teigen(1,3)/) 
       Eh=ellip(i,10); vh=ellip(i,11)
       call CMat(Em,vm,Cm); call CMat(Eh,vh,Ch); dC=Cm-Ch
       call SolveSix(Cm,stresst,straint)
       call SolveSix(Cm-matmul(dC,S2),matmul(dC,straint)+matmul(Ch,eigent),    &
          eigent)
       call Vec2Mat(eigent(:,1),Teigen)
       do j=1,nobs
          vert(:,1)=ocoord(j,:)-ellip(i,:3) ! Relative coordinate
          vert=matmul(matmul(R_init,Rb),vert)
          call EshD4(vm,a,vert(:,1),D4,fderphi,tderpsi) 
          call EshDisp(vm,eigent(:,1),fderphi,tderpsi,disp)
          dispt(:,1)=disp
          ! Rotate back
          dispt=matmul(matmul(R,Rb_init),dispt)
          !dispt=matmul(transpose(matmul(R_init,Rb)),dispt)
          ! Record displacement
          sol(j,:3)=sol(j,:3)+dispt(:,1)
          if (vert(1,1)**2/a(1)**2+vert(2,1)**2/a(2)**2+vert(3,1)**2/a(3)**2   &
             <=1) then ! J-th obs interior to i-th inclusion 
             ! Total elastic stress
             stresst=stresst+matmul(Cm,(matmul(S2,eigent)-eigent))
          else ! J-th obs exterior to i-th inclusion
             Ttmp=f0 
             do k=1,3  
                do l=1,3 
                   do m=1,3 
                      do n=1,3 
                         Ttmp(k,l)=Ttmp(k,l)+D4(k,l,m,n)*Teigen(m,n) 
                      end do
                   end do
                end do
             end do
             Vtmp(:,1)=(/Ttmp(1,1),Ttmp(2,2),Ttmp(3,3),Ttmp(1,2),Ttmp(2,3),    &
                       Ttmp(1,3)/)
             ! Total elastic stress
             stresst=stresst+matmul(Cm,Vtmp)
          end if
          ! Rotate back to origin coordinate
          call Vec2Mat(stresst(:,1),Tstress)
          Tstress=matmul(matmul(matmul(R,Rb_init),Tstress),                    &
                  transpose(matmul(R,Rb_init)))
          stresst(:,1)=(/Tstress(1,1),Tstress(2,2),Tstress(3,3),               &
                       Tstress(1,2),Tstress(2,3),Tstress(1,3)/) 
          ! Record stress         
          sol(j,4:9)=sol(j,4:9)+stresst(:,1)
          stresst(:,1)=stress ! Reset to background stress 
       end do ! nobs
    end do ! nellip
  end subroutine EshSol    

  ! Translation 3x3x3x3 tensor to 6x6, FE order 1-11,2-22,3-33,4-12,5-23,6-13
  subroutine T4T2(T4,T2)
    implicit none
    real(8) :: T4(3,3,3,3),T2(6,6)
    T2(1,:)=(/T4(1,1,1,1),T4(1,1,2,2),T4(1,1,3,3),T4(1,1,1,2)+T4(1,1,2,1),     &
              T4(1,1,2,3)+T4(1,1,3,2),T4(1,1,1,3)+T4(1,1,3,1)/)
    T2(2,:)=(/T4(2,2,1,1),T4(2,2,2,2),T4(2,2,3,3),T4(2,2,1,2)+T4(2,2,2,1),     &
              T4(2,2,2,3)+T4(2,2,3,2),T4(2,2,1,3)+T4(2,2,3,1)/)
    T2(3,:)=(/T4(3,3,1,1),T4(3,3,2,2),T4(3,3,3,3),T4(3,3,1,2)+T4(3,3,2,1),     &
              T4(3,3,2,3)+T4(3,3,3,2),T4(3,3,1,3)+T4(3,3,3,1)/)
    T2(4,:)=(/T4(1,2,1,1),T4(1,2,2,2),T4(1,2,3,3),T4(1,2,1,2)+T4(1,2,2,1),     &
              T4(1,2,2,3)+T4(1,2,3,2),T4(1,2,1,3)+T4(1,2,3,1)/)
    T2(5,:)=(/T4(2,3,1,1),T4(2,3,2,2),T4(2,3,3,3),T4(2,3,1,2)+T4(2,3,2,1),     &
              T4(2,3,2,3)+T4(2,3,3,2),T4(2,3,1,3)+T4(2,3,3,1)/)
    T2(6,:)=(/T4(1,3,1,1),T4(1,3,2,2),T4(1,3,3,3),T4(1,3,1,2)+T4(1,3,2,1),     & 
              T4(1,3,2,3)+T4(1,3,3,2),T4(1,3,1,3)+T4(1,3,3,1)/)   
  end subroutine T4T2    

  ! Rotate T2(6,6) tensor by R(3,3), FE order  
  subroutine T2Rot(T2,R,T2R)
    implicit none  
    real(8) :: T2(6,6),T2R(6,6),R(3,3),matR(6,6),matRIv(6,6)
    matR(:,1)=(/R(1,1)**2,R(1,2)**2,R(1,3)**2,R(1,1)*R(1,2),R(1,2)*R(1,3),     &
                R(1,1)*R(1,3)/)
    matR(:,2)=(/R(2,1)**2,R(2,2)**2,R(2,3)**2,R(2,1)*R(2,2),R(2,2)*R(2,3),     &
                R(2,1)*R(2,3)/)
    matR(:,3)=(/R(3,1)**2,R(3,2)**2,R(3,3)**2,R(3,1)*R(3,2),R(3,2)*R(3,3),     &
                R(3,1)*R(3,3)/)
    matR(:,4)=(/f2*R(1,1)*R(2,1),f2*R(1,2)*R(2,2),f2*R(1,3)*R(2,3),            & 
                R(1,1)*R(2,2)+R(1,2)*R(2,1),R(1,2)*R(2,3)+R(1,3)*R(2,2),       &
                R(1,1)*R(2,3)+R(1,3)*R(2,1)/) 
    matR(:,5)=(/f2*R(2,1)*R(3,1),f2*R(2,2)*R(3,2),f2*R(2,3)*R(3,3),            & 
                R(2,1)*R(3,2)+R(2,2)*R(3,1),R(2,2)*R(3,3)+R(2,3)*R(3,2),       &
                R(2,1)*R(3,3)+R(2,3)*R(3,1)/)
    matR(:,6)=(/f2*R(1,1)*R(3,1),f2*R(1,2)*R(3,2),f2*R(1,3)*R(3,3),            & 
                R(1,1)*R(3,2)+R(1,2)*R(3,1),R(1,2)*R(3,3)+R(1,3)*R(3,2),       &
                R(1,1)*R(3,3)+R(1,3)*R(3,1)/)
    call MatInv(matR,matRIv)
    T2R=matmul(matmul(matRIv,T2),matR)
  end subroutine T2Rot    

  subroutine OkSol(E,v,rects,ocoord,crest,sol) 
    ! Okada wrapper for multiple rectangle faults
    ! E, nu: elastic moduli;
    ! rects: (nrect,9) 1~3 fault center (x,y,z); 4,5 fault length in strike    
    ! and dip; 6 dip angle (rad); 7~9 strike (left lateral+) dip
    ! (up+) and tensile (open) dislocations;
    ! ocoord: (nobs_loc,3) observation locations;
    ! crest: topography height;
    ! sol: (nobs_loc,9) ux ... yz, sigma_xx .. sigma_xz 
    implicit none
    integer :: nobs,nrect,i,j,iret
    real(8) :: E,v,nu,lbd,alpha,rects(:,:),ocoord(:,:),crest,sol(:,:),x,y,z,   &
       deep, ux,uy,uz,exx,eyy,ezz,exy,eyz,exz,eyx,ezy,ezx,strain(6,1),         &
       stress(6,1),Cm(6,6),vectmp(3,1),R(3,3),tenstmp(3,3),ang(3)
    call CMat(E,v,Cm)
    nu=E/f2/(1+v)
    lbd=E*v/(f1+v)/(f1-f2*v)
    alpha=(lbd+nu)/(lbd+f2*nu)
    nrect=size(rects,1); nobs=size(ocoord,1)
    ang=(/f0,f0,pi/f2/)
    call Ang2Mat(ang,R,f1)
    do i=1,nrect
       do j=1,nobs 
          x=ocoord(j,1)-rects(i,1) 
          y=ocoord(j,2)-rects(i,2)
          z=ocoord(j,3)-crest         ! Adjust height against crest 
          deep=-(rects(i,3)-crest)    ! deep positive
          ! Translate to Okada coordinate
          vectmp(:,1)=(/x,y,z/) 
          vectmp=matmul(R,vectmp) 
          x=vectmp(1,1); y=vectmp(2,1); z=vectmp(3,1)
          call dc3d(alpha,x,y,z,deep,rects(i,6)*1.8D2/pi,rects(i,4)/f2,        &
             rects(i,4)/f2,rects(i,5)/f2,rects(i,5)/f2,rects(i,7),rects(i,8),  &
             rects(i,9),ux,uy,uz,exx,eyx,ezx,exy,eyy,ezy,exz,eyz,ezz,iret)
          ! Translate back to Eshelby coordinate
          vectmp(:,1)=(/ux,uy,uz/) 
          vectmp=matmul(transpose(R),vectmp)  
          ux=vectmp(1,1); uy=vectmp(2,1); uz=vectmp(3,1)
          tenstmp=reshape((/exx,eyx,ezx,exy,eyy,ezy,exz,eyz,ezz/),(/3,3/),     &
                          (/f0,f0/),(/2,1/))  
          tenstmp=matmul(matmul(transpose(R),tenstmp),R)
          strain(:,1)=(/tenstmp(1,1),tenstmp(2,2),tenstmp(3,3),(tenstmp(1,2)+  &
                        tenstmp(2,1))/f2,(tenstmp(2,3)+tenstmp(3,2))/f2,       &
                        (tenstmp(1,3)+tenstmp(3,1))/f2/)
          stress=matmul(Cm,strain)
          sol(j,:3)=sol(j,:3)+(/ux,uy,uz/)
          sol(j,4:9)=sol(j,4:9)+stress(:,1)
       end do
    end do
  end subroutine OkSol  

end module esh3d
