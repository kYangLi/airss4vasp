!==================================================================================!
!                                niggli_reduce                                     !
!==================================================================================!
!                                                                                  !
! This file is part of the AIRSS structure prediction package.                     !
!                                                                                  !
! AIRSS is free software; you can redistribute it and/or modify it under the terms !
! of the GNU General Public License version 2 as published by the Free Software    !
! Foundation.                                                                      !
!                                                                                  !
! This program is distributed in the hope that it will be useful, but WITHOUT ANY  !
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A  !
! PARTICULAR PURPOSE.  See the GNU General Public License for more details.        !           
!                                                                                  !
! You should have received a copy of the GNU General Public License along with this!
! program; if not, write to the Free Software Foundation, Inc., 51 Franklin Street,!                   
! Fifth Floor, Boston, MA  02110-1301, USA.                                        !
!                                                                                  !
!----------------------------------------------------------------------------------!
! This subroutine performs a niggli reduction                                      !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

subroutine niggli_reduce(Lm,Cm)
  
  use constants

  implicit none

  real(kind=dp), dimension(3,3), intent(inout) :: Lm
  integer,       dimension(3,3), intent(out)   :: Cm

  ! -----

  integer, parameter :: maxit=1000
  real(kind=dp),parameter :: eps=1e-10_dp

  integer       :: Cc(3,3)

  real(kind=dp) :: A,B,C,xi,eta,zeta

  integer :: it
  integer :: l,m,n
  integer :: i,j,k
  
  logical :: change

  ! ** Following Atsushi Togoto's implementation of
  ! Numerically stable algorithms for the computation of reduced unit cells,
  ! R. W. Grosse-Kunstleve, N. K. Sauter and P. D. Adams, Acta Cryst., A60, 1-6 (2004)
  ! http://atztogo.github.io/niggli/#algorithm 05/01/16

  Cm=nint(ident)
  Cc=nint(ident)
  change=step0()
  !write (*,*) nint(A),nint(B),nint(C),nint(xi),nint(eta),nint(zeta)

  do it=1,maxit

     change=step1()
     
     change=step2()
     if(change) cycle

     change=step3()

     change=step4()

     change=step5()
     if(change) cycle

     change=step6()
     if(change) cycle

     change=step7()
     if(change) cycle

     change=step8()
     if(change) cycle

     exit
     
  end do

contains

  subroutine reset()

    Lm=matmul(Lm,transpose(Cc))
    Cm=matmul(Cm,transpose(Cc))
    change=step0()
    !write (*,*) nint(A),nint(B),nint(C),nint(xi),nint(eta),nint(zeta)

  end subroutine reset
  
  function step0()

    logical :: step0

    A=   dot_product(Lm(:,1),Lm(:,1))
    B=   dot_product(Lm(:,2),Lm(:,2))
    C=   dot_product(Lm(:,3),Lm(:,3))
    xi=  dot_product(Lm(:,2),Lm(:,3))*2.0_dp
    eta= dot_product(Lm(:,3),Lm(:,1))*2.0_dp
    zeta=dot_product(Lm(:,1),Lm(:,2))*2.0_dp

    l=0 ; m=0 ; n=0 

    if(xi.lt.-eps) l=-1
    if(xi.gt.eps) l=1
    if(eta.lt.-eps) m=-1
    if(eta.gt.eps) m=1
    if(zeta.lt.-eps) n=-1
    if(zeta.gt.eps) n=1
    
    step0=.true.

  end function step0

  function step1()

    logical :: step1

    if((A.gt.B+eps).or.((.not.(abs(A-B).gt.eps)).and.(abs(xi).gt.abs(eta)+eps))) then
       Cc=0
       Cc(2,1)=-1
       Cc(1,2)=-1
       Cc(3,3)=-1
       call reset()
       step1=.true.
    else
       step1=.false.
    end if

  end function step1

  function step2()

    logical :: step2

    if((B.gt.C+eps).or.((.not.(abs(B-C).gt.eps)).and.(abs(eta).gt.abs(zeta)+eps))) then
       Cc=0
       Cc(1,1)=-1
       Cc(3,2)=-1
       Cc(2,3)=-1
       call reset()
       step2=.true.
    else
       step2=.false.
    end if

  end function step2

  function step3()

    logical :: step3

    if(l*m*n.eq.1) then

       if(l.eq.-1) then
          i=-1
       else
          i=1
       end if
       if(m.eq.-1) then
          j=-1
       else
          j=1
       end if
       if(n.eq.-1) then
          k=-1
       else
          k=1
       end if

       Cc=0
       Cc(1,1)=i
       Cc(2,2)=j
       Cc(3,3)=k

       call reset()
       step3=.true.
    else
       step3=.false.
    end if

  end function step3

  function step4()

    logical :: step4
    integer :: r

    if((l.eq.-1).and.(m.eq.-1).and.(n.eq.-1)) then
       step4=.false.
       return
    end if
    
    
    ! ** this needs to be updated to the new step4, r etc
    
    if((l*m*n.eq.0).or.(l*m*n.eq.-1)) then

       i=1 ; j=1 ; k=1
       
       if(l.eq.1) i=-1
       if(l.eq.0) r=1
       if(m.eq.1) j=-1
       if(m.eq.0) r=2
       if(n.eq.1) k=-1
       if(n.eq.0) r=3
        
       if(i*j*k.eq.-1) then
          if(r.eq.1) i=-1
          if(r.eq.2) j=-1
          if(r.eq.3) k=-1
       end if

       Cc=0
       Cc(1,1)=i
       Cc(2,2)=j
       Cc(3,3)=k

       call reset()
       step4=.true.
    else
       step4=.false.
    end if

  end function step4
  
  function step5()

    logical :: step5

    if ((abs(xi).gt.B+eps).or.(.not.(abs(B-xi).gt.eps).and.(2.0_dp*eta.lt.zeta-eps))&
         .or.(.not.(abs(B+xi).gt.eps).and.(zeta.lt.-eps))) then

       Cc=0
       Cc(1,1)=1
       Cc(2,2)=1
       Cc(3,3)=1
       if(xi.gt.0.0_dp) Cc(3,2)=-1  
       if(xi.lt.0.0_dp) Cc(3,2)=1  

       call reset()
       step5=.true.
    else
       step5=.false.
    end if

  end function step5

  function step6()

    logical :: step6

    if ((abs(eta).gt.A+eps).or.(.not.(abs(A-eta).gt.eps).and.(2.0_dp*xi.lt.zeta-eps))&
         .or.(.not.(abs(A+eta).gt.eps).and.(zeta.lt.-eps))) then

       Cc=0
       Cc(1,1)=1
       Cc(2,2)=1
       Cc(3,3)=1
       if(eta.gt.0.0_dp) Cc(3,1)=-1  
       if(eta.lt.0.0_dp) Cc(3,1)=1
       
       call reset()
       step6=.true.
    else
       step6=.false.
    end if

  end function step6

  function step7()

    logical :: step7

    if ((abs(zeta).gt.A+eps).or.(.not.(abs(A-zeta).gt.eps).and.(2.0_dp*xi.lt.eta-eps))&
         .or.(.not.(abs(A+zeta).gt.eps).and.(eta.lt.-eps))) then

       Cc=0
       Cc(1,1)=1
       Cc(2,2)=1
       Cc(3,3)=1
       if(zeta.gt.0.0_dp) Cc(2,1)=-1  
       if(zeta.lt.0.0_dp) Cc(2,1)=1
       
       call reset()
       step7=.true.
    else
       step7=.false.
    end if

  end function step7

  function step8()

    logical :: step8

    if((xi+eta+zeta+A+B.lt.-eps).or.(.not.(abs(xi+eta+zeta+A+B).gt.eps).and.(2.0_dp*(A+eta)+zeta.gt.eps))) then

       Cc=0
       Cc(1,1)=1
       Cc(2,2)=1
       Cc(3,3)=1
       Cc(3,1)=1
       Cc(3,2)=1
       
       call reset()
       step8=.true.
    else
       step8=.false.
    endif
    
  end function step8

end subroutine niggli_reduce
