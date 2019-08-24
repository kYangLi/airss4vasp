!==================================================================================!
!                                     algor                                        !
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
! This module provides useful numerical algorithms                                 !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module algor

  use constants

  implicit none

contains

  function sle(A,y)

    real(kind=dp), dimension(:,:), intent(in) :: A
    real(kind=dp), dimension(:),   intent(in) :: y

    real(kind=dp) :: sle(size(A(:,1))),AA(size(A(:,1)),size(A(:,2)))
   
    integer :: n,info
    
    n=size(sle)

    sle=y
    AA=A
    call dpotrf('L', n, AA, n, info)
    if (info /= 0) stop 'sle: matrix is not positive definite, or illegal value in dpotrf!'
    call dpotrs('L', n, 1, AA, n, sle, n, info)
    if (info /= 0) stop 'sle: illegal value in dpotrf!'

  end function sle
    

  function inv_svd(A,det) result(Ainv)
    
    real(kind=dp), dimension(:,:),    intent(in)    :: A
    real(kind=dp), optional,          intent(out)   :: det
    
    real(kind=dp), dimension(size(A,1),size(A,1))   :: Ainv

    real(kind=dp) :: work(size(A,1)*10),u(size(A,1),size(A,1)),vt(size(A,1),size(A,1))
    real(kind=dp) :: mxv,sing(size(A,1)),s(size(A,1),size(A,1))
    real(kind=dp) :: alpha=1e-12_dp

    integer :: i,n,info,nc
    
    Ainv = A
    n = size(A,1)

    call dgesvd('A','A',n,n,Ainv,n,sing,u,n,vt,n,work,10*n,info)

    mxv=maxval(sing)
    
    s=0.0_dp
    nc=0
    do i=1,n
       if(sing(i).gt.alpha*mxv) then
          s(i,i)=1.0_dp/sing(i)
       else
          nc=nc+1
       end if
    end do
    
    !write (*,'(a,i5)') 'Number of singular values: ',nc
    
    Ainv=matmul(matmul(transpose(vt),s),transpose(u))

    if(present(det)) then
       det=sing(1)
       do i=2,n
          if(abs(sing(i)).gt.0.0_dp) det=det*sing(i)
       end do
    end if
    
  end function inv_svd

  
end module algor
