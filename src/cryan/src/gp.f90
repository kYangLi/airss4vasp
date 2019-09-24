! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!==================================================================================!
!                                       gp                                         !
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
! This module performs gaussian process regression                                 !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module gp

  use constants
  use algor

  implicit none

  private

  integer :: ntrain, nsets

  real(kind=dp), allocatable, dimension(:,:)     :: C,Cinv,ytrain,alpha
  real(kind=dp), allocatable, dimension(:)       :: Cs,xtrain

  real(kind=dp) :: sigma_n,l_scale,determinant,likelihood

  public :: gp_initialise
  public :: gp_predict
  
contains

  subroutine gp_initialise(xt,yt,ls,sn)

    real(kind=dp), dimension(:),     intent(in) :: xt
    real(kind=dp), dimension(:,:),   intent(in) :: yt
    real(kind=dp),                   intent(in) :: ls
    real(kind=dp),                   intent(in) :: sn

    integer :: ns
    
    ntrain  = size(xt(:))
    nsets   = size(yt(1,:))
    
    l_scale=ls
    
    sigma_n = sn

    if(allocated(xtrain)) deallocate(xtrain,ytrain,alpha,C,Cinv,Cs)

    allocate(xtrain(ntrain),ytrain(ntrain,nsets),alpha(ntrain,nsets))
    allocate(C(ntrain,ntrain),Cinv(ntrain,ntrain),Cs(ntrain))

    xtrain=xt
    ytrain=yt

    ! ** Construct the covariance matrix 
       
    C=covariance(Cinv,det=determinant)
   
    do ns=1,nsets

       ! ** Store alpha=Cinv * Yt

       alpha(:,ns)=matmul(Cinv,ytrain(:,ns))

    end do
    
    likelihood=gp_likelihood()

    write (stderr,'(a,g20.5)') 'Likelihood: ',likelihood  

  end subroutine gp_initialise

  function gp_likelihood()

    real(kind=dp) :: gp_likelihood

    integer :: ns

    gp_likelihood=0.0_dp
    do ns=1,nsets
       gp_likelihood=gp_likelihood+&
            (dot_product(ytrain(:,ns),alpha(:,ns))+log(determinant)+ntrain*log(tpi))/2.0_dp
    end do

  end function gp_likelihood
  
  function gp_predict(x,n)

    real(kind=dp),               intent(in)  :: x
    integer,                     intent(in)  :: n
    
    real(kind=dp) :: gp_predict
    
    real(kind=dp), allocatable, dimension(:) :: Cstar
    
    integer :: nt

    allocate(Cstar(ntrain))
    
    do nt=1,ntrain
       Cstar(nt) = kernel(x,xtrain(nt))
    end do
    
    gp_predict=dot_product(Cstar(:),alpha(:,n))

  end function gp_predict
  
  function covariance(covinv,det)

    real(kind=dp), optional, dimension(ntrain,ntrain),        intent(out) :: covinv
    real(kind=dp), optional,                                  intent(out) :: det
    
    real(kind=dp), dimension(ntrain,ntrain) :: covariance

    integer :: nt1,nt2

    do nt1=1,ntrain
       do nt2=1,ntrain
          covariance(nt1,nt2)=kernel(xtrain(nt1),xtrain(nt2),addn=.true.)
       end do
    end do

    if(present(covinv)) then
       if(present(det)) then
          covinv=inv_svd(covariance,det)
       else
          covinv=inv_svd(covariance)
       end if
    endif
    
  end function covariance

  function kernel(x,xp,addn)

    real(kind=dp),            intent(in)  :: x
    real(kind=dp),            intent(in)  :: xp
    logical,       optional,  intent(in)  :: addn

    real(kind=dp) :: kernel

    logical :: xeqxp

    xeqxp=((x-xp)**2.lt.epsilon(1.0_dp))
    if(present(addn).and.(.not.addn)) xeqxp=.false.
    if(.not.present(addn)) xeqxp=.false.

    if(l_scale.lt.-delta) then

       ! ** RBF

       kernel=exp(-(x-xp)**2/abs(l_scale)**2/2.0_dp)

    else if(l_scale.gt.delta) then

       ! ** Matern32

       kernel=(1.0_dp+abs(x-xp)/abs(l_scale))*exp(-abs(x-xp)/abs(l_scale))
       
    else

       ! ** EXP

       kernel=exp(-abs(x-xp)/10.0_dp)
       
    end if

    if(xeqxp) kernel=kernel+sigma_n**2

  end function kernel

end module gp
