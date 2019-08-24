!==================================================================================!
!                                       rng                                        !
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
! This module handles random number generation                                     !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module rng

  use constants

  implicit none

!!$  interface ! Needed for ifort
!!$     integer function getpid()
!!$     end function getpid
!!$  end interface

  private

  public :: init_pseudorandom
  public :: random_single
  public :: random_double
  public :: random_triple
  public :: random_integer
  
  integer, public, allocatable, dimension(:) :: seed

contains

  subroutine init_pseudorandom(sstring)

    integer, parameter                 :: maxr = HUGE(0) - 1
    integer(kind=8)                    :: count,countrate,countmax,pid
    integer                            :: i,j,iters,size
    real(kind=dp)                      :: rr
    character(len=*), optional         :: sstring
 
    pid=getpid()

    ! ** Seed the random number generator

    call random_seed(size=size)
    
    if(allocated(seed)) deallocate(seed)
    
    allocate(seed(size))

    if(present(sstring).and.(len_trim(sstring).gt.0)) then
       seed=0
       read(sstring,*,end=99,err=99) seed
       call random_seed(put=seed)
       return
99     stop 'init_pseudorandom : error reading seed'
    end if

    call random_seed(get=seed)

    do j=1,size
       call system_clock(count,countrate,countmax)
       iters = int(mod(count,pid),4)
       do i=1,iters+1
          call random_number(rr)
       end do       
       seed(j) = int(2*(rr-0.5_dp)*maxr)
    end do

    call random_seed(put=seed)
    
  end subroutine init_pseudorandom

  function random_single()
   
    real(kind=dp) :: random_single
   
    call random_number(random_single)

  end function random_single

  function random_double()
   
    real(kind=dp) :: random_double(2)
   
    call random_number(random_double)

  end function random_double

  function random_triple()
   
    real(kind=dp) :: random_triple(3)
   
    call random_number(random_triple)

  end function random_triple

  function random_integer(n)

    integer, intent(in) :: n

    integer :: random_integer
    
    real(kind=dp) :: r
    
    call random_number(r)
    
    random_integer=int(r*n)+1
    
  end function random_integer
  
end module rng
