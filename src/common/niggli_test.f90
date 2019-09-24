!==================================================================================!
!                                   niggli_test                                    !
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
! This program tests niggli_reduce                                                 !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program niggli_test

  use constants

  implicit none

  real(kind=dp) :: Lmat(3,3),abc(6)
  integer       :: Cmat(3,3)

  ! ** The test case from Acta Cryst. (1976). A32, 297-298
  
  abc(1)=3.0_dp
  abc(2)=5.19615242270663188058233902452_dp
  abc(3)=2.0_dp 
  abc(4)=1.81374177014849137201968687323_dp/dgrd
  abc(5)=1.91063323624901855632771420503_dp/dgrd
  abc(6)=2.35413687121593310838372729124_dp/dgrd
  
  write (stdout,'(3f10.5)') abc
  write (stdout,*)

  Lmat(:,1) = (/abc(1),0.0_dp,0.0_dp/)
  Lmat(:,2) = (/abc(2)*cos(dgrd*abc(6)),abc(2)*sin(dgrd*abc(6)),0.0_dp/)
  Lmat(1,3) = abc(3)*cos(dgrd*abc(5))
  Lmat(2,3) = abc(3)*(cos(dgrd*abc(4))-cos(dgrd*abc(5))*cos(dgrd*abc(6)))/sin(dgrd*abc(6))
  Lmat(3,3) = sqrt(abc(3)**2-Lmat(1,3)**2-Lmat(2,3)**2)
  
  !write (stdout,'(3f10.5)') transpose(Lmat)
  write (stdout,'(3f10.5)') (Lmat)
  write (stdout,*)

  call niggli_reduce(Lmat,Cmat)
  write (stdout,*)

  !write (stdout,'(3f10.5)') transpose(Lmat)
  write (stdout,'(3f10.5)') (Lmat)
  write (stdout,*)
  
  abc(1) = sqrt(Lmat(1,1)**2+Lmat(2,1)**2+Lmat(3,1)**2)
  abc(2) = sqrt(Lmat(1,2)**2+Lmat(2,2)**2+Lmat(3,2)**2)
  abc(3) = sqrt(Lmat(1,3)**2+Lmat(2,3)**2+Lmat(3,3)**2)
  abc(4) = acos(dot_product(Lmat(:,2),Lmat(:,3))/abc(2)/abc(3))/dgrd
  abc(5) = acos(dot_product(Lmat(:,1),Lmat(:,3))/abc(1)/abc(3))/dgrd
  abc(6) = acos(dot_product(Lmat(:,1),Lmat(:,2))/abc(1)/abc(2))/dgrd
  
  write (stdout,'(3f10.5)') abc

end program niggli_test

