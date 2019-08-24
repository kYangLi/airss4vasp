!==================================================================================!
!                                      opt                                         !
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
! This module performs the geometry optimisation                                   !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module opt

  use constants
  use cell
  use pp

  implicit none

  private

  public :: opt_tpsd
  public :: opt_cd

  logical, public :: optimise=.true.
  integer, public :: maxsteps=9999

  !---------------------------------------------------!

  ! ** Private data

  real(kind=dp) :: step=0.001_dp ! ** Set the intial step size
  real(kind=dp) :: etol=1e-10_dp ! ** Set the energy/enthalpy tolerance

contains

  subroutine opt_tpsd(converged)

    !-------------------------------------------------------------!
    ! Two-point Step Size Gradient Methods - Barzilai and Borwein !
    ! IMA Journal of Numerical Analysis (1988) 8, 141-148         !
    !-------------------------------------------------------------!

    logical, intent(out) :: converged

    real(kind=dp) :: e,h,h0,dh
    real(kind=dp) :: ion_positions0(3,num_ions),lattice_car0(3,3)
    real(kind=dp) :: g(3,num_ions),g0(3,num_ions),sig(3,3),s(3,3),s0(3,3)
    real(kind=dp) :: xg,gg,xx

    integer :: i,j,n,nn,steps

    real(kind=dp) :: gamma


    dh = huge(1.0_dp)
    h  = 1.0_dp
    g  = 1.0_dp
    s  = 1.0_dp

    ion_positions0 = ion_positions
    lattice_car0   = lattice_car
    
    steps     = 0
    converged = .true.

    do while(abs(dh).gt.etol)
       
       steps=steps+1

       if(steps.gt.maxsteps) then 
          converged=.false.
          exit
       end if

       h0=h
       g0=g
       s0=s

       if(.not.fix_cell) then
          call eval_pp(e,g,sig)
          sig=sig
       else
          call eval_pp(e,g)
          sig=0.0_dp
       end if

       s = matmul(sig-external_pressure,lattice_car)

!!$       h = e + (external_pressure(1,1)+external_pressure(2,2)+external_pressure(3,3))/3.0_dp*volume
       h = e + external_pressure(3,3)*volume
       
       if(.not.optimise) exit

       xg=0.0_dp
       gg=0.0_dp
!!$       xx=0.0_dp

       do n=1,num_ions
          do i=1,3
             xg = xg + (ion_positions(i,n)-ion_positions0(i,n))*(g(i,n)-g0(i,n))
             gg = gg + (g(i,n)-g0(i,n))*(g(i,n)-g0(i,n))
!!$             xx = xx + (ion_positions(i,n)-ion_positions0(i,n))*(ion_positions(i,n)-ion_positions0(i,n))
          end do
       end do
       do i=1,3
          do j=1,3
             xg = xg + (lattice_car(i,j)-lattice_car0(i,j))*(s(i,j)-s0(i,j))
             gg = gg + (s(i,j)-s0(i,j))*(s(i,j)-s0(i,j))
!!$             xx = xx + (lattice_car(i,j)-lattice_car0(i,j))*(lattice_car(i,j)-lattice_car0(i,j))
          end do
       end do

       if(abs(xg).gt.0.0_dp) step = abs(xg/gg)

       ion_positions0 = ion_positions
       ion_positions  = ion_positions + step*g

       ! * Convert to fractional

       do n=1,num_ions
          ion_positions0(:,n) = matmul(lattice_rec,ion_positions0(:,n))
          ion_positions(:,n)  = matmul(lattice_rec,ion_positions(:,n))
          ion_positions0(:,n) = ion_positions0(:,n) - floor(ion_positions(:,n))
          ion_positions(:,n)  = ion_positions(:,n)  - floor(ion_positions(:,n))
       end do

       if(.not.fix_cell) then
          lattice_car0   = lattice_car
          lattice_car    = lattice_car + step*s
          call update_cell()
       end if

       ! * Convert to absolute

       do n=1,num_ions
          ion_positions0(:,n) = ion_positions0(1,n)*lattice_car(:,1) + ion_positions0(2,n)*lattice_car(:,2) &
               + ion_positions0(3,n)*lattice_car(:,3)
          ion_positions(:,n) = ion_positions(1,n)*lattice_car(:,1) + ion_positions(2,n)*lattice_car(:,2) &
               + ion_positions(3,n)*lattice_car(:,3)
       end do

       if(.not.quiet) then
          write (unit_conv,'(f35.25)') h
          call flush(unit_conv)
       end if
       dh=h-h0
       
    end do

    write (*,'(a11,f20.10)') 'Volume:    ',  volume
    write (*,'(a11,f20.10)') 'Pressure:  ',  (sig(1,1)+sig(2,2)+sig(3,3))/3.0_dp
    write (*,'(a11,f20.10)') 'Energy:    ',  e
    write (*,'(a11,f20.10)') 'Enthalpy:  ',  h
    write (*,'(a11,6f10.5)') 'Stress:    ',  sig(1,1),sig(2,2),sig(3,3),sig(1,2),sig(1,3),sig(2,3)

  end subroutine opt_tpsd

  subroutine opt_cd(converged)

    logical, intent(out) :: converged

    integer, parameter :: nsteps=100

    integer :: n,i,j

    real(kind=dp) :: pos0,g(3,num_ions),h,e,emin,pmin,ebest,sig(3,3)=0.0_dp
    converged=.true.
    do 

       do n=1,num_ions
          do i=1,3
             
             pos0=ion_positions(i,n)
             
             emin=huge(1.0_dp)

             do j=1,nsteps
             
                ion_positions(i,n)=pos0+real(j-nsteps/2,dp)/real(nsteps/2,dp)*10.0_dp
                call eval_pp(e)
                if(e.lt.emin) then
                   emin=e
                   pmin=ion_positions(i,n)
                end if

             end do
             ion_positions(i,n)=pmin
          end do
       end do
       call eval_pp(e)
       if(e.lt.ebest) then
          ebest=e
       else
          exit
       end if
    end do

    h=e

    write (*,'(a11,f20.10)') 'Volume:    ',  volume
    write (*,'(a11,f20.10)') 'Pressure:  ',  (sig(1,1)+sig(2,2)+sig(3,3))/3.0_dp
    write (*,'(a11,f20.10)') 'Energy:    ',  e
    write (*,'(a11,f20.10)') 'Enthalpy:  ',  h
    write (*,'(a11,6f10.5)') 'Stress:    ',  sig(1,1),sig(2,2),sig(3,3),sig(1,2),sig(1,3),sig(2,3)

  end subroutine opt_cd

end module opt
