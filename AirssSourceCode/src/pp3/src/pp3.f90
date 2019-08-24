!==================================================================================!
!                                      cell                                        !
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
! Thisis the main program of pp3 - a simple 3D pair potential optimiser            !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program pp3

  use constants
  use cell
  use pp
  use opt

  implicit none

  logical :: converged

  ! ** Get the command line argument

  call get_arguments(seedname)

  ! ** Read the cell data

  call read_cell()

  if (optimise.and.(.not.fix_cell)) call compact_cell()

  ! ** Read the pair potential

  call read_pp()
  !call write_pp()

  ! ** Initialise the pair potential

  call init_pp()

  ! ** Perform the geometry optimisation

  call opt_tpsd(converged)
  if (optimise.and.(.not.fix_cell)) call compact_cell()
  call opt_tpsd(converged)

  if(.not.converged) write(*,'(a)') "FAILED"

  if(optimise.and.(.not.fix_cell)) call compact_cell()

  ! ** Write the final structure

  if(optimise) call write_cell()

contains

  ! *******************************************************************************************

  subroutine get_arguments(name)

    character(len=80), intent(out) :: name

    character(len=80), allocatable,dimension(:) :: cargs

    character(len=2) :: flags(10)=(/"-f","-c","-l","-d","-n","-e","-q","-m","-v","-p"/)

    integer :: i,iargc,num_args

    num_args = iargc()

    allocate(cargs(num_args))

    do i=1,num_args
       call getarg(i,cargs(i))
    end do

    name = trim(cargs(num_args))

    do i=1,num_args-1
       if(cargs(i)(1:1)=="-") then
          if(.not.any(cargs(i).eq.flags)) goto 99
       end if
    end do

    if(any(cargs.eq."-n")) optimise=.false.   ! * Do not relax the structure
    if(any(cargs.eq."-e")) forceshift=.false. ! * Energy shift only
    if(any(cargs.eq."-q")) quiet=.true.       ! * Quiet, minimal output
    if(any(cargs.eq."-m")) then               ! * Maximum number of optimisation steps
       do i=1,num_args-1
          if(cargs(i).eq."-m") exit
       end do
       read(cargs(i+1),*) maxsteps
    end if
    if(any(cargs.eq."-v")) then               ! * Volume override
       do i=1,num_args-1
          if(cargs(i).eq."-v") exit
       end do
       read(cargs(i+1),*) voveride
    end if
    if(any(cargs.eq."-p")) then               ! * Pressure override
       do i=1,num_args-1
          if(cargs(i).eq."-p") exit
       end do
       read(cargs(i+1),*) poveride
    end if
    
    deallocate(cargs)

    return

99  write (*,*) 'Usage: pp3 [-n] [-e] [-q] [-m maxsteps] [-v voveride] [-p poveride] <seedname>'
    write (*,*) '-n : Do not relax the structure'
    write (*,*) '-e : Energy shift only, no force shift'
    write (*,*) '-q : Quiet - minimal output'
    write (*,*) '-m : Maximum number of optimisation steps'
    write (*,*) '-v : Volume overide'
    write (*,*) '-p : Pressure overide'

    stop

  end subroutine get_arguments

  ! *******************************************************************************************

end program pp3
