! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!==================================================================================!
!                                  buildcell                                       !
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
! This is the main program of buildcell, a cell builder                            !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program buildcell

  use constants
  use cell
  use symmetry
  use rng
  use build
  use pp
  use opt

  implicit none

  integer       :: status,nreject,ntotal,nattempts,symmstat,nstar

  real(kind=dp) :: temp,time0,start_time,timed,time,over

  logical :: conv=.false.,had_success=.false.

  ntotal    = 0
  
  ! ** Initialise the timer

  call init_timer()

  ! ** Get the command line argument

  call get_arguments()

  ! ** Read the cell data

  call read_cell()

  ! ** Seed the random number generator

  call get_seed()
  
  call init_pseudorandom(seed_string)
  
  ! ** Reload the cell buffer

99 call reset_cellbuff()

  ! ** Eliminate commented out sections

  call comments_cell()

  ! ** Eliminate blank lines

  call blank_cell()
  
  ! ** Expansion

  call expand_cell()

  ! ** Atoms supplied by species

  call by_species_cell()

  ! ** Adjust the cell, allowing for variable sets

  call adjust_cell()

  ! ** Initialise the unit cell

  call init_cell(status)

  if(status.eq.888) then
     initialised=.false.
     cell_consistent=.false.
     goto 99
  end if
  
  if(autoslack) slack=slack+slack_step

  if(slack.gt.0.0_dp) write (stderr,'(a,f6.3)') 'Applying slack of:',slack 
  
  ! ** Identify the sets of atoms

  call gen_sets()
  
  ! ** Build the crystal cell

  nreject   = 0
  nattempts = 0
  nstar     = 0

  status  = 1

  had_success=.false.

  do while (status/=0)

     status=0
     
     ! ** Get the new lattice vectors
     
     call gen_lattice(status)
     
     if(status.eq.1) then
        initialised=.false.
        cell_consistent=.false.
        write(stderr,*)
        goto 99
     end if

     nreject=nreject+1
     ntotal=ntotal+1     
     
     ! ** Combined position generation and incell separation check

     call gen_position_check(status)

     if(status.eq.999) then
        nstar=nstar+1
     else
        nstar=0
        had_success=.true.
     end if
     
     ! ** Full separation check, push atoms if needed

     call push_opt(status)

     ! ** Check the co-ordination

     !call check_coordination(status)

     ! ** Assign the spin to atoms, if required

     call assign_spin(status)

     if(status.eq.888) then
        nudgen=1 ! Make more use of special positions
        initialised=.false.
        cell_consistent=.false.
        write(stderr,*)
        goto 99
     end if

     ! ** Check if we are spending too long in this symmetry

     call cpu_time(temp)
     if(((temp-time0).gt.maxtime).or.((nstar.ge.10).and..not.had_success)) then
        call cpu_time(time0)
        if(push) write(stderr,*)
        write (stderr,'(a,i8)') 'Number of rejects:    ',nreject
        if(status==0) exit
        nreject=0
        nstar=0
        had_success=.false.
        nattempts=nattempts+1
        if(nattempts.gt.10) then
           initialised=.false.
           cell_consistent=.false.
           nudgen=nudgen+1
           if(nudgen.gt.num_ions-adjgen) nudgen=0
           goto 99
        end if
        call select_symmetry(symmstat)
        if(symmstat.gt.0) then
           initialised=.false.
           cell_consistent=.false.
           goto 99
        end if
        write (stderr,'(a,a,a,i3)') 'Symm: ',symmetry_name,' Nops: ',num_symm
        initialised=.false.
        cell_consistent=.false.
     end if
     
  end do

  if(push) write(stderr,*)

  ! ** Now adjust using the pair-potentials, relax and shake if required

  if((overlap.gt.-999.9_dp)) then

     ! ** Decide whether to vary the cell shape or not

     if((fix_cell).or.(abs(cellamp).lt.delta)) celladapt=.false.

     if(celladapt) write (stderr,'(a)') 'Varying cell, volume fixed'

     write (stderr,'(a)') 'Using hard sphere potentials'

     if(dorash) then

        write (stderr,'(a)') 'Relax and shake'

        call rash(conv,over)
     else
        call opt_tpsd(conv,over)
     end if

     if(.not.conv)  write (stderr,'(a)') 'Did not converge'

     conv=.false.

     write (stderr,'(a,f16.5/)') 'Overlap:',over/real(num_ions_total,dp)
     
     if(overlap+delta.gt.over/real(num_ions_total,dp)) conv=.true.

     if(.not.conv) then
        initialised=.false.
        cell_consistent=.false.
        goto 99
     end if

     call compact_cell()

     call update_cell()

  end if

  if(ntotal.gt.1) write (stderr,'(a,i8)') 'Total rejects:     ',ntotal-1

  if(sum(rejected).gt.0) then
     write (stderr,*)
     if(rejected(1)>0) write (stderr,'(a,f6.2,a)') &
          'Initial separation:       ', real(rejected(1),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(2)>0) write (stderr,'(a,f6.2,a)') &
          'Star contraction:         ', real(rejected(2),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(3)>0) write (stderr,'(a,f6.2,a)') &
          'Secondary separation:     ', real(rejected(3),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(4)>0) write (stderr,'(a,f6.2,a)') &
          'Tertiary separation:      ', real(rejected(4),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(5)>0) write (stderr,'(a,f6.2,a)') &
          'Clearing gap:             ', real(rejected(5),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(6)>0) write (stderr,'(a,f6.2,a)') &
          'Too many neighbours:      ', real(rejected(6),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(7)>0) write (stderr,'(a,f6.2,a)') &
          'Wrong type of neighbours: ', real(rejected(7),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(8)>0) write (stderr,'(a,f6.2,a)') &
          'Wrong coordination:       ', real(rejected(8),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(9)>0) write (stderr,'(a,f6.2,a)') &
          'Bond angle too small:     ', real(rejected(9),kind=dp)/(sum(rejected))*100.0_dp,'%'
     if(rejected(10)>0) write (stderr,'(a,f6.2,a)') &
          'Bond angle too large:     ', real(rejected(10),kind=dp)/(sum(rejected))*100.0_dp,'%'
     write (stderr,*)
  end if

  ! ** Write the lattice vectors

  call write_cell()

  ! ** Report timing info

  call report_timer()

contains

  ! *******************************************************************************************

  subroutine get_arguments()

    integer :: iargc,num_args

    num_args = iargc()

    if(num_args.ne.0) then
       write (stderr,'(a)') 'Usage: buildcell < in.cell > out.cell'
       stop
    end if

  end subroutine get_arguments

  subroutine init_timer()

    call cpu_time(start_time)
    call cpu_time(time0)
    timed=0.0_dp

  end subroutine init_timer

  subroutine start_timer()

    call cpu_time(time)

  end subroutine start_timer

  subroutine stop_timer()

    call cpu_time(temp)

    timed=timed+temp-time

  end subroutine stop_timer

  subroutine report_timer()

    write (stderr,'(a,f15.5)') 'Timed:',timed
    call cpu_time(temp)
    temp=temp-start_time
    write (stderr,'(a,f15.5)') 'Total:',temp
    write (stderr,'(a,f11.1,a)') 'Perc: ',timed/temp*100.0_dp,'%'

  end subroutine report_timer

  ! *******************************************************************************************

end program buildcell
