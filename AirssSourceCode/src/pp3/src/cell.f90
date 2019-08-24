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
! This module read, knows and writes the unit cell                                 !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module cell

  use constants

  implicit none

  private

  public :: read_cell
  public :: write_cell
  public :: update_cell
  public :: compact_cell

  ! ** Public data defining the system
  
  real(kind=dp), public :: lattice_abc(6)
  real(kind=dp), public :: lattice_car(3,3)
  real(kind=dp), public :: lattice_rec(3,3)

  real(kind=dp), public :: volume,voveride=0.0_dp
  real(kind=dp), public :: external_pressure(3,3),poveride=-9999999.99_dp

  integer,       public                                 :: num_ions
  integer,       public                                 :: num_spec
  integer,       public                                 :: num_symm

  integer,       public, dimension(:), allocatable      :: ion_species
  integer,       public, dimension(:,:), allocatable    :: ion_equiv

  real(kind=dp), public, dimension(:,:),   allocatable  :: ion_positions
  real(kind=dp), public, dimension(:,:),   allocatable  :: external_forces

  logical, public, dimension(:), allocatable  :: ion_cons

  
  real(kind=dp), public :: symm_ops(3,4,48)

  character(len=80), public :: seedname

  character(len=6), public, dimension(:),     allocatable :: ion_names, ion_names_spec

  logical, public :: fix_cell=.false.
  logical, public :: cluster=.false.
  logical, public :: ionabs
  logical, public :: quiet=.false.

  integer, public :: unit_conv=21

  !---------------------------------------------------!

  ! ** Private IO data

  integer, parameter                           :: unit_cell=20
  integer, parameter                           :: unit_out_cell=21
  character(len=80)                            :: cellfile,convfile

  ! ** Counters and temp variables

  integer :: i,j,k,m,n,ni,num_pot,ns,na,nb,nc

  real(kind=dp) :: ang,temp,dist,scale,v(3),v1(3),x,y,z

  character(len=120)                           :: ctemp,ctemp2 
  character(len=80)                            :: kpspacing
  character(len=80), dimension(:), allocatable :: spec_pot

contains

  subroutine read_cell
    
    ! ** Open the cell file

    write(cellfile,'(a)') trim(seedname)//".cell"

    open(unit=unit_cell,file=cellfile,status="old",err=999)

    lattice_abc = 0.0_dp

    ! ** Read the lattice vectors

    do
       read(unit_cell,'(a)',end=100) ctemp

       if((index(ctemp,'%BLOCK LATTICE_ABC')>0).or.(index(ctemp,'%BLOCK lattice_abc')>0)) then

          read(unit_cell,*) lattice_abc(1:3)
          read(unit_cell,*) lattice_abc(4:6)

          ang = 1.0_dp-cos(dgrd*lattice_abc(4))**2-cos(dgrd*lattice_abc(5))**2-cos(dgrd*lattice_abc(6))**2+&
               2.0_dp*cos(dgrd*lattice_abc(4))*cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6))

          volume = abs(lattice_abc(1)*lattice_abc(2)*lattice_abc(3)*sqrt(abs(ang)))

       end if
       if(trim(ctemp)=='FIX_ALL_CELL : true') fix_cell=.true.
       if(trim(ctemp)=='#CLUSTER') then
          cluster=.true.
          fix_cell=.true.
       end if
    end do
    
100 rewind(unit_cell)

    do
       read(unit_cell,'(a)',end=101) ctemp
       if((index(ctemp,'%BLOCK LATTICE_CART')>0).or.(index(ctemp,'%BLOCK lattice_cart')>0)) then

          read(unit_cell,*) lattice_car(1:3,1)
          read(unit_cell,*) lattice_car(1:3,2)
          read(unit_cell,*) lattice_car(1:3,3)

          volume = abs(lattice_car(1,1)*(lattice_car(2,2)*lattice_car(3,3)-lattice_car(3,2)*lattice_car(2,3))+&
               lattice_car(2,1)*(lattice_car(3,2)*lattice_car(1,3)-lattice_car(1,2)*lattice_car(3,3))+&
               lattice_car(3,1)*(lattice_car(1,2)*lattice_car(2,3)-lattice_car(2,2)*lattice_car(1,3)))
          read(unit_cell,*) ctemp

       end if

    end do

101 rewind(unit_cell)    


    ! ** Read the cell contents

    num_ions=0
    ionabs = .false.

    do
       read(unit_cell,'(a)',end=103) ctemp
       if((index(ctemp,'%BLOCK POSITIONS_FRAC')>0).or.(index(ctemp,'%BLOCK positions_frac')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK POSITIONS_FRAC')>0).or.(index(ctemp,'%ENDBLOCK positions_frac')>0)) exit
             num_ions=num_ions+1
          end do

          exit
       end if

    end do

103 rewind(unit_cell)

    if(num_ions==0) then

       do
          read(unit_cell,'(a)',end=999) ctemp
          if((index(ctemp,'%BLOCK POSITIONS_ABS')>0).or.(index(ctemp,'%BLOCK positions_abs')>0)) then
             do
                read(unit_cell,'(a)') ctemp
                if((index(ctemp,'%ENDBLOCK POSITIONS_ABS')>0).or.(index(ctemp,'%ENDBLOCK positions_abs')>0)) exit
                num_ions=num_ions+1
             end do
             
             exit
          end if
          
       end do
       
       ionabs = .true.
       rewind(unit_cell)
    
    end if

    if(allocated(ion_positions)) deallocate(ion_positions,ion_names,ion_species,ion_equiv,ion_names_spec,external_forces)
    allocate(ion_positions(3,num_ions),external_forces(3,num_ions),ion_cons(num_ions))
    allocate(ion_names(num_ions),ion_names_spec(num_ions))
    allocate(ion_species(num_ions))

    external_forces=0.0_dp
    do
       read(unit_cell,'(a)',end=110) ctemp
       if((index(ctemp,'%BLOCK POSITIONS_')>0).or.(index(ctemp,'%BLOCK positions_')>0)) then
          
          do i=1,num_ions
             read(unit_cell,'(a)') ctemp
             if(index(ctemp,'#')>0) then
                read(ctemp(1:index(ctemp,'#')),*) ion_names(i),ion_positions(:,i)
                if(index(ctemp,'F=')>0) then
                   read(ctemp(index(ctemp,'F=')+2:),*) external_forces(:,i)
                end if
             else
                read(ctemp,*) ion_names(i),ion_positions(:,i)
             endif
          end do

          exit
       end if

    end do
    
110    rewind(unit=unit_cell)

    ! ** Read pseudopotentials

    num_pot = 0

    do
       read(unit_cell,'(a)',end=106) ctemp
       if((index(ctemp,'%BLOCK SPECIES_POT')>0).or.(index(ctemp,'%BLOCK species_pot')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK SPECIES_POT')>0).or.(index(ctemp,'%ENDBLOCK species_pot')>0)) exit
             num_pot=num_pot+1
          end do
          exit
       end if

    end do
    
106 rewind(unit=unit_cell)    

    if(allocated(spec_pot)) deallocate(spec_pot)
    allocate(spec_pot(num_pot))
    
    do
       read(unit_cell,'(a)',end=107) ctemp
       if((index(ctemp,'%BLOCK SPECIES_POT')>0).or.(index(ctemp,'%BLOCK species_pot')>0)) then
          do i=1,num_pot
             read(unit_cell,'(a)',end=999) spec_pot(i)
          end do
          exit
       end if
    end do
107 rewind(unit=unit_cell)    

    ! ** Read other info
    
    kpspacing = ' '
    
    do
       read(unit_cell,'(a)',end=108) ctemp
       if(index(ctemp,'KPOINTS_MP_SPACING')>0) kpspacing=trim(ctemp)
       
    end do
108 rewind(unit=unit_cell)   

    external_pressure=0.0_dp
    do
       read(unit_cell,'(a)',end=109) ctemp
       if((index(ctemp,'%BLOCK EXTERNAL_PRESSURE')>0).or.(index(ctemp,'%BLOCK external_pressure')>0)) then
          read(unit_cell,*,end=999) external_pressure(1,1),external_pressure(1,2),external_pressure(1,3)
          read(unit_cell,*,end=999) external_pressure(2,2),external_pressure(2,3)
          read(unit_cell,*,end=999) external_pressure(3,3)
          external_pressure(2,1)=external_pressure(1,2)
          external_pressure(3,1)=external_pressure(1,3)
          external_pressure(3,2)=external_pressure(2,3)
          exit
       end if
    end do
109 rewind(unit=unit_cell)    

    num_symm=0
    symm_ops=0.0_dp
    do
       read(unit_cell,'(a)',end=111) ctemp
       if((index(ctemp,'%BLOCK SYMMETRY_OPS')>0).or.(index(ctemp,'%BLOCK symmetry_ops')>0)) then
          do
             num_symm=num_symm+1
             read(unit_cell,*,err=99,end=999) symm_ops(1,1,num_symm),symm_ops(2,1,num_symm),symm_ops(3,1,num_symm)
             read(unit_cell,*,err=99,end=999) symm_ops(1,2,num_symm),symm_ops(2,2,num_symm),symm_ops(3,2,num_symm)
             read(unit_cell,*,err=99,end=999) symm_ops(1,3,num_symm),symm_ops(2,3,num_symm),symm_ops(3,3,num_symm)
             read(unit_cell,*,err=99,end=999) symm_ops(1,4,num_symm),symm_ops(2,4,num_symm),symm_ops(3,4,num_symm)
          end do
99        num_symm=num_symm-1
          exit
       end if
    end do
111 rewind(unit=unit_cell)   


    ! ** Read the ionic constraints
    
    ion_cons=.false.
    do
       read(unit_cell,'(a)',end=112) ctemp
       if((index(ctemp,'%BLOCK IONIC_CONSTRAINTS')>0).or.(index(ctemp,'%BLOCK ionic_constraints')>0)) then
          do
             read(unit_cell,'(a)') ctemp
             if((index(ctemp,'%ENDBLOCK IONIC_CONSTRAINTS')>0).or.(index(ctemp,'%ENDBLOCK ionic_constraints')>0)) exit
             read(ctemp,*) i,ctemp2,j,x,y,z
             n=0
             do ni=1,num_ions
                if(ctemp2.eq.ion_names(ni)) n=n+1
                if(n.eq.j) exit
             end do
             ion_cons(ni)=.true.
                             
          end do
          exit
       end if

    end do
112 rewind(unit=unit_cell)
    
    allocate(ion_equiv(num_symm,num_ions))

    if(poveride.gt.-1000.00_dp) then
       external_pressure=0.0_dp
       external_pressure(1,1)=poveride
       external_pressure(2,2)=poveride
       external_pressure(3,3)=poveride
    end if

    if(.not.all(lattice_abc==0.0_dp)) then

       ! ** Convert to cartesian

       lattice_car(:,1) = (/lattice_abc(1),0.0_dp,0.0_dp/)
       lattice_car(:,2) = (/lattice_abc(2)*cos(dgrd*lattice_abc(6)),lattice_abc(2)*sin(dgrd*lattice_abc(6)),0.0_dp/)
       lattice_car(1,3) = lattice_abc(3)*cos(dgrd*lattice_abc(5))
       lattice_car(2,3) = lattice_abc(3)*(cos(dgrd*lattice_abc(4))-cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6)))&
            /sin(dgrd*lattice_abc(6))
       lattice_car(3,3) = sqrt(lattice_abc(3)**2-lattice_car(1,3)**2-lattice_car(2,3)**2)

    else

       lattice_abc(1) = sqrt(lattice_car(1,1)**2+lattice_car(2,1)**2+lattice_car(3,1)**2)
       lattice_abc(2) = sqrt(lattice_car(1,2)**2+lattice_car(2,2)**2+lattice_car(3,2)**2)
       lattice_abc(3) = sqrt(lattice_car(1,3)**2+lattice_car(2,3)**2+lattice_car(3,3)**2)
       lattice_abc(4) = acos(dot_product(lattice_car(:,2),lattice_car(:,3))/lattice_abc(2)/lattice_abc(3))/dgrd
       lattice_abc(5) = acos(dot_product(lattice_car(:,1),lattice_car(:,3))/lattice_abc(1)/lattice_abc(3))/dgrd
       lattice_abc(6) = acos(dot_product(lattice_car(:,1),lattice_car(:,2))/lattice_abc(1)/lattice_abc(2))/dgrd

    end if

    ! ** Calculate the reciprocal lattice vectors

    lattice_rec(1,1)=lattice_car(2,2)*lattice_car(3,3)-lattice_car(3,2)*lattice_car(2,3)
    lattice_rec(2,1)=lattice_car(2,3)*lattice_car(3,1)-lattice_car(3,3)*lattice_car(2,1)
    lattice_rec(3,1)=lattice_car(2,1)*lattice_car(3,2)-lattice_car(3,1)*lattice_car(2,2)
    lattice_rec(1,2)=lattice_car(3,2)*lattice_car(1,3)-lattice_car(1,2)*lattice_car(3,3)
    lattice_rec(2,2)=lattice_car(3,3)*lattice_car(1,1)-lattice_car(1,3)*lattice_car(3,1)
    lattice_rec(3,2)=lattice_car(3,1)*lattice_car(1,2)-lattice_car(1,1)*lattice_car(3,2)
    lattice_rec(1,3)=lattice_car(1,2)*lattice_car(2,3)-lattice_car(2,2)*lattice_car(1,3)
    lattice_rec(2,3)=lattice_car(1,3)*lattice_car(2,1)-lattice_car(2,3)*lattice_car(1,1)
    lattice_rec(3,3)=lattice_car(1,1)*lattice_car(2,2)-lattice_car(2,1)*lattice_car(1,2)

    lattice_rec(:,:)=lattice_rec(:,:)/volume

    ! ** Convert to absolute coordinates

    if(.not.ionabs) then
       do n=1,num_ions
          ion_positions(:,n) = ion_positions(:,n) - floor(ion_positions(:,n))
          ion_positions(:,n) = ion_positions(1,n)*lattice_car(:,1) + ion_positions(2,n)*lattice_car(:,2) &
               + ion_positions(3,n)*lattice_car(:,3)
       end do
    else
       do n=1,num_ions
          ion_positions(:,n) = matmul(lattice_rec,ion_positions(:,n)) 
          ion_positions(:,n) = ion_positions(:,n) - floor(ion_positions(:,n))
          ion_positions(:,n) = ion_positions(1,n)*lattice_car(:,1) + ion_positions(2,n)*lattice_car(:,2) &
               + ion_positions(3,n)*lattice_car(:,3)
       end do
    end if

    ! ** Scale to a new volume specified on command line

    if(voveride>0.0_dp) then
       
       scale = (voveride/volume)**(1.0_dp/3.0_dp)
       
       lattice_abc(1:3) = lattice_abc(1:3)*scale
       lattice_car      = lattice_car*scale
       lattice_rec      = lattice_rec/scale
       ion_positions    = ion_positions*scale
       volume           = volume*scale**3

    end if

    ! ** How many distinct "species"

    num_spec = 0
    do n=1,num_ions
       if(.not.any(ion_names(1:n-1)==ion_names(n))) then
          num_spec = num_spec+1
          ion_names_spec(num_spec) = ion_names(n)
       end if

    end do

    ! ** Set the ion_species array correctly

    do n=1,num_ions
       do m=1,num_spec
          if(ion_names(n)==ion_names_spec(m)) then
             ion_species(n) = m
             exit
          end if
       end do
    end do

    ! ** Identify symmetry equivalent ions

    if(num_symm.gt.1) then

       do n=1,num_ions
          do m=1,num_ions
             
             do ns=1,num_symm
                
                
                do na=-1,1
                   do nb=-1,1
                      do nc=-1,1
                         
                         v(:) = ion_positions(:,n) - matmul(symm_ops(1:3,1:3,ns),ion_positions(:,m)) - &
                              symm_ops(1,4,ns)*lattice_car(1:3,1) - &
                              symm_ops(2,4,ns)*lattice_car(1:3,2) - &
                              symm_ops(3,4,ns)*lattice_car(1:3,3) - &
                              na*lattice_car(:,1)-nb*lattice_car(:,2)-nc*lattice_car(:,3)

                         if(dot_product(v,v).lt.delta) ion_equiv(ns,n)=m

                      end do
                   end do
                end do
                
                
             end do
             
       
          end do
       end do
       
    end if

    ! ** Open a file for the convergence data

    if(.not.quiet) then

       write(convfile,'(a)') trim(seedname)//".conv"       
       open(unit=unit_conv,file=convfile,status="unknown",err=998)

    end if

    close(unit_cell)
    
    return

998 stop 'There is a problem opening the conv file. Stopping.'
999 stop 'There is a problem reading the cell information. Stopping.'

  end subroutine read_cell  

  subroutine write_cell()

    write(cellfile,'(a)') trim(seedname)//"-out.cell"

    open(unit=unit_out_cell,file=cellfile,status="unknown",err=998)

    ! ** Write the lattice vectors

    write (unit_out_cell,'(a)') "%BLOCK LATTICE_ABC"
    write (unit_out_cell,'(3f20.15)') lattice_abc(1:3)
    write (unit_out_cell,'(3f20.15)') lattice_abc(4:6)
    write (unit_out_cell,'(a)') "%ENDBLOCK LATTICE_ABC"
    write (unit_out_cell,*)

!!$    write (unit_out_cell,'(a)') "%BLOCK LATTICE_CART"
!!$    write (unit_out_cell,'(3f20.15)') lattice_car(1:3,1)
!!$    write (unit_out_cell,'(3f20.15)') lattice_car(1:3,2)
!!$    write (unit_out_cell,'(3f20.15)') lattice_car(1:3,3)
!!$    write (unit_out_cell,'(a)') "%ENDBLOCK LATTICE_CART"
!!$    write (unit_out_cell,*)

    ! ** Write the cell contents

!!$    write (unit_out_cell,'(a)') "%BLOCK POSITIONS_ABS"
!!$    do i=1,num_ions
!!$       if((trim(ion_names(i)).ne.'Z')) &
!!$            write (unit_out_cell,'(a4,3f20.15)') trim(ion_names(i)),ion_positions(:,i)
!!$    end do
!!$    write (unit_out_cell,'(a)') "%ENDBLOCK POSITIONS_ABS"
!!$    write (unit_out_cell,*)

    write (unit_out_cell,'(a)') "%BLOCK POSITIONS_FRAC"
    do i=1,num_ions
       v(:) = matmul(lattice_rec,ion_positions(:,i)) 
       if((trim(ion_names(i)).ne.'Z')) &
            write (unit_out_cell,'(a4,3f20.15)') trim(ion_names(i)),v(:)
    end do
    write (unit_out_cell,'(a)') "%ENDBLOCK POSITIONS_FRAC"
    write (unit_out_cell,*)

    if(fix_cell) then
       write (unit_out_cell,'(a)') "FIX_ALL_CELL = .true."
       write (unit_out_cell,*)
    end if

    ! ** Write the kpoint spacing, if required

    if(len_trim(kpspacing)>0) then
       write (unit_out_cell,'(a)') kpspacing
       write (unit_out_cell,*)
    end if


    if(num_symm.eq.0) then

       ! ** Let CASTEP find the symmetry

       write (unit_out_cell,'(a)') 'SYMMETRY_GENERATE'
       write (unit_out_cell,'(a)') 'SNAP_TO_SYMMETRY'
       write (unit_out_cell,*)

    else

       ! ** Use the supplied symmetry operations

       write (unit_out_cell,'(a)') '%BLOCK SYMMETRY_OPS'
       do i=1,num_symm
          write (unit_out_cell,'(3f20.15)') symm_ops(:,1,i)
          write (unit_out_cell,'(3f20.15)') symm_ops(:,2,i)
          write (unit_out_cell,'(3f20.15)') symm_ops(:,3,i)
          write (unit_out_cell,'(3f20.15)') symm_ops(:,4,i)
       end do
       write (unit_out_cell,'(a)') '%ENDBLOCK SYMMETRY_OPS'
       write (unit_out_cell,*)

    end if

    ! ** Write out the pseudopotential information

    if(num_pot>0) then
       write (unit_out_cell,'(a)') "%BLOCK SPECIES_POT"
       do i=1,num_pot
          write (unit_out_cell,'(a)') spec_pot(i)
       end do
       write (unit_out_cell,'(a)') "%ENDBLOCK SPECIES_POT"
       write (unit_out_cell,*)
    end if

    ! ** Write out the applied external pressure

    write (unit_out_cell,'(a)') "%BLOCK EXTERNAL_PRESSURE"
    write (unit_out_cell,'(3f10.5)') external_pressure(1,1),external_pressure(1,2),external_pressure(1,3)
    write (unit_out_cell,'(2F10.5)') external_pressure(2,2),external_pressure(2,3)
    write (unit_out_cell,'(1F10.5)') external_pressure(3,3)
    write (unit_out_cell,'(a)') "%ENDBLOCK EXTERNAL_PRESSURE"
    write (unit_out_cell,*) 

    close(unit=unit_out_cell)

    return
    
998 stop 'There is a problem writing the cell information. Stopping.'

  end subroutine write_cell

  subroutine update_cell()

    ! ** Calculate the volume

    volume = lattice_car(1,1)*(lattice_car(2,2)*lattice_car(3,3)-lattice_car(3,2)*lattice_car(2,3))+&
         lattice_car(2,1)*(lattice_car(3,2)*lattice_car(1,3)-lattice_car(1,2)*lattice_car(3,3))+&
         lattice_car(3,1)*(lattice_car(1,2)*lattice_car(2,3)-lattice_car(2,2)*lattice_car(1,3))
    
    ! ** Calculate the ABC lattice

    lattice_abc(1) = sqrt(lattice_car(1,1)**2+lattice_car(2,1)**2+lattice_car(3,1)**2)
    lattice_abc(2) = sqrt(lattice_car(1,2)**2+lattice_car(2,2)**2+lattice_car(3,2)**2)
    lattice_abc(3) = sqrt(lattice_car(1,3)**2+lattice_car(2,3)**2+lattice_car(3,3)**2)
    lattice_abc(4) = acos(dot_product(lattice_car(:,2),lattice_car(:,3))/lattice_abc(2)/lattice_abc(3))/dgrd
    lattice_abc(5) = acos(dot_product(lattice_car(:,1),lattice_car(:,3))/lattice_abc(1)/lattice_abc(3))/dgrd
    lattice_abc(6) = acos(dot_product(lattice_car(:,1),lattice_car(:,2))/lattice_abc(1)/lattice_abc(2))/dgrd
    
    ! ** Calculate the reciprocal lattice vectors
    
    lattice_rec(1,1)=lattice_car(2,2)*lattice_car(3,3)-lattice_car(3,2)*lattice_car(2,3)
    lattice_rec(2,1)=lattice_car(2,3)*lattice_car(3,1)-lattice_car(3,3)*lattice_car(2,1)
    lattice_rec(3,1)=lattice_car(2,1)*lattice_car(3,2)-lattice_car(3,1)*lattice_car(2,2)
    lattice_rec(1,2)=lattice_car(3,2)*lattice_car(1,3)-lattice_car(1,2)*lattice_car(3,3)
    lattice_rec(2,2)=lattice_car(3,3)*lattice_car(1,1)-lattice_car(1,3)*lattice_car(3,1)
    lattice_rec(3,2)=lattice_car(3,1)*lattice_car(1,2)-lattice_car(1,1)*lattice_car(3,2)
    lattice_rec(1,3)=lattice_car(1,2)*lattice_car(2,3)-lattice_car(2,2)*lattice_car(1,3)
    lattice_rec(2,3)=lattice_car(1,3)*lattice_car(2,1)-lattice_car(2,3)*lattice_car(1,1)
    lattice_rec(3,3)=lattice_car(1,1)*lattice_car(2,2)-lattice_car(2,1)*lattice_car(1,2)
    
    lattice_rec(:,:)=lattice_rec(:,:)/volume

    volume=abs(volume)

  end subroutine update_cell

  subroutine compact_cell

    real(kind=dp) :: lc0(3,3),prod(3,3)

    integer :: i,j

    logical :: changed

    changed=.true.

    do 

       lc0=lattice_car

       changed = .false.

       do i=1,3
          do j=1,3
             prod(i,j) = dot_product(lc0(:,i),lc0(:,j))
          end do
       end do

       do i=1,3
          do j=1,3
             if(i.ne.j) then

                if((abs(prod(i,j)).gt.(prod(j,j)/2.0_dp+1e-13_dp)).and.(.not.changed)) then
                   lattice_car(:,i) = lc0(:,i) - lc0(:,j)*prod(i,j)/abs(prod(i,j))
                   changed=.true.
                end if

             end if
          end do
       end do

       if(.not.changed) exit

    end do

    call update_cell()

  end subroutine compact_cell

end module cell
