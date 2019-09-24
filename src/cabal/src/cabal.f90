!==================================================================================!
!                                     cabal                                        !
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
! This program converts structure formats                                          !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program cabal

  use constants

  implicit none

  integer, parameter :: max_words=120
  integer, parameter :: max_species=20
  
  integer :: num_lines=1000000 ! ** STDIN can only be read once, so need to estimate a maximum input length

  integer :: num_ions,num_species

  integer, allocatable, dimension(:) :: num_words,ion_index
  
  real(kind=dp)                                   :: lattice_abc(6)=0.0_dp,lattice_car(3,3)=0.0_dp,lattice_rec(3,3)=0.0_dp
  real(kind=dp)                                   :: lattice_volume=0.0_dp
  real(kind=dp)                                   :: pad=5.0_dp
  real(kind=dp),      allocatable, dimension(:)   :: ion_occ,ion_spin
  real(kind=dp),      allocatable, dimension(:,:) :: ion_fractional,ion_absolute
  

  character(len=10)                                       :: inext,outext,species_names(max_species)=''
  character(len=10),          allocatable, dimension(:)   :: ion_names
  character(len=max_words*2), allocatable, dimension(:)   :: buff
  character(len=max_words*2), allocatable, dimension(:,:) :: buff_words

  logical :: x2x

  ! * =======================================================================
    
  call get_arguments()

  call read_buff()

  call read_input()
  
  call consolidate()

  if(x2x) call niggli()
  
  call write_output()
  
contains

  subroutine get_arguments()
    
    integer :: iargc

    character(len=20) :: ctemp

    if((iargc().lt.2).or.(iargc().gt.3)) then
       write (stderr,'(a)') 'Usage: cabal in out [x] < seed.in > seed.out'
       write (stderr,'(a)') '  in==out : Niggli reduce'
       write (stderr,'(a)') '  supports castep+,cell,shx,res,gulp*,cif*,psi4*,xtl,xyz(e),config+,poscar'
       write (stderr,'(a)') '  *output only +input only'
       write (stderr,'(a)')
       write (stderr,'(a)') '  x : for xyz input, pad unit cell with x'
       stop
    end if
    
    call getarg(1,inext)
    call getarg(2,outext)

    if(iargc().eq.3) then
       call getarg(3,ctemp)
       read (ctemp,*) pad
    end if

    x2x=inext.eq.outext
    
  end subroutine get_arguments

  subroutine read_buff()

    integer :: n,nw,indx

    character(len=max_words*2) :: line

    allocate(buff(num_lines))
    n=0
    do 
       n=n+1
       if(n.gt.num_lines) stop 'cabal: increase num_lines'
       read(5,'(a)',end=99,err=100) buff(n)
    end do
100 stop 'cabal: problem reading input data'
99  continue
    
    num_lines=n-1

    allocate(buff_words(max_words,num_lines),num_words(num_lines))

    do n=1,num_lines

       nw=0
       line=trim(adjustl(detab(buff(n))))
       do while (len_trim(line).gt.0)
          nw=nw+1
          indx=index(line,' ')
          buff_words(nw,n)=trim(line(:indx))
          line=adjustl(line(indx+1:))
       end do
       num_words(n)=nw

    end do
    
  end subroutine read_buff

  subroutine read_input()

    select case (inext)
    case('castep')
       call read_castep()
    case('cell')
       call read_cell()
    case('res','shx')
       call read_res()
    case('xyz')
       call read_xyz()
    case('xyze')
       call read_xyze()
    case('gulp')
       call read_gulp()
    case('xtl')
       call read_xtl()
    case('cif')
       call read_cif()
    case('psi4')
       call read_psi4()
    case('config')
       call read_config()
    case('poscar')
       call read_poscar()
    case default
       stop 'cabal: input format not recognised'
    end select
    
  end subroutine read_input

  subroutine write_output()

    select case (outext)
    case('cell')
       call write_cell()
    case('castep')
       call write_castep()
    case('res','shx')
       call write_res()
    case('xyz')
       call write_xyz()
    case('xyze')
       call write_xyze()
    case('gulp')
       call write_gulp()
    case('xtl')
       call write_xtl()
    case('cif')
       call write_cif()
    case('psi4')
       call write_psi4()
    case('config')
       call write_config()
    case('poscar')
       call write_poscar()
    case default
       stop 'cabal: output format not recognised'
    end select
    
  end subroutine write_output

  subroutine consolidate()

    integer :: n,ni

    if(all(lattice_car.eq.0.0_dp)) call abc2car(lattice_abc,lattice_car)

    if(all(lattice_abc.eq.0.0_dp)) call car2abc(lattice_car,lattice_abc)

    call car2rec(lattice_car,lattice_rec)

    if(all(ion_fractional.eq.0.0_dp)) then
       do ni=1,num_ions
          ion_fractional(:,ni)=matmul(lattice_rec,ion_absolute(:,ni))
       end do
    end if

    if(all(ion_absolute.eq.0.0_dp)) then
       do ni=1,num_ions
          ion_absolute(:,ni)=matmul(lattice_car,ion_fractional(:,ni))
       end do
    end if

    if(all(ion_index.eq.0)) then
       num_species=0
       species_names=''
       do ni=1,num_ions
          if(.not.any(species_names==ion_names(ni))) then
             num_species=num_species+1
             species_names(num_species)=ion_names(ni)
          end if
       end do

       do ni=1,num_ions
          do n=1,num_species
             if(ion_names(ni).eq.species_names(n)) exit
          end do
          ion_index(ni)=n
       end do
       
    end if

    lattice_volume=car_volume(lattice_car)

  end subroutine consolidate

  subroutine niggli()

    real(kind=dp) :: lattice_car0(3,3)
    integer       :: Cmat(3,3)

    integer :: ni

    interface
       subroutine niggli_reduce(Lm,Cm)
         use constants
         real(kind=dp), dimension(3,3), intent(inout) :: Lm
         integer,       dimension(3,3), intent(out)   :: Cm
       end subroutine niggli_reduce
    end interface
    
    lattice_car0=lattice_car
    call niggli_reduce(lattice_car,Cmat)
    call car2rec(lattice_car,lattice_rec)    
    call car2abc(lattice_car,lattice_abc)

    do ni=1,num_ions
       ion_fractional(:,ni)=matmul(lattice_rec,matmul(lattice_car0,ion_fractional(:,ni)))
       ion_absolute(:,ni)=matmul(lattice_car,ion_fractional(:,ni))
    end do
    
  end subroutine niggli
  

  ! * =======================================================================

  ! ** Castep output file

  subroutine read_castep()

    integer :: n,i,ni

    do n=1,num_lines
       if(index(buff(n),'Real Lattice(A)').gt.0) then
          do i=1,3
             read(buff_words(1:3,n+i),*) lattice_car(1:3,i)
          end do
       end if
    end do
    
    do n=1,num_lines
       if(index(buff(n),'Total number of ions in cell =').gt.0) read(buff_words(8,n),*) num_ions
    end do

    call init_ions()

    do n=1,num_lines
       if(index(buff(n),'Number           u          v          w').gt.0) then
          do ni=1,num_ions
             read(buff_words(2,n+ni+1),*) ion_names(ni)
             ion_names(ni)=trim(strip(ion_names(ni)))
             read(buff_words(4:6,n+ni+1),*) ion_fractional(1:3,ni)
          end do
       end if
       if(index(buff(n),'Species          Ion Spin      s      p      d      f').gt.0) then
          do ni=1,num_ions
             read(buff_words(10,n+2*ni),*) ion_spin(ni)
          end do
       end if
    end do
    
  end subroutine read_castep

  subroutine write_castep()
    stop 'write_castep: not implemented'
  end subroutine write_castep
  
  ! ** Castep cell file format
  
  subroutine read_cell()

    integer :: n,i,ni,nstart,nfinish,indx

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_CART'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_CART'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.eq.4) then
          do i=1,3
             read(buff_words(1:3,nstart+i),*) lattice_car(1:3,i)
          end do
       else if(nfinish-nstart.eq.5) then
          if(up(buff_words(1,nstart+1)).ne.'ANG') stop 'cabal: LATTICE_CART units not recognised'
          do i=1,3
             read(buff_words(1:3,nstart+i+1),*) lattice_car(1:3,i)
          end do
       else
          stop 'cabal: problem reading LATTICE_CART in cell data'
       end if
    end if

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_ABC'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'LATTICE_ABC'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.ne.3) stop 'cabal: problem reading LATTICE_ABC in cell data'
       read(buff_words(1:3,nstart+1),*) lattice_abc(1:3)
       read(buff_words(1:3,nstart+2),*) lattice_abc(4:6)
    end if

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_FRAC'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_FRAC'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.lt.2) stop 'cabal: problem reading POSITIONS_FRAC in cell data'
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          ion_names(ni)=trim(buff_words(1,nstart+ni))
          read(buff_words(2:4,nstart+ni),*) ion_fractional(1:3,ni)
          indx=index(buff(nstart+ni),'SPIN=')
          if(indx.gt.0) then
             read(buff(nstart+ni)(indx+5:),*) ion_spin(ni)
          end if
       end do
    end if

    nstart=0
    nfinish=0
    do n=1,num_lines
       if((up(buff_words(1,n)).eq.'%BLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_ABS'))) then
          nstart=n
       end if
       if((up(buff_words(1,n)).eq.'%ENDBLOCK').and.((up(buff_words(2,n)).eq.'POSITIONS_ABS'))) then
          nfinish=n
          exit
       end if
    end do
    if(nstart.lt.nfinish) then
       if(nfinish-nstart.lt.2) stop 'cabal: problem reading POSITIONS_ABS in cell data'
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          ion_names(ni)=trim(buff_words(1,nstart+ni))
          read(buff_words(2:4,nstart+ni),*) ion_absolute(1:3,ni)
          indx=index(buff(nstart+ni),'SPIN=')
          if(indx.gt.0) then
             read(buff(nstart+ni)(indx+5:),*) ion_spin(ni)
          end if
       end do
    end if
    
    return
    
  end subroutine read_cell

  subroutine write_cell()

    integer :: ni

    write (stdout,'(a)') '%BLOCK LATTICE_CART'
    write (stdout,'(3f20.14)') lattice_car
    write (stdout,'(a)') '%ENDBLOCK LATTICE_CART'
    write (stdout,*)
    write (stdout,'(a)') '%BLOCK POSITIONS_FRAC'
    if(any(ion_spin.ne.0.0_dp)) then
       do ni=1,num_ions
          write (stdout,'(a,3f20.14," SPIN=",f7.3)') trim(adjustl(species_names(ion_index(ni)))),ion_fractional(:,ni),ion_spin(ni)
       end do
    else
       do ni=1,num_ions
          write (stdout,'(a,3f20.14)') trim(adjustl(species_names(ion_index(ni)))),ion_fractional(:,ni)
       end do
    end if

    write (stdout,'(a)') '%ENDBLOCK POSITIONS_FRAC'

  end subroutine write_cell

  ! ** SHLX res (result) file format
  
  subroutine read_res()

    integer :: n,i,ns,ni
        
    if(.not.(buff_words(1,1).eq.'TITL')) stop 'cabal: first line of res/shx data should start with TITL'
    if(.not.(buff_words(1,num_lines).eq.'END')) stop 'cabal: last line of res/shx data should start be END'

    do n=1,num_lines
       if(buff_words(1,n).eq.'CELL') then
          do i=1,6
             read(buff_words(i+2,n),*) lattice_abc(i)
          end do
       end if
       if(buff_words(1,n).eq.'SFAC') then
          num_species=num_words(n)-1
          do ns=1,num_species
             read(buff_words(ns+1,n),*) species_names(ns)
          end do
          exit
       end if
    end do

    num_ions=num_lines-n-1

    call init_ions()

    do ni=1,num_ions
       ion_names(ni)=trim(buff_words(1,n+ni))
       read(buff_words(2,n+ni),*) ion_index(ni)
       read(buff_words(3:5,n+ni),*) ion_fractional(:,ni)
       read(buff_words(6,n+ni),*) ion_occ(ni)
       if(num_words(n+ni).gt.6) read(buff_words(7,n+ni),*) ion_spin(ni)
    end do
    
    return
    
  end subroutine read_res

  subroutine write_res()

    integer :: n

    character(len=40) :: fmt,ctemp,ctemp2
    
    if((buff_words(1,1).eq.'TITL').and.(buff_words(10,1).eq.'n')) then
       write(stdout,'(a)') trim(adjustl(buff(1)))
    else
       write(ctemp,'(f15.6)') lattice_volume
       write(ctemp2,*) num_ions
       write (stdout,'(a,a,a,a,f4.1,a,f4.1,2i2,1x,a,a)') 'TITL cabal-',trim(inext),'2',&
            trim(outext),0.0_dp,' '//trim(adjustl(ctemp)),0.0_dp,0,0,trim(adjustl(ctemp2)),' (n/a) n - 1'
    end if
!!$    write (stdout,'(a)')        'REM'
!!$    write (stdout,'(a,a)')      'REM Converted by cabal ',date()
!!$    write (stdout,'(a)')        'REM'
    write (stdout,'(a,6f11.5)') 'CELL 1.54180',lattice_abc
    write (stdout,'(a)') 'LATT -1'
    write (ctemp,*) num_species
    write (fmt,*) "(a,"//trim(adjustl(ctemp))//"a3)"
    write (stdout,trim(adjustl(fmt))) 'SFAC ',species_names(1:num_species)
    if(any(ion_spin.ne.0.0_dp)) then
       do n=1,num_ions
          write (stdout,'(a4,i4,3f17.13,f4.1,f7.2)') &
               species_names(ion_index(n)),ion_index(n),ion_fractional(:,n),ion_occ(n),ion_spin(n)
       end do
    else
       do n=1,num_ions
          write (stdout,'(a4,i4,3f17.13,f4.1)') &
               species_names(ion_index(n)),ion_index(n),ion_fractional(:,n),ion_occ(n)
       end do
    end if
    write (stdout,'(a)') 'END'
    
  end subroutine write_res

  ! ** XYZ file format

  subroutine read_xyz()

    integer :: ni
    real(kind=dp) :: lat(3),com(3),moi(3,3),eig(3),dii
    
    read(buff_words(1,1),*) num_ions

    call init_ions()

    do ni=1,num_ions
       ion_names(ni)=trim(buff_words(1,ni+2))
       read(buff_words(2:4,ni+2),*) ion_absolute(:,ni)
    end do

    com=sum(ion_absolute,2)/real(num_ions)

    moi=0.0_dp
    do ni=1,num_ions
       dii=dot_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
       moi(1,1)=moi(1,1)+dii
       moi(2,2)=moi(2,2)+dii
       moi(3,3)=moi(3,3)+dii
       moi=moi-outer_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
    end do
    
    call diag(moi,eig)
    
    do ni=1,num_ions
       ion_absolute(:,ni)=matmul(moi,ion_absolute(:,ni)-com)
    end do
    
    lat(1)=maxval(ion_absolute(1,:))-minval(ion_absolute(1,:))+pad
    lat(2)=maxval(ion_absolute(2,:))-minval(ion_absolute(2,:))+pad
    lat(3)=maxval(ion_absolute(3,:))-minval(ion_absolute(3,:))+pad
    
    lattice_abc(1:6) =(/lat(1),lat(2),lat(3),90.0_dp,90.0_dp,90.0_dp/)

    do ni=1,num_ions
       ion_absolute(:,ni)=ion_absolute(:,ni)+lat/2.0_dp
    end do
    
  end subroutine read_xyz
  
  subroutine write_xyz()

    integer :: ni
    character(len=10) :: ctemp
    real(kind=dp) :: com(3),moi(3,3),eig(3),dii

    com=sum(ion_absolute,2)/real(num_ions)

    moi=0.0_dp
    do ni=1,num_ions
       dii=dot_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
       moi(1,1)=moi(1,1)+dii
       moi(2,2)=moi(2,2)+dii
       moi(3,3)=moi(3,3)+dii
       moi=moi-outer_product(ion_absolute(:,ni)-com,ion_absolute(:,ni)-com)
    end do

!!$    write (stderr,'(3f10.5)') moi
!!$    write (stderr,*)

    call diag(moi,eig)

!!$    write (stderr,'(3f10.5)') eig
!!$    write (stderr,*)
!!$    
!!$    write (stderr,'(3f10.5)') moi
!!$    write (stderr,*)

    moi=inv(moi)

    do ni=1,num_ions
       ion_absolute(:,ni)=matmul(moi,ion_absolute(:,ni)-com)
    end do
    
    write (ctemp,'(i10)') num_ions
    write (stdout,'(a/)') trim(adjustl(ctemp))
    do ni=1,num_ions
       write (stdout,'(a,3f20.13)') trim(adjustl(species_names(ion_index(ni)))),ion_absolute(:,ni)
    end do
    
  end subroutine write_xyz

  ! ** Extended XYZ file format

  subroutine read_xyze()
    
    integer :: ni,indx1,indx2
    real(kind=dp) :: lat(9),com(3)
    
    read(buff_words(1,1),*) num_ions

    indx1=index(buff(2),'"')
    indx2=index(buff(2),'"',back=.true.)
    
    read(buff(2)(indx1+1:indx2-1),*) lat

    lattice_car=reshape(lat,(/3,3/))
    
    call init_ions()

    com=0.0_dp
    do ni=1,num_ions
       ion_names(ni)=trim(buff_words(1,ni+2))
       read(buff_words(2:4,ni+2),*) ion_absolute(:,ni)
       com=com+ion_absolute(:,ni)/real(num_ions,dp)
    end do

    do ni=1,num_ions
       ion_absolute(:,ni)=ion_absolute(:,ni)-com(:)
    end do
    
  end subroutine read_xyze
  
  subroutine write_xyze()

    integer :: ni
    character(len=10) :: ctemp
    
    write (ctemp,'(i10)') num_ions
    write (stdout,'(a)') trim(adjustl(ctemp))
    write (stdout,'(a,9f20.6,a)') 'Lattice="',reshape(lattice_car,(/9/)),'" Properties=species:S:1:pos:R:3'
    do ni=1,num_ions
       write (stdout,'(a,3f20.13)') trim(adjustl(species_names(ion_index(ni)))),ion_absolute(:,ni)
    end do
    
  end subroutine write_xyze

  ! ** gulp file format

  subroutine read_gulp()
    stop 'read_gulp: not implemented'
  end subroutine read_gulp

  subroutine write_gulp()

    integer :: ni    

    write (stdout,'(a/a/a/a/a)') 'opti prop','title','Output','end','cell'
    write (stdout,'(6f8.3)') lattice_abc
    write (stdout,'(a)') 'frac'
    do ni=1,num_ions
       write (stdout,'(a,a6,3f18.13)') trim(adjustl(species_names(ion_index(ni)))),'core',ion_fractional(:,ni)
    end do
    
  end subroutine write_gulp

  ! ** xtl file format

  subroutine read_xtl()

    integer :: n,nstart,nfinish,ni

    do n=1,num_lines
       if(up(buff_words(1,n)).eq.'CELL') then
          read(buff_words(1:6,n+1),*) lattice_abc
          exit
       end if
    end do

    nstart=0
    nfinish=0
    do n=1,num_lines
       if(up(buff_words(1,n)).eq.'EOF') then
          nfinish=n
          exit
       end if
       if((up(buff_words(1,n)).eq.'ATOMS').and.((up(buff_words(1,n+1)).eq.'NAME'))) nstart=n+1
    end do

    if(nstart.lt.nfinish) then
       if(nfinish-nstart.lt.2) stop 'cabal: problem reading positions in xtl data'
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          read(buff_words(1:4,nstart+ni),*) ion_names(ni),ion_fractional(1:3,ni)
          ion_names(ni)=trim(strip(ion_names(ni)))
       end do
       
    end if
    
  end subroutine read_xtl

  subroutine write_xtl()

    integer :: n
    
    write (stdout,'(4a)') 'TITLE cabal-',trim(inext),'2',trim(outext)
    write (stdout,'(a/,6f11.5)') 'CELL',lattice_abc
    write (stdout,'(a/a)') 'ATOMS','NAME         X           Y           Z'
    
    do n=1,num_ions
       write (stdout,'(a4,3f17.13)') &
            species_names(ion_index(n)),ion_fractional(:,n)
    end do
    write (stdout,'(a)') 'EOF'
       
  end subroutine write_xtl

  ! ** cif file format

  subroutine read_cif()

    stop 'read_cif: not implemented'
    
  end subroutine read_cif

  subroutine write_cif()

    integer :: ni

    write (stdout,'(a/)')       "data_cif"
    write (stdout,'(a/)')       "_audit_creation_method             'generated by cabal'"
    write (stdout,'(a)')        "_symmetry_space_group_name_H-M     'P1'"
    write (stdout,'(a)')        "_symmetry_Int_Tables_number        1"
    write (stdout,'(a/)')       "_symmetry_cell_setting             triclinic"
    write (stdout,'(a,f10.4)')  "_cell_length_a                     ",lattice_abc(1) 
    write (stdout,'(a,f10.4)')  "_cell_length_b                     ",lattice_abc(2) 
    write (stdout,'(a,f10.4)')  "_cell_length_c                     ",lattice_abc(3) 
    write (stdout,'(a,f10.4)')  "_cell_angle_alpha                  ",lattice_abc(4) 
    write (stdout,'(a,f10.4)')  "_cell_angle_beta                   ",lattice_abc(5) 
    write (stdout,'(a,f10.4)')  "_cell_angle_gamma                  ",lattice_abc(6) 
    write (stdout,'(a)')        "loop_"
    write (stdout,'(a)')        "_atom_site_label"
    write (stdout,'(a)')        "_atom_site_fract_x"
    write (stdout,'(a)')        "_atom_site_fract_y"
    write (stdout,'(a)')        "_atom_site_fract_z"
    write (stdout,'(a)')        "_atom_site_occupancy"
    do ni=1,num_ions
       write (stdout,'(a5,4f9.5)') trim(species_names(ion_index(ni))),ion_fractional(:,ni),ion_occ(ni)
    end do
    
  end subroutine write_cif

  ! ** psi4 file format

  subroutine read_psi4()

    real(kind=dp) :: com(3)
    
    integer :: ni,nstart,nfinish,n
    
    nstart=0
    nfinish=0
    do n=1,num_lines
       if(index(buff(n),'Cartesian Geometry (in Angstrom)').gt.0) nstart=n
       if(index(buff(n),'Saving final (previous) structure.').gt.0) nfinish=n
    end do

    if(nfinish-nstart.lt.2) then
       stop 'cabal: problem reading atomic positions from psi4 data'
    else
       num_ions=nfinish-nstart-1
       call init_ions()
       do ni=1,num_ions
          ion_names(ni)=trim(buff_words(1,nstart+ni))
          read(buff_words(2:4,nstart+ni),*) ion_absolute(1:3,ni)
       end do
       
    end if

    com=sum(ion_absolute,2)/real(num_ions,dp)
    
    lattice_car=ident*maxval(abs(ion_absolute))*8.0_dp

    do ni=1,num_ions
       ion_absolute(:,ni)=ion_absolute(:,ni)-com+sum(lattice_car,2)/2.0_dp
    end do
    
  end subroutine read_psi4

  subroutine write_psi4()

    integer :: ni

    write (stdout,'(a)') 'molecule {'
    do ni=1,num_ions
       write (stdout,'(a,3f18.13)') trim(species_names(ion_index(ni))),ion_absolute(:,ni)
    end do
    write (stdout,'(a)') '}'
    
  end subroutine write_psi4

  ! ** DL_POLY config file format

  subroutine read_config()

    integer :: i,n,indx

    character(len=10) :: ctemp
    
    read(buff_words(3,2),*) num_ions

    do i=1,3
       read(buff(2+i),*) lattice_car(:,i)
    end do
    
    call init_ions()
    
    do n=1,num_ions
       ctemp=buff_words(1,4*n+2)
       indx=index(ctemp,'_')
       if(indx.gt.1) ctemp=ctemp(:indx-1)
       ctemp(1:1)=up(ctemp(1:1))
       ion_names(n)=ctemp
       read(buff(4*n+3),*) ion_absolute(:,n)
    end do
    
  end subroutine read_config

  subroutine write_config()

    stop 'write_config not implemented'
    
  end subroutine write_config

  ! ** VASP POSCAR format

  subroutine read_poscar()

    integer :: i,j,n,indx,nisp

    real(kind=dp) :: scale

    character(len=10) :: ctemp
    

    read(buff(2),*) scale
    do i=1,3
       read(buff(2+i),*) lattice_car(:,i)
    end do
    lattice_car=lattice_car*scale
    
    num_species=num_words(6)

    num_ions=0
    do i=1,num_species
       read(buff_words(i,7),*) n
       num_ions=num_ions+n
    end do
    
    call init_ions()

    n=0
    do i=1,num_species
       read(buff_words(i,7),*) nisp
       do j=1,nisp
          n=n+1
          ion_names(n)=buff_words(i,6)
       end do
    end do

    if(trim(adjustl(buff(8))).ne."Direct") stop 'cabal: read_poscar - direct only'
    
    do i=1,num_ions
       read(buff(8+i),*) ion_fractional(:,i)
    end do
           
  end subroutine read_poscar

  subroutine write_poscar()

    integer :: n,ns,num_type(num_species)
    
    write (stdout,'(a)') 'cabal: POSCAR file converted at '//trim(date())
    write (stdout,*) 1.0_dp
    write (stdout,'(3f15.10)') lattice_car
    write (stdout,*) species_names(1:num_species)
    do n=1,num_species
       num_type(n)=count(ion_names.eq.species_names(n))
    end do
    write (stdout,*) num_type
    write (stdout,'(a)') ' Direct'
    do ns=1,num_species
       do n=1,num_ions
          if(ion_names(n).eq.species_names(ns)) write (stdout,'(3f15.10)') ion_fractional(:,n)
       end do
    end do
    
  end subroutine write_poscar
  
  !* =======================================================================

  subroutine init_ions()

    if(allocated(ion_names)) stop 'cabal: multiple ion definitions'
    
    allocate(ion_names(num_ions),ion_index(num_ions),ion_occ(num_ions),ion_spin(num_ions))
    allocate(ion_fractional(3,num_ions),ion_absolute(3,num_ions))

    ion_names=''
    ion_index=0
    ion_occ=1.0_dp
    ion_spin=0.0_dp
    ion_fractional=0.0_dp
    ion_absolute=0.0_dp
    
  end subroutine init_ions
  
  subroutine car2abc(car,abc)

    real(kind=dp), intent(in)  :: car(3,3)
    real(kind=dp), intent(out) :: abc(6)

    abc(1) = sqrt(car(1,1)**2+car(2,1)**2+car(3,1)**2)
    abc(2) = sqrt(car(1,2)**2+car(2,2)**2+car(3,2)**2)
    abc(3) = sqrt(car(1,3)**2+car(2,3)**2+car(3,3)**2)
    abc(4) = acos(dot_product(car(:,2),car(:,3))/abc(2)/abc(3))/dgrd
    abc(5) = acos(dot_product(car(:,1),car(:,3))/abc(1)/abc(3))/dgrd
    abc(6) = acos(dot_product(car(:,1),car(:,2))/abc(1)/abc(2))/dgrd

  end subroutine car2abc

  subroutine abc2car(abc,car)

    real(kind=dp), intent(in) :: abc(6)
    real(kind=dp), intent(out):: car(3,3)
    
    car(:,1) = (/abc(1),0.0_dp,0.0_dp/)
    car(:,2) = (/abc(2)*cos(dgrd*abc(6)),abc(2)*sin(dgrd*abc(6)),0.0_dp/)
    car(1,3) = abc(3)*cos(dgrd*abc(5))
    car(2,3) = abc(3)*(cos(dgrd*abc(4))-cos(dgrd*abc(5))*cos(dgrd*abc(6)))/sin(dgrd*abc(6))
    car(3,3) = sqrt(abc(3)**2-car(1,3)**2-car(2,3)**2)
        
  end subroutine abc2car

  subroutine car2rec(car,rec)

    real(kind=dp), intent(in)  :: car(3,3)
    real(kind=dp), intent(out) :: rec(3,3)

    real(kind=dp) :: volume

    volume = car(1,1)*(car(2,2)*car(3,3)-car(3,2)*car(2,3))+&
         car(2,1)*(car(3,2)*car(1,3)-car(1,2)*car(3,3))+&
         car(3,1)*(car(1,2)*car(2,3)-car(2,2)*car(1,3))

    if(abs(volume).lt.epsilon(1.0_dp)) stop 'cabal: zero volume cell detected'
    
    ! ** Calculate the reciprocal lattice vectors

    rec(1,1)=car(2,2)*car(3,3)-car(3,2)*car(2,3)
    rec(2,1)=car(2,3)*car(3,1)-car(3,3)*car(2,1)
    rec(3,1)=car(2,1)*car(3,2)-car(3,1)*car(2,2)
    rec(1,2)=car(3,2)*car(1,3)-car(1,2)*car(3,3)
    rec(2,2)=car(3,3)*car(1,1)-car(1,3)*car(3,1)
    rec(3,2)=car(3,1)*car(1,2)-car(1,1)*car(3,2)
    rec(1,3)=car(1,2)*car(2,3)-car(2,2)*car(1,3)
    rec(2,3)=car(1,3)*car(2,1)-car(2,3)*car(1,1)
    rec(3,3)=car(1,1)*car(2,2)-car(2,1)*car(1,2)

    rec(:,:)=rec(:,:)/volume

  end subroutine car2rec 
  
  function car_volume(car)

    real(kind=dp), dimension(3,3), intent(in) :: car

    real(kind=dp) :: car_volume
    
    car_volume = car(1,1)*(car(2,2)*car(3,3)-car(3,2)*car(2,3))+&
         car(2,1)*(car(3,2)*car(1,3)-car(1,2)*car(3,3))+&
         car(3,1)*(car(1,2)*car(2,3)-car(2,2)*car(1,3))
    
    if(abs(car_volume).lt.epsilon(1.0_dp)) stop 'cabal: zero volume cell detected'

  end function car_volume

  function up(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=240)            :: up

    integer :: i,ic,ia,iz,ishift

    ! ** Initialise system dependant bounds in collating sequence
    
    ia     = ichar('a')
    iz     = ichar('z') 
    ishift = ichar('A')-ia

   ! ** Raise all lowercase to upper case

    up=string

    do i=1,200
       ic = ichar(up(i:i))
       if((ic.ge.ia).and.(ic.le.iz)) up(i:i) = char(ishift+ic) 
    end do

  end function up

  function strip(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=400)           :: strip

    integer                      :: i,ia,iz,iAA,iZZ,itest

    ! ** Initialise system dependant bounds in collating sequence
    
    ia=ichar('a')
    iz=ichar('z')
    iAA=ichar('A')
    iZZ=ichar('Z')

   ! ** Strip all characters after the first non-alphabetical character

    do i=1,len(string)
       itest=ichar(string(i:i))
       if(((itest.lt.ia).or.(itest.gt.iz)).and.((itest.lt.iAA).or.(itest.gt.iZZ))) exit
    end do
    
    strip=string(1:i-1)

  end function strip

  function detab(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=400)           :: detab

    integer                      :: i,indx

    detab=string
    do while (scan(detab,char(9)).gt.0)
       i=scan(detab,char(9))
       detab=trim(detab(1:i-1))//'      '//trim(detab(i+1:))
    end do
    
  end function detab
  
  function date()
    
    implicit none
    
    character(len=80) :: date
    character(len=80) :: time,month,zone
    character(len=2)  :: suffix  

    call date_and_time(date,time,zone)
    
    ! ** Reformat into human readable form
    
    ! * First, the time (including zone information)
    
    if(zone(2:2).eq.'0') then
       zone = '(GMT'//zone(1:1)//zone(3:3)//'.'//zone(4:4)//')'
    else
       zone = '(GMT'//zone(1:1)//zone(2:3)//'.'//zone(4:4)//')'
       end if

       time = time(1:2)//':'//time(3:4)//':'//time(5:6)//' '//trim(zone)

       ! * Now, the date

       ! Find the month

       select case(date(5:6))
       case ('01') ; month='January'
       case ('02') ; month='February'
       case ('03') ; month='March'
       case ('04') ; month='April'
       case ('05') ; month='May'
       case ('06') ; month='June'
       case ('07') ; month='July'
       case ('08') ; month='August'
       case ('09') ; month='September'
       case ('10') ; month='October'
       case ('11') ; month='November'
       case ('12') ; month='December'
       case default
          month = 'Month Unknown'
       end select

       ! Date suffix

       select case(date(7:8))
       case ('01') ; suffix='st'
       case ('02') ; suffix='nd'
       case ('03') ; suffix='rd'
       case ('21') ; suffix='st'
       case ('22') ; suffix='nd'
       case ('23') ; suffix='rd'
       case ('31') ; suffix='st'
       case default
          suffix='th'
       end select

       ! * Finally, combine into a single string

       if(date(7:7).eq.'0') then
          date = trim(time)//' '//date(8:8)//suffix//' '//&
               trim(month)//' '//date(1:4)
       else
          date = trim(time)//' '//date(7:8)//suffix//' '//&
               trim(month)//' '//date(1:4)
       end if

  end function date

  function outer_product(a,b) result(c)

    real(kind=dp), intent(in) :: a(3),b(3)

    real(kind=dp) :: c(3,3)
    
    c=spread(a,dim=2,ncopies=3)*spread(b,dim=1,ncopies=3)
   
  end function outer_product

  subroutine diag(a,b) 

    real(kind=dp), intent(inout) :: a(3,3)

    real(kind=dp), intent(out) :: b(3)

    integer :: info

    real(kind=dp) :: work(12)
    
    call dsyev('V','U',3,a,3,b,work,12,info)

    if(info.ne.0) stop 'error in diag'
    
  end subroutine diag

  function inv(A) result(Ainv)

    real(kind=dp), dimension(:,:), intent(in) :: A
    real(kind=dp), dimension(size(A,1),size(A,2)) :: Ainv

    real(kind=dp), dimension(size(A,1)) :: work 
    integer,       dimension(size(A,1)) :: ipiv
    integer :: n, info

    external dgetrf
    external dgetri

    Ainv=A
    n=size(A,1)

    call dgetrf(n,n,Ainv,n,ipiv,info)

    if (info.ne.0)  stop 'inv: matrix is numerically singular!'

    call dgetri(n,Ainv,n,ipiv,work,n,info)

    if (info.ne.0) stop 'inv: matrix inversion failed!'

  end function inv
  
end program cabal

