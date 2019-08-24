! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!==================================================================================!
!                                     cell                                         !
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
! This module reads, knows and writes the unit cell                                !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module cell

  use iso_fortran_env
  use constants
  use symmetry
  use rng

  implicit none

!!$  interface ! Needed for ifort
!!$     integer function getpid()
!!$     end function getpid
!!$  end interface

  private

  public :: read_cell
  public :: reset_cellbuff
  public :: get_seed
  public :: comments_cell
  public :: blank_cell
  public :: expand_cell
  public :: by_species_cell
  public :: adjust_cell
  public :: init_cell
  public :: check_chem
  public :: write_cell
  public :: gen_sets
  public :: select_symmetry
  public :: symm_cell
  public :: volume_cell
  public :: compact_cell
  public :: update_cell
  public :: up
  public :: strip
  public :: heap_sort_index

  ! ** Public data defining the system
  
  real(kind=dp), public :: lattice_abc(6)
  real(kind=dp), public :: lattice_car(3,3)
  real(kind=dp), public :: lattice_car_orig(3,3)
  real(kind=dp), public :: lattice_rec(3,3)
  real(kind=dp), public :: lattice_rec_orig(3,3)
  real(kind=dp), public :: setting_perm(3,3)
  real(kind=dp), public :: volume

  integer,       public                                :: num_ions
  integer,       public                                :: num_ions_total
  integer,       public                                :: num_spec
  integer,       public                                :: num_sets
  integer,       public                                :: num_cells
  integer,       public                                :: num_form
  integer,       public                                :: num_vacancies
  integer,       public                                :: num_pairs
  integer,       public                                :: num_fails
  integer,       public                                :: num_atom
  integer,       public                                :: adjgen
  integer,       public                                :: nudgen=0
  integer,       public                                :: nrank
  integer,       public                                :: focus_species
  integer,       public                                :: pushmax
  integer,       public, dimension(:,:),   allocatable :: ion_coord
  integer,       public, dimension(:,:,:), allocatable :: ion_set_coord
  integer,       public, dimension(:),     allocatable :: ion_ion_set
  integer,       public, dimension(:),     allocatable :: ninset
  integer,       public, dimension(:),     allocatable :: indx_nearestn
  integer,       public, dimension(:),     allocatable :: ion_mult
  integer,       public, dimension(:),     allocatable :: ion_fails
  integer,       public, dimension(:),     allocatable :: ion_species
  integer,       public, dimension(:,:),   allocatable :: ion_equiv
  integer,       public, dimension(:,:),   allocatable :: ion_set_mult
  integer,       public, dimension(3,3)                :: supercell

  real(kind=dp), public, dimension(:,:),   allocatable :: ion_positions,ion_new_positions,ion_positions_temp
  real(kind=dp), public, dimension(:,:),   allocatable :: ion_temp_positions,ion_hold_positions,ion_push
  real(kind=dp), public, dimension(:,:,:), allocatable :: ion_set_positions,ion_ion_vec
  real(kind=dp), public, dimension(:,:),   allocatable :: nearestn
  real(kind=dp), public, dimension(:),     allocatable :: posamp,posamp_temp,minamp,angamp
  real(kind=dp), public, dimension(:),     allocatable :: xamp,yamp,zamp,xamp_temp,yamp_temp,zamp_temp
  real(kind=dp), public, dimension(:),     allocatable :: ion_rad
  real(kind=dp), public, dimension(:),     allocatable :: ion_min
  real(kind=dp), public, dimension(:),     allocatable :: ion_occ
  real(kind=dp), public, dimension(:),     allocatable :: posamp_set,minamp_set,angamp_set
  real(kind=dp), public, dimension(:),     allocatable :: xamp_set,yamp_set,zamp_set
  real(kind=dp), public, dimension(:),     allocatable :: spinvec,afmvec
  real(kind=dp), public, dimension(:,:),   allocatable :: ion_set_rad
  real(kind=dp), public, dimension(:,:),   allocatable :: ion_set_min
  real(kind=dp), public, dimension(:,:),   allocatable :: ion_set_occ
  real(kind=dp), public, dimension(:,:),   allocatable :: ion_set_centre
  real(kind=dp), public, dimension(:),     allocatable :: ion_set_sphere
  
  character(len=16), public, dimension(:),     allocatable :: species_names
  character(len=16), public, dimension(:),     allocatable :: ion_names
  character(len=16), public, dimension(:),     allocatable :: ion_set
  character(len=16), public, dimension(:,:),   allocatable :: ion_nn
  character(len=16), public, dimension(:,:,:), allocatable :: ion_set_nn
  character(len=16), public, dimension(:,:),   allocatable :: ion_set_names

  character(len=16),  public, dimension(2,100)              :: pair_names
  character(len=40),  public                                :: perm_names
  character(len=40),  public                                :: spinlist
  character(len=4),   public                                :: crystalsystem
  character(len=20),  public                                :: csystem
  character(len=400), public                                :: seed_string=''

  logical,           public, dimension(:),     allocatable :: ion_fix,ion_fix_temp
  logical,           public, dimension(:),     allocatable :: ion_bfix,ion_bfix_temp
  logical,           public, dimension(:),     allocatable :: ion_perm
  logical,           public, dimension(:),     allocatable :: ion_adatom,ion_adatom_temp
  logical,           public, dimension(:),     allocatable :: ion_athole,ion_athole_temp
  logical,           public, dimension(:),     allocatable :: ion_moved
  logical,           public, dimension(:,:),   allocatable :: ion_set_fix
  logical,           public, dimension(:,:),   allocatable :: ion_set_bfix
  
  real(kind=dp), public :: cons
  real(kind=dp), public :: acons
  real(kind=dp), public :: minsep,minmin,minmax
  real(kind=dp), public :: targvol
  real(kind=dp), public :: cellconvec(6)  
  real(kind=dp), public :: cellamp
  real(kind=dp), public :: breakamp
  real(kind=dp), public :: rash_angamp
  real(kind=dp), public :: rash_posamp
  real(kind=dp), public :: permfrac
  real(kind=dp), public :: minbangle
  real(kind=dp), public :: maxbangle
  real(kind=dp), public :: slack=0.0_dp, slack_step=0.01_dp
  real(kind=dp), public :: maxtime
  real(kind=dp), public :: overlap
  real(kind=dp), public :: width
  real(kind=dp), public :: shift
  real(kind=dp), public :: vacuum
  real(kind=dp), public :: radius
  real(kind=dp), public :: coreradius
  real(kind=dp), public :: ellipsvec(3)
  real(kind=dp), public :: three
  real(kind=dp), public :: pushstep  
  real(kind=dp), public :: spin,spinmod
  real(kind=dp), public :: len_min_use,len_max_use
  real(kind=dp), public :: lookup(3,27),vtemp(3)
  real(kind=dp), public :: pair_minsep(100)
  real(kind=dp), public :: hole_rad,hole_pos(3)

  logical, public :: symmgen
  logical, public :: fix_cell
  logical, public :: fix_all_cell
  logical, public :: fix_cell_vol
  logical, public :: fix_caxis
  logical, public :: fix_abaxis
  logical, public :: con_cell
  logical, public :: internalsymm
  logical, public :: ionabs
  logical, public :: tight
  logical, public :: fix_vol
  logical, public :: free
  logical, public :: molecules
  logical, public :: permute
  logical, public :: supcell
  logical, public :: findsup
  logical, public :: compact
  logical, public :: slab
  logical, public :: surface
  logical, public :: sphere
  logical, public :: core
  logical, public :: ellipsoid
  logical, public :: cylinder
  logical, public :: cluster
  logical, public :: symmapprox
  logical, public :: autoslack=.false.
  logical, public :: push
  logical, public :: translate=.true.
  logical, public :: octet
  logical, public :: dorash
  logical, public :: celladapt
  logical, public :: cell_consistent=.false.
  logical, public :: havespin
  logical, public :: remove
  logical, public :: flip

  integer, public :: rejected(100)=0

  !------------------------------------------------------------------------!
  ! **                       Generation parameters                      ** !

  real(kind=dp), public, parameter :: ang_min=0.0_dp, ang_max=180.0_dp
  real(kind=dp), public, parameter :: len_min=0.0_dp, len_max=1.0_dp
  real(kind=dp), public, parameter :: vol_min=0.95_dp, vol_max=1.05_dp

  !------------------------------------------------------------------------!

  ! ** Private IO data
  
  integer, parameter :: unit_minsep=31
  integer                                       :: norecords,nrecords,orecords,icluster
  character(len=400), dimension(10000)          :: origbuff
  character(len=400), dimension(:), allocatable :: cellbuff,cellbuffup,cellbuffcopy,cellbuffsave

  ! ** Private structure data

  real(kind=dp) :: scale_vol

  logical                            :: symmorphic,notsymm,randsymm
  logical, dimension(:), allocatable :: flag

  integer                            :: num_form_min,adjgen_min
  integer                            :: num_form_max,adjgen_max

  ! ** Counters and temp variables

  integer :: i,j,k,m,n,nn,ni,nj,ns,nops,num_pot,num_hub,sel(230),nstart,indx,num,den
  integer :: symmno,symmin,symmax,nsel,itmp(9),lowso,highso,nummin,nummax
  integer :: num_ions_orig

  integer, dimension(:), allocatable :: mult

  integer, dimension(6) :: cellcon=(/1,2,3,4,5,6/)
  
  integer, parameter, dimension(10) :: numsops=(/1,2,3,4,6,8,12,16,24,48/)
  integer, parameter, dimension(13) :: numsclops=(/1,2,3,5,4,6,7,8,9,10,11,12,24/)

  real(kind=dp) :: ang,rn(3),temp,v(3),vmin,vmax,targvolauto
  
  character(len=400)                            :: ctemp,ctemp2,ctemp4,ctemp6,composition
  character(len=40)                             :: ctemp5,label,label2,tempfile
  character(len=10)                             :: ctemp3,vac_name,carray(100)
  character(len=400)                            :: kpspacing(20)
  character(len=400)                            :: symmetrytol
  character(len=400)                            :: autotargvol,autominsep
  character(len=400), dimension(:), allocatable :: spec_pot,hub_u
  character(len=16),  dimension(:), allocatable :: ion_names_temp, ion_set_temp, set_label

  logical :: inblock,open,auto=.false.

contains

  subroutine read_cell

    ! ** Determine the length of the file, and read it in from STDIN

    norecords=0
    do 
       read(stdin,'(a)',end=997) origbuff(norecords+1)
       norecords=norecords+1
    end do
997 return

  end subroutine read_cell

  subroutine reset_cellbuff()

    if(allocated(cellbuff)) deallocate(cellbuff,cellbuffup,cellbuffcopy)
    allocate(cellbuff(norecords),cellbuffup(norecords),cellbuffcopy(norecords))
 
    ! ** Read in the cell file to cellbuff

    cellbuff(1:norecords)=origbuff(1:norecords)
    nrecords=norecords

    do n=1,nrecords
       cellbuff(n) = trim(adjustl(cellbuff(n)))
       cellbuffup(n) = up(cellbuff(n))
    end do

  end subroutine reset_cellbuff

  subroutine comments_cell()

    ! ** Eliminate portions of the input commented out using '##'

    i=0 ; n=1
    do

       if(n+i.gt.nrecords) then
          nrecords=nrecords-i
          exit
       end if

       indx=index(cellbuff(n+i),'##')
       
       if(indx.eq.1) then
          i=i+1
       else
          if(indx.gt.1) then 
             cellbuff(n)=cellbuff(n+i)(:indx-1)
          else 
             cellbuff(n)=cellbuff(n+i)
          end if
          n=n+1
       end if

    end do    

  end subroutine comments_cell

  subroutine blank_cell()

    ! ** Eliminate blank lines

    i=0 ; n=1
    do

       if(n+i.gt.nrecords) then
          nrecords=nrecords-i
          exit
       end if
       
       if(len_trim(cellbuff(n+i)).eq.0) then
          i=i+1
       else
          cellbuff(n)=cellbuff(n+i)
          n=n+1
       end if

    end do

  end subroutine blank_cell
  
  subroutine expand_cell()

    ! ** Expand {a,b,c}

    do n=1,nrecords
       i=index(cellbuff(n),"{")
       j=index(cellbuff(n),"}")
       if(i.gt.0) then
          ctemp=cellbuff(n)(i+1:j-1) 
          m=0
          do
             m=m+1
             nn=index(ctemp,",")
             if(nn>0) then
                carray(m)=ctemp(1:nn-1)
                ctemp=ctemp(nn+1:)
             else
                carray(m)=trim(ctemp)
                exit
             end if
          end do

          rn = random_triple()
          nn = 1+int(rn(1)*m)

          cellbuff(n)=cellbuff(n)(1:i-1)//trim(carray(nn))//cellbuff(n)(j+1:)

       end if

    end do
    
    return

  end subroutine expand_cell

  subroutine get_seed()

    ! ** Note - because get_seed is called before comments_cell, #SEED cannot be commented out
    
    seed_string = ''
    do n=1,norecords
       ctemp=origbuff(n)
       if(index(up(ctemp),'#SEED=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#SEED=')+6:))
          read(ctemp2,'(a)') seed_string
       end if
    end do

  end subroutine get_seed
  
  subroutine by_species_cell()

    integer :: na,nspec,indx,natom,natom_min,natom_max
    character(len=400), allocatable, dimension(:) :: cspec,cmod
    integer, allocatable, dimension(:) :: atom2spec

    natom=0
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#NATOM=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#NATOM=')+7:))
          if(index(ctemp2,'-')>0) then
             read(ctemp2(:index(ctemp2,'-')-1),*) natom_min
             read(ctemp2(index(ctemp2,'-')+1:),*) natom_max
             rn = random_triple()       
             natom = natom_min+int(rn(1)*(natom_max-natom_min+1))
          else
             read(ctemp2,*) natom
          endif
       end if
    end do

    num_atom=natom
    
    nspec=0
    do n=1,nrecords

       if(index(up(cellbuff(n)),'#SPECIES=')>0) then

          do m=1,nrecords
             if(index(up(cellbuff(m)),'%BLOCK POSITIONS_')>0) then
                stop 'Both SPECIES and POSITIONS_ block defined'
             end if
          end do

          indx=index(up(cellbuff(n)),'#SPECIES=')
          ctemp=adjustl(cellbuff(n)(indx+9:))
          nspec=1
          do m=1,len(ctemp)
             if(ctemp(m:m).eq.',') nspec=nspec+1
          end do

          allocate(cspec(nspec),cmod(nspec))

          do ns=1,nspec
             indx=index(ctemp,',')
             if(indx.gt.0) then
                ctemp2=ctemp(:indx-1)
             else
                ctemp2=ctemp
             end if
             ctemp=ctemp(indx+1:)

             indx=index(ctemp2,'%')
             if(indx.gt.0) then
                cspec(ns)=ctemp2(:indx-1)
                cmod(ns)=ctemp2(indx+1:)
             else
                cspec(ns)=ctemp2
                cmod(ns)=''
             end if

          end do

       end if

    end do

    scale_vol=1.0_dp

    if(nspec.gt.0) then

       if(natom.gt.0) then

          if(allocated(cellbuffcopy)) deallocate(cellbuffcopy)
          allocate(cellbuffcopy(nrecords))
          cellbuffcopy=cellbuff(1:nrecords)
          deallocate(cellbuff)
          nrecords=nrecords+natom+2
          allocate(cellbuff(nrecords))
          cellbuff=''
          cellbuff(1:size(cellbuffcopy))=cellbuffcopy(:)

          allocate(atom2spec(natom))

          atom2spec=1

          do ns=2,nspec
             rn = random_triple()
             na = int(rn(1)*(natom+1))
             if(na.lt.natom) atom2spec(na+1:)=atom2spec(na+1:)+1
          end do

          n=size(cellbuffcopy)
          cellbuff(n+1)='%BLOCK POSITIONS_FRAC'
          do na=1,natom
             write(ctemp2,*) na
             cellbuff(n+1+na)=trim(cspec(atom2spec(na)))//' 0.0 0.0 0.0 # '//trim(adjustl(ctemp2))
          end do
          cellbuff(n+natom+2)='%ENDBLOCK POSITIONS_FRAC'

          scale_vol=real(natom,dp)

       else

          if(allocated(cellbuffcopy)) deallocate(cellbuffcopy)
          allocate(cellbuffcopy(nrecords))
          cellbuffcopy=cellbuff(1:nrecords)
          deallocate(cellbuff)
          nrecords=nrecords+nspec+2
          allocate(cellbuff(nrecords))
          cellbuff=''
          cellbuff(1:size(cellbuffcopy))=cellbuffcopy(:)

          n=size(cellbuffcopy)
          cellbuff(n+1)='%BLOCK POSITIONS_FRAC'
          do ns=1,nspec
             write(ctemp2,*) ns
             cellbuff(n+1+ns)=trim(cspec(ns))//' 0.0 0.0 0.0 # '//trim(adjustl(ctemp2))//' % '//trim(cmod(ns))
          end do
          cellbuff(n+nspec+2)='%ENDBLOCK POSITIONS_FRAC'

       end if

       deallocate(cspec)

    end if

  end subroutine by_species_cell

  subroutine adjust_cell()

    integer :: count_attempt=0
    
    ! ** Adjust the cellbuff to allow for variable numbers of sets

    if(allocated(cellbuffcopy)) deallocate(cellbuffcopy)
    if(allocated(cellbuffsave)) deallocate(cellbuffsave)
    allocate(cellbuffcopy(nrecords),cellbuffsave(nrecords))
    cellbuffsave(1:nrecords) = cellbuff(1:nrecords)
    deallocate(cellbuff)
    orecords=nrecords
    nrecords=nrecords+10000 ! *** FIX THIS ***
    allocate(cellbuff(nrecords))

996 num_ions_orig=0
    num_ions=0
    n=0
    nn=0
    cellbuffcopy=cellbuffsave

!!$    write (stderr,*) '---',orecords
!!$    do i=1,orecords
!!$       write (stderr,*) trim(cellbuffcopy(i))
!!$    end do
!!$    write (stderr,*) '---'
        
    do
       n=n+1 ; if(n.gt.orecords) exit
       nn=nn+1
       ctemp=cellbuffcopy(n)
       cellbuff(nn)=ctemp
       if(index(up(ctemp),'%BLOCK POSITIONS_')>0) then
          nstart=n
          do
             n=n+1 ; if(n.gt.orecords) exit
             ctemp=cellbuffcopy(n)
             if(index(up(ctemp),'%ENDBLOCK POSITIONS_')>0) exit
             num_ions_orig=num_ions_orig+1
             if(index(ctemp,'NUM=')>0) then
                ctemp2 = ctemp(index(ctemp,'NUM=')+4:)
                ctemp2 = ctemp2(:index(ctemp2," "))

                if(index(ctemp2,'-')>0) then
                   read(ctemp2(:index(ctemp2,'-')-1),*) nummin
                   read(ctemp2(index(ctemp2,'-')+1:),*) nummax
                else
                   read(ctemp2,*) nummin
                   nummax=nummin
                endif

                ! * Make sure all in the same unit are repeated the same number of times
                                
                rn = random_triple()
                i = nummin+int(rn(1)*(nummax-nummin+1))

                ctemp4=cellbuffcopy(n)(index(cellbuffcopy(n),'NUM='):)
                ctemp4=ctemp4(index(ctemp4,' '):)
                write(ctemp2,*) i

                cellbuffcopy(n)=&
                     trim(cellbuffcopy(n)(:index(cellbuffcopy(n),'NUM=')-1))//' NUM='//trim(adjustl(ctemp2))//trim(ctemp4)
                
                nummin=i
                nummax=i

             else

                nummin=1
                nummax=1

                ! * Check if a repeat number has already been given for this label

                if(index(cellbuffcopy(n),'#')>0) then
                   read(cellbuffcopy(n)(index(cellbuffcopy(n),'#')+1:),*) label
                else
                   label=''
                end if

                do i=nstart+1,n-1
                   if(index(cellbuffcopy(i),'#')>0) then
                      read(cellbuffcopy(i)(index(cellbuffcopy(i),'#')+1:),*) label2
                   else
                      label2=''
                   end if

                   if(label2.eq.label) then
                      ctemp4=cellbuffcopy(i)
                      if(index(ctemp4,'NUM=')>0) then
                         ctemp2 = ctemp4(index(ctemp4,'NUM=')+4:)
                         ctemp2 = ctemp2(:index(ctemp2," "))

                         if(index(ctemp2,'-')>0) then
                            read(ctemp2(:index(ctemp2,'-')-1),*) nummin
                            read(ctemp2(index(ctemp2,'-')+1:),*) nummax
                         else
                            read(ctemp2,*) nummin
                            nummax=nummin                
                         endif
                         exit
                      end if
                   end if
                end do

             end if

             if(index(ctemp,"%").eq.0) ctemp=trim(ctemp)//" %"

             rn = random_triple()
             i = nummin+int(rn(1)*(nummax-nummin+1))

             do j=1,i
                nn=nn+1
                if((index(ctemp,"#")>0).and.(index(ctemp,"%")>0).and.(i>1)) then
                   write (ctemp3,'(i5)') j
                   ctemp2=trim(ctemp(index(ctemp,"#")+1:index(ctemp,"%")-1))//"-"//trim(adjustl(ctemp3))
                   cellbuff(nn)=trim(ctemp(:index(ctemp,"#")))//trim(ctemp2)//" "//trim(ctemp(index(ctemp,"%"):))
                else
                   cellbuff(nn)=trim(ctemp)
                endif
                num_ions=num_ions+1
             end do

          end do

          exit
       end if

    end do

    cellbuff(nn+1:nn+orecords-n+1)=cellbuffcopy(n:orecords)
    cellbuff(nn+orecords-n+2:)=''
    nrecords=nn+orecords-n+1
    
    ! ** Scale the target volume by the number of actual atoms chosen

    if(num_ions > 0) then
       scale_vol = scale_vol*real(num_ions,dp)/real(num_ions_orig,dp)
    else
       count_attempt=count_attempt+1
       if(count_attempt.lt.1000) then
          goto 996 ! * Try again
       else
          stop 'adjust_cell : num_ions <= 0'
       end if
    end if
    
    do n=1,nrecords
      
       if(index(cellbuff(n),"UNTIE")>0) then
          read(cellbuff(n)(index(cellbuff(n),'#')+1:),*) label
          
          nn=0
          do j=1,nrecords
             if(index(cellbuff(j),'#')>0) then
                read(cellbuff(j)(index(cellbuff(j),'#')+1:),*) label2
                if(label.eq.label2) then
                   nn=nn+1
                   write (ctemp3,'(i5)') nn
                   ctemp=cellbuff(j)
                   m=index(ctemp,trim(label))+len_trim(label)
                   cellbuff(j)=trim(ctemp(:m))//"-"//trim(adjustl(ctemp3))//" "//trim(ctemp(m+1:))
                end if
             end if
          end do
          
       end if

    end do

    return

  end subroutine adjust_cell

  subroutine init_cell(stat)

    integer, intent(inout) :: stat

    integer :: symmstat

    real(kind=dp) :: rmin,ellipar,x0

    stat=999

    ! ** Symmorphic groups only? 

    symmorphic=any(cellbuffup.eq.'#SYMMORPHIC')

    fix_cell    = .false.
    fix_caxis   = .false.
    fix_abaxis  = .false.
    lattice_abc = 0.0_dp
    lattice_car = 0.0_dp

    ! ** Read the lattice vectors

    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK LATTICE_ABC')>0) then

          n=n+1 ; if(n.gt.nrecords) exit
          read(cellbuff(n),*) lattice_abc(1:3)
          n=n+1 ; if(n.gt.nrecords) exit
          read(cellbuff(n),*) lattice_abc(4:6)

          ang = 1.0_dp-cos(dgrd*lattice_abc(4))**2-cos(dgrd*lattice_abc(5))**2-cos(dgrd*lattice_abc(6))**2+&
               2.0_dp*cos(dgrd*lattice_abc(4))*cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6))

          volume = lattice_abc(1)*lattice_abc(2)*lattice_abc(3)*sqrt(abs(ang))

          lattice_car(:,1) = (/lattice_abc(1),0.0_dp,0.0_dp/)
          lattice_car(:,2) = (/lattice_abc(2)*cos(dgrd*lattice_abc(6)),lattice_abc(2)*sin(dgrd*lattice_abc(6)),0.0_dp/)
          lattice_car(1,3) = lattice_abc(3)*cos(dgrd*lattice_abc(5))
          lattice_car(2,3) = lattice_abc(3)*(cos(dgrd*lattice_abc(4))-cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6)))&
               /sin(dgrd*lattice_abc(6))
          lattice_car(3,3) = sqrt(lattice_abc(3)**2-lattice_car(1,3)**2-lattice_car(2,3)**2)

          n=n+1 ; if(n.gt.nrecords) exit
          ctemp=cellbuff(n)
          if(trim(up(ctemp))=='#FIX')   fix_cell=.true.
          if(trim(up(ctemp))=='#CFIX')  fix_caxis=.true.
          if(trim(up(ctemp))=='#ABFIX') fix_abaxis=.true.
          exit
       end if

    end do

    ! ** Fix the unit cell?

    symmgen=any(cellbuffup.eq."SYMMETRY_GENERATE")

    ! ** Fix the unit cell?

    fix_all_cell=any(cellbuffup.eq."FIX_ALL_CELL : TRUE")

    ! ** Fix the volume

    fix_cell_vol=any(cellbuffup.eq."FIX_VOL : TRUE")

    ! ** Cell constraints

    con_cell=.false.
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#CELLCON=')>0) then
          con_cell=.true.
          ctemp2=trim(ctemp(index(up(ctemp),'#CELLCON=')+9:))
          read(ctemp2,*) cellconvec
       end if
    end do
    
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SYSTEM=')>0) then
          if(con_cell) stop 'SYSTEM and CELLCON defined'
          con_cell=.true.
          ctemp2=trim(ctemp(index(up(ctemp),'#SYSTEM=')+8:))
          read(ctemp2,*) crystalsystem
          select case(crystalsystem)
          case('Tric')
             cellconvec=(/0,0,0,0,0,0/)
             csystem='Triclinic'
          case('Mono')
             cellconvec=(/0,0,0,90,0,90/)
             csystem='Monoclinic'
          case('Hexa')
             cellconvec=(/-1,-1,0,90,90,120/)
             csystem='Hexagonal'
          case('Rhom','Trig')
             cellconvec=(/-1,-1,-1,-1,-1,-1/)
             csystem='Rhombohedral'
          case('Orth')
             cellconvec=(/0,0,0,90,90,90/)
             csystem='Orthorhombic'
          case('Tetr')
             cellconvec=(/-1,-1,0,90,90,90/)
             csystem='Tetragonal'
          case('Cubi')
             cellconvec=(/-1,-1,-1,90,90,90/)
             csystem='Cubic'
          case default
             stop 'Crystal system not recognised'
          end select
          write (stderr,'(a,a)') 'Crystal system: ',csystem
          ! *** PERMUTE?
       end if
    end do

    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK LATTICE_CART')>0) then

          do i=1,3
             n=n+1 ; if(n.gt.nrecords) exit
             read(cellbuff(n),*) lattice_car(1:3,i)
          end do

          volume = lattice_car(1,1)*(lattice_car(2,2)*lattice_car(3,3)-lattice_car(3,2)*lattice_car(2,3))+&
               lattice_car(2,1)*(lattice_car(3,2)*lattice_car(1,3)-lattice_car(1,2)*lattice_car(3,3))+&
               lattice_car(3,1)*(lattice_car(1,2)*lattice_car(2,3)-lattice_car(2,2)*lattice_car(1,3))

          lattice_abc(1) = sqrt(lattice_car(1,1)**2+lattice_car(2,1)**2+lattice_car(3,1)**2)
          lattice_abc(2) = sqrt(lattice_car(1,2)**2+lattice_car(2,2)**2+lattice_car(3,2)**2)
          lattice_abc(3) = sqrt(lattice_car(1,3)**2+lattice_car(2,3)**2+lattice_car(3,3)**2)
          lattice_abc(4) = acos(dot_product(lattice_car(:,2),lattice_car(:,3))/lattice_abc(2)/lattice_abc(3))/dgrd
          lattice_abc(5) = acos(dot_product(lattice_car(:,1),lattice_car(:,3))/lattice_abc(1)/lattice_abc(3))/dgrd
          lattice_abc(6) = acos(dot_product(lattice_car(:,1),lattice_car(:,2))/lattice_abc(1)/lattice_abc(2))/dgrd

          n=n+1 ; if(n.gt.nrecords) exit
          ctemp=cellbuff(n)
          if(trim(up(ctemp))=='#FIX')   fix_cell   = .true.
          if(trim(up(ctemp))=='#CFIX')  fix_caxis  = .true.
          if(trim(up(ctemp))=='#ABFIX') fix_abaxis = .true.
          exit
       end if

    end do

    ! ** Prepare a supercell

    supercell = 0
    supercell(1,1) = 1
    supercell(2,2) = 1
    supercell(3,3) = 1
    supcell = .false.
    findsup = .false.
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SUPERCELL')>0) then
          supcell = .true.
          ctemp2=trim(adjustl(ctemp(index(up(ctemp),'#SUPERCELL=')+11:)))
          itmp=0
          m=0
          do while(len(trim(ctemp2))>0)
             m=m+1
             read (ctemp2(:index(ctemp2,' ')),*) itmp(m)
             ctemp2 = trim(ctemp2(index(ctemp2,' ')+1:))
          end do

          select case(m)
          case(1)
             supercell(1,1)=itmp(1)
             findsup=.true.
          case(3)
             supercell(1,1)=itmp(1)          
             supercell(2,2)=itmp(2)     
             supercell(3,3)=itmp(3)
          case(9)
             supercell(1,1)=itmp(1)          
             supercell(2,1)=itmp(2)     
             supercell(3,1)=itmp(3)
             supercell(1,2)=itmp(4)          
             supercell(2,2)=itmp(5)     
             supercell(3,2)=itmp(6)
             supercell(1,3)=itmp(7)          
             supercell(2,3)=itmp(8)     
             supercell(3,3)=itmp(9)
          case default
             write (stderr,*) 'Requested supercell not understood'
             stop
          end select

       end if
    end do

    ! ** Choose a slab geometry for the supercell? 

    slab=any(cellbuffup.eq.'#SLAB')

    ! ** Are we trying to build a suface model? 

    surface=any(cellbuffup.eq.'#SURFACE')

    if(surface) slab=.true.

    ! ** Vacuum padding

    vacuum = 0.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#VACUUM=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#VACUUM=')+8:))
          read(ctemp2,*) vacuum
       end if
    end do

    ! ** Are we trying to build a cluster? 

    icluster=1

    cluster=any(cellbuffup.eq.'#CLUSTER')

    if(cluster) icluster=-1

    ! ** Confining sphere potential?

    radius=-1.0_dp

    sphere=.false.
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SPHERE=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#SPHERE=')+8:))
          sphere=.true.
          read(ctemp2,*) radius
          if(radius.lt.0.0_dp) radius=abs(radius)*real(num_ions,dp)**(1.0_dp/3.0_dp)
          if(radius.gt.0.0_dp) then
             write (stderr,'(a,f5.2)') 'Using a confining sphere of radius ',radius
          else
             write (stderr,'(a,f5.2)') 'Using an attractive point potential with strength ',abs(radius)
          end if

          ellipsvec=1.0_dp

       end if
    end do

    core=.false.
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#CORE=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#CORE=')+6:))
          core=.true.
          read(ctemp2,*) coreradius
          if(coreradius.gt.0.0_dp) then
             write (stderr,'(a,f7.3)') 'Using a repulsive core of radius ',coreradius
          end if
       end if
    end do

    ! ** Confining ellipsoid - input values controls the aspect ration (0 -> sphere, 1 max variation)

    ellipsoid=.false.

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#ELLIPSOID=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#ELLIPSOID=')+11:))
          ellipsoid=.true.
          sphere=.true.
          read(ctemp2,*) rmin,ellipar
          if(rmin.lt.0.0_dp) then
             rmin=abs(rmin)
             radius=rmin*real(num_ions,dp)**(1.0_dp/3.0_dp)
          else
             radius=rmin
             rmin=radius/real(num_ions,dp)**(1.0_dp/3.0_dp)
          end if

          if(ellipar.lt.0.0_dp) stop 'ellipsoid out of range - 0 for spherical, larger for more variance'

          ellipsvec=random_triple()**(ellipar*1.5_dp)
          ellipsvec=ellipsvec/(ellipsvec(1)*ellipsvec(2)*ellipsvec(3))**(1.0_dp/3.0_dp)

       end if
    end do

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#PANCAKE=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#PANCAKE=')+9:))
          ellipsoid=.true.
          sphere=.true.
          read(ctemp2,*) rmin,ellipar
          if(rmin.lt.0.0_dp) then
             rmin=abs(rmin)
             radius=rmin*real(num_ions,dp)**(1.0_dp/3.0_dp)
          else
             radius=rmin
             rmin=radius/real(num_ions,dp)**(1.0_dp/3.0_dp)
          end if

          if((ellipar.lt.0.0_dp).or.(ellipar.gt.1.0_dp)) stop 'pancake out of range - 0 for flat, 1 for full variance'

          x0=min(1.0_dp,(rmin/radius/2.0_dp)**1.5_dp)
          ellipsvec(3)=(x0+ellipar*random_single()*(1.0_dp-x0))**(2.0_dp/3.0_dp)

          ellipsvec(1)=1.0_dp/sqrt(ellipsvec(3))
          ellipsvec(2)=1.0_dp/sqrt(ellipsvec(3))

       end if
    end do

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#CIGAR=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#CIGAR=')+7:))
          ellipsoid=.true.
          sphere=.true.
          read(ctemp2,*) rmin,ellipar
          if(rmin.lt.0.0_dp) then
             rmin=abs(rmin)
             radius=rmin*real(num_ions,dp)**(1.0_dp/3.0_dp)
          else
             radius=rmin
             rmin=radius/real(num_ions,dp)**(1.0_dp/3.0_dp)
          end if

          if((ellipar.lt.0.0_dp).or.(ellipar.gt.1.0_dp)) stop 'cigar out of range - 0 for needle, 1 for full variance'

          x0=min(1.0_dp,(rmin/radius*2.0_dp)**3.0_dp)
          ellipsvec(1)=(x0+ellipar*random_single()*(1.0_dp-x0))**(2.0_dp/3.0_dp)

          ellipsvec(2)=ellipsvec(1)

          ellipsvec(3)=1.0_dp/ellipsvec(1)**2

       end if
    end do

    if(cluster.and.(vacuum.gt.0.0_dp).and.sphere) then

       if(ellipsoid) write (stderr,'(a,3f12.3)') 'Using a confining ellipsoid: ',radius*ellipsvec

       fix_cell=.true.
       if((all(lattice_abc.eq.0.0_dp)).and.(all(lattice_car.eq.0.0_dp))) then
          lattice_abc(1)=2*radius*ellipsvec(1)+vacuum
          lattice_abc(2)=2*radius*ellipsvec(2)+vacuum
          lattice_abc(3)=2*radius*ellipsvec(3)+vacuum
          lattice_abc(4)=90.0_dp
          lattice_abc(5)=90.0_dp
          lattice_abc(6)=90.0_dp
       end if
    end if

    ! ** Confining cylinder?

    cylinder=.false.
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#CYLINDER=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#CYLINDER=')+10:))
          cylinder=.true.
          read(ctemp2,*) radius
          if(radius.gt.0.0_dp) then
             write (stderr,'(a,f5.2)') 'Using a confining cylinder of radius ',radius
          else
             write (stderr,'(a,f5.2)') 'Using an attractive line potential with strength ',abs(radius)
          end if
       end if
    end do

    ! ** Find a supercell consistent with the number of cells requested

    if(findsup) call num_cells_to_super()

    ! ** If a supercell calculation, how many cells?

    if(supcell) then
       num_cells = num_super_cell(supercell)
       if(num_cells.le.0) then
          write (stderr,*) 'Zero or negative volume supercell'
          stop
       end if
    else
       num_cells = 1
    end if

    ! ** Set P1 by default

    call symm_by_number(1)  

    ! * Set the lowest ranked space group acceptable (all by default)

    nrank = 230
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SGRANK=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#SGRANK=')+8:))
          read(ctemp2,*) nrank
       end if
    end do

    ! ** Read the symmetry information

    symm_ops = 0.0_dp
    setting_perm = ident
    num_symm = 0
    internalsymm = .true.

    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK SYMMETRY_OPS')>0) then
          internalsymm = .false.
          symmetry_name="external"
          do
             n=n+1 ; if(n.gt.nrecords) exit
             ctemp=cellbuff(n)
             if(index(up(ctemp),'%ENDBLOCK SYMMETRY_OPS')>0) exit
             if(len_trim(ctemp).gt.0) num_symm=num_symm+1
          end do
          exit
       end if

    end do
    num_symm=num_symm/4

    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK SYMMETRY_OPS')>0) then
          do i=1,num_symm
             n=n+1 ; if(n.gt.nrecords) exit
             read(cellbuff(n),*) symm_ops(:,1,i)
             n=n+1 ; if(n.gt.nrecords) exit
             read(cellbuff(n),*) symm_ops(:,2,i)
             n=n+1 ; if(n.gt.nrecords) exit
             read(cellbuff(n),*) symm_ops(:,3,i)
             n=n+1 ; if(n.gt.nrecords) exit
             read(cellbuff(n),*) symm_ops(:,4,i)
          end do
          exit
       end if

    end do

    ! * Check if the symmetry is defined by number of operations

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SYMMOPS=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#SYMMOPS=')+9:))

          if(index(ctemp2,'~')>0) then
             symmapprox=.true.
             ctemp2=ctemp2(index(ctemp2,'~')+1:)
          endif

          ! * Possibly choose from a range of symmops

          if(index(ctemp2,'-')>0) then
             read(ctemp2(:index(ctemp2,'-')-1),*) lowso
             read(ctemp2(index(ctemp2,'-')+1:),*) highso

             if(cluster) then
                if(.not.any(lowso.eq.numsclops)) stop 'Number of symmetry operations not allowed for a cluster'
                if(.not.any(highso.eq.numsclops)) stop 'Number of symmetry operations not allowed for a cluster'
             else                
                if(.not.any(lowso.eq.numsops)) stop 'Number of symmetry operations not allowed for a crystal'
                if(.not.any(highso.eq.numsops)) stop 'Number of symmetry operations not allowed for a crystal'
             end if

             do
                rn = random_triple()
                if(cluster) then
                   nops = numsclops(1+int(rn(1)*13))
                else
                   nops = numsops(1+int(rn(1)*10))
                end if
                if((nops.ge.min(lowso,highso)).and.(nops.le.max(lowso,highso))) exit
             end do

          else             
             read(ctemp2,*) nops
          endif

          if (cluster) then

             select case(nops)

             case(1)
                ns=1
                sel(1) = 1
             case(2)
                ns=3
                sel(1) = 2
                sel(2) = 3
                sel(3) = 4
             case(3)
                ns=1
                sel(1) = 5
             case(4)
                ns=5
                sel(1) = 6
                sel(2) = 15                
                sel(3) = 20 
                sel(4) = 25               
                sel(5) = 34
             case(5)
                ns=1
                sel(1) = 7
             case(6)
                ns=5
                sel(1) = 8
                sel(2) = 16
                sel(3) = 21
                sel(4) = 26
                sel(5) = 35
             case(7)
                ns=1
                sel(1) = 9
             case(8)
                ns=7
                sel(1) = 10
                sel(2) = 17
                sel(3) = 22
                sel(4) = 27
                sel(5) = 30
                sel(6) = 32
                sel(7) = 36
             case(9)
                ns=1
                sel(1) = 11
             case(10)
                ns=5
                sel(1) = 12
                sel(2) = 18
                sel(3) = 23
                sel(4) = 28
                sel(5) = 37
             case(11)
                ns=1
                sel(1) = 13
             case(12)
                ns=8
                sel(1) = 14
                sel(2) = 19
                sel(3) = 24
                sel(4) = 29
                sel(5) = 31
                sel(6) = 33
                sel(7) = 38
                sel(8) = 39
             case(24)
                ns=2
                sel(1) = 40
                sel(2) = 41
             case default
                write (stderr,*) 'Requested number of operations not available'
                stop
             end select

          else

             select case(nops)

             case(1)
                ns=1
                sel(1) = 1
             case(2)
                ns=0
                do i=2,9
                   ns=ns+1
                   sel(ns)=i
                end do
             case(3)
                ns=0
                do i=143,146
                   ns=ns+1
                   sel(ns)=i
                end do
             case(4)
                ns=0
                do i=10,46
                   ns=ns+1
                   sel(ns)=i
                end do
                do i=75,82
                   ns=ns+1
                   sel(ns)=i
                end do
             case(6)
                ns=0
                do i=147,161
                   ns=ns+1
                   sel(ns)=i
                end do
                do i=168,174
                   ns=ns+1
                   sel(ns)=i
                end do
             case(8)
                ns=0
                do i=47,74
                   ns=ns+1
                   sel(ns)=i
                end do
                do i=83,122
                   ns=ns+1
                   sel(ns)=i
                end do
             case(12)
                ns=0
                do i=162,167
                   ns=ns+1
                   sel(ns)=i
                end do
                do i=175,190
                   ns=ns+1
                   sel(ns)=i
                end do
                do i=195,199
                   ns=ns+1
                   sel(ns)=i
                end do
             case(16)
                ns=0
                do i=123,142
                   ns=ns+1
                   sel(ns)=i
                end do
             case(24)
                ns=0
                do i=200,220
                   ns=ns+1
                   sel(ns)=i
                end do
             case(48)
                ns=0
                do i=221,230
                   ns=ns+1
                   sel(ns)=i
                end do
             case default
                write (stderr,*) 'Requested number of operations not available'
                stop
             end select

          end if

          randsymm=.true.
          nsel=ns

          call select_symmetry(symmstat)

          if(symmstat.gt.0) then
             stat=888
             return
          endif

       end if
    end do

    ! * Check if the symmetry is defined by name

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SYMM=')>0) then
          ctemp2=ctemp(index(up(ctemp),'#SYMM=')+6:)
          if(index(ctemp2,'~')>0) then
             symmapprox=.true.
             ctemp2=ctemp2(index(ctemp2,'~')+1:)
          endif
          ctemp=ctemp2
          call symm_by_name(ctemp)

          call random_setting()

       end if
    end do

    ! * Check if the symmetry is defined by number

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SYMMNO=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#SYMMNO=')+8:))
          if(index(ctemp2,'~')>0) then
             symmapprox=.true.
             ctemp2=ctemp2(index(ctemp2,'~')+1:)
          endif
          if(index(ctemp2,',')>0) then

             m=0
             do i=1,len(ctemp2)
                if(ctemp2(i:i)==',') m=m+1
             end do

             ns = 0
             do i=1,m+1

                if(index(ctemp2,',').gt.0) then
                   ctemp = ctemp2(:index(ctemp2,',')-1)
                else
                   ctemp=ctemp2
                end if

                if(index(ctemp,'-')>0) then

                   read(ctemp(:index(ctemp,'-')-1),*) symmin
                   read(ctemp(index(ctemp,'-')+1:),*) symmax

                   if(symmin.gt.symmax) goto 999

                   do ni=symmin,symmax
                      ns=ns+1
                      sel(ns) = ni
                   end do

                else
                   ns=ns+1
                   read(ctemp,*) symmno
                   sel(ns) = symmno
                end if

                ctemp2=ctemp2(index(ctemp2,',')+1:)

             end do

             rn = random_triple()
             symmno = 1+int(rn(1)*(ns))
             call symm_by_number(icluster*sel(symmno))  

          else if(index(ctemp2,'-')>0) then

             read(ctemp2(:index(ctemp2,'-')-1),*) symmin
             read(ctemp2(index(ctemp2,'-')+1:),*) symmax

             rn = random_triple()
             symmno = symmin+int(rn(1)*(symmax-symmin+1))
             call symm_by_number(icluster*symmno)  
          else
             read(ctemp2,*) symmno
             call symm_by_number(icluster*symmno)
          end if

          call random_setting()

       end if
    end do

    ! * We at least want the identity

    if(num_symm==0) then
       if(.not.cluster)then
          symmetry_name="P1"
       else
          symmetry_name="C1"
       end if
       num_symm=1
       symm_ops(1,1,1) = 1.0_dp
       symm_ops(2,2,1) = 1.0_dp
       symm_ops(3,3,1) = 1.0_dp
    end if

    ! ** Check symmetry is consistent with requested geometry

    notsymm=check_symm()

    ! ** Check symmetry is consistent with cell constraints

    if((.not.cluster).and.(fix_cell.and.(.not.notsymm))) notsymm=.not.check_lattsymm()

    if(notsymm) then
       stat=888
       return
    endif

    ! * Report some information about the symmetry to be used

    write (stderr,'(a,a,a,i3)') 'Symm: ',symmetry_name,' Nops: ',num_symm

    ! * How many formula units (or rather, copies of input cell)

    num_form=-1
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#NFORM=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#NFORM=')+7:))
          if(index(ctemp2,'-')>0) then
             read(ctemp2(:index(ctemp2,'-')-1),*) num_form_min
             read(ctemp2(index(ctemp2,'-')+1:),*) num_form_max
             rn = random_triple()       
             num_form = num_form_min+int(rn(1)*(num_form_max-num_form_min+1))
          else
             read(ctemp2,*) num_form
          endif
       end if
    end do

    if(num_form.gt.0) write (stderr,'(a,i4)') "Formula units: ",num_form

    ! * Adjust the number of general positions

    adjgen=0
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#ADJGEN=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#ADJGEN=')+8:))
          if(index(ctemp2,'-')>0) then
             read(ctemp2(:index(ctemp2,'-')-1),*) adjgen_min
             read(ctemp2(index(ctemp2,'-')+1:),*) adjgen_max
             rn = random_triple()       
             adjgen = num_form_min+int(rn(1)*(adjgen_max-adjgen_min+1))
          else
             read(ctemp2,*) adjgen
          endif
       end if
    end do

    adjgen=adjgen+nudgen

    if(adjgen.gt.0) write (stderr,'(a,i5)') "Adjusting number of general positions: ",adjgen

    ! ** Read the cell contents

    num_ions=0
    ionabs = .false.

    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK POSITIONS_FRAC')>0) then
          do
             n=n+1 ; if(n.gt.nrecords) exit
             ctemp=cellbuff(n)
             if(index(up(ctemp),'%ENDBLOCK POSITIONS_FRAC')>0) exit
             num_ions=num_ions+1
          end do

          exit
       end if

    end do

    if(num_ions==0) then

       n=0
       do
          n=n+1 ; if(n.gt.nrecords) exit
          ctemp=cellbuff(n)
          if(index(up(ctemp),'%BLOCK POSITIONS_ABS')>0) then
             do
                n=n+1 ; if(n.gt.nrecords) exit
                ctemp=cellbuff(n)
                if(index(up(ctemp),'%ENDBLOCK POSITIONS_ABS')>0) exit
                num_ions=num_ions+1
             end do

             exit
          end if

       end do

       ionabs = .true.

    end if

    ! ** Allocate the cell data arrays

    if(allocated(ion_new_positions))  deallocate(ion_new_positions,ion_temp_positions,ion_hold_positions,ion_nn,ion_fix,ion_perm,&
         ion_push,ion_adatom,ion_moved,ion_min,ion_occ,ion_mult,ion_fails,ion_names,ion_names_temp,ion_rad,ion_coord,spinvec,&
         ion_bfix,ion_positions,ion_set,ion_equiv,posamp,minamp,angamp,xamp,yamp,zamp,ion_ion_vec,ion_species,species_names,&
         ion_set_temp,ion_positions_temp,ion_fix_temp,ion_bfix_temp,afmvec,ion_adatom_temp,posamp_temp,xamp_temp,yamp_temp,&
         zamp_temp,ion_athole,ion_athole_temp)

    n=max(num_ions*num_symm*num_cells,num_ions*num_form*num_symm)

    allocate(ion_new_positions(3,n))
    allocate(ion_hold_positions(3,n))
    allocate(ion_push(3,n))
    allocate(ion_temp_positions(3,n))
    allocate(ion_nn(2,n))
    allocate(ion_fix(n))
    allocate(ion_fix_temp(n))    
    allocate(ion_bfix(n))
    allocate(ion_bfix_temp(n))
    allocate(ion_perm(n))
    allocate(ion_adatom(n))
    allocate(ion_adatom_temp(n))
    allocate(ion_athole(n))
    allocate(ion_athole_temp(n))
    allocate(ion_moved(n))
    allocate(ion_min(n))
    allocate(ion_occ(n))
    allocate(ion_mult(n))
    allocate(ion_fails(n))
    allocate(ion_equiv(num_symm,n))
    allocate(ion_names(n))
    allocate(ion_names_temp(n))
    allocate(ion_rad(n))
    allocate(ion_coord(2,n))
    allocate(spinvec(n))
    allocate(afmvec(n))
    allocate(ion_positions(3,n))
    allocate(ion_positions_temp(3,n))
    allocate(ion_set(n))
    allocate(ion_set_temp(n))    
    allocate(posamp(n))
    allocate(posamp_temp(n))
    allocate(minamp(n))
    allocate(angamp(n))
    allocate(xamp(n))
    allocate(yamp(n))
    allocate(zamp(n))
    allocate(xamp_temp(n))
    allocate(yamp_temp(n))
    allocate(zamp_temp(n))
    allocate(ion_ion_vec(3,n,n))
    allocate(ion_species(n))
    allocate(species_names(n))

    ! ** Overall position amplitude

    if(radius.lt.0.0_dp) then
       posamp = -1.0_dp
    else
       posamp = radius
    end if
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#POSAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#POSAMP=')+8:))
          read(ctemp2,*) posamp(1)
          posamp = posamp(1)
       end if
    end do

    ! ** Overall minimum position amplitude

    minamp = 0.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#MINAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#MINAMP=')+8:))
          read(ctemp2,*) minamp(1)
          minamp = minamp(1)
       end if
    end do

    ! ** Overall position x-amplitude

    xamp = -1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#XAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#XAMP=')+6:))
          read(ctemp2,*) xamp(1)
          xamp = xamp(1)
       end if
    end do

    ! ** Overall position y-amplitude

    yamp = -1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#YAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#YAMP=')+6:))
          read(ctemp2,*) yamp(1)
          yamp = yamp(1)
       end if
    end do

    ! ** Overall position z-amplitude

    zamp = -1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#ZAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#ZAMP=')+6:))
          read(ctemp2,*) zamp(1)
          zamp = zamp(1)
       end if
    end do

    ! ** Overall angular amplitude

    angamp = -1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#ANGAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#ANGAMP=')+8:))
          read(ctemp2,*) angamp(1)
          angamp = angamp(1)
       end if
    end do

    ! ** Symmetry breaking amplitude

    breakamp = -1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#BREAKAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#BREAKAMP=')+10:))
          read(ctemp2,*) breakamp
       end if
    end do

    ! ** Overall coordination

    ion_coord = -1
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#COORD=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#COORD=')+7:))
          if(index(ctemp2,'-').gt.0) then
             read(ctemp2(:index(ctemp2,'-')-1),*) ion_coord(1,1)
             read(ctemp2(index(ctemp2,'-')+1:),*) ion_coord(2,1)
             ion_coord(1,:) = ion_coord(1,1)
             ion_coord(2,:) = ion_coord(2,1)
          else
             read(ctemp2,*) ion_coord(1,1)
             ion_coord = ion_coord(1,1)
          end if
       end if
    end do

    ! ** Number of constraint failures to tolerate

    num_fails = 0
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#NFAILS=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#NFAILS=')+8:))
          read(ctemp2,*) num_fails
       end if
    end do

    ! ** Overall radius

    ion_rad = 0.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#RAD=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#RAD=')+5:))
          read(ctemp2,*) ion_rad(1)
          ion_rad = ion_rad(1)
       end if
    end do

    ion_nn=' '
    ion_set=' '
    ion_names=' '

    ion_occ=1.0_dp
    ion_positions=0.0_dp
    ion_mult=-1
    ion_fix=.false.
    ion_bfix=.false.
    ion_perm=.false.
    ion_adatom=.false.
    ion_athole=.false.
    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK POSITIONS_')>0) then

          do i=1,num_ions
             n=n+1 ; if(n.gt.nrecords) exit
             ctemp=cellbuff(n)
             if(index(ctemp,'#')>0) then
                read(ctemp(1:index(ctemp,'#')),*) ion_names(i),ion_positions(:,i)
                if(index(ctemp,'%')>0) then
                   read(ctemp(index(ctemp,'#')+1:index(ctemp,'%')-1),*) ion_set(i)
                   if(index(up(ctemp),'POSAMP')>0) read(ctemp(index(up(ctemp),'POSAMP')+7:),*) posamp(i)
                   if(index(up(ctemp),'MINAMP')>0) read(ctemp(index(up(ctemp),'MINAMP')+7:),*) minamp(i)
                   if(index(up(ctemp),'ANGAMP')>0) read(ctemp(index(up(ctemp),'ANGAMP')+7:),*) angamp(i)
                   if(index(up(ctemp),'XAMP')>0)   read(ctemp(index(up(ctemp),'XAMP')+5:),*) xamp(i)
                   if(index(up(ctemp),'YAMP')>0)   read(ctemp(index(up(ctemp),'YAMP')+5:),*) yamp(i)
                   if(index(up(ctemp),'ZAMP')>0)   read(ctemp(index(up(ctemp),'ZAMP')+5:),*) zamp(i)
                   if(index(up(ctemp),'RAD')>0)    read(ctemp(index(up(ctemp),'RAD')+4:),*) ion_rad(i)
                   if(index(up(ctemp),'OCC')>0) then
                      ctemp2=ctemp(index(up(ctemp),'OCC')+4:)
                      if(index(ctemp2,'/').gt.0) then
                         read(ctemp2(:index(ctemp2,'/')-1),*) num
                         read(ctemp2(index(ctemp2,'/')+1:),*) den
                         ion_occ(i)=real(num,dp)/real(den,dp)
                      else   
                         read(ctemp2,*) ion_occ(i)
                      endif
                   end if
                   if(index(up(ctemp),'MULT')>0)   read(ctemp(index(up(ctemp),'MULT')+5:),*) ion_mult(i)
                   if(ion_mult(i) > 0) ion_occ(i)=real(ion_mult(i),dp)/real(num_symm,dp)
                   if(index(up(ctemp),'FIX')>0) then
                      ion_fix(i)=.true.
                      ion_bfix(i)=.true.
                   end if
                   if(index(up(ctemp),'NOMOVE')>0) ion_bfix(i)=.true.
                   if(index(up(ctemp),'PERM')>0) ion_perm(i)=.true.
                   if(index(up(ctemp),'ADATOM')>0) ion_adatom(i)=.true.
                   if(index(up(ctemp),'ATHOLE')>0) ion_athole(i)=.true.

                   if(index(up(ctemp),'COORD')>0) then 
                      ctemp2=ctemp(index(up(ctemp),'COORD')+6:)
                      if(index(ctemp2(:index(ctemp2,' ')),'-').gt.0) then
                         read(ctemp2(:index(ctemp2,'-')-1),*) ion_coord(1,i)
                         read(ctemp2(index(ctemp2,'-')+1:),*) ion_coord(2,i)
                      else
                         read(ctemp2,*) ion_coord(1,i)
                         ion_coord(2,i) = ion_coord(1,i)
                      end if
                   end if
                   if(index(up(ctemp),'NN')>0) then
                      ctemp2=ctemp(index(up(ctemp),'NN')+3:)
                      if(ctemp2(1:1)=='-') then
                         read(ctemp2(2:),*) ion_nn(1,i)
                      else if(ctemp2(1:1)=='+') then
                         read(ctemp2(2:),*) ion_nn(2,i)
                      else 
                         goto 999
                      end if
                   end if
                else
                   read(ctemp(index(ctemp,'#')+1:),*) ion_set(i)
                end if
             else
                read(ctemp,*) ion_names(i),ion_positions(:,i)
                write(ctemp2,*) i
                write (ion_set(i),'(a,a)') "Ion-",trim(adjustl(ctemp2))
             endif
          end do

          exit
       end if

    end do

    ! ** Fill in for the symmetry related ions

    do ns=2,num_symm
       ion_coord(:,1+(ns-1)*num_ions:num_ions+(ns-1)*num_ions)=ion_coord(:,1:num_ions)
    end do

    do ns=2,num_symm
       ion_nn(:,1+(ns-1)*num_ions:num_ions+(ns-1)*num_ions)=ion_nn(:,1:num_ions)
    end do

    do ns=2,num_symm
       ion_names(1+(ns-1)*num_ions:num_ions+(ns-1)*num_ions)=ion_names(1:num_ions)
    end do

    if(allocated(nearestn)) deallocate(nearestn,indx_nearestn)
    allocate(nearestn(3,20),indx_nearestn(20))

    ! ** Read pseudopotentials

    num_pot = 0
    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK SPECIES_POT')>0) then
          do
             n=n+1 ; if(n.gt.nrecords) exit
             ctemp=cellbuff(n)
             if(index(up(ctemp),'%ENDBLOCK SPECIES_POT')>0) exit
             num_pot=num_pot+1
          end do
          exit
       end if

    end do

    if(allocated(spec_pot)) deallocate(spec_pot)
    allocate(spec_pot(num_pot))
    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK SPECIES_POT')>0) then
          do i=1,num_pot
             n=n+1 ; if(n.gt.nrecords) exit
             spec_pot(i)=trim(cellbuff(n))
          end do
          exit
       end if
    end do

    ! ** Read Hubbard U information

    num_hub = 0
    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK HUBBARD_U')>0) then
          do
             n=n+1 ; if(n.gt.nrecords) exit
             ctemp=cellbuff(n)
             if(index(up(ctemp),'%ENDBLOCK HUBBARD_U')>0) exit
             num_hub=num_hub+1
          end do
          exit
       end if

    end do

    if(allocated(hub_u)) deallocate(hub_u)
    allocate(hub_u(num_hub))
    n=0
    do
       n=n+1 ; if(n.gt.nrecords) exit
       ctemp=cellbuff(n)
       if(index(up(ctemp),'%BLOCK HUBBARD_U')>0) then
          do i=1,num_hub
             n=n+1 ; if(n.gt.nrecords) exit
             hub_u(i)=trim(cellbuff(n))
          end do
          exit
       end if
    end do

    ! ** Read other info

    kpspacing = ''
    i=0
    do n=1,nrecords
       ctemp=cellbuff(n)
       if((index(up(ctemp),'KPOINTS_MP_')>0).or.(index(up(ctemp),'KPOINT_MP_')>0)) then
          i=i+1
          kpspacing(i)=trim(ctemp)       
       end if
    end do

    symmetrytol = ''

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'SYMMETRY_TOL')>0) then
          symmetrytol=trim(ctemp)       
       end if
    end do

    ! ** Determine the constraint on the cell generation
    !    - 0 means total freedom, 1 means totally fixed (to cubic), 0.4 default (aspect ratio of 4)

    cons = 0.4_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#CONS=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#CONS=')+6:))
          read(ctemp2,*) cons
       end if
    end do

    ! ** Angle constraint

    acons = 0.5_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#ACONS=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#ACONS=')+7:))
          read(ctemp2,*) acons
       end if
    end do

    ! ** Focus on a particular composition

    focus_species=0
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#FOCUS=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#FOCUS=')+7:))
          read(ctemp2,*) focus_species
       end if
    end do

    ctemp=formula(ns)
    if((focus_species.gt.0).and.(focus_species.ne.ns)) then
       stat=888
       return
    end if

    ! ** Check chemistry

    octet=any(cellbuffup.eq.'#OCTET')

    if(octet) then
       call check_chem(stat)
       if(stat.eq.888) return
    end if

    ! ** Randomly flip (mirror image) fragments

    flip=any(cellbuffup.eq.'#FLIP')

    ! ** Remove atoms that end up on top of each other

    remove=any(cellbuffup.eq.'#REMOVE')

    ! ** Cut a hole?

    hole_rad=-1.0_dp
    hole_pos=random_triple()
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#HOLE=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#HOLE=')+6:))
          read(ctemp2,*) hole_rad
       end if
       ! ** In fractional coordinates
       if(index(up(ctemp),'#HOLEPOS=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#HOLEPOS=')+9:))
          read(ctemp2,*) hole_pos
       end if
    end do

    do i=1,num_ions
       if(ion_athole(i)) ion_positions(:,i)=hole_pos
    end do

    ! ** Determine the minimum floating atom separation 

    minsep = -999.9_dp
    pair_names = ''
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#MINSEP=')>0) then

          if(index(up(ctemp),'AUTO')>0) then

             ctemp=ctemp(:index(up(ctemp),'AUTO')-1)

             auto=.false.
             composition = formula(ns)

             tempfile=trim(adjustl(composition))//'.minsep'
             ctemp6='comp2minsep '//trim(composition)//' >> '//trim(tempfile) ; call system(ctemp6)
             open(unit=unit_minsep,file=tempfile,status='UNKNOWN')
             autotargvol=''
             autominsep=''
             do
                read(unit_minsep,'(a)',err=456,end=456) autotargvol
                read(unit_minsep,'(a)',err=456,end=456) autominsep
             end do
456          close(unit=unit_minsep)

             if(len_trim(autotargvol(index(up(autotargvol),'#TARGVOL=')+9:)).gt.0) then
                ctemp6=autotargvol(index(up(autotargvol),'#TARGVOL=')+9:)
                read(ctemp6,*) targvolauto
                targvolauto=targvolauto/real(ns,dp)*real(num_ions,dp)
                auto=.true.
             end if
             if(len_trim(autominsep(index(up(autominsep),'#MINSEP=')+8:)).gt.0) then
                ctemp=autominsep
                auto=.true.
             end if

             write (stderr,'(a)') 'Note - automatic minsep/targvol for '//trim(adjustl(composition))
             write (stderr,'(a)') trim(autotargvol)
             write (stderr,'(a)') trim(autominsep)

          end if


          ctemp4=trim(ctemp(index(up(ctemp),'#MINSEP=')+8:))
          ctemp5=adjustl(trim(ctemp4(:index(up(ctemp4),' '))))
          if (index(ctemp5,"-").gt.0) then
             read(ctemp5(:index(ctemp5,"-")-1),*) minmin
             read(ctemp5(index(ctemp5,"-")+1:),*) minmax
             minsep=minmin
          else
             read(ctemp5,*) minmin
             minmax=minmin
             minsep=minmin
          end if

          ctemp4=adjustl(trim(ctemp4(index(up(ctemp4),' ')+1:)))
          num_pairs=0
          pair_minsep=minsep
          do while (len_trim(adjustl(ctemp4)).gt.0)
             num_pairs=num_pairs+1
             ctemp5=adjustl(trim(ctemp4(:index(up(ctemp4),' '))))
             read(ctemp5(:index(ctemp5,"-")-1),*) pair_names(1,num_pairs)
             read(ctemp5(index(ctemp5,"-")+1:index(ctemp5,"=")-1),*) pair_names(2,num_pairs)
             ctemp5=adjustl(trim(ctemp5(index(ctemp5,"=")+1:)))
             if (ctemp5(1:1)=='-') then
                read(ctemp5,*) pair_minsep(num_pairs)
             else if(index(ctemp5,"-").gt.0) then
                read(ctemp5(:index(ctemp5,"-")-1),*) minmin
                read(ctemp5(index(ctemp5,"-")+1:),*) minmax
                rn = random_triple()
                pair_minsep(num_pairs)=minmin+rn(1)*(minmax-minmin)
             else
                read(ctemp5,*) pair_minsep(num_pairs)
             end if
             ctemp4=adjustl(trim(ctemp4(index(up(ctemp4),' ')+1:)))
          end do

       end if
    end do
    
    ! ** Determine the target volume 

    targvol = volume
    fix_vol = .false.
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#TARGVOL=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#TARGVOL=')+9:))
          if(index(ctemp2,"-").gt.0) then
             read(ctemp2(:index(ctemp2,"-")-1),*) vmin
             read(ctemp2(index(ctemp2,"-")+1:),*) vmax
             rn = random_triple()
             targvol=vmin+rn(1)*(vmax-vmin)
          else
             read(ctemp2,*) targvol
          endif
          fix_vol=.true.
       end if
    end do

    ! ** Determine the cell shaking amplitude

    cellamp = -1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#CELLAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#CELLAMP=')+9:))
          read(ctemp2,*) cellamp
       end if
    end do

    ! ** How many vacancies in the final cell

    num_vacancies = 0
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#VACANCIES=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#VACANCIES=')+11:))
          if(index(ctemp2,"@")>0) then
             read(ctemp2(:index(ctemp2,"@")-1),*) num_vacancies
             read(ctemp2(index(ctemp2,"@")+1:),*) vac_name             
          else
             read(ctemp2,*) num_vacancies
          end if
       end if
    end do

    ! ** Tight packing? 

    tight=any(cellbuffup.eq.'#TIGHT')

    ! ** PUSH - PUt and SHake.

    push=.not.any(cellbuff.eq.'#NOPUSH') ! * Later on, deactivate for units

    ! ** Stepsize for PUSH

    pushstep = 0.25_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#PUSHSTEP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#PUSHSTEP=')+10:))
          read(ctemp2,*) pushstep
       end if
    end do

    ! ** Max PUSH steps

    pushmax = 100
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#PUSHMAX=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#PUSHMAX=')+9:))
          read(ctemp2,*) pushmax
       end if
    end do

    ! ** What overlap to tolerate - using hard-sphere potentials

    overlap = -999.9_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#OVERLAP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#OVERLAP=')+9:))
          read(ctemp2,*) overlap
       end if
    end do

    ! ** Use a three body hard sphere potential

    three = -999.9_dp
    do n=1,nrecords
       ctemp=cellbuff(n)

       if(index(up(ctemp),'#THREE=')>0) then
          ctemp4=trim(ctemp(index(up(ctemp),'#THREE=')+7:))
          ctemp5=adjustl(trim(ctemp4(:index(up(ctemp4),' '))))

          read(ctemp5,*) three
          stop 'Three body term not currently activated'
       end if

    end do

    ! ** Use a confining slab potential?

    width = -999.9_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#WIDTH=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#WIDTH=')+7:))
          read(ctemp2,*) width
          write (stderr,'(a,f5.2)') 'Using a slab spacer of width ',width
       end if
    end do

    shift = 0.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SHIFT=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#SHIFT=')+7:))
          read(ctemp2,*) shift
          if(width.gt.0.0_dp) write (stderr,'(a,f5.2)') 'Shifted by ',shift
       end if
    end do

    ! ** Use relax and shake algorithm (RASH) 

    dorash=any(cellbuffup.eq.'#RASH')

    rash_posamp = 1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#RASH_POSAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#RASH_POSAMP=')+13:))
          read(ctemp2,*) rash_posamp
       end if
    end do

    rash_angamp = 30.0_dp 
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#RASH_ANGAMP=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#RASH_ANGAMP=')+13:))
          read(ctemp2,*) rash_angamp
       end if
    end do


    ! ** Vary the cell shape during optimisation 

    celladapt=any(cellbuffup.eq.'#CELLADAPT')

    ! ** What is the minimum bond angle

    minbangle = 0.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#MINBANGLE=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#MINBANGLE=')+11:))
          read(ctemp2,*) minbangle
       end if
    end do

    ! ** What is the maximum bond angle

    maxbangle = 180.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#MAXBANGLE=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#MAXBANGLE=')+11:))
          read(ctemp2,*) maxbangle
       end if
    end do

    ! ** What is the attempt time

    maxtime = 1.0_dp ! * 1 second max by default
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#MAXTIME=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#MAXTIME=')+9:))
          read(ctemp2,*) maxtime
       end if
    end do

    ! ** Fractional slack on the bonding

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SLACK=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#SLACK=')+7:))
          read(ctemp2,*) slack
       end if
    end do

    ! ** Automatically increase slack

    if(.not.autoslack) then
       do n=1,nrecords
          ctemp=cellbuff(n)
          if(index(up(ctemp),'#AUTOSLACK=')>0) then
             ctemp2=trim(ctemp(index(up(ctemp),'#AUTOSLACK=')+11:))
             autoslack=.true.
             read(ctemp2,*) slack
             slack=slack-slack_step
          end if
       end do
    end if

    ! ** Determine the variable target volume 

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#VARVOL=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#VARVOL=')+8:))
          read(ctemp2,*) targvol
       end if
    end do

    ! ** Applied spin

    spin    = 0.0_dp
    spinmod = 0.0_dp
    spinlist = ''
    havespin = .false.
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#SPIN=')>0) then
          havespin=.true.
          ctemp2=trim(ctemp(index(up(ctemp),'#SPIN=')+6:))
          read(ctemp2,*,end=88) spin,spinmod,spinlist
88        read(ctemp2,*) spin,spinmod
       end if
    end do

    ! ** Scale the target volume by the number of actual ions in the cells

    targvol = targvol * scale_vol

    if(auto) targvol=targvolauto

    ! ** If a supercell calculation - scale the target volume

    targvol = targvol * num_cells

    ! ** Molecular units? 

    molecules=any(cellbuffup.eq.'#MOLECULES')

    ! ** What fraction to permute?

    permfrac = 1.0_dp
    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#PERMFRAC=')>0) then
          ctemp2=trim(ctemp(index(up(ctemp),'#PERMFRAC=')+10:))
          read(ctemp2,*) permfrac
       end if
    end do

    ! ** Permute the ions? 

    perm_names=''
    permute=.false.

    ! ** Specify by species?

    do n=1,nrecords
       ctemp=cellbuff(n)
       if(index(up(ctemp),'#PERMUTE')>0) permute=.true.
       if(index(up(ctemp),'#PERMUTE=')>0) then
          perm_names=trim(ctemp(index(up(ctemp),'#PERMUTE=')+9:))
       end if
    end do

    do i=1,num_ions
       if(index(perm_names,trim(ion_names(i))).gt.0) ion_perm(i)=.true.
    end do

    if(permute.and.(.not.any(ion_perm))) ion_perm=.true.

    ! ** Set the permute flags randomly

    do n=1,num_ions
       if(random_single().ge.permfrac) ion_perm(n)=.false.
    end do

    if (.not.any(ion_perm)) permute=.false.

    ! ** Compact the final cell? Not if the cell is supposed to be fixed in any way, not if a cluster

    compact=((.not.cluster).and.(.not.(fix_cell.or.fix_caxis.or.fix_abaxis)).and.(cellamp.lt.0.0_dp))
    
    ! ** Override to force compaction
    
    compact=compact.or.any(cellbuffup.eq.'#COMPACT')
    
    compact=compact.and.(.not.(any(cellbuffup.eq.'#NOCOMPACT'))) ! NOT CORRECT
    
    if(any(cellbuffup.eq.'#COMPACT').and.any(cellbuffup.eq.'#NOCOMPACT')) then
       write (stderr,'(a)') 'Warning: COMPACT and NOCOMPACT'
    end if

    if(con_cell.and.compact) then
       write (stderr,'(a)') 'Cell contraints - setting compact=.false.'
       compact=.false.
    end if
 
    if(compact) then
       write (stderr,'(a)') 'Compacting the unit cell using a Niggli reduction'
    else
       write (stderr,'(a)') 'Cell not compacted'
    end if
    
    temp = cons*(len_max-len_min)
    len_min_use = len_min+temp/2.0_dp
    len_max_use = len_max-temp/2.0_dp

    if(.not.all(lattice_abc==0.0_dp)) then

       ! ** Convert to cartesian

       lattice_car(:,1) = (/lattice_abc(1),0.0_dp,0.0_dp/)
       lattice_car(:,2) = (/lattice_abc(2)*cos(dgrd*lattice_abc(6)),lattice_abc(2)*sin(dgrd*lattice_abc(6)),0.0_dp/)
       lattice_car(1,3) = lattice_abc(3)*cos(dgrd*lattice_abc(5))
       lattice_car(2,3) = lattice_abc(3)*(cos(dgrd*lattice_abc(4))-cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6)))&
            /sin(dgrd*lattice_abc(6))
       lattice_car(3,3) = sqrt(lattice_abc(3)**2-lattice_car(1,3)**2-lattice_car(2,3)**2)

    end if

    ! ** Calculate the reciprocal lattice vectors, and recalculate volume etc

    call update_cell()

!!$    if(.not.internalsymm) then
!!$
!!$       stop 'coding in progress'
!!$
!!$       ! ** Identify the equivalent ions, and cull, if external symmetries provided
!!$
!!$       call cell_symm_equiv()
!!$
!!$       do ni=1,num_ions
!!$          do ns=1,num_symm
!!$             write (stderr,*) ns,ni,ion_equiv(ns,ni)
!!$          end do
!!$       end do
!!$
!!$       ion_mult=0
!!$       do ni=1,num_ions
!!$          if(ion_occ(ni).eq.0.0_dp) cycle
!!$          do ns=1,num_symm
!!$             if(ion_occ(ion_equiv(ns,ni)).eq.0.0_dp) cycle
!!$             ion_occ(ion_equiv(ns,ni))=0.0_dp
!!$             ion_mult(ni)=ion_mult(ni)+1
!!$          end do
!!$       end do
!!$       write (stderr,*)
!!$       ion_occ=0.0_dp
!!$       do ni=1,num_ions
!!$
!!$          if(ion_mult(ni).eq.0.0_dp) cycle
!!$          ion_occ(ni)=real(ion_mult(ni),dp)/real(num_symm,dp)
!!$          write (stderr,*) ni,ion_mult(ni),ion_occ(ni)
!!$
!!$       end do
!!$
!!$    end if

    ! ** If we need to use special positions and MULT have not been specified

    !if((num_form.gt.0).and.(num_form.ne.num_symm)) then

!!$    do i=1,num_ions
!!$       write (stderr,*) i,ion_names(i),ion_fix(i),ion_bfix(i),ion_positions(:,i)
!!$    end do
!!$    write (stderr,*)

!!$    if((num_form.gt.0).and.internalsymm) then
    if(num_form.gt.0) then

       targvol=targvol*real(num_form,dp)/real(num_symm,dp)

       m=0
       ion_names_temp=ion_names
       ion_set_temp=ion_set
       ion_positions_temp=ion_positions
       ion_fix_temp=ion_fix
       ion_bfix_temp=ion_bfix
       ion_adatom_temp=ion_adatom
       ion_athole_temp=ion_athole
       posamp_temp=posamp
       xamp_temp=xamp
       yamp_temp=yamp
       zamp_temp=zamp
       do j=1,num_form
          do i=1,num_ions
             m=m+1
             ion_names(m)=ion_names_temp(i)
             ion_fix(m)=ion_fix_temp(i)
             ion_bfix(m)=ion_bfix_temp(i)
             ion_adatom(m)=ion_adatom_temp(i)
             ion_athole(m)=ion_athole_temp(i)
             ion_positions(:,m)=ion_positions_temp(:,i)
             posamp(m)=posamp_temp(i)
             xamp(m)=xamp_temp(i)
             yamp(m)=yamp_temp(i)
             zamp(m)=zamp_temp(i)
             write(ctemp,*) j
             ion_set(m)=trim(ion_set_temp(i))//"-"//trim(adjustl(ctemp))
          end do
       end do

       num_ions=m

       ion_names_temp=ion_names
       ion_set_temp=ion_set
       ion_positions_temp=ion_positions
       ion_fix_temp=ion_fix
       ion_bfix_temp=ion_bfix
       ion_adatom_temp=ion_adatom
       ion_athole_temp=ion_athole
       posamp_temp=posamp
       xamp_temp=xamp
       yamp_temp=yamp
       zamp_temp=zamp

       allocate(flag(num_ions),mult(num_ions))

       ! ** Put like species (and fix states) together

       flag=.true.
       m=0
       do i=1,num_ions

          do j=1,num_ions

             if((ion_names_temp(i).eq.ion_names_temp(j)).and.flag(j)&
                  .and.(ion_bfix_temp(i).eqv.ion_bfix_temp(j))&
                  .and.(ion_fix_temp(i).eqv.ion_fix_temp(j))) then
                m=m+1
                ion_names(m)=ion_names_temp(j)
                ion_set(m)=ion_set_temp(j)
                ion_fix(m)=ion_fix_temp(j)
                ion_bfix(m)=ion_bfix_temp(j)
                ion_adatom(m)=ion_adatom_temp(j)
                ion_athole(m)=ion_athole_temp(j)
                posamp(m)=posamp_temp(j)
                xamp(m)=xamp_temp(j)
                yamp(m)=yamp_temp(j)
                zamp(m)=zamp_temp(j)
                ion_positions(:,m)=ion_positions_temp(:,j)
                flag(j)=.false.
             end if

          end do

       end do

       ! ** Work out the size of the blocks, and possible multiplicities

       ion_names_temp=ion_names
       ion_positions_temp=ion_positions
       ion_set_temp=ion_set
       ion_fix_temp=ion_fix
       ion_bfix_temp=ion_bfix
       ion_adatom_temp=ion_adatom
       ion_athole_temp=ion_athole
       posamp_temp=posamp

       flag=.true.
       m=0
       do i=1,num_ions

          if(flag(i)) then

             ! ** Identify a block of ions

             ni=0
             inblock=.false.
             do j=i,num_ions
                if(ion_bfix_temp(j).and.(flag(j))) then
                   if(inblock) exit
                   ni=1
                   exit
                else
                   if((ion_names_temp(j).eq.ion_names_temp(i)).and.flag(j)) then 
                      ni=ni+1
                      if(i.ne.j) flag(j)=.false.
                      inblock=.true.
                   else
                      if(inblock) exit
                   end if
                end if
             end do

             ! ** Evaluate a possible multiplicity

             mult=0
             call choose_mult(ni,num_symm,mult,n)
             if(any(mult.gt.1)) write (stderr,'(a6,a3,800i3)') 'Mult: ',trim(ctemp),mult(1:n)
             
             ! ** Set the ionic data

             do j=1,n
                m=m+1
                ion_mult(m)=mult(j)
                ion_names(m)=ion_names_temp(i-1+j)
                write (ctemp,*) ion_mult(m)
                ion_set(m)=trim(adjustl(ion_set_temp(i-1+j)))//"-"//trim(adjustl(ctemp))
                ion_positions(:,m)=ion_positions_temp(:,i-1+j)
                ion_fix(m)=ion_fix_temp(i-1+j)
                ion_bfix(m)=ion_bfix_temp(i-1+j)
                ion_adatom(m)=ion_adatom_temp(i-1+j)
                ion_athole(m)=ion_athole_temp(i-1+j)
                posamp(m)=posamp_temp(i-1+j)
                xamp(m)=xamp_temp(i-1+j)
                yamp(m)=yamp_temp(i-1+j)
                zamp(m)=zamp_temp(i-1+j)
             end do


          end if

       end do

       num_ions=m
       do i=1,num_ions
          if(ion_mult(i) > 0) then
             ion_occ(i)=real(ion_mult(i),dp)/real(num_symm,dp)
          elseif (ion_mult(i) < 0) then
             ion_occ(i) = real(ion_mult(i),dp)
          endif
       end do

       deallocate(flag,mult)

       ! ** Unset the multiplicities for the fixed ions

       do ni=1,num_ions
          if(ion_bfix(ni)) then
             ion_mult(ni)=-1
             ion_occ(ni)=1.0_dp
          end if
       end do

    end if

!!$    do i=1,num_ions
!!$       write (stderr,*) i,ion_names(i),ion_occ(i),ion_mult(i),ion_positions(:,i)
!!$    end do
!!$    write (stderr,*)
    ! ** Generate a super cell

    call super_cell()

    ! ** Cut hole

    call cut_hole()

    ! ** Check number of atoms -- symmetry?

    if((num_atom.gt.0).and.(num_symm.eq.1)) then

       n=0
       do ni=1,num_ions
          if(ion_occ(ni).gt.1.0_dp-delta) n=n+1
       end do

       if(n.ne.num_atom) then
          write (stderr,*) num_atom,n,num_ions,num_form,num_symm
          stat=888
          return
       end if
       
    end if

    ! ** Save the original cell for shaking

    lattice_car_orig=lattice_car
    lattice_rec_orig=lattice_rec

    stat=0

    return

999 write (stderr,*) 'There is a problem reading the cell information. Stopping.'
    stop

  end subroutine init_cell

  subroutine gen_sets()

    real(kind=dp) :: dist
    
    ! ** Number of sets of atoms

    num_sets = 0
    do n=1,num_ions
       if(.not.any(ion_set(1:n-1)==ion_set(n))) num_sets = num_sets+1
    end do
    
    if(allocated(set_label)) then
       deallocate(set_label,ion_set_centre,ion_set_sphere,ion_set_positions)
       deallocate(ninset,ion_set_min,posamp_set,minamp_set,angamp_set)
       deallocate(ion_set_names,ion_set_rad,ion_set_coord,ion_set_nn,ion_set_fix,ion_set_bfix)
       deallocate(xamp_set,yamp_set,zamp_set,ion_ion_set,ion_set_occ,ion_set_mult)
    end if

    allocate(set_label(num_sets))
    allocate(ion_set_centre(3,num_sets))
    allocate(ion_set_sphere(num_sets))
    allocate(ion_set_positions(3,num_ions,num_sets))
    allocate(ninset(num_sets))
    allocate(ion_set_min(num_ions,num_sets))
    allocate(posamp_set(num_sets))
    allocate(minamp_set(num_sets))
    allocate(angamp_set(num_sets))
    allocate(ion_set_names(num_ions,num_sets))
    allocate(ion_set_rad(num_ions,num_sets))
    allocate(ion_set_coord(2,num_ions,num_sets))
    allocate(ion_set_nn(2,num_ions,num_sets))
    allocate(xamp_set(num_sets))
    allocate(yamp_set(num_sets))
    allocate(zamp_set(num_sets))
    allocate(ion_ion_set(num_ions*num_symm))
    allocate(ion_set_occ(num_ions,num_sets))
    allocate(ion_set_mult(num_ions,num_sets))
    allocate(ion_set_fix(num_ions,num_sets))
    allocate(ion_set_bfix(num_ions,num_sets))

    ! ** Find out the labels for the sets

    m = 0

    do n=1,num_ions

       if(.not.any(ion_set(1:n-1)==ion_set(n))) then 
          m = m+1
          set_label(m)  = ion_set(n)
          posamp_set(m) = posamp(n) ! ** NOTE: take from the first in each set, ignore otherwise
          minamp_set(m) = minamp(n)
          angamp_set(m) = angamp(n)
          xamp_set(m)   = xamp(n)
          yamp_set(m)   = yamp(n)
          zamp_set(m)   = zamp(n)
       end if

    end do

    ! ** Divide up into sets

    ninset = 0

    do n=1,num_ions

       do m=1,num_sets

          if(ion_set(n)==set_label(m)) then
             ninset(m) = ninset(m)+1
             ion_set_positions(:,ninset(m),m) = ion_positions(:,n)
             ion_set_names(ninset(m),m) = ion_names(n)
             ion_set_coord(:,ninset(m),m) = ion_coord(:,n)
             ion_set_nn(:,ninset(m),m) = ion_nn(:,n)
             ion_set_rad(ninset(m),m) = ion_rad(n)
             ion_set_occ(ninset(m),m) = ion_occ(n)
             ion_set_mult(ninset(m),m) = ion_mult(n)
             ion_set_fix(ninset(m),m) = ion_fix(n)
             ion_set_bfix(ninset(m),m) = ion_bfix(n)
          end if

       end do

    end do
    
    ! ** Which ion is in which set

    ni=0
    do n=1,num_sets
       do ns=1,num_symm
          do m=1,ninset(n)
             ni=ni+1
             ion_ion_set(ni)=n+(ns-1)*10000
          end do
       end do
    end do

    ! ** Convert to absolute coordinates

    if(.not.ionabs) then
       do n=1,num_sets
          do m=1,ninset(n)
             ion_set_positions(:,m,n) = ion_set_positions(1,m,n)*lattice_car(:,1) &
                  + ion_set_positions(2,m,n)*lattice_car(:,2) &
                  + ion_set_positions(3,m,n)*lattice_car(:,3)
          end do
       end do
    end if
    
    ! ** Store relative the geometrical centre

    do n=1,num_sets
       do i=1,3

          ion_set_centre(i,n) = sum(ion_set_positions(i,1:ninset(n),n))/real(ninset(n),dp)

          ion_set_positions(i,1:ninset(n),n) = ion_set_positions(i,1:ninset(n),n) - ion_set_centre(i,n)

       end do
    end do

    ! ** Calculate the radius of a sphere that holds the set

    do n=1,num_sets

       ion_set_sphere=0.0_dp
       
       do i=1,ninset(n)
          dist=sqrt(dot_product(ion_set_positions(:,i,n),ion_set_positions(:,i,n)))
          if(dist.gt.ion_set_sphere(n)) ion_set_sphere(n)=dist
       end do
       
    end do
    
    ! ** Work out the minimum bond length within a set
    
    ion_set_min = 0.0_dp

    do n=1,num_sets

       if(ninset(n).gt.1) then

          do m=1,ninset(n)

             ion_set_min(m,n) = huge(1.0_dp)
             
             do i=1,ninset(n)
                
                v = ion_set_positions(:,m,n)-ion_set_positions(:,i,n)
                
                do j=1,27
                   
                   if(m.ne.i) then
                      
                      dist = sqrt(dot_product(lookup(:,j)+v,lookup(:,j)+v))
                      
                      if(dist<ion_set_min(m,n)) ion_set_min(m,n) = dist
                      
                   end if
                   
                end do
                
             end do

          end do

       end if
       
    end do
    
    ! * Scale by 125% to give a bit of space, and shrink a bit if "tight"

    if(.not.tight) then
       ion_set_min = ion_set_min*1.25_dp
    else
       ion_set_min = ion_set_min*0.90_dp
    end if

    ! * Override using minsep

    if(minsep.ge.0.0_dp) then
       ion_set_min = minsep
    end if

    ! * Deactivate translation to cell for clusters

    if(cluster) translate=.false.

    ! * Add vacuum to cell
    
    if((abs(lattice_abc(4)-90.0_dp).lt.1e-5_dp).and.(abs(lattice_abc(5)-90.0_dp).lt.1e-5_dp)&
         .and.(vacuum.ne.0.0_dp).and.(cellamp.lt.0.0_dp).and.(.not.cluster)) then

       write (stderr,'(a,f6.2)') 'Padding vacuum:',vacuum

       lattice_car(3,3)=lattice_car(3,3)+vacuum
       cell_consistent=.false.
       
       call update_cell()

    end if
    
  end subroutine gen_sets

  subroutine write_cell()

    real(kind=dp) :: shift(3)
    character(len=255) :: hostname,pwd,user
    
    ! ** Write header

    call get_environment_variable("PWD",pwd)
    call hostnm(hostname)
    call get_environment_variable("USER",user)
        
    write (stdout,'(a)')      '#              Generated by Buildcell'
    write (stdout,'(a)')      '#'
    write (stdout,'(a,a)')    '# at ',trim(date())
    write (stdout,'(a,a)')    '# in ',trim(pwd)
    write (stdout,'(a,a)')    '# on ',trim(hostname)
    write (stdout,'(a,a)')    '# by ',trim(user)
    write (stdout,'(a)')      '#'
#ifdef COMPAT    
    write (stdout,'(a,a)')    '# compiler ','unknown'
    write (stdout,'(a,a)')    '# options ','unknown'
#else
    write (stdout,'(a,a)')    '# compiler ',trim(compiler_version())
    write (stdout,'(a,a)')    '# options ',trim(compiler_options())
#endif
    write (stdout,'(a)')      '#'
    write (stdout,'(a,50i12)')'# seed',seed
    write (stdout,'(a)')      '#'
    write (stdout,'(a)')      '#      Author: C. J. Pickard (cjp20@cam.ac.uk)'
    write (stdout,'(a,a)')    '#                Space group: ',trim(symmetry_name)
    write (stdout,'(a)')      '#'

    write (stdout,*)

    ! ** Write the lattice vectors

!!$    write (stdout,'(a)') "%BLOCK LATTICE_CART"
!!$    write (stdout,'(3f15.10)') lattice_car(:,1)
!!$    write (stdout,'(3f15.10)') lattice_car(:,2)
!!$    write (stdout,'(3f15.10)') lattice_car(:,3)
!!$    write (stdout,'(a,/)') "%ENDBLOCK LATTICE_CART"

    write (stdout,'(a)') "%BLOCK LATTICE_ABC"
    write (stdout,'(3f13.5)') lattice_abc(1:3)
    write (stdout,'(3f13.5)') lattice_abc(4:6)
    write (stdout,'(a,/)') "%ENDBLOCK LATTICE_ABC"

    ! ** Write the cell contents

    if(permute) then

       if(any(ion_perm)) then
          do i=1,num_ions*num_symm
             if (ion_perm(i)) then
                ion_names_temp(i) = 'ABCD'
             else
                ion_names_temp(i) = ion_names(i)
             end if
          end do
       else
          ion_names_temp = 'ABCD'
       end if
       
       do i=1,num_ions*num_symm
          if(.not.ion_perm(i)) cycle

          do
             rn = random_triple()
             n = 1+int(rn(1)*(num_ions*num_symm))
             
             if(ion_names_temp(n)=='ABCD') then
                ion_names_temp(n)=ion_names(i)
                exit
             end if             
          end do

       end do

       ion_names=ion_names_temp

    end if
    
    ! ** Add some vacancies

    if(num_vacancies.gt.0) then
       if(num_vacancies.ge.num_ions*num_symm) stop 'num_vacancies too large - no ions left'
       i=0
       do 
          rn = random_triple()       
          n = 1+int(rn(1)*(num_ions*num_symm))
          if((ion_names(n).ne.'Z').and.(ion_names(n).eq.vac_name)) then
             i=i+1
             ion_names(n)='Z'
          end if
          if(i.eq.num_vacancies) exit
       end do

    end if

    ! ** Write out the ions

    if(cluster) then
       !if(.true.) then
       ! * Move the centre of mass to the centre of the cell
       
       shift=0.0_dp

       do i=1,num_ions*num_symm
          shift(:)=shift(:)+matmul(lattice_rec,ion_new_positions(:,i)) 
       end do

       shift=0.5_dp-shift/real(num_ions*num_symm,dp)

    else

       shift=0.0_dp

    end if
    
    temp=0.0_dp
    write (stdout,'(a)') "%BLOCK POSITIONS_FRAC"
    do i=1,num_ions*num_symm
       v(:) = matmul(lattice_rec,ion_new_positions(:,i)) + shift(:)
       if((trim(ion_names(i)).ne.'Z').and.(ion_occ(i)==1.0_dp)) then
          ctemp=strip(ion_names(i))
          if(index(ctemp,'-').gt.0) ctemp=ctemp(1:index(ctemp,'-')-1)
          if(havespin) then
             write (stdout,'(a4,3f12.7,a6,f10.5)') trim(ctemp),v(:),"SPIN=",spinvec(i)
             temp=temp+spinvec(i)
          else
             write (stdout,'(a4,3f17.7)') trim(ctemp),v(:)
          end if
       end if
    end do
    write (stdout,'(a,/)') "%ENDBLOCK POSITIONS_FRAC"

    if(havespin) write(stdout,'(a,f7.2,/)') '# Note: total spin =',temp
    
    ! ** Write the kpoint spacing, if required

    do i=1,20
       if(len_trim(kpspacing(i))>0) then
          write (stdout,'(a,/)') trim(kpspacing(i))
       end if
    end do

    if(len_trim(symmetrytol)>0) then
       write (stdout,'(a,/)') trim(symmetrytol)
    end if

    ! ** Fix cell if required

    if(fix_all_cell) then
       write (stdout,'(a,/)') "FIX_ALL_CELL : true"
    end if    

    ! ** Fix cell volume if required

    if(fix_cell_vol) then
       write (stdout,'(a,/)') "FIX_VOL : true"
    end if

    ! ** Fix lattice vectors

    if(fix_caxis) cellcon=(/1,2,0,0,0,3/)
    
    if(fix_abaxis) cellcon=(/0,0,1,2,3,0/)       
    
    if(((sum(abs(cellcon-(/1,2,3,4,5,6/))).gt.0).and.(.not.compact).and.(.not.fix_all_cell))&
         .or.fix_caxis.or.fix_abaxis) then
       write (stdout,'(a)') "%BLOCK CELL_CONSTRAINTS"
       write (stdout,'(3i2)') cellcon
       write (stdout,'(a,/)') "%ENDBLOCK CELL_CONSTRAINTS"       
    end if
    
    ! ** Carry the cluster flag through

    if(cluster) then
       write (stdout,'(a,/)') '#CLUSTER'
    end if

    ! ** Write out the pseudopotential information

    if(num_pot>0) then
       write (stdout,'(a)') "%BLOCK SPECIES_POT"
       do i=1,num_pot
          write (stdout,'(a)') trim(spec_pot(i))
       end do
       write (stdout,'(a,/)') "%ENDBLOCK SPECIES_POT"
    end if

    ! ** Write out the Hubbard U information

    if(num_hub>0) then
       write (stdout,'(a)') "%BLOCK HUBBARD_U"
       do i=1,num_hub
          write (stdout,'(a)') trim(hub_u(i))
       end do
       write (stdout,'(a,/)') "%ENDBLOCK HUBBARD_U"
    end if

    ! ** Write the symmetry operations (if more than just identity, and a crystal)

    if((.not.cluster).and.(num_symm>1).and.(breakamp<0.0_dp).or.symmgen) then
       
       write (stdout,'(a)') "SNAP_TO_SYMMETRY"
       write (stdout,'(a,/)') "SYMMETRY_GENERATE"
       
!!$       write (stdout,'(a)') "%BLOCK SYMMETRY_OPS"
!!$       do i=1,num_symm
!!$          
!!$          if(internalsymm) then
!!$             Sym(1:3,1:3,i)=matmul(lattice_car,matmul(symm_ops(1:3,1:3,i),lattice_rec))
!!$          else
!!$             Sym(1:3,1:3,i)=symm_ops(1:3,1:3,i)
!!$          endif
!!$
!!$          write (stdout,'(3f20.15)') Sym(1:3,1,i)
!!$          write (stdout,'(3f20.15)') Sym(1:3,2,i)
!!$          write (stdout,'(3f20.15)') Sym(1:3,3,i)
!!$          write (stdout,'(3f20.15)') symm_ops(1:3,4,i)
!!$ 
!!$       end do
!!$       write (stdout,'(a,/)') "%ENDBLOCK SYMMETRY_OPS"
    end if

    ! ** Write constraints

    if(any(ion_fix)) then

       write (stdout,'(a)') "FIX_COM : false"
       write (stdout,'(a)') "%BLOCK IONIC_CONSTRAINTS"
       n=0
       do i=1,num_ions*num_symm
          if(ion_fix(i).and.(ion_occ(i).gt.delta)) then
             m=0
             do j=1,i
                if(ion_occ(j).lt.delta) cycle
                if(ion_names(i).eq.ion_names(j)) m=m+1
             end do
             n=n+1
             write (stdout,'(i5,5x,a8,i5,3f12.8)') n,strip(ion_names(i)),m,(/1.0_dp,0.0_dp,0.0_dp/)
             n=n+1
             write (stdout,'(i5,5x,a8,i5,3f12.8)') n,strip(ion_names(i)),m,(/0.0_dp,1.0_dp,0.0_dp/)
             n=n+1
             write (stdout,'(i5,5x,a8,i5,3f12.8)') n,strip(ion_names(i)),m,(/0.0_dp,0.0_dp,1.0_dp/)
          end if
       end do
       write (stdout,'(a,/)') "%ENDBLOCK IONIC_CONSTRAINTS"

    end if

  end subroutine write_cell

  subroutine select_symmetry(stat)

    integer, intent(out) :: stat
    
    integer :: count
    
    if(randsymm) then

       count=0
       notsymm=.true.
       do while(notsymm.and.(count.lt.10000))

          count=count+1
          rn = random_triple()       
          symmno = 1+int(rn(1)*(nsel))
          call symm_by_number(icluster*sel(symmno)) 

          call random_setting()

          ! ** Check if this space group is consistent with requested geometry
          
          notsymm=check_symm()
          
          ! ** Check if this space group is consistent with fixed lattice constraints
          
          if((.not.cluster).and.(fix_cell.and.(.not.notsymm))) notsymm=.not.check_lattsymm()

          ! ** Check if space group common enough

          if((.not.cluster).and.(.not.notsymm)) notsymm=.not.any(sel(symmno).eq.sgrank(1:nrank))
                    
       end do
       
    end if

    if(notsymm) then
       stat=111
    else
       stat=0
    end if
    
  end subroutine select_symmetry

  function check_symm()

    logical :: check_symm
    
    if(cluster) then
       check_symm=.false.
    else if(slab) then
       check_symm=.false.
       if(any(symm_ops(3,1,1:num_symm).gt.delta)) check_symm=.true.
       if(any(symm_ops(3,2,1:num_symm).gt.delta)) check_symm=.true.
       if(any(symm_ops(1,3,1:num_symm).gt.delta)) check_symm=.true.
       if(any(symm_ops(2,3,1:num_symm).gt.delta)) check_symm=.true.
       if(any(abs(symm_ops(3,3,1:num_symm)).lt.1.0_dp-delta)) check_symm=.true.
       if(any(symm_ops(3,4,1:num_symm).gt.delta)) check_symm=.true.
       if (surface) then
          if(any(symm_ops(3,3,1:num_symm).lt.1.0_dp-delta)) check_symm=.true.
       end if
    else
       check_symm=(symmorphic.and.any(abs(symm_ops(:,4,1:num_symm)).gt.delta)) 
    end if

  end function check_symm
  
  subroutine random_setting()

    integer :: perm(3,3)

    !return ! ** Not currently working
    
    perm=0
    do while(any(sum(abs(perm),2).eq.0))
       perm=0
       perm(int(1+random_single()*3),1)=1
       perm(int(1+random_single()*3),2)=1
       perm(int(1+random_single()*3),3)=1
    end do
    
    do ns=1,num_symm

       symm_ops(1:3,1:3,ns)=matmul(matmul(perm,symm_ops(1:3,1:3,ns)),transpose(perm))
       symm_ops(1:3,4,ns)=matmul(perm,symm_ops(1:3,4,ns))

    end do
        
    cang=matmul(perm,cang)
    clen=matmul(perm,clen)
    
    ! ** Store

    setting_perm=perm
    
  end subroutine random_setting

  function check_lattsymm()

    logical :: check_lattsymm

    real(kind=dp) :: small=1e-5_dp

    check_lattsymm=.true.

    do i=1,3
       if(cang(i).gt.0.0_dp) then
          if(abs(cang(i)-lattice_abc(3+i)).gt.small) check_lattsymm=.false.
       end if
    end do

    do i=1,3
       do j=1,3

          if((cang(i).lt.0.0_dp).and.(cang(j).lt.0.0_dp)) then
             if(.not.(abs(lattice_abc(3+i)-lattice_abc(3+j)).lt.small)) then
                check_lattsymm=.false.
             end if
          end if

          if((clen(i).lt.0.0_dp).and.(clen(j).lt.0.0_dp)) then
             if(.not.(abs(lattice_abc(i)-lattice_abc(j)).lt.small)) then
                check_lattsymm=.false.
             end if
          end if

       end do
    end do

  end function check_lattsymm

  subroutine update_cell()

    integer :: n,n1,n2,n3

    if(cell_consistent) return

    if(all(lattice_car==0.0_dp)) return
    
    ! ** Calculate the volume

    volume = volume_cell(lattice_car)
    
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

    ! ** Calculate the lookup of lattice vectors

    n=0
    do n1=-1,1
       do n2=-1,1
          do n3=-1,1
             n=n+1
             lookup(:,n) = n1*lattice_car(:,1)+n2*lattice_car(:,2)+n3*lattice_car(:,3)
          end do
       end do
    end do
    
    ! ** Update symmetry operators
    
    do ns=1,num_symm
    
       if(internalsymm) then
          Sym(1:3,1:3,ns)=matmul(lattice_car,matmul((symm_ops(1:3,1:3,ns)),lattice_rec))
       else
          Sym(1:3,1:3,ns)=symm_ops(1:3,1:3,ns)
       endif

    end do

    ! ** Store the cell constraints
    
    call gen_cellcon()
    
    ! ** We now have a consistent cell
    
    cell_consistent=.true.

  end subroutine update_cell

  function volume_cell(lc)

    real(kind=dp) :: lc(3,3)

    real(kind=dp) :: volume_cell
    
    volume_cell = lc(1,1)*(lc(2,2)*lc(3,3)-lc(3,2)*lc(2,3))+lc(2,1)*(lc(3,2)*lc(1,3)-lc(1,2)*lc(3,3))+&
         lc(3,1)*(lc(1,2)*lc(2,3)-lc(2,2)*lc(1,3))
    
  end function volume_cell

  subroutine compact_cell()
    
    integer :: Cmat(3,3)
    
    if(compact) then
       call niggli_reduce(lattice_car,Cmat)    
       if(.not.all(Cmat.eq.nint(ident))) then
          do ns=1,num_symm
             symm_ops(1:3,1:3,ns)=matmul(invert3x3(real(Cmat,dp)),matmul(symm_ops(1:3,1:3,ns),real(Cmat,dp)))
             symm_ops(1:3,4,ns)=matmul(invert3x3(real(Cmat,dp)),symm_ops(1:3,4,ns))
          end do
          cell_consistent=.false.
       end if
    end if    
    
  end subroutine compact_cell
  
  subroutine num_cells_to_super()

    integer       :: nmax,stemp(3,3),isgn,k

    nmax=supercell(1,1)
    
    do
       do i=1,3
          do j=1,3
             rn = random_triple()       
             supercell(i,j) = int(rn(1)*(2*nmax+1))-nmax
          end do
       end do
       if(slab) then
          supercell(3,1:2) = (/0,0/)
          supercell(:,3) = (/0,0,1/)
       end if
       if(num_super_cell(supercell).eq.nmax) exit
    end do
    
    call compact_super_cell(supercell)

    if(slab) then
       stemp=supercell
       do i=1,3
          if((stemp(1,i).eq.0).and.((stemp(2,i).eq.0))) exit
       end do
       
       supercell(:,3) = (/0,0,1/)

       k=0
       do j=1,3
          if(j.eq.i) cycle
          k=k+1
          supercell(:,k)=stemp(:,j)
       end do

       isgn=num_super_cell(supercell)/abs(num_super_cell(supercell))
       
       supercell(:,1)=supercell(:,1)*isgn
       
    end if

    write (stderr,'(a,9i4)') "Using supercell: ",supercell
    
  end subroutine num_cells_to_super

  subroutine compact_super_cell(scell)

    integer, dimension(3,3), intent(inout) :: scell

    real(kind=dp) :: sc0(3,3)
    integer       :: Cmat(3,3)

    sc0=real(scell)
    call niggli_reduce(sc0,Cmat)
    scell=nint(sc0)
   
  end subroutine compact_super_cell
  
  subroutine super_cell

    integer                                    :: nc,ni,nii,i,j,k,nuse,nzero
    character(len=20)                          :: ctemp
    real(kind=dp), allocatable, dimension(:,:) :: supercell_vec,old_ion_positions
    real(kind=dp), dimension(3,3)              :: old_lattice_car,mtemp
    real(kind=dp), dimension(3)                :: vfrac

    if(.not.supcell) return

    ! ** Convert ion_positions to absolute coordinates

    if(.not.ionabs) then

       do ni=1,num_ions
          if(ion_adatom(ni)) cycle
          ion_positions(:,ni) = ion_positions(1,ni)*lattice_car(:,1)+&
               ion_positions(2,ni)*lattice_car(:,2)+&
               ion_positions(3,ni)*lattice_car(:,3)
       end do

    end if

    ! ** Store the old lattice vectors

    old_lattice_car = lattice_car

    ! ** Make new cell

    lattice_car = matmul(lattice_car,supercell)

    ! ** Update the cell data for consistency

    cell_consistent=.false.

    call update_cell()

    ! ** Determine the supercell sub-cell vectors

    allocate(supercell_vec(3,num_cells),old_ion_positions(3,num_ions))

    nuse=1

555 supercell_vec = 0.0_dp
    mtemp = matmul(lattice_rec,old_lattice_car)
    nzero=0
    nc=1
    supercell_vec(:,nc)=0.0_dp
    do i=-nuse,nuse 
       do j=-nuse,nuse
          do k=-nuse,nuse
             if((i==0).and.(j==0).and.(k==0)) cycle

             vfrac(:) = i*mtemp(:,1)+j*mtemp(:,2)+k*mtemp(:,3)   

             if(all(vfrac.gt.-delta).and.all(vfrac.lt.1.0_dp-delta)) then
                nc=nc+1
                if(nc.gt.num_cells) then
                   stop 'Supercell subvector generation failed - nc > num_cells'
                end if
                supercell_vec(:,nc) = i*old_lattice_car(:,1)+j*old_lattice_car(:,2)+k*old_lattice_car(:,3)
             end if

          end do
       end do
    end do

    if(nc.lt.num_cells) then
       !write (stderr,*) nc,num_cells
       !stop 'Supercell subvector generation failed - nc < num_cells'
       nuse=nuse*2
       goto 555
    end if

    ! ** Make copies of the cell contents    

    old_ion_positions(:,1:num_ions) = ion_positions(:,1:num_ions)
    nii=0
    do nc=1,num_cells
       do ni=1,num_ions

          if((.not.ion_adatom(ni)).or.(nc.eq.1)) then

             nii=nii+1

             ! * New ionic positons

             ion_positions(:,nii) = old_ion_positions(:,ni)+supercell_vec(:,nc)

             ! * Fill in the rest of the atomic data

             ion_nn(:,nii)    = ion_nn(:,ni)
             ion_min(nii)     = ion_min(ni)
             ion_occ(nii)     = ion_occ(ni)
             ion_mult(nii)    = ion_mult(ni)
             ion_names(nii)   = ion_names(ni)
             ion_rad(nii)     = ion_rad(ni)
             ion_fix(nii)     = ion_fix(ni)
             ion_bfix(nii)    = ion_bfix(ni)
             ion_perm(nii)    = ion_perm(ni)
             ion_adatom(nii)  = ion_adatom(ni)
             ion_athole(nii)  = ion_athole(ni)
             ion_coord(:,nii) = ion_coord(:,ni)
             write (ctemp,*) nc
             ion_set(nii)    = trim(ion_set(ni))//'-'//trim(adjustl(ctemp))
             posamp(nii)     = posamp(ni)
             minamp(nii)     = minamp(ni)
             angamp(nii)     = angamp(ni)
             xamp(nii)       = xamp(ni)
             yamp(nii)       = yamp(ni)
             zamp(nii)       = zamp(ni)

          endif

       end do
    end do
    
    ! ** Update num_ions to account for copies

    num_ions=nii

    ! ** Convert ionic positions to fractional, in the new supercell

    if(.not.ionabs) then

       do ni=1,num_ions
          if(ion_adatom(ni)) cycle
          ion_positions(:,ni)=matmul(lattice_rec,ion_positions(:,ni))
       end do

    end if

    deallocate(supercell_vec,old_ion_positions)

  end subroutine super_cell

  integer function num_super_cell(mat)

    implicit none
    
    integer, dimension(3,3), intent(inout) :: mat

    integer :: det,i,j,k

    det=0
    do i=1,3
       j=mod(i,3)+1
       k=mod(j,3)+1
       det=det+mat(1,i)*(mat(2,j)*mat(3,k)-mat(2,k)*mat(3,j))
    end do

    num_super_cell = det

  end function num_super_cell
  
  subroutine choose_mult(ntot,nsymm,mult,nmult)

    integer, intent(in) :: ntot
    integer, intent(in) :: nsymm
    integer, intent(out), dimension(:) :: mult
    integer, intent(out) :: nmult

    integer :: factors(10),numf(10),nfact,nf,ncount,adj_inequiv,ninequiv,nleft,minequiv

    if(ntot.le.0) stop 'choose_mult : positive number of ions required'

    ! ** Select the possible multiplicities

    select case(nsymm)
    case(1)
       nfact=1
       factors(1) = 1
    case(2)
       nfact=2
       factors(1) = 1
       factors(2) = 2
    case(3)
       nfact=2
       factors(1) = 1
       factors(2) = 3
    case(4)
       nfact=3
       factors(1) = 1
       factors(2) = 2
       factors(3) = 4
    case(5)
       nfact=2
       factors(1) = 1
       factors(2) = 5
    case(6)
       nfact=4
       factors(1) = 1
       factors(2) = 2
       factors(3) = 3
       factors(4) = 6
    case(7)
       nfact=2
       factors(1) = 1
       factors(2) = 7
    case(8)
       nfact=4
       factors(1) = 1
       factors(2) = 2
       factors(3) = 4
       factors(4) = 8
    case(9)
       nfact=3
       factors(1) = 1
       factors(2) = 3
       factors(3) = 9
    case(10)
       nfact=4
       factors(1) = 1
       factors(2) = 2
       factors(3) = 5
       factors(4) = 10
    case(11)
       nfact=2
       factors(1) = 1
       factors(2) = 11
    case(12)
       nfact=6
       factors(1) = 1
       factors(2) = 2
       factors(3) = 3
       factors(4) = 4
       factors(5) = 6
       factors(6) = 12
    case(16)
       nfact=5
       factors(1) = 1
       factors(2) = 2
       factors(3) = 4
       factors(4) = 8
       factors(5) = 16
    case(20)
       nfact=6
       factors(1) = 1
       factors(2) = 2
       factors(3) = 4       
       factors(4) = 5
       factors(5) = 10
       factors(6) = 20
    case(24)
       nfact=8
       factors(1) = 1
       factors(2) = 2
       factors(3) = 3
       factors(4) = 4
       factors(5) = 6
       factors(6) = 8
       factors(7) = 12
       factors(8) = 24
    case(48)
       nfact=10
       factors(1) = 1
       factors(2) = 2
       factors(3) = 3
       factors(4) = 4
       factors(5) = 6
       factors(6) = 8
       factors(7) = 12
       factors(8) = 16
       factors(9) = 24
       factors(10) = 48
    case default
       write (stderr,*) 'Requested nsymm not possible'
       stop
    end select

    ! ** If symmetry is only be approximately applied (filling general positions only)

    if(symmapprox) then       
       nmult=max(ntot/maxval(factors(1:nfact))-adjgen,0)
       mult(1:nmult) = maxval(factors(1:nfact))       
       mult(nmult+1:nmult+ntot-sum(mult(1:nmult)))=-1
       nmult=nmult+ntot-sum(mult(1:nmult))
       return
    end if

    ! ** Assume maximum use of most general positions

    numf=0
    nleft=ntot
    do nf=nfact,1,-1
       numf(nf)=nleft/factors(nf)
       nleft=nleft-numf(nf)*factors(nf)
    end do

    minequiv=sum(numf)

    if(adjgen.gt.0) then

       adj_inequiv=min(minequiv+adjgen,ntot)
       ninequiv=1+int(random_single()*adj_inequiv)
       ncount=0
       do
          ncount=ncount+1

          if(ncount>1000) then
             adj_inequiv=min(adj_inequiv+1,ntot)
             ninequiv=1+int(random_single()*adj_inequiv)
             ncount=0
          end if

          numf=0
          do n=1,ninequiv
             nf=1+int(random_single()*(nfact))
             numf(nf)=numf(nf)+1
          end do

          if(dot_product(numf(1:nfact),factors(1:nfact)).eq.ntot) exit

       end do
    end if

    ! ** Set the multiplicities

    nmult=0
    do nf=nfact,1,-1
       mult(nmult+1:nmult+numf(nf))=factors(nf)
       nmult=nmult+numf(nf)
    end do

  end subroutine choose_mult

  subroutine cell_symm_equiv()

    real(kind=dp) :: VO(3),posi(3),posj(3)
    real(kind=dp) :: symm_thresh=1e-6_dp

    ion_equiv=0

    if(num_symm==1) then
       do ni=1,num_ions
          ion_equiv(1,ni)=ni
       end do
       return 
    end if

    ! ** Which atoms are symmetry related to which?    
    
    do ni=1,num_ions
       if(ion_occ(ni).lt.1.0_dp-delta) cycle
       posi(:)=ion_positions(1,ni)*lattice_car(:,1)+ion_positions(2,ni)*lattice_car(:,2)+ion_positions(3,ni)*lattice_car(:,3)
       do nj=1,num_ions
          if(ion_occ(nj).lt.1.0_dp-delta) cycle
          posj(:)=ion_positions(1,nj)*lattice_car(:,1)+ion_positions(2,nj)*lattice_car(:,2)+ion_positions(3,nj)*lattice_car(:,3)

          do ns=1,num_symm
             
             VO=posi(:)-matmul(Sym(1:3,1:3,ns),posj(:))-&
                  symm_ops(1,4,ns)*lattice_car(1:3,1) - &
                  symm_ops(2,4,ns)*lattice_car(1:3,2) - &
                  symm_ops(3,4,ns)*lattice_car(1:3,3)
             
             VO = matmul(lattice_rec,VO)
             VO = VO - nint(VO)
             VO = matmul(lattice_car,VO)
             
             if(dot_product(VO,VO).lt.symm_thresh) ion_equiv(ns,ni)=nj

          end do

       end do
    end do
    
  end subroutine cell_symm_equiv

  subroutine symm_cell()

    real(kind=dp) :: lwork(3,3),lc_symm(3,3)
    integer :: symm_count(3),i,j
    
    lc_symm = 0.0_dp
    symm_count = 0
    do ns = 1, num_symm
       lwork = matmul(lattice_car, transpose(symm_ops(1:3,1:3,ns)))

       do i=1,3
          do j=1,3
             if (all(abs(lwork(i,:)-lattice_car(j,:)) < delta )) then
                symm_count(j) = symm_count(j) + 1
                lc_symm(j,:) = lc_symm(j,:) + lwork(i,:)
             end if
             !
             ! Include inversion symmetry even if not in space group - all Bravais lattices have this
             !
             if (all(abs(lwork(i,:)+lattice_car(j,:)) < delta )) then
                symm_count(j) = symm_count(j) + 1
                lc_symm(j,:) = lc_symm(j,:) - lwork(i,:)
             end if
          end do
       end do
    end do

    do j=1,3
       lc_symm(j, :) = lc_symm(j, :) / real(symm_count(j),dp)
    end do

    lattice_car=lc_symm
    
    cell_consistent=.false.

  end subroutine symm_cell
  
  function formula(nspec)

    integer, intent(out) :: nspec

    character(len=400) :: formula

    ! *---

    integer :: n,ni,ns,num(100),nform,element_index(117)
    real(kind=dp) :: element_order(117)
    character(len=16) :: name(100),cnum(100)
    logical :: have_c

    nspec=0
    name=''
    do ni=1,num_ions
       if(.not.any(name==ion_names(ni))) then
          nspec=nspec+1
          name(nspec)=ion_names(ni)
       end if
    end do

    do ns=1,nspec
       num(ns)=count(ion_names(1:num_ions).eq.name(ns))
    end do

    nform = num(1)
    do ns = 2,nspec
       nform = gcd_rec(nform,num(ns))
    end do
    num=num/nform

    do ns=1,nspec
       if(num(ns).gt.1) then
          write (cnum(ns),'(i10)') num(ns)
       else
          cnum(ns)=''
       end if
    end do

    if(any(name(:).eq."C  ")) have_c=.true.
    
    do ns=1,nspec
       do n=1,117
          if(name(ns).eq.elements_alpha(n)) exit
       end do
       if(name(ns).eq."O  ") then
          element_order(ns)=huge(1.0_dp)
       else if(name(ns).eq."C  ") then
          element_order(ns)=0.0_dp
       else if((name(ns).eq."H  ").and.have_c) then
          element_order(ns)=0.1_dp
       else
          element_order(ns)=real(n,dp)
       end if
       element_index(ns)=ns
    end do
    
    call heap_sort_index(nspec,element_order,element_index)

    write(formula,*) &
         (trim(name(element_index(ns))),&
         trim(adjustl(cnum(element_index(ns)))),ns=1,nspec)
    
  end function formula

  subroutine check_chem(stat)

    integer, intent(inout) :: stat

    integer :: ni,n,count_valence

    count_valence=0
    do ni=1,num_ions
       if(ion_occ(ni).eq.0.0_dp) cycle
       do n=1,117
          if(ion_names(ni).eq.elements_alpha(n)) exit
       end do
       if(elements_valence(n).eq.0) stop 'valence not known'
       count_valence=count_valence+elements_valence(n)
    end do

    if(mod(count_valence,8).ne.0) stat=888
    
  end subroutine check_chem

  subroutine gen_cellcon()

    integer :: n,m,inc    
    
    cellcon=0
    inc=1
    cellcon(1)=inc
    do n=2,3
       if(any(abs(clen(n)-clen(1:n-1)).lt.delta).and.(abs(clen(n)).gt.delta)) then
          do m=1,n-1
             if(abs(clen(n)-clen(m)).lt.delta) then
                cellcon(n)=cellcon(m)
                exit
             end if
          end do
       else
          inc=inc+1
          cellcon(n)=inc
       end if
       
    end do

    if(cang(1).lt.delta) then
       inc=inc+1
       cellcon(4)=inc
    else
       cellcon(4)=0
    end if
    do n=5,6
       if(cang(n-3).lt.delta) then
          if(any(abs(cang(n-3)-cang(1:n-4)).lt.delta).and.(abs(cang(n-3)).gt.delta)) then
             
             do m=1,n-4
                if(abs(cang(n-3)-cang(m)).lt.delta) then
                   cellcon(n)=cellcon(3+m)
                   exit
                end if
             end do
             
          else
             
             inc=inc+1
             cellcon(n)=inc
          end if
       else
          cellcon(n)=0
       endif
       
    end do    
    
  end subroutine gen_cellcon

  subroutine cut_hole()

    real(kind=dp) :: v(3),dist
    
    if(hole_rad.lt.0.0_dp) return

    write (stderr,'(a,f7.3,a,3f7.3)') 'Cutting a hole, radius ',hole_rad,' and position ',hole_pos
    
    do ni=1,num_ions
       if(ion_adatom(ni)) cycle
       v(:)=(ion_positions(1,ni)-hole_pos(1))*lattice_car(:,1)+&
            (ion_positions(2,ni)-hole_pos(2))*lattice_car(:,2)+&
            (ion_positions(3,ni)-hole_pos(3))*lattice_car(:,3)
       do j=1,27         
          dist = sqrt(dot_product(lookup(:,j)+v,lookup(:,j)+v))
          if(dist<hole_rad) ion_occ(ni)=0.0_dp
       end do
          
    end do

  end subroutine cut_hole
  
  function up(string)

    implicit none

    character(len=*), intent(in) :: string
    character(len=400)           :: up

    integer                      :: i,ic,ia,iz,ishift

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

  recursive function gcd_rec(u, v) result(gcd)
    integer             :: gcd
    integer, intent(in) :: u, v

    if (mod(u, v) /= 0) then
       gcd = gcd_rec(v, mod(u, v))
    else
       gcd = v
    end if
  end function gcd_rec

  subroutine heap_sort_index(num_items,weight,indx)

    !=========================================================================!
    ! This subroutine sorts the list of weights into descending order.        !
    ! The weights are unchanged, with index returning the result              !
    !                                                                         !
    ! This is a heap sort                                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   num_items (input) :: The number of items to sort                      !
    !   weight (in)       :: The weights of each item.                        !
    !   indx (inout)      :: The indices of the weights in decending order    !
    !-------------------------------------------------------------------------!
    ! Parent module variables used:                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Modules used:                                                           !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Key Internal Variables:                                                 !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    ! Necessary conditions:                                                   !
    !   None                                                                  !
    !-------------------------------------------------------------------------!
    !                                  !
    !=========================================================================!

    implicit none

    ! Arguments

    integer, intent(in) :: num_items
    real(kind=dp), dimension(num_items), intent(inout) :: weight
    integer, dimension(num_items), intent(inout) :: indx

    ! Local variables

    integer :: i,ir,j,l,indxa ! Loop counters
    real(kind=dp) :: wta
    real(kind=dp), dimension(num_items) :: wgtemp

    if(num_items.lt.2) return

    wgtemp=weight
!!$    do i=1,num_items
!!$       indx(i)=i
!!$    end do

    l=num_items/2+1
    ir=num_items

    do
       if(l.gt.1) then
          l=l-1
          wta=wgtemp(l)
          indxa=indx(l)
       else
          wta=wgtemp(ir)
          indxa=indx(ir)
          wgtemp(ir)=wgtemp(1)
          indx(ir)=indx(1)
          ir=ir-1
          if(ir.eq.1) then
             wgtemp(1)=wta
             indx(1)=indxa
             return
          end if
       end if
       i=l
       j=l+l
20     if(j.le.ir) then
          if(j.lt.ir) then
             if(wgtemp(j).lt.wgtemp(j+1)) j=j+1
          end if
          if(wta.lt.wgtemp(j)) then
             wgtemp(i)=wgtemp(j)
             indx(i)=indx(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          goto 20
       end if
       wgtemp(i)=wta
       indx(i)=indxa
    end do

  end subroutine heap_sort_index

  function invert3x3(mat)

    real(kind=dp), intent(in) :: mat(3,3)
    real(kind=dp)             :: invert3x3(3,3)

    integer :: i,j
    real(kind=dp) :: cof(3,3),det

    cof(1,1)=mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2)
    cof(1,2)=mat(2,3)*mat(3,1)-mat(2,1)*mat(3,3)
    cof(1,3)=mat(2,1)*mat(3,2)-mat(2,2)*mat(3,1)

    cof(2,1)=mat(1,3)*mat(3,2)-mat(1,2)*mat(3,3)
    cof(2,2)=mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1)
    cof(2,3)=mat(1,2)*mat(3,1)-mat(1,1)*mat(3,2)

    cof(3,1)=mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2)
    cof(3,2)=mat(1,3)*mat(2,1)-mat(1,1)*mat(2,3)
    cof(3,3)=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1)

    det=mat(1,1)*cof(1,1)+mat(1,2)*cof(1,2)+mat(1,3)*cof(1,3)

    if(abs(det).lt.tiny(1.0_dp)) stop 'invert3x3 : det=0'

    do i=1,3
       do j=1,3
          invert3x3(i,j)=cof(j,i)/det
       end do
    end do

  end function invert3x3
  
end module cell
