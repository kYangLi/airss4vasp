! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!==================================================================================!
!                                     cryan                                        !
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
! This is the main program of cryan, the crystal structure analyser                !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

program cryan

  use constants
  use rng
  use gp
  use modularity

  implicit none

  integer, parameter                              :: unit_qhull=30,unit_elim_sed=31,unit_elim_rm=32,unit_unit=33,unit_adj=34
  integer                                         :: num_lines=10000000
  integer                                         :: num_interp=1000,ngrid=10000
  integer                                         :: num_structures
  integer                                         :: num_compositions
  integer                                         :: num_total
  integer                                         :: num_tasks
  integer                                         :: num_top=huge(1)
  integer                                         :: num_tri  
  integer                                         :: num_hull  
  integer                                         :: num_normal
  integer                                         :: num_units=-1
  integer                                         :: hull_dim
  integer                                         :: num_species
  integer                                         :: num_points
  integer                                         :: element_index(118)
  integer                                         :: ns,i,j,k,na,nb,nc,nt,n,seed_indx,stat,l,nin,ninc,ne,natmpt,nas,nbs
  integer                                         :: n1,n2,nrep,dof,indx1,indx2,num_penth,num_pestruc,ntrain,indx
  integer                                         :: nref(1),num_struc
  integer,            allocatable, dimension(:)   :: sort_index
  integer,            allocatable, dimension(:)   :: comp_index
  integer,            allocatable, dimension(:)   :: comp_index_inc  
  integer,            allocatable, dimension(:)   :: comp_num
  integer,            allocatable, dimension(:)   :: hull_index
  integer,            allocatable, dimension(:)   :: similar
  integer,            allocatable, dimension(:)   :: orig_copies
  integer,            allocatable, dimension(:)   :: numform
  integer,            allocatable, dimension(:)   :: iwork,iwork2
  integer,            allocatable, dimension(:)   :: in
  integer,            allocatable, dimension(:)   :: clustocc
  integer,            allocatable, dimension(:)   :: dcom
  integer,            allocatable, dimension(:)   :: comp_factor
  integer,            allocatable, dimension(:,:) :: tri_index
  integer,            allocatable, dimension(:,:) :: adjacancy
  integer,            allocatable, dimension(:,:) :: comp_vec
  integer,            allocatable, dimension(:,:) :: clustsim
  integer,            allocatable, dimension(:,:) :: dia
  integer,            allocatable, dimension(:,:) :: eab
  integer,            allocatable, dimension(:,:) :: cmat
  integer,            allocatable, dimension(:,:) :: cvec

  real(kind=dp)                                   :: rmax=20.0_dp,bondlength=0.0_dp,bondscale=0.0_dp,bl2,dist2,R0(3),thresh
  real(kind=dp)                                   :: dmod=1.0_dp,alpha=-1.0_dp,vscale,lstsq(2),pmin,pmax,press
  real(kind=dp)                                   :: comp_threshold,elim_threshold,enth,temp,maxq,lscale,dimensionality=-1.0_dp
  real(kind=dp)                                   :: kern3(-1:1,-1:1,-1:1),emat3(-1:1,-1:1,-1:1),dos_smear=0.1_dp
  real(kind=dp)                                   :: eref,eform,ehull,etemp,sc_threshold,kern(-1:1,-1:1),emat(-1:1,-1:1)
  real(kind=dp)                                   :: emin,emax,press_add=0.0_dp,minvol,maxvol,scale_a,scale_b,minda,mindb
  real(kind=dp)                                   :: diff,newdiff,delta_e=huge(1.0_dp),element_order(118),last,vtemp(3)
  real(kind=dp),    allocatable, dimension(:)     :: sort_array,vec,comp_enth,E,work,S,So,comp_energy,point_energy,pressure,enthalpy
  real(kind=dp),    allocatable, dimension(:)     :: comp_efromhull,energy,disclap,ptrain,sdos,egrid,vol
  real(kind=dp),    allocatable, dimension(:,:)   :: maxdiff,LL,MM,SS,B,clustdist,unit_save,comp_frac,hull_normal,etrain,enref
  real(kind=dp),    allocatable, dimension(:,:,:) :: peinterp

  character(len=20)                               :: task(100)
  character(len=10)                               :: species_names(100)
  character(len=40)                               :: formula='',formula_convert='',formula_ref=''
  character(len=40)                               :: comp_1='',comp_2='',comp_3='',comp_4=''
  character(len=240)                              :: char_copies,char_new_copies  
  character(len=40)                               :: comp_structure=''
  character(len=40)                               :: ctemp4
  character(len=240)                              :: ctemp,ctemp2,ctemp3,unitfile,fmt,label
  character(len=40),  allocatable, dimension(:)   :: composition,struct_name,pestruct
  character(len=240), allocatable, dimension(:)   :: resbuff,ctext
  character(len=240), allocatable, dimension(:)   :: sumbuff

  logical, allocatable, dimension(:) :: taken
  logical, allocatable, dimension(:) :: moved
  logical, allocatable, dimension(:) :: comp_include
  logical, allocatable, dimension(:) :: comp_onhull

  logical :: delete=.false.
  logical :: cluster=.false.
  logical :: have_c=.false.
  logical :: xmgrace=.false.
  logical :: have_spin=.false.
  logical :: have_dos=.false.
  logical :: long=.false.
  logical :: first
  logical :: modular=.false.
  logical :: connected
  logical :: notsymm=.false.
  logical :: xyz=.false.
  logical :: off=.false.
  logical :: comout=.false.
  logical :: adj=.false.
  logical :: weight=.false.
  logical :: newcomp
  logical :: sloppy=.false.

  ! ** Define the crystal type

  type :: crystal
     integer                                          :: num_ions
     integer                                          :: num_copies
     integer                                          :: num_dist
     integer                                          :: num_species
     integer                                          :: num_form
     integer                                          :: num_units
     integer                                          :: community
     integer                                          :: composition
     integer,           allocatable, dimension(:)     :: species_num
     integer,           allocatable, dimension(:)     :: dist_indx
     integer,           allocatable, dimension(:)     :: in_unit
     integer,           allocatable, dimension(:)     :: num_unit
     integer,           allocatable, dimension(:)     :: unit_num_asym
     integer,           allocatable, dimension(:,:)   :: dist_ions
     integer,           allocatable, dimension(:,:)   :: adjacancy
     integer,           allocatable, dimension(:)     :: ion_names_indx
     integer,           allocatable, dimension(:,:)   :: unit_names_indx
     real(kind=dp)                                    :: pressure
     real(kind=dp)                                    :: volume
     real(kind=dp)                                    :: enthalpy
     real(kind=dp)                                    :: relative_enthalpy
     real(kind=dp)                                    :: spin
     real(kind=dp)                                    :: modspin
     real(kind=dp)                                    :: del
     real(kind=dp)                                    :: vol_calc
     real(kind=dp)                                    :: e_from_hull
     real(kind=dp)                                    :: dimensionality
     real(kind=dp)                                    :: lattice_abc(6)
     real(kind=dp)                                    :: lattice_car(3,3)
     real(kind=dp)                                    :: lattice_rec(3,3)
     real(kind=dp),     allocatable, dimension(:)     :: species_frac
     real(kind=dp),     allocatable, dimension(:)     :: distall
     real(kind=dp),     allocatable, dimension(:,:)   :: positions_frac
     real(kind=dp),     allocatable, dimension(:,:)   :: positions_cart
     real(kind=dp),     allocatable, dimension(:,:,:) :: units_cart
     real(kind=dp),     allocatable, dimension(:,:)   :: dist
     real(kind=dp),     allocatable, dimension(:,:)   :: laplacian
     character(len=240)                               :: structure_label
     character(len=20)                                :: structure_label_short
     character(len=40)                                :: symmetry
     character(len=40)                                :: formula
     character(len=10), allocatable, dimension(:)     :: ion_names
     character(len=10), allocatable, dimension(:,:)   :: unit_names
     character(len=10), allocatable, dimension(:)     :: unit_pointgp
     character(len=10), allocatable, dimension(:)     :: species_names
     character(len=10), allocatable, dimension(:,:)   :: dist_names
     logical,           allocatable, dimension(:)     :: dist_incell
     logical                                          :: on_hull
  end type crystal

  type(crystal), allocatable, dimension(:) :: structure

  !------------------------------------!

  ! ** Get arguments and determine run type

  call get_arguments()

  ! ** Initialise random number generator

  call init_pseudorandom()

  ! ** Read the structure data into a buffer

  call read_buff()

  ! ** Read the crystal structure data, find out the actual number of structures

  call read_res()

  ! ** Calculate the enthalpy relative to the most stable to each stoichiometry

  call compositions()

  ! ** Choose what to do with the data

  do nt=1,num_tasks

     select case (task(nt))

     case('formula') ; call formula_pick()

     case ('formula_convert') ; call formula_conv()

     case ('enthalpy') ; call enthalpy_plot()

     case('connect') ; call connect()

     case('geometry') ; call geometry()

     case('rank') ; call rank()

     case('summary') ; call summary()

     case('compare') ; call compare()

     case('eliminate') ; call eliminate()

     case('struct_comm') ; call struct_comm()

     case('maxwell') ; call maxwell()

     case('phull') ; call phull()

     case('dos') ; call dos()

     case default ; stop 'Task unknown'
        
     end select
     
  end do

contains

  subroutine read_buff()

    integer :: n

    allocate(resbuff(num_lines))
    n=0
    do 
       n=n+1
       if(n.gt.num_lines) stop 'increase num_lines'
       read(5,'(a)',end=99,err=100) resbuff(n)
    end do
100 stop 'Problem reading structure data'
99  continue
    
    num_lines=n-1
    
    num_structures=count(resbuff(1:num_lines).eq."END")

    allocate(structure(num_structures))

  end subroutine read_buff

  subroutine read_res()

    integer            :: ns,nl,ni,nstart,nend,nrem,indx1,indx2,indx3,nlen,ntokens
    character(len=240) :: ctemp,tokens(100)
    character(len=10)  :: name(100)

    !---------------------------------

    nend=1
    ns=0
    do while (nend.lt.num_lines)
       ns=ns+1

       ! ** Find the start of the next record

       do nstart=nend,num_lines
          if(index(resbuff(nstart),"TITL").gt.0) exit
       end do

       ! ** Find the end of the record

       do nend=nstart,num_lines
          if(resbuff(nend) .eq. "END") exit
       end do

       ! ** Back track in case of partial records

       nrem=0
       do nstart=nend,1,-1
          if(index(resbuff(nstart),"TITL").gt.0) exit
          if(resbuff(nstart)(1:3).eq."REM") nrem=nrem+1
       end do

       ! ** Read the header data

       indx1=index(resbuff(nstart),"(")
       indx2=index(resbuff(nstart),")")

       if(indx1.eq.0) indx1=index(resbuff(nstart)," n -")

       structure(ns)%num_ions=-1
       structure(ns)%spin=0.0_dp
       structure(ns)%modspin=0.0_dp
       structure(ns)%del=-1.0_dp
       structure(ns)%community=-1
       structure(ns)%num_units=-1
       structure(ns)%dimensionality=-1.0_dp

       ntokens=0
       ctemp=resbuff(nstart)(:indx1-1)
       do while (len_trim(ctemp).gt.0)
          ntokens=ntokens+1
          indx3=index(ctemp," ")
          tokens(ntokens)=trim(ctemp(:indx3))
          ctemp=adjustl(ctemp(indx3+1:))
       end do

       structure(ns)%structure_label="Unknown"
       structure(ns)%pressure=0.0_dp
       structure(ns)%volume=-1.0_dp
       structure(ns)%enthalpy=huge(1.0_dp)
       structure(ns)%spin=0.0_dp
       structure(ns)%modspin=0.0_dp
       structure(ns)%del=-1.0_dp
       structure(ns)%num_ions=-1

       select case(ntokens) 
       case(9)
          read(tokens(2),*,err=334) structure(ns)%structure_label
          read(tokens(3),*,err=334) structure(ns)%pressure
          read(tokens(4),*,err=334) structure(ns)%volume
          read(tokens(5),*,err=334) structure(ns)%enthalpy
          read(tokens(6),*,err=334) structure(ns)%spin
          read(tokens(7),*,err=334) structure(ns)%modspin
          read(tokens(8),*,err=334) structure(ns)%del
          read(tokens(9),*,err=334) structure(ns)%num_ions
334       continue
       case(8)
          read(tokens(2),*,err=335) structure(ns)%structure_label
          read(tokens(3),*,err=335) structure(ns)%pressure
          read(tokens(4),*,err=335) structure(ns)%volume
          read(tokens(5),*,err=335) structure(ns)%enthalpy
          read(tokens(6),*,err=335) structure(ns)%spin
          read(tokens(7),*,err=335) structure(ns)%modspin
          read(tokens(8),*,err=335) structure(ns)%num_ions
335       continue
       case(7)
          read(tokens(2),*,err=336) structure(ns)%structure_label
          read(tokens(3),*,err=336) structure(ns)%pressure
          read(tokens(4),*,err=336) structure(ns)%volume
          read(tokens(5),*,err=336) structure(ns)%enthalpy
          read(tokens(6),*,err=336) structure(ns)%spin
          read(tokens(7),*,err=336) structure(ns)%modspin
336       continue
       case(6)
          read(tokens(2),*,err=337) structure(ns)%structure_label
          read(tokens(3),*,err=337) structure(ns)%pressure
          read(tokens(4),*,err=337) structure(ns)%volume
          read(tokens(5),*,err=337) structure(ns)%enthalpy
337       continue      
       case(5)
          read(tokens(2),*,err=338) structure(ns)%structure_label
          read(tokens(3),*,err=338) structure(ns)%pressure
          read(tokens(4),*,err=338) structure(ns)%volume
          read(tokens(5),*,err=338) structure(ns)%enthalpy
338       continue
       case(4)
          read(tokens(2),*,err=339) structure(ns)%structure_label
          read(tokens(3),*,err=339) structure(ns)%volume
          read(tokens(4),*,err=339) structure(ns)%enthalpy
339       continue
       case default
          continue
       end select

       if(structure(ns)%structure_label.ne."Unknown") then
          structure(ns)%symmetry=resbuff(nstart)(indx1+1:indx2-1)
          if(len_trim(structure(ns)%symmetry).eq.0) structure(ns)%symmetry="-"
          indx1=index(resbuff(nstart),"n - ")
          read(resbuff(nstart)(indx1+3:),*,err=444,end=444) structure(ns)%num_copies
       else
          structure(ns)%symmetry="-"
          structure(ns)%num_copies=1
       end if

       if(structure(ns)%num_ions.ne.nend-nstart-nrem-4) then
          if(structure(ns)%num_ions<0) then
             structure(ns)%num_ions=nend-nstart-nrem-4
          else
             write (*,*) structure(ns)%num_ions,nend,nstart,nrem,num_lines
             do nl=nstart,nend
                write (*,*) trim(resbuff(nl))
             end do
             stop 'num_ions .ne. nend-nstart-4'
          end if
       end if

       allocate(structure(ns)%positions_frac(3,structure(ns)%num_ions))
       allocate(structure(ns)%positions_cart(3,structure(ns)%num_ions))
       allocate(structure(ns)%ion_names(structure(ns)%num_ions))
       allocate(structure(ns)%ion_names_indx(structure(ns)%num_ions))

       ! ** Read the structure data

       read(resbuff(nstart+nrem+1),*,err=340,end=340) ctemp,ctemp,structure(ns)%lattice_abc(1:6)
340    continue
       structure(ns)%lattice_abc(4:6) = structure(ns)%lattice_abc(4:6)*asin(1.0_dp)/90.0_dp

       ni=0
       do nl=nstart+nrem+4,nstart+nrem+4+structure(ns)%num_ions-1
          ni=ni+1
          read(resbuff(nl),*,err=444,end=444) structure(ns)%ion_names(ni),ctemp,structure(ns)%positions_frac(:,ni)

          ! ** Remove numbers from ion name

          if(scan(structure(ns)%ion_names(ni),'1234567890').gt.0) then
             structure(ns)%ion_names(ni)=structure(ns)%ion_names(ni)(:scan(structure(ns)%ion_names(ni),'1234567890')-1)
          end if
          ! ** IS THIS NEEDED? Seems so
          do i=1,3
             structure(ns)%positions_frac(i,ni)=structure(ns)%positions_frac(i,ni)-floor(structure(ns)%positions_frac(i,ni))
          end do
       end do

       ! ** Index the ion_names

       do ni=1,structure(ns)%num_ions
          do i=1,118
             if(structure(ns)%ion_names(ni).eq.elements_alpha(i)) then
                structure(ns)%ion_names_indx(ni)=i
                exit
             end if
          end do
       end do

       ! ** Analyse the formula

       structure(ns)%num_species=0
       name=''
       do ni=1,structure(ns)%num_ions
          if(.not.any(name==structure(ns)%ion_names(ni))) then
             structure(ns)%num_species=structure(ns)%num_species+1
             name(structure(ns)%num_species)=structure(ns)%ion_names(ni)
          end if
       end do

       allocate(structure(ns)%species_names(structure(ns)%num_species))
       allocate(structure(ns)%species_num(structure(ns)%num_species))
       allocate(structure(ns)%species_frac(structure(ns)%num_species))

       structure(ns)%species_names=name(1:structure(ns)%num_species)

       do nl=1,structure(ns)%num_species
          structure(ns)%species_num(nl)=count(structure(ns)%ion_names.eq.structure(ns)%species_names(nl))
       end do

       structure(ns)%num_form = structure(ns)%species_num(1)
       do ni = 2,size(structure(ns)%species_num)
          structure(ns)%num_form = gcd_rec(structure(ns)%num_form,structure(ns)%species_num(ni))
       end do

       structure(ns)%species_num=structure(ns)%species_num/structure(ns)%num_form

       do ni=1,structure(ns)%num_species
          if(structure(ns)%species_num(ni).gt.1) then
             write (name(ni),'(i10)') structure(ns)%species_num(ni)
          else
             name(ni)=''
          end if
       end do

       ! Construct a chemical formula according to the Hill system (with oxygen last)

       have_c=.false.

       if(any(structure(ns)%species_names(:).eq."C  ")) have_c=.true.

       do ni=1,structure(ns)%num_species
          do n=1,118
             if(structure(ns)%species_names(ni).eq.elements_alpha(n)) exit
          end do
          if(structure(ns)%species_names(ni).eq."O  ") then
             element_order(ni)=huge(1.0_dp)
          else if(structure(ns)%species_names(ni).eq."C  ") then
             element_order(ni)=0.0_dp
          else if((structure(ns)%species_names(ni).eq."H  ").and.have_c) then
             element_order(ni)=0.1_dp
          else
             element_order(ni)=real(n,dp)
          end if
          element_index(ni)=ni
       end do

       call heap_sort_index(structure(ns)%num_species,element_order,element_index)

       write(structure(ns)%formula,*) &
            (trim(structure(ns)%species_names(element_index(ni))),&
            trim(adjustl(name(element_index(ni)))),ni=1,structure(ns)%num_species)

       structure(ns)%formula = adjustl(structure(ns)%formula)

       do ni=1,structure(ns)%num_species
          structure(ns)%species_frac(ni)=real(structure(ns)%species_num(ni),dp)/&
               real(sum(structure(ns)%species_num(1:structure(ns)%num_species)),dp)
       end do

       ! ** Add small random number to enthalpy to ensure no exact degeneracy

       structure(ns)%enthalpy=structure(ns)%enthalpy+random_single()*delta
       
       ! ** Add a pressure, update enthalpy

       if(press_add.ne.0.0_dp) then
          structure(ns)%enthalpy=structure(ns)%enthalpy+press_add*structure(ns)%volume/evbyang3
       end if

       ! ** Convert to cartesian

       structure(ns)%lattice_car(:,1) = (/structure(ns)%lattice_abc(1),0.0_dp,0.0_dp/)
       structure(ns)%lattice_car(:,2) = (/structure(ns)%lattice_abc(2)*cos(structure(ns)%lattice_abc(6)),&
            structure(ns)%lattice_abc(2)*sin(structure(ns)%lattice_abc(6)),0.0_dp/)
       structure(ns)%lattice_car(1,3) = structure(ns)%lattice_abc(3)*cos(structure(ns)%lattice_abc(5))
       structure(ns)%lattice_car(2,3) = structure(ns)%lattice_abc(3)*(cos(structure(ns)%lattice_abc(4))-&
            cos(structure(ns)%lattice_abc(5))*cos(structure(ns)%lattice_abc(6)))/sin(structure(ns)%lattice_abc(6))
       structure(ns)%lattice_car(3,3) = sqrt(structure(ns)%lattice_abc(3)**2-structure(ns)%lattice_car(1,3)**2-&
            structure(ns)%lattice_car(2,3)**2)

       ! ** Calculate the volume

       structure(ns)%vol_calc = structure(ns)%lattice_car(1,1)*(structure(ns)%lattice_car(2,2)*structure(ns)%lattice_car(3,3)-&
            structure(ns)%lattice_car(3,2)*structure(ns)%lattice_car(2,3))+&
            structure(ns)%lattice_car(2,1)*(structure(ns)%lattice_car(3,2)*structure(ns)%lattice_car(1,3)-&
            structure(ns)%lattice_car(1,2)*structure(ns)%lattice_car(3,3))+&
            structure(ns)%lattice_car(3,1)*(structure(ns)%lattice_car(1,2)*structure(ns)%lattice_car(2,3)-&
            structure(ns)%lattice_car(2,2)*structure(ns)%lattice_car(1,3))

       if(structure(ns)%volume<0.0_dp) structure(ns)%volume=structure(ns)%vol_calc 

       ! ** Reciprocal lattice

       structure(ns)%lattice_rec(1,1)=structure(ns)%lattice_car(2,2)*structure(ns)%lattice_car(3,3)-&
            structure(ns)%lattice_car(3,2)*structure(ns)%lattice_car(2,3)
       structure(ns)%lattice_rec(2,1)=structure(ns)%lattice_car(2,3)*structure(ns)%lattice_car(3,1)-&
            structure(ns)%lattice_car(3,3)*structure(ns)%lattice_car(2,1)
       structure(ns)%lattice_rec(3,1)=structure(ns)%lattice_car(2,1)*structure(ns)%lattice_car(3,2)-&
            structure(ns)%lattice_car(3,1)*structure(ns)%lattice_car(2,2)
       structure(ns)%lattice_rec(1,2)=structure(ns)%lattice_car(3,2)*structure(ns)%lattice_car(1,3)-&
            structure(ns)%lattice_car(1,2)*structure(ns)%lattice_car(3,3)
       structure(ns)%lattice_rec(2,2)=structure(ns)%lattice_car(3,3)*structure(ns)%lattice_car(1,1)-&
            structure(ns)%lattice_car(1,3)*structure(ns)%lattice_car(3,1)
       structure(ns)%lattice_rec(3,2)=structure(ns)%lattice_car(3,1)*structure(ns)%lattice_car(1,2)-&
            structure(ns)%lattice_car(1,1)*structure(ns)%lattice_car(3,2)
       structure(ns)%lattice_rec(1,3)=structure(ns)%lattice_car(1,2)*structure(ns)%lattice_car(2,3)-&
            structure(ns)%lattice_car(2,2)*structure(ns)%lattice_car(1,3)
       structure(ns)%lattice_rec(2,3)=structure(ns)%lattice_car(1,3)*structure(ns)%lattice_car(2,1)-&
            structure(ns)%lattice_car(2,3)*structure(ns)%lattice_car(1,1)
       structure(ns)%lattice_rec(3,3)=structure(ns)%lattice_car(1,1)*structure(ns)%lattice_car(2,2)-&
            structure(ns)%lattice_car(2,1)*structure(ns)%lattice_car(1,2)

       structure(ns)%lattice_rec(:,:)=structure(ns)%lattice_rec(:,:)/structure(ns)%vol_calc

       ! ** Construct the short structure label

       nlen=len(trim(structure(ns)%structure_label))

       if(nlen.gt.20) then
          structure(ns)%structure_label_short=structure(ns)%structure_label(1:8)//"*"//structure(ns)%structure_label(nlen-8:nlen)
       else
          structure(ns)%structure_label_short=structure(ns)%structure_label(1:20)
       end if

       cycle

444    continue

       if(allocated(structure(ns)%positions_frac)) then
          deallocate(structure(ns)%positions_frac)
          deallocate(structure(ns)%positions_cart)
          deallocate(structure(ns)%ion_names)
          deallocate(structure(ns)%ion_names_indx)
       end if

       ns=ns-1

    end do

    ! ** Set actual number of complete records read

    num_structures=ns

    if(any(structure(:)%modspin.gt.0.0_dp)) have_spin=.true.
    if(any(structure(:)%del.gt.-1.0_dp)) have_dos=.true.

    return

  end subroutine read_res

  subroutine formula_pick()
    
    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       if (structure(ns)%formula.ne.formula) structure(ns)%num_copies=0
    end do
    
  end subroutine formula_pick
  
  subroutine formula_conv()

    ! ** Try to convert the composition to a new one specified in formula_convert

    ! ** Identify the original composition space

    num_species=0
    species_names=''

    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       do n=1,structure(ns)%num_species
          if(all(structure(ns)%species_names(n).ne.species_names(1:num_species))) then
             num_species=num_species+1
             species_names(num_species)=structure(ns)%species_names(n)
          end if
       end do
    end do

    allocate(cmat(num_species,1))

    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       !write (stderr,*) structure(ns)%formula,comp2vec(structure(ns)%formula),species_names(1:num_species)

       cmat(:,1) = comp2vec(structure(ns)%formula)
       do n=1,num_species
          if(cmat(n,1).eq.2) then
             do i=1,structure(ns)%num_species
                if(species_names(n).eq.(structure(ns)%species_names(i))) then
                   structure(ns)%species_names(i)='O'
                   do j=1,structure(ns)%num_ions
                      if(structure(ns)%ion_names(j).eq.species_names(n)) structure(ns)%ion_names(j)='O' 
                   end do
                end if
             end do
          end if

          if(cmat(n,1).eq.1) then
             do i=1,structure(ns)%num_species
                if(species_names(n).eq.(structure(ns)%species_names(i))) then
                   structure(ns)%species_names(i)='Ti'
                   do j=1,structure(ns)%num_ions
                      if(structure(ns)%ion_names(j).eq.species_names(n)) structure(ns)%ion_names(j)='Ti' 
                   end do
                end if
             end do
          end if

       end do


    end do


  end subroutine formula_conv

  subroutine enthalpy_plot()

    if(allocated(composition)) deallocate(composition,numform)
    allocate(composition(num_structures),numform(num_structures),sumbuff(num_structures),energy(num_structures))

    composition=''
    numform=0
    num_compositions=0
    num_total=0
    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       if(cluster) then

          newcomp=.true.
          do nc=1,num_compositions
             if((structure(ns)%formula.eq.composition(nc)).and.(structure(ns)%num_form.eq.numform(nc))) newcomp=.false.
          end do

          if(newcomp) then
             num_compositions=num_compositions+1
             composition(num_compositions) = structure(ns)%formula
             numform(num_compositions) = structure(ns)%num_form
          end if
          num_total = num_total + structure(ns)%num_copies

       else
          if(.not.any(structure(ns)%formula.eq.composition(1:num_compositions))) then
             num_compositions=num_compositions+1
             composition(num_compositions) = structure(ns)%formula
          end if
          num_total = num_total + structure(ns)%num_copies
       end if
    end do

    write (stderr,'(a)') 'Constructing enthalpy plots'

    allocate(struct_name(num_structures),pestruct(num_structures),pressure(num_structures),&
         enthalpy(num_structures))

    do nc=1,num_compositions

       write (stdout,'(2a)') 'Composition: ',composition(nc)

       num_penth=1
       num_pestruc=0

       do ns=1,num_structures
          if(composition(nc).ne.structure(ns)%formula) cycle

          indx1=index(structure(ns)%structure_label,'_')
          indx2=index(structure(ns)%structure_label,'-')
          struct_name(num_penth)=structure(ns)%structure_label(:indx1-1)//structure(ns)%structure_label(indx2:)
          read(structure(ns)%structure_label(indx1+1:indx2-1),*,err=333) ctemp
          if(index(ctemp,'m').gt.0) ctemp(index(ctemp,'m'):index(ctemp,'m'))='-'
          if(index(ctemp,'p').gt.0) ctemp(index(ctemp,'p'):index(ctemp,'p'))='.'
          read(ctemp,*) pressure(num_penth)            
          enthalpy(num_penth)=structure(ns)%enthalpy/structure(ns)%num_form



          if(.not.any(struct_name(num_penth).eq.(pestruct(1:num_pestruc)))) then
             num_pestruc=num_pestruc+1
             pestruct(num_pestruc)=struct_name(num_penth)
          end if

          num_penth=num_penth+1

333       continue
       end do
       num_penth=num_penth-1

       lstsq=least_square(num_penth,pressure,enthalpy)

       do n=1,num_penth
          enthalpy(n)=enthalpy(n)-lstsq(1)*pressure(n)-lstsq(2)
       end do

       sort_array=pressure
       do n=1,num_penth
          sort_index(n)=n
       end do
       call heap_sort_index(num_penth,sort_array,sort_index)           

       allocate(peinterp(num_interp,2,num_pestruc),enref(num_structures,num_pestruc))

       pmin=minval(pressure(1:num_penth))
       pmax=maxval(pressure(1:num_penth))

       do na=1,num_pestruc

          ntrain=0

          do n=1,num_penth
             ns=sort_index(n) 
             if(struct_name(ns).eq.pestruct(na)) ntrain=ntrain+1
          end do

          allocate(ptrain(ntrain),etrain(ntrain,1))

          ntrain=0
          do n=1,num_penth
             ns=sort_index(n)
             if(struct_name(ns).eq.pestruct(na)) then
                ntrain=ntrain+1
                ptrain(ntrain)=pressure(ns)
                etrain(ntrain,1)=enthalpy(ns)
             end if
          end do

          call gp_initialise(ptrain,etrain,ls=lscale*(pmax-pmin),sn=1e-3_dp)

          deallocate(ptrain,etrain)

          do n=1,num_interp
             press=pmin+(pmax-pmin)*real(n-1,dp)/real(num_interp-1,dp)

             peinterp(n,1,na)=press
             peinterp(n,2,na)=gp_predict(press,1)

          end do

          do n=1,num_penth
             ns=sort_index(n) 
             enref(ns,na)=gp_predict(pressure(ns),1)
          end do

       end do

       nref=minloc(peinterp(num_interp/2,2,1:num_pestruc))

       do n=1,num_pestruc
          if(formula_ref.eq.pestruct(n)) then
             nref(1)=n
             exit
          end if
       end do
       
       call write_enthalpy()

       deallocate(peinterp,enref)

    end do

    deallocate(struct_name,pressure,enthalpy)

  end subroutine enthalpy_plot
  
  subroutine compositions()

    allocate(composition(num_structures),comp_enth(num_structures),comp_index(num_structures),comp_num(num_structures),&
         numform(num_structures),sort_index(num_structures),sort_array(num_structures))

    composition=''
    numform=0
    num_compositions=0
    num_total=0
    do ns=1,num_structures
       structure(ns)%composition=0
       do nc=1,num_compositions
          if(cluster) then
             if((structure(ns)%formula.eq.composition(nc)).and.structure(ns)%num_form.eq.numform(nc)) then
                structure(ns)%composition=nc
                exit
             end if
          else
             if(structure(ns)%formula.eq.composition(nc)) then
                structure(ns)%composition=nc
                exit
             end if
          end if
       end do
       if(structure(ns)%composition.eq.0) then
          num_compositions=num_compositions+1
          composition(num_compositions)=structure(ns)%formula
          numform(num_compositions)=structure(ns)%num_form
          structure(ns)%composition=num_compositions
       end if
       num_total = num_total + structure(ns)%num_copies
    end do

    comp_enth=huge(1.0_dp)
    comp_index=-1
    comp_num=0
    do ns=1,num_structures
       if (structure(ns)%num_copies.eq.0) cycle
       if(structure(ns)%enthalpy/structure(ns)%num_form.lt.comp_enth(structure(ns)%composition)) then
          comp_enth(structure(ns)%composition)=structure(ns)%enthalpy/structure(ns)%num_form
          comp_index(structure(ns)%composition)=ns
       end if
       comp_num(structure(ns)%composition)=comp_num(structure(ns)%composition)+structure(ns)%num_copies
    end do

    do ns=1,num_structures
       structure(ns)%relative_enthalpy=structure(ns)%enthalpy/structure(ns)%num_form-comp_enth(structure(ns)%composition)
       if(.not.any(task.eq.'phull').and.&
            (structure(ns)%relative_enthalpy*structure(ns)%num_form/structure(ns)%num_ions.gt.delta_e)) structure(ns)%num_copies=0
    end do

  end subroutine compositions

  subroutine connect()

    ! **************************************************
    ! ** Calculate the connectivity of the structures **
    ! **************************************************

    do ns=1,num_structures

       if(structure(ns)%num_copies>0) then

          if(max(bondlength,bondscale).gt.0.0_dp) then

             n=structure(ns)%num_ions
             l=n*(3+n/2)
             allocate(LL(n,n),E(n),work(l))

             call calculate_distances(ns,rmax,nosort=.true.,laplacian=max(bondlength,bondscale))

             LL=structure(ns)%laplacian

             call dsyev('V','U',n,LL,n,E,work,l,stat)

             structure(ns)%num_units=count(E.lt.delta)

             deallocate(LL,E,work)

             if((num_units.gt.0).and.(num_units.ne.structure(ns)%num_units)) structure(ns)%num_copies=0

             ! ** Calculate Laplacian for 2x2x2 supercell

             n=structure(ns)%num_ions*2**3
             l=n*(3+n/2)
             allocate(LL(n,n),E(n),work(l))

             call calculate_distances(ns,rmax,nosort=.true.,laplacian=max(bondlength,bondscale),superlap=LL)

             call dsyev('V','U',n,LL,n,E,work,l,stat)

             structure(ns)%dimensionality=log(8.0_dp/(real(count(E.lt.delta),dp)/real(structure(ns)%num_units,dp)))/log(2.0_dp)

             deallocate(LL,E,work)

             if((dimensionality.ge.0.0_dp).and.((dimensionality.gt.structure(ns)%dimensionality+delta)&
                  .or.(dimensionality.lt.structure(ns)%dimensionality-delta))) structure(ns)%num_copies=0


          end if
       end if

    end do

  end subroutine connect

  subroutine geometry()

    ! ******************************************************
    ! ** Calculate the atomic geometry for the structures **
    ! ******************************************************

    do ns=1,num_structures

       allocate(structure(ns)%in_unit(structure(ns)%num_ions))
       structure(ns)%in_unit=0

       if(structure(ns)%num_copies>0) then

          if(max(bondlength,bondscale).gt.0.0_dp) then

             n=structure(ns)%num_ions
             l=n*(3+n/2)
             allocate(LL(n,n),MM(n,n),B(n,n),SS(n,n),S(n),So(n),E(n),work(l),iwork(n),iwork2(n))
             allocate(taken(n),moved(n),in(n),dcom(n),dia(n,n),eab(n,n))
             call calculate_distances(ns,rmax,nosort=.true.,laplacian=max(bondlength,bondscale))

             if (modular) then

                structure(ns)%num_units=1
                structure(ns)%in_unit=1

                maxq=-huge(1.0_dp)
                i=0
                nrep=0
                do while((nrep.lt.3).and.(i.lt.10)) ! Can change this.
                   i=i+1
                   call init_modularity(structure(ns)%adjacancy,dmod) 
                   write (stderr,'(a,i3,i3,f10.5,i3,i3,i5)') &
                        "Iteration: ",i,nrep!,dmod,rasym,rhigh,int((1.0_dp-dmod)*rasym)
                   call optimise_rash()
                   if(all(com.eq.structure(ns)%in_unit)) nrep=nrep+1
                   if(q.gt.maxq) then
                      maxq=q
                      structure(ns)%in_unit=com
                      structure(ns)%num_units=nunit
                      write (stderr,'(a,i3,f8.5)') 'Maximum q: ',nunit,maxq
                      i=0
                      nrep=0
                   end if
                end do

                if(comout) then

                   write(ctemp,'(f5.3)') dmod
                   ctemp=trim(adjustl(structure(ns)%structure_label))//"-"//trim(adjustl(ctemp))//".com"
                   open(unit=unit_unit,form="formatted",file=ctemp)

                   write (unit_unit,*) structure(ns)%num_ions,structure(ns)%num_units,dmod,maxq
                   do i=1,structure(ns)%num_ions
                      write (unit_unit,*) structure(ns)%in_unit(i)
                   end do

                   close(unit=unit_unit)

                end if

                if(adj) then

                   ctemp=trim(adjustl(structure(ns)%structure_label))//".adj"
                   open(unit=unit_unit,form="unformatted",file=ctemp)
                   write(unit_unit) structure(ns)%adjacancy
                   close(unit=unit_unit)

                end if

             else

                LL=structure(ns)%laplacian

                call dsyev('V','U',n,LL,n,E,work,l,stat)

                ! ** Identify units

                structure(ns)%num_units=0
                taken=.false.
                do i=1,structure(ns)%num_ions
                   if(E(i).lt.delta) then
                      structure(ns)%num_units=structure(ns)%num_units+1
                      do j=1,structure(ns)%num_ions
                         iwork(j)=j
                      end do
                      call heap_sort_index(structure(ns)%num_ions,LL(:,i),iwork)
                      first=.true.
                      do j=1,structure(ns)%num_ions
                         if(abs(LL(iwork(j),i)).gt.delta) then
                            if(first) then
                               first=.false.
                               temp=0.0_dp
                            else
                               temp=abs(LL(iwork(j),i)-last)
                            end if
                            last=LL(iwork(j),i)

                            if(temp.lt.delta) then
                               if(.not.taken(iwork(j))) then
                                  structure(ns)%in_unit(iwork(j))=structure(ns)%num_units
                                  taken(iwork(j))=.true.
                               else
                                  first=.true.
                               end if
                            else
                               goto 101
                            end if
                         end if
                      end do
                   end if
101                continue
                end do

             end if

             ! ** Store units

             allocate(structure(ns)%num_unit(structure(ns)%num_units))

             do nc=1,structure(ns)%num_units
                structure(ns)%num_unit(nc)=0
                do i=1,structure(ns)%num_ions
                   if(structure(ns)%in_unit(i).eq.nc) structure(ns)%num_unit(nc)= structure(ns)%num_unit(nc)+1
                end do
             end do

             allocate(structure(ns)%units_cart(3,maxval(structure(ns)%num_unit),structure(ns)%num_units))
             allocate(structure(ns)%unit_names(maxval(structure(ns)%num_unit),structure(ns)%num_units))
             allocate(structure(ns)%unit_names_indx(maxval(structure(ns)%num_unit),structure(ns)%num_units))
             allocate(structure(ns)%unit_pointgp(structure(ns)%num_units))
             allocate(structure(ns)%unit_num_asym(structure(ns)%num_units))
             allocate(unit_save(3,maxval(structure(ns)%num_unit)))

             do nc=1,structure(ns)%num_units
                nin=0
                do i=1,structure(ns)%num_ions
                   if(structure(ns)%in_unit(i).eq.nc) then
                      nin=nin+1
                      structure(ns)%units_cart(:,nin,nc)=structure(ns)%positions_cart(:,i)
                      structure(ns)%unit_names(nin,nc)=structure(ns)%ion_names(i)
                      do j=1,118
                         if(structure(ns)%unit_names(nin,nc).eq.elements_alpha(j)) then
                            structure(ns)%unit_names_indx(nin,nc)=j
                            exit
                         end if
                      end do
                   end if
                end do
             end do

             ! ** Join units

             do nc=1,structure(ns)%num_units

                unit_save(:,:)=structure(ns)%units_cart(:,:,nc)
                connected=.false.
                natmpt=0
                do while(.not.connected)
                   if((natmpt.gt.10).or.(structure(ns)%num_unit(nc).le.1)) exit

                   structure(ns)%units_cart(:,:,nc)=unit_save(:,:)
                   call join_unit(structure(ns)%units_cart(:,:,nc),structure(ns)%num_unit(nc),&
                        structure(ns)%unit_names_indx(:,nc))
                   call compact_unit(structure(ns)%units_cart(:,:,nc),structure(ns)%num_unit(nc),&
                        structure(ns)%unit_names_indx(:,nc),ne)

                   ! ** Final check - Is the unit actually connected

                   LL=0.0_dp

                   do i=1,structure(ns)%num_unit(nc)
                      do j=i+1,structure(ns)%num_unit(nc)

                         if(bondscale.gt.0.0_dp) then
                            bl2=bondscale*(elements_radius(structure(ns)%unit_names_indx(i,nc))+&
                                 elements_radius(structure(ns)%unit_names_indx(j,nc)))
                            bl2=bl2**2
                         else
                            bl2=bondlength**2
                         end if

                         vtemp(:) = structure(ns)%units_cart(:,j,nc)-structure(ns)%units_cart(:,i,nc)

                         dist2=dot_product(vtemp,vtemp)

                         if(dist2.lt.bl2) then
                            LL(i,j)=-1.0
                            LL(j,i)=-1.0
                         end if

                      end do
                   end do

                   do i=1,structure(ns)%num_unit(nc)
                      LL(i,i)=-sum(LL(i,1:structure(ns)%num_unit(nc)))
                   end do

                   n=structure(ns)%num_ions
                   l=n*(3+n/2)
                   call dsyev('N','L',structure(ns)%num_unit(nc),LL,n,E,work,l,stat)
                   do i=1,structure(ns)%num_unit(nc)
                      if(E(i).gt.delta) exit
                   end do

                   if(i.le.2) then 
                      connected=.true.
                   else
                      write(stderr,'(a)',advance='no') '.'
                      flush(stderr)
                      natmpt=natmpt+1
                   end if

                end do

                ! ** If not connected - give up, and go back to the original

                if(.not.connected) structure(ns)%units_cart(:,:,nc)=unit_save(:,:)

                ! ** Put centre of mass into unit cell

                R0=0.0_dp
                do i=1,structure(ns)%num_unit(nc)
                   R0(:)=R0(:)+structure(ns)%units_cart(:,i,nc)
                end do

                R0(:) = R0(:)/real(structure(ns)%num_unit(nc),dp)
                R0(:) = floor(matmul(structure(ns)%lattice_rec,R0))
                R0(:) = R0(1)*structure(ns)%lattice_car(:,1) + &
                     R0(2)*structure(ns)%lattice_car(:,2) + R0(3)*structure(ns)%lattice_car(:,3)
                do i=1,structure(ns)%num_unit(nc)
                   structure(ns)%units_cart(:,i,nc) = structure(ns)%units_cart(:,i,nc) - R0(:)
                end do

                if(structure(ns)%num_unit(nc).gt.0) then

                   ! ** Output unit as OFF

                   write(ctemp,*) nc
                   unitfile=trim(adjustl(structure(ns)%structure_label))//"-"//trim(adjustl(ctemp))//".dat"
                   open(unit=unit_unit,form="formatted",file=unitfile)

                   write(ctemp,*) structure(ns)%num_unit(nc),0,0
                   if(alpha.eq.0.0_dp) then
                      write (unit_unit,'(a)') "OFF"
                      write (unit_unit,'(a)') trim(adjustl(ctemp))
                   end if
                   do i=1,structure(ns)%num_unit(nc)
                      write (unit_unit,'(3f15.7)') structure(ns)%units_cart(:,i,nc)
                   end do

                   close(unit=unit_unit)

                   ! ** Generate the convex hulls

                   ctemp4=adjustl(unitfile(:index(unitfile,".")-1))
                   if(alpha.gt.0.0_dp) then
                      write (ctemp2,*) alpha
                      ctemp="hull -aa"//trim(adjustl(ctemp2))//" -A -afoff -oN -oF"&
                           //trim(ctemp4)&
                           //" < "//trim(adjustl(unitfile))//" 2> /dev/null"
                      call system(ctemp)
                   else if(alpha.lt.0.0_dp) then
                      ctemp="hull -A -afoff -oN -oF"//trim(ctemp4)&
                           //" < "//trim(adjustl(unitfile))//" 2> /dev/null"
                      call system(ctemp)
                   else
                      ctemp="conv_hull -Q A0.9999 "//trim(adjustl(unitfile))//" > "//trim(adjustl(unitfile))//".ch 2> /dev/null"
                      call system(ctemp)
                      ctemp="mv "//trim(adjustl(unitfile))//".ch "//trim(adjustl(unitfile(:index(unitfile,".dat")-1)))//".off"
                      call system(ctemp)
                   end if

                   ! ** Write the XYZ files

                   write(ctemp,*) nc
                   unitfile=trim(adjustl(structure(ns)%structure_label))//"-"//trim(adjustl(ctemp))//".xyz"
                   open(unit=unit_unit,form="formatted",file=unitfile)

                   write(ctemp,*) structure(ns)%num_unit(nc)
                   write (unit_unit,'(a)') trim(adjustl(ctemp))
                   write (unit_unit,'(a)') "XYZ"

                   vtemp=vtemp/real(structure(ns)%num_unit(nc),dp)
                   do i=1,structure(ns)%num_unit(nc)
                      write (unit_unit,'(a,3f15.7)') trim(structure(ns)%unit_names(i,nc)),structure(ns)%units_cart(:,i,nc)!-vtemp(:)
                   end do

                   close(unit=unit_unit)                    

                   ! ** Use external program to find point group, and count number in asymmetric unit

                   if(notsymm) then
                      structure(ns)%unit_pointgp(nc)="Unit"
                      structure(ns)%unit_num_asym(nc)=structure(ns)%num_unit(nc)                    
                   else
                      if(structure(ns)%num_unit(nc).gt.2) then
                         ctemp="csymm "//trim(adjustl(unitfile))//" > "//trim(adjustl(unitfile))//".pg"
                         call system(ctemp)
                         open(unit=unit_unit,form="formatted",file=trim(adjustl(unitfile))//".pg")
                         read(unit_unit,*) structure(ns)%unit_pointgp(nc)
                         if(structure(ns)%unit_pointgp(nc).eq."C1") then
                            structure(ns)%unit_num_asym(nc)=structure(ns)%num_unit(nc)
                         else
                            read(unit_unit,*) structure(ns)%unit_num_asym(nc)
                         end if
                         close(unit=unit_unit)
                         ctemp="rm "//trim(adjustl(unitfile))//".pg"
                         call system(ctemp)
                      else if(structure(ns)%num_unit(nc).eq.2) then
                         if(structure(ns)%unit_names(1,nc).eq.structure(ns)%unit_names(2,nc)) then
                            structure(ns)%unit_pointgp(nc)='D(inf)h'
                            structure(ns)%unit_num_asym(nc)=1
                         else
                            structure(ns)%unit_pointgp(nc)='C(inf)v'
                            structure(ns)%unit_num_asym(nc)=2
                         end if
                      else if(structure(ns)%num_unit(nc).eq.1) then
                         structure(ns)%unit_pointgp(nc)="R3"
                         structure(ns)%unit_num_asym(nc)=1
                      end if
                   end if

                   ! ** Tidy up

                   write(ctemp,*) nc
                   if(.not.xyz) then
                      unitfile=trim(adjustl(structure(ns)%structure_label))//"-"//trim(adjustl(ctemp))//".{xyz,dat}"
                   else
                      unitfile=trim(adjustl(structure(ns)%structure_label))//"-"//trim(adjustl(ctemp))//".dat"
                   end if
                   ctemp="rm -f "//trim(adjustl(unitfile))
                   call system(ctemp)

                end if

             end do

             ctemp="off_util "//trim(adjustl(structure(ns)%structure_label))//"-*.off 2>/dev/null | off_color -f K 2>/dev/null > "&
                  //trim(adjustl(structure(ns)%structure_label))//".off"
             call system(ctemp)
             ctemp="rm -f "//trim(adjustl(structure(ns)%structure_label))//"-*.off"
             call system(ctemp)
             if(.not.off) then
                ctemp="rm -f "//trim(adjustl(structure(ns)%structure_label))//".off"
                call system(ctemp)
             end if

             deallocate(moved,unit_save)

          else

             ! ** Bondlength zero case

             structure(ns)%num_units=structure(ns)%num_ions

             allocate(structure(ns)%num_unit(structure(ns)%num_units))

             structure(ns)%num_unit=1

             allocate(structure(ns)%units_cart(3,maxval(structure(ns)%num_unit),structure(ns)%num_units))
             allocate(structure(ns)%unit_names(maxval(structure(ns)%num_unit),structure(ns)%num_units))
             allocate(structure(ns)%unit_pointgp(structure(ns)%num_units))
             allocate(structure(ns)%unit_num_asym(structure(ns)%num_units))

             structure(ns)%unit_pointgp='R3'
             structure(ns)%unit_num_asym=1

             do nc=1,structure(ns)%num_units
                structure(ns)%in_unit(nc)=nc
             end do

             do na=1,structure(ns)%num_ions

                structure(ns)%positions_cart(:,na) = &
                     structure(ns)%positions_frac(1,na)*structure(ns)%lattice_car(:,1)+&
                     structure(ns)%positions_frac(2,na)*structure(ns)%lattice_car(:,2)+&
                     structure(ns)%positions_frac(3,na)*structure(ns)%lattice_car(:,3)

             end do

             do nc=1,structure(ns)%num_units
                nin=0
                do i=1,structure(ns)%num_ions
                   if(structure(ns)%in_unit(i).eq.nc) then
                      nin=nin+1
                      structure(ns)%units_cart(:,nin,nc)=structure(ns)%positions_cart(:,i)
                      structure(ns)%unit_names(nin,nc)=structure(ns)%ion_names(i) 
                   end if
                end do
             end do

             do nc=1,structure(ns)%num_units
                R0=0.0_dp
                do i=1,structure(ns)%num_unit(nc)
                   R0(:)=R0(:)+structure(ns)%units_cart(:,i,nc)
                end do

                R0(:) = R0(:)/real(structure(ns)%num_unit(nc),dp)
                R0(:) = floor(matmul(structure(ns)%lattice_rec,R0))
                R0(:) = R0(1)*structure(ns)%lattice_car(:,1) + &
                     R0(2)*structure(ns)%lattice_car(:,2) + R0(3)*structure(ns)%lattice_car(:,3)
                do i=1,structure(ns)%num_unit(nc)
                   structure(ns)%units_cart(:,i,nc) = structure(ns)%units_cart(:,i,nc) - R0(:)
                end do
             end do

          end if

          ! ** Compare units

          allocate(clustdist(maxval(structure(ns)%num_unit)**2,structure(ns)%num_units))
          allocate(clustsim(structure(ns)%num_units,structure(ns)%num_units))
          allocate(clustocc(structure(ns)%num_units))

          ! * Calculate distance fingerprints

          clustdist=huge(1.0_dp)

          do na=1,structure(ns)%num_units

             n=0
             do i=1,structure(ns)%num_unit(na)
                do j=1,structure(ns)%num_unit(na)
                   n=n+1
                   vtemp(:)=structure(ns)%units_cart(:,i,na)-structure(ns)%units_cart(:,j,na)

                   clustdist(n,na)=sqrt(dot_product(vtemp,vtemp))

                end do
             end do

             call heap_sort(n,clustdist(:,na))

          end do

          ! * Compare fingerprints

          clustsim=0

          do na=1,structure(ns)%num_units
             clustsim(na,na)=1
             do nb=1,structure(ns)%num_units                    

                if(structure(ns)%num_unit(na).ne.structure(ns)%num_unit(nb)) cycle
                n=structure(ns)%num_unit(na)**2
                if(maxval(abs(clustdist(1:n,na)-clustdist(1:n,nb))).lt.thresh) clustsim(na,nb)=1

             end do
          end do

          ! * Check composition

          do na=1,structure(ns)%num_units
             do nb=1,structure(ns)%num_units
                if(structure(ns)%num_unit(na).ne.structure(ns)%num_unit(nb)) then
                   clustsim(na,nb)=0
                   cycle
                end if

                if(structure(ns)%num_unit(na).eq.1) then

                   if(structure(ns)%unit_names(1,na).ne.structure(ns)%unit_names(1,nb)) then
                      clustsim(na,nb)=0
                      cycle  
                   end if

                else

                   if(structure(ns)%unit_pointgp(na).ne.structure(ns)%unit_pointgp(nb)) then
                      clustsim(na,nb)=0
                      cycle  
                   end if

                   ! ** largely ignoring this possibility for the moment

                end if

             end do
          end do

          ! * Reduce the cell

          clustocc=0
          do na=1,structure(ns)%num_units
             clustocc(na) = sum(clustsim(:,na))                
             do nb=1,structure(ns)%num_units
                if(clustsim(na,nb).gt.0) clustsim(:,nb)=0
             end do
          end do

          n = clustocc(1)
          do na = 2,structure(ns)%num_units
             if(clustocc(na)==0) cycle
             n = gcd_rec(n,clustocc(na))
          end do
          clustocc=clustocc/n

          nb=0
          do nc=1,structure(ns)%num_units
             if(clustocc(nc)==0) cycle
             nb=nb+structure(ns)%num_unit(nc)                 
          end do

          vscale= real(nb,dp)/real(structure(ns)%num_ions,dp)

          deallocate(clustdist,clustsim)

          ! ** Calculation of the degrees of freedom

          dof=-6
          if(structure(ns)%num_ions.eq.1) dof=dof+3
          if(structure(ns)%num_ions.eq.2) dof=dof+1

          do nc=1,structure(ns)%num_units
             if(structure(ns)%num_unit(nc).gt.2) then
                dof=dof+6
             else if(structure(ns)%num_unit(nc).eq.2) then
                dof=dof+5
             else if(structure(ns)%num_unit(nc).eq.1) then
                dof=dof+3
             end if
             if(clustocc(nc)==0) cycle
             dof=dof+3*structure(ns)%unit_num_asym(nc)-6
             if(structure(ns)%unit_num_asym(nc).eq.1) then 
                dof=dof+3
             else if(structure(ns)%unit_num_asym(nc).eq.2) then 
                dof=dof+1
             end if
          end do



          write (stderr,*) 
          write (stderr,'(a,i5,f7.3)') "Degrees of freedom: ",dof,real(dof,dp)/real(3*structure(ns)%num_ions-5,dp)
          write (stderr,*) 
          write (stderr,'(a,i5,f7.3)') "Number of units:    ",structure(ns)%num_units
          write (stderr,*)  


          ! ** Output units

          write (stdout,'(a)') '%BLOCK LATTICE_CART'
          write (stdout,'(3f10.5)') structure(ns)%lattice_car(:,1)
          write (stdout,'(3f10.5)') structure(ns)%lattice_car(:,2)
          write (stdout,'(3f10.5)') structure(ns)%lattice_car(:,3)
          write (stdout,'(a)') '%ENDBLOCK LATTICE_CART'
          write (stdout,*) 
          write (ctemp,'(f10.2)') structure(ns)%volume*vscale
          write (stdout,'(a)') '#TARGVOL='//trim(adjustl(ctemp))
          write (stdout,*) 
          write (stdout,'(a)') '%BLOCK POSITIONS_ABS'
          nb=0
          do nc=1,structure(ns)%num_units
             if(clustocc(nc)==0) cycle
             nb=nb+1
             write(ctemp,'(i5)') clustocc(nc)
             ctemp=' % NUM='//trim(adjustl(ctemp))
             write(ctemp2,'(i5)') nb
             ctemp2=trim(adjustl(ctemp2))//"-"//trim(adjustl(structure(ns)%unit_pointgp(nc)))
             do i=1,structure(ns)%num_unit(nc)
                write (stdout,'(a5,3f10.5,a,a,a)') &
                     trim(structure(ns)%unit_names(i,nc)),structure(ns)%units_cart(:,i,nc),' # ',&
                     trim(ctemp2),trim(ctemp)
                ctemp=''
             end do
          end do
          write (stdout,'(a)') '%ENDBLOCK POSITIONS_ABS'
          write (stdout,*) 
          write (stdout,'(a)') '#SYMMOPS=1'
          write (stdout,'(a)') '##SGRANK=20'
          write(ctemp,'(i3)') n
          write (stdout,'(a)') '#NFORM='//trim(adjustl(ctemp))
          write (stdout,'(a)') '#SLACK=0.25'
          write (stdout,'(a)') '#OVERLAP=0.1'
          write (stdout,'(a)') '#COMPACT'

          if(max(bondlength,bondscale).gt.0.0_dp) then

             deallocate(LL,MM,SS,B,E,S,So,work,iwork,iwork2,taken,in,dcom,dia,eab)

          end if

          deallocate(clustocc)


       end if

       ! ** Update ionic positions

       do nc=1,structure(ns)%num_units
          nin=0
          do i=1,structure(ns)%num_ions
             if(structure(ns)%in_unit(i).eq.nc) then
                nin=nin+1
                structure(ns)%positions_cart(:,i)=structure(ns)%units_cart(:,nin,nc)
                structure(ns)%positions_frac(:,i)=matmul(structure(ns)%lattice_rec,structure(ns)%units_cart(:,nin,nc))
             end if
          end do
       end do


    end do

    ! ** Identify the composition space

    num_species=0
    species_names=''

    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       do n=1,structure(ns)%num_species
          if(all(structure(ns)%species_names(n).ne.species_names(1:num_species))) then
             num_species=num_species+1
             species_names(num_species)=structure(ns)%species_names(n)
          end if
       end do
    end do

    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle

       call calculate_distances(ns,rmax,nosort=.true.)

       write (stdout,'("#MINSEP=1.0 ")',advance='no')

       do na=1,num_species
          do nb=na,num_species

             do n=1,structure(ns)%num_dist
                if(structure(ns)%distall(structure(ns)%dist_indx(n)).lt.delta) cycle
                n1=structure(ns)%in_unit(structure(ns)%dist_ions(1,structure(ns)%dist_indx(n)))
                n2=structure(ns)%in_unit(structure(ns)%dist_ions(2,structure(ns)%dist_indx(n)))                    
                if((n1.eq.n2).and.structure(ns)%dist_incell(structure(ns)%dist_indx(n))) cycle
                if((species_names(na).eq.structure(ns)%dist_names(1,structure(ns)%dist_indx(n)))&
                     .and.(species_names(nb).eq.structure(ns)%dist_names(2,structure(ns)%dist_indx(n)))) exit
             end do

             if(structure(ns)%distall(structure(ns)%dist_indx(n)).ge.10.0_dp) then
                write (stdout,'(a,"-",a,"=",f4.1," ")',advance="no") &
                     trim(species_names(na)),trim(species_names(nb)),structure(ns)%distall(structure(ns)%dist_indx(n))
             else
                write (stdout,'(a,"-",a,"=",f4.2," ")',advance="no") &
                     trim(species_names(na)),trim(species_names(nb)),structure(ns)%distall(structure(ns)%dist_indx(n))
             end if

          end do
       end do
       write (stdout,*)
    end do

  end subroutine geometry

  subroutine rank()

    ! ***********************************************
    ! ** Rank the structures according to enthalpy **
    ! ***********************************************

    n=0
    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       n=n+1
       sort_array(n)=structure(ns)%enthalpy/structure(ns)%num_form
       sort_index(n)=ns
    end do

    call heap_sort_index(n,sort_array,sort_index)

    do ns=1,min(num_top,n)

       if(structure(sort_index(ns))%num_copies>0) then
          if(ns.eq.1) then
             enth=structure(sort_index(ns))%enthalpy/structure(sort_index(ns))%num_form
          else
             enth=structure(sort_index(ns))%relative_enthalpy
          end if
          ctemp=''
          if(have_spin) then
             write(ctemp,'(2f6.2)') structure(sort_index(ns))%spin/structure(sort_index(ns))%num_form,&
                  structure(sort_index(ns))%modspin/structure(sort_index(ns))%num_form
          end if
          ctemp2=''
          if(have_dos) then
             if(structure(sort_index(ns))%del.ge.0.0_dp) then
                write(ctemp2,'(f6.3)') structure(sort_index(ns))%del/structure(sort_index(ns))%num_form
             else
                ctemp2='      '
             end if
          end if

          if(long) then
             fmt='(a40,f9.2,f10.3,f12.3,2a,i4,1x,a12,1x,a7,i5)'
             label=structure(sort_index(ns))%structure_label
          else
             fmt='(a20,f9.2,f10.3,f12.3,2a,i4,1x,a12,1x,a7,i5)'
             label=structure(sort_index(ns))%structure_label_short
          end if

          write (stdout,fmt,advance='NO') &
               label,&
               structure(sort_index(ns))%pressure,&
               structure(sort_index(ns))%volume/structure(sort_index(ns))%num_form,&
               enth,&
               trim(ctemp),&
               trim(ctemp2),&
               structure(sort_index(ns))%num_form,&
               structure(sort_index(ns))%formula,&
               structure(sort_index(ns))%symmetry,&
               structure(sort_index(ns))%num_copies

          if(structure(sort_index(ns))%community.gt.0) then
             write (stdout,'(i5)',advance='NO') structure(sort_index(ns))%community
          end if

          if(structure(sort_index(ns))%dimensionality.ge.0.0_dp) then
             write (stdout,'(i5,f5.2)') structure(sort_index(ns))%num_units,structure(sort_index(ns))%dimensionality
          else
             write (stdout,*)
          end if

       end if
    end do


  end subroutine rank

  subroutine summary()

    ! *********************************************************************************
    ! ** Summarise the results - print the most stable structure at each composition **
    ! *********************************************************************************

    allocate(sumbuff(num_structures),energy(num_structures))

    ! ** Identify the composition space - code duplication with geometry

    num_species=0
    species_names=''

    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       do n=1,structure(ns)%num_species
          if(all(structure(ns)%species_names(n).ne.species_names(1:num_species))) then
             num_species=num_species+1
             species_names(num_species)=structure(ns)%species_names(n)
          end if
       end do
    end do

    allocate(cvec(num_species,num_compositions),disclap(num_compositions))

    disclap=huge(1.0_dp)

    ! ** Construct the output string for each composition

    do nc=1,num_compositions


       indx=comp_index(nc)

       ctemp=''
       if(have_spin) then
          write(ctemp,'(2f6.2)') structure(indx)%spin/structure(indx)%num_form,&
               structure(indx)%modspin/structure(indx)%num_form
       end if
       ctemp2=''
       if(have_dos) then
          if(structure(indx)%del.ge.0.0_dp) then
             write(ctemp2,'(f6.3)') structure(indx)%del/structure(indx)%num_form
          else
             ctemp2='      '
          end if
       end if

       energy(nc)=structure(indx)%enthalpy

       if(cluster) then

          ! ** Problem here with multiple species - FIX
          write (ctemp3,'(10i4)') comp2vec(composition(nc))*numform(nc)

          cvec(:,nc)=comp2vec(composition(nc))*numform(nc)

          if(long) then
             fmt='(a40,f12.3,1x,3a,1x,a12,1x,a7,i5,i6)'
             label=structure(indx)%structure_label
          else
             fmt='(a20,f12.3,1x,3a,1x,a12,1x,a7,i5,i6)'
             label=structure(indx)%structure_label_short
          end if

          write (sumbuff(nc),fmt) &
               label,&
               structure(indx)%enthalpy,&
               trim(ctemp),&
               trim(ctemp2),&
               trim(ctemp3),&
               structure(indx)%formula,&
               structure(indx)%symmetry,&
               structure(indx)%num_copies,&
               comp_num(nc)

       else

          if(long) then
             fmt='(a40,f9.2,f10.3,f12.3,2a,i4,1x,a12,1x,a7,i5,i6)'
             label=structure(indx)%structure_label
          else
             fmt='(a20,f9.2,f10.3,f12.3,2a,i4,1x,a12,1x,a7,i5,i6)'
             label=structure(indx)%structure_label_short
          end if

          write (sumbuff(nc),fmt) &
               label,&
               structure(indx)%pressure,&
               structure(indx)%volume/structure(indx)%num_form,&
               structure(indx)%enthalpy/structure(indx)%num_form,&
               trim(ctemp),&
               trim(ctemp2),&
               structure(indx)%num_form,&
               structure(indx)%formula,&
               structure(indx)%symmetry,&
               structure(indx)%num_copies,&
               comp_num(nc)

       end if
    end do

    ! ** Calculate the discrete Laplacian for clusters

    if(cluster) then
       select case(num_species)
       case(1)
          disclap=huge(1.0_dp)

          do nc=1,num_compositions

             emat=huge(1.0_dp)
             do na=1,num_compositions

                if(cvec(1,na).eq.cvec(1,nc)-1) emat(-1,0)=energy(na)
                if(cvec(1,na).eq.cvec(1,nc)) emat(0,0)=energy(na)
                if(cvec(1,na).eq.cvec(1,nc)+1) emat(1,0)=energy(na)

             end do

             if(cvec(1,nc).eq.1) emat(-1,0)=0.0_dp

             if(all(emat(:,0).lt.huge(1.0_dp))) then

                kern(-1,0)=1.0_dp
                kern(0,0)=-2.0_dp
                kern(1,0)=1.0_dp

                disclap(nc)=0.0_dp
                do i=-1,1
                   disclap(nc)=disclap(nc)+emat(i,0)*kern(i,0)
                end do

             end if

          end do

       case(2)
          disclap=huge(1.0_dp)
          do nc=1,num_compositions

             emat=huge(1.0_dp)
             do na=1,num_compositions

                if((cvec(1,na).eq.cvec(1,nc)-1).and.(cvec(2,na).eq.cvec(2,nc)-1)) emat(-1,-1)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc)-1)) emat(0,-1)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)+1).and.(cvec(2,na).eq.cvec(2,nc)-1)) emat(1,-1)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)-1).and.(cvec(2,na).eq.cvec(2,nc))) emat(-1,0)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc))) emat(0,0)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)+1).and.(cvec(2,na).eq.cvec(2,nc))) emat(1,0)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)-1).and.(cvec(2,na).eq.cvec(2,nc)+1)) emat(-1,+1)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc)+1)) emat(0,+1)=energy(na)
                if((cvec(1,na).eq.cvec(1,nc)+1).and.(cvec(2,na).eq.cvec(2,nc)+1)) emat(1,+1)=energy(na)

             end do

             emat(1,1)=0.0_dp
             emat(1,-1)=0.0_dp
             emat(-1,1)=0.0_dp
             emat(-1,-1)=0.0_dp

             if((cvec(1,nc).eq.0).or.(cvec(2,nc).eq.0)) then
                if(cvec(1,nc).eq.0) then
                   emat(-1,0)=2*emat(0,0)-emat(1,0)
                   if(cvec(2,nc).eq.1) emat(0,-1)=0.0_dp
                endif
                if(cvec(2,nc).eq.0) then
                   emat(0,-1)=2*emat(0,0)-emat(0,1)
                   if(cvec(1,nc).eq.1) emat(-1,0)=0.0_dp                       
                end if
             end if

             if(all(emat.lt.huge(1.0_dp))) then

                kern(-1,-1)=0.0_dp
                kern(-1,0)=1.0_dp
                kern(-1,1)=0.0_dp
                kern(0,-1)=1.0_dp
                kern(0,0)=-4.0_dp
                kern(0,1)=1.0_dp
                kern(1,-1)=0.0_dp
                kern(1,0)=1.0_dp
                kern(1,1)=0.0_dp

                disclap(nc)=0.0_dp
                do i=-1,1
                   do j=-1,1
                      disclap(nc)=disclap(nc)+emat(i,j)*kern(i,j)
                   end do
                end do

             end if

          end do

       case(3)

          disclap=huge(1.0_dp)
          do nc=1,num_compositions

             emat3=0.0_dp

             emat3(0,-1,0)=huge(1.0_dp)
             emat3(-1,0,0)=huge(1.0_dp)
             emat3(0,0,0)=huge(1.0_dp)
             emat3(1,0,0)=huge(1.0_dp)
             emat3(0,1,0)=huge(1.0_dp)
             emat3(0,0,-1)=huge(1.0_dp)
             emat3(0,0,1)=huge(1.0_dp)

             do na=1,num_compositions

                if(cvec(3,na).eq.cvec(3,nc)) then
                   if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc)-1)) emat3(0,-1,0)=energy(na)
                   if((cvec(1,na).eq.cvec(1,nc)-1).and.(cvec(2,na).eq.cvec(2,nc))) emat3(-1,0,0)=energy(na)
                   if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc))) emat3(0,0,0)=energy(na)
                   if((cvec(1,na).eq.cvec(1,nc)+1).and.(cvec(2,na).eq.cvec(2,nc))) emat3(1,0,0)=energy(na)
                   if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc)+1)) emat3(0,1,0)=energy(na)
                end if

                if(cvec(3,na).eq.cvec(3,nc)-1) then
                   if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc))) emat3(0,0,-1)=energy(na)
                end if

                if(cvec(3,na).eq.cvec(3,nc)+1) then
                   if((cvec(1,na).eq.cvec(1,nc)).and.(cvec(2,na).eq.cvec(2,nc))) emat3(0,0,1)=energy(na)
                end if

             end do

!!$                 if((cvec(1,nc).eq.0).or.(cvec(2,nc).eq.0)) then
!!$                    if(cvec(1,nc).eq.0) then
!!$                       emat(-1,0)=2*emat(0,0)-emat(1,0)
!!$                       if(cvec(2,nc).eq.1) emat(0,-1)=0.0_dp
!!$                    endif
!!$                    if(cvec(2,nc).eq.0) then
!!$                       emat(0,-1)=2*emat(0,0)-emat(0,1)
!!$                       if(cvec(1,nc).eq.1) emat(-1,0)=0.0_dp                       
!!$                    end if
!!$                 end if

             if(all(emat3.lt.huge(1.0_dp))) then

                kern3=0.0_dp

                kern3(-1,0,0)=1.0_dp
                kern3(0,-1,0)=1.0_dp
                kern3(0,0,0)=-6.0_dp
                kern3(0,1,0)=1.0_dp
                kern3(1,0,0)=1.0_dp

                kern3(0,0,-1)=1.0_dp

                kern3(0,0,1)=1.0_dp


                disclap(nc)=0.0_dp
                do i=-1,1
                   do j=-1,1
                      do k=-1,1
                         disclap(nc)=disclap(nc)+emat3(i,j,k)*kern3(i,j,k)
                      end do
                   end do
                end do

             end if

          end do
       case default
          disclap=huge(1.0_dp)
       end select

    end if

    ! ** Sort the output

    do nc=1,num_compositions
       sort_array(nc)=-energy(nc)
       sort_index(nc)=nc
    end do
    call heap_sort_index(num_compositions,sort_array,sort_index)

    do nc=1,num_compositions
       if(disclap(sort_index(nc)).lt.huge(1.0_dp)) then
          write (stdout,'(a,f10.3)') trim(sumbuff(sort_index(nc))),disclap(sort_index(nc))
       else
          write (stdout,'(a,a7)') trim(sumbuff(sort_index(nc))),'      ~'             
       end if
    end do

    deallocate(composition,numform,sumbuff,energy,disclap)

    write (stderr,'(a,i7)') 'Number of structures   :',num_total
    write (stderr,'(a,i7)') 'Number of compositions :',num_compositions

  end subroutine summary

  subroutine compare()

    ! *********************************************************
    ! ** Compare the specified structure with all the others **
    ! *********************************************************

    do seed_indx=1,num_structures
       if(comp_structure .eq. structure(seed_indx)%structure_label) exit
    end do

    if(seed_indx.gt.num_structures) stop 'compare: structure not found'

    call calculate_distances(seed_indx,rmax)

    minda=minval(structure(seed_indx)%distall,mask=structure(seed_indx)%distall.gt.delta)


    allocate(maxdiff(num_structures,num_structures))
    maxdiff = 0.0_dp

    do nb=1,num_structures

       if(sloppy) then
          if(structure(seed_indx)%num_species.ne.structure(nb)%num_species) then
             maxdiff(nb,seed_indx)=huge(1.0_dp)
             cycle
          end if
          if(product(structure(seed_indx)%species_num).ne.product(structure(nb)%species_num)) then
             maxdiff(nb,seed_indx)=huge(1.0_dp)
             cycle
          end if
          if(sum(structure(seed_indx)%species_num).ne.sum(structure(nb)%species_num)) then
             maxdiff(nb,seed_indx)=huge(1.0_dp)
             cycle
          end if

       else

          if(structure(nb)%formula.ne.structure(seed_indx)%formula) then
             maxdiff(nb,seed_indx)=huge(1.0_dp)
             cycle
          end if

       end if

       if(nb.eq.seed_indx) then
          maxdiff(nb,seed_indx)=huge(1.0_dp)
          cycle
       end if

       call calculate_distances(nb,rmax)

       mindb=minval(structure(nb)%distall,mask=structure(nb)%distall.gt.delta)

       scale_a=((structure(nb)%volume/structure(nb)%num_form+structure(seed_indx)%volume/structure(seed_indx)%num_form)/&
            (2.0_dp*structure(seed_indx)%volume/structure(seed_indx)%num_form))**(1.0_dp/3.0_dp)
       scale_b=((structure(nb)%volume/structure(nb)%num_form+structure(seed_indx)%volume/structure(seed_indx)%num_form)/&
            (2.0_dp*structure(nb)%volume/structure(nb)%num_form))**(1.0_dp/3.0_dp)


       do n=1,min(structure(nb)%num_dist*structure(seed_indx)%num_form,structure(seed_indx)%num_dist*structure(nb)%num_form)
          if(structure(seed_indx)%distall((n-1)/structure(nb)%num_form+1)<rmax/1.75_dp) then
             diff=abs(structure(nb)%distall((n-1)/structure(seed_indx)%num_form+1)*scale_b-&
                  structure(seed_indx)%distall((n-1)/structure(nb)%num_form+1)*scale_a)
             if(diff>maxdiff(nb,seed_indx)) maxdiff(nb,seed_indx) = diff
          end if
       end do

    end do

    do nb=1,num_structures
       if(maxdiff(nb,seed_indx).lt.comp_threshold) then

          enth=structure(nb)%enthalpy/structure(nb)%num_form
          ctemp=''
          if(have_spin) then
             write(ctemp,'(2f6.2)') structure(nb)%spin/structure(nb)%num_form,&
                  structure(nb)%modspin/structure(nb)%num_form
          end if
          ctemp2=''
          if(have_dos) then
             if(structure(nb)%del.ge.0.0_dp) then
                write(ctemp2,'(f6.3)') structure(nb)%del/structure(nb)%num_form
             else
                ctemp2='      '
             end if
          end if

          if(long) then
             fmt='(a40,f9.2,f10.3,f12.3,2a,i4,1x,a12,1x,a7,i5,f10.4)'
             label=structure(nb)%structure_label
          else
             fmt='(a20,f9.2,f10.3,f12.3,2a,i4,1x,a12,1x,a7,i5,f10.4)'
             label=structure(nb)%structure_label_short
          end if

          write (stdout,fmt) &
               label,&
               structure(nb)%pressure,&
               structure(nb)%volume/structure(nb)%num_form,&
               enth,&
               trim(ctemp),&
               trim(ctemp2),&
               structure(nb)%num_form,&
               structure(nb)%formula,&
               structure(nb)%symmetry,&
               structure(nb)%num_copies,&
               maxdiff(nb,seed_indx)

       end if

    end do

    write (stderr,'(a,i6)') "Matches: ",count(maxdiff(:,seed_indx).lt.comp_threshold)

    deallocate(maxdiff)

  end subroutine compare

  subroutine eliminate()

    ! **********************************************
    ! ** Identify and remove duplicate structures **
    ! **********************************************

    allocate(similar(num_structures),orig_copies(num_structures))

    do ns=1,num_structures
       sort_array(ns)=structure(ns)%enthalpy/structure(ns)%num_form
       sort_index(ns)=ns
    end do

    call heap_sort_index(num_structures,sort_array,sort_index)

    do na=1,num_structures
       nas=sort_index(na)
       if((structure(nas)%num_copies.gt.0).and.(elim_threshold.gt.0.0_dp)) then
          call calculate_distances(nas,rmax)           
          orig_copies(nas)=structure(nas)%num_copies
       else
          orig_copies(nas)=1
       end if
    end do

    ! ** SORT FIRST ***

    if(elim_threshold.gt.0.0_dp) then

       do na=1,num_structures
          nas=sort_index(na)
          if(structure(nas)%num_copies==0) cycle
          minda=minval(structure(nas)%distall,mask=structure(nas)%distall.gt.delta)

          similar = 1

          do nb=na+1,num_structures
             nbs=sort_index(nb)
             if(structure(nbs)%num_copies==0) cycle
             mindb=minval(structure(nbs)%distall,mask=structure(nbs)%distall.gt.delta)

             if(sloppy) then
                if(structure(nas)%num_species.ne.structure(nbs)%num_species) then
                   similar(nbs)=0
                   cycle
                end if
                if(product(structure(nas)%species_num).ne.product(structure(nbs)%species_num)) then
                   similar(nbs)=0
                   cycle
                end if
                if(sum(structure(nas)%species_num).ne.sum(structure(nbs)%species_num)) then
                   similar(nbs)=0
                   cycle
                end if

             else
                if(structure(nas)%formula.ne.structure(nbs)%formula) then
                   similar(nbs)=0
                   cycle
                end if
             end if

             scale_a=((structure(nbs)%volume/structure(nbs)%num_form+structure(nas)%volume/structure(nas)%num_form)/&
                  (2.0_dp*structure(nas)%volume/structure(nas)%num_form))**(1.0_dp/3.0_dp)
             scale_b=((structure(nbs)%volume/structure(nbs)%num_form+structure(nas)%volume/structure(nas)%num_form)/&
                  (2.0_dp*structure(nbs)%volume/structure(nbs)%num_form))**(1.0_dp/3.0_dp)

             diff = 0.0_dp

             do n=1,min(structure(nas)%num_dist*structure(nbs)%num_form,structure(nbs)%num_dist*structure(nas)%num_form)
                if(structure(nas)%distall(n/structure(nbs)%num_form+1)<rmax/1.75_dp) then
                   newdiff=abs(structure(nas)%distall((n-1)/structure(nbs)%num_form+1)*scale_a-&
                        structure(nbs)%distall((n-1)/structure(nas)%num_form+1)*scale_b)
                   if(newdiff>elim_threshold*(minda*scale_a+mindb*scale_b)/2.0_dp) then
                      similar(nbs)=0         
                      exit
                   end if

                end if
             end do

          end do

          do nb=na+1,num_structures
             nbs=sort_index(nb)
             if(structure(nas)%num_copies==0) cycle

             if(similar(nbs).gt.0) then
                structure(nas)%num_copies=structure(nas)%num_copies+structure(nbs)%num_copies
                structure(nbs)%num_copies=0
             end if
          end do

       end do

    end if

    deallocate(similar)

    if(delete) then

       open(unit=unit_elim_sed,form="formatted",file="elim_sed.sh")
       open(unit=unit_elim_rm,form="formatted",file="elim_rm.sh")

       do na=1,num_structures
          nas=sort_index(na)
          if(orig_copies(nas).ne.structure(nas)%num_copies) then

             if(structure(nas)%num_copies.gt.0) then
                write(char_copies,*) orig_copies(nas)
                write(char_new_copies,*) structure(nas)%num_copies
                ctemp="sed 's/n - "//trim(adjustl(char_copies))//"/n - "&
                     //trim(adjustl(char_new_copies))//"/g' "&
                     //trim(adjustl(structure(nas)%structure_label))//".res > "&
                     //trim(adjustl(structure(nas)%structure_label))//".res-new ; mv "&
                     //trim(adjustl(structure(nas)%structure_label))//".res-new "&
                     //trim(adjustl(structure(nas)%structure_label))//".res"
                write (unit_elim_sed,'(a)') trim(adjustl(ctemp))

             else
                ctemp=""//trim(structure(nas)%structure_label)//".res"
                write (unit_elim_rm,'(a)') trim(adjustl(ctemp))
             end if

          end if

       end do

       flush(unit_elim_sed)
       flush(unit_elim_rm)

       if(elim_threshold.gt.0.0_dp) then
          stat=0
          call system('sh elim_sed.sh')
          if (stat.ne.0) then
             ctemp='Problem executing external command :: '//trim(ctemp)
             write (stderr,'(a)') trim(ctemp)
             stop 
          end if
       end if

       stat=0
       call system('cat elim_rm.sh | xargs rm -f')
       if (stat.ne.0) then
          ctemp='Problem executing external command :: '//trim(ctemp)
          write (stderr,'(a)') trim(ctemp)
          stop 
       end if
       close(unit=unit_elim_sed,status="delete")
       close(unit=unit_elim_rm,status="delete")

    end if

  end subroutine eliminate

  subroutine struct_comm()

    real(kind=dp) :: dist(num_structures,num_structures)
    
    ! ****************************************
    ! ** Identify communities of structures **
    ! ****************************************

    allocate(adjacancy(num_structures,num_structures))

    do na=1,num_structures
       if((structure(na)%num_copies.gt.0).and.(sc_threshold.gt.0.0_dp)) then
          call calculate_distances(na,rmax)           
       end if
    end do

    adjacancy=0
    dist=0.0_dp

    if(sc_threshold.gt.0.0_dp) then

       do na=1,num_structures
          if(structure(na)%num_copies==0) cycle

          adjacancy(:,na) = 1

          do nb=1,num_structures
             if(structure(nb)%num_copies==0) cycle

             if(structure(na)%formula.ne.structure(nb)%formula) then
                adjacancy(nb,na)=0
                cycle
             end if
             diff = 0.0_dp

             do n=1,min(structure(na)%num_dist,structure(nb)%num_dist)
                if(structure(na)%distall(n)<rmax/1.75_dp) then
                   newdiff=abs(structure(na)%distall(n)-structure(nb)%distall(n))
                   if(newdiff>sc_threshold) then 
                      adjacancy(nb,na)=0         
                      cycle
                   end if
                   dist(na,nb)=dist(na,nb)+(structure(na)%distall(n)-structure(nb)%distall(n))**2
                end if
             end do

          end do

       end do

    end if

    call init_modularity(adjacancy,dmod) 

    call optimise_rash()

    write (stderr,'(a,i3,f8.5)') 'Maximum q: ',nunit,q

    do na=1,num_structures
       structure(na)%community=com(na)
    end do

    if(adj) then

       open(unit=unit_adj,form="formatted",file="struc_comm.adj")

       write(unit_adj,*) num_structures
       write(unit_adj,*)
       do na=1,num_structures
          write (unit_adj,*) 'C',0.0_dp,0.0_dp,0.0_dp,structure(na)%enthalpy
       end do
       
       write(unit_adj,*) (count(adjacancy.gt.0)-num_structures)/2
       do na=1,num_structures
          do nb=na+1,num_structures
             if(adjacancy(na,nb).gt.0) write (unit_adj,*) na,nb,sqrt(dist(na,nb))
          end do
       end do

       close(unit_adj)
       
    end if

    deallocate(adjacancy)

  end subroutine struct_comm

  subroutine maxwell()

    ! *******************************************************************
    ! ** Calculate the stable compositions using a Maxwell contruction **
    ! *******************************************************************

    open(unit=unit_qhull,form="formatted",file="qconvex.in")

    allocate(comp_index_inc(num_structures))

    
    ! ** Identify the composition space

    num_species=0
    species_names=''

    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       do n=1,structure(ns)%num_species
          if(all(structure(ns)%species_names(n).ne.species_names(1:num_species))) then
             num_species=num_species+1
             species_names(num_species)=structure(ns)%species_names(n)
          end if
       end do
    end do

    ! ** Check inclusion of structure

    if((len_trim(comp_1).gt.0).and.(len_trim(comp_2).gt.0)) then

       if  (len_trim(comp_4).gt.0) then

          allocate(cmat(num_species,4),comp_vec(4,num_compositions))
          allocate(comp_frac(4,num_compositions),comp_energy(num_compositions))
          allocate(point_energy(4))

          cmat(:,1) = comp2vec(comp_1)
          cmat(:,2) = comp2vec(comp_2)
          cmat(:,3) = comp2vec(comp_3)
          cmat(:,4) = comp2vec(comp_4)

          num_points=4

       elseif (len_trim(comp_3).gt.0) then

          allocate(cmat(num_species,3),comp_vec(3,num_compositions))
          allocate(comp_frac(3,num_compositions),comp_energy(num_compositions))
          allocate(point_energy(3))

          cmat(:,1) = comp2vec(comp_1)
          cmat(:,2) = comp2vec(comp_2)
          cmat(:,3) = comp2vec(comp_3)

          num_points=3

       else

          allocate(cmat(num_species,2),comp_vec(2,num_compositions))
          allocate(comp_frac(2,num_compositions),comp_energy(num_compositions))
          allocate(point_energy(2))

          cmat(:,1) = comp2vec(comp_1)
          cmat(:,2) = comp2vec(comp_2)

          num_points=2

       end if

    else

       allocate(cmat(num_species,num_species),comp_vec(num_species,num_compositions))
       allocate(comp_frac(num_species,num_compositions),comp_energy(num_compositions))
       allocate(point_energy(num_species))

       num_points=num_species

       select case(num_points)
       case(1)
          comp_1=species_names(1)
       case(2)
          comp_1=species_names(1)
          comp_2=species_names(2)
       case(3)
          comp_1=species_names(1)
          comp_2=species_names(2)
          comp_3=species_names(3)
       case(4)
          comp_1=species_names(1)
          comp_2=species_names(2)
          comp_3=species_names(3)
          comp_4=species_names(4)
       end select

       cmat=0
       do ns=1,num_points
          cmat(ns,ns)=1
       end do

    end if

    allocate(comp_include(num_compositions),comp_factor(num_compositions),&
         comp_onhull(num_compositions),comp_efromhull(num_compositions))

    comp_include=.true.

    do nc=1,num_compositions

       comp_vec(:,nc)=vec2decomp(cmat,comp2vec(composition(nc)),comp_factor(nc))

       if(sum(abs(comp_vec(:,nc))).eq.0) comp_include(nc)=.false.

    end do

    ! ** Calculate the convex hull

    write (unit_qhull,*) num_points,count(comp_include)+num_points
    allocate(vec(num_points-1))
    do ns=1,num_points
       vec=0.0_dp
       if(ns.gt.1) vec(ns-1)=1.0_dp
       write (unit_qhull,'(8f17.8)') vec,1000000.00_dp
    end do

    ninc=0
    comp_energy=0.0_dp
    comp_frac=0.0_dp
    point_energy=0.0_dp
    do nc=1,num_compositions
       if(.not.comp_include(nc)) cycle
       ninc=ninc+1

       indx=comp_index(nc)
       comp_index_inc(ninc)=indx

       do i=1,num_points
          comp_frac(i,ninc)=real(comp_vec(i,nc),dp)/real(sum(comp_vec(:,nc)),dp)
       end do

       comp_energy(ninc)=comp_factor(nc)*structure(indx)%enthalpy/&
            real(structure(indx)%num_form*sum(comp_vec(:,nc)),dp)

       if(count(comp_vec(:,nc).gt.0).eq.1) then
          point_energy(maxloc(comp_vec(:,nc)))=comp_energy(ninc)
       end if

    end do

    ninc=0
    do nc=1,num_compositions
       if(comp_include(nc)) then
          ninc=ninc+1
          do n=1,num_points
             comp_energy(ninc)=comp_energy(ninc)-comp_frac(n,ninc)*point_energy(n)
          end do
          write (unit_qhull,'(8f17.8)') comp_frac(1:num_points-1,ninc),comp_energy(ninc)
       end if
    end do

    flush(unit_qhull)

    ! ** It is recommended the qconvex is called as an external program - use system

    if(num_points.gt.2) then
       ctemp="qconvex i Qt Fx n < qconvex.in > qconvex.out "
    else
       ctemp="qconvex Fx n < qconvex.in > qconvex.out "
    end if
    stat=0
    call system(ctemp)
    ctemp='Problem executing external command :: '//trim(ctemp)
    if (stat.ne.0) then
       write (stderr,'(a)') trim(ctemp)
       stop 
    end if

    close(unit_qhull,status="delete")

    ! * Read the external data back in

    open(unit=unit_qhull,form="formatted",file="qconvex.out")
    if(num_points.gt.2) then
       read(unit_qhull,*) num_tri
       allocate(tri_index(3,num_tri))
       do n=1,num_tri
          read(unit_qhull,*) tri_index(:,n)
       end do
    end if
    read(unit_qhull,*) num_hull
    allocate(hull_index(num_hull))
    do n=1,num_hull
       read(unit_qhull,*) hull_index(n)
    end do
    hull_index=hull_index+1-num_points
    read(unit_qhull,*) hull_dim
    read(unit_qhull,*) num_normal
    allocate(hull_normal(hull_dim,num_normal))
    hull_normal=0.0_dp
    do n=1,num_normal
       read(unit_qhull,*) hull_normal(:,n)
    end do
    close(unit_qhull,status="delete")

    ! * Eliminate structures too far from the hull

    ninc=0
    do nc=1,num_compositions
       if(.not.comp_include(nc)) cycle
       ninc=ninc+1
       comp_onhull(ninc)=.false.
       comp_efromhull(ninc)=huge(1.0_dp)
       do ns=1,num_structures
          if(structure(ns)%num_copies==0) cycle
          if(structure(ns)%formula.ne.composition(nc)) cycle

          eref=0
          do n=1,num_points
             eref=eref+comp_frac(n,ninc)*point_energy(n)
          end do

          eform=comp_factor(nc)*structure(ns)%enthalpy/real(structure(ns)%num_form*sum(comp_vec(:,nc)))-eref

          ehull=-huge(1.0_dp)

          do n=2,num_normal

             if(abs(hull_normal(num_points,n)).gt.delta) then

                etemp = -hull_normal(hull_dim,n)/hull_normal(num_points,n)-&
                     dot_product(hull_normal(1:num_points-1,n),comp_frac(1:num_points-1,ninc))/hull_normal(num_points,n)
                if(etemp.gt.ehull) ehull=etemp

             end if

          end do

          structure(ns)%e_from_hull=eform-ehull

          if((structure(ns)%e_from_hull).gt.(delta_e+1e-6_dp)) structure(ns)%num_copies=0

          if((structure(ns)%e_from_hull).lt.(1e-6_dp)) comp_onhull(ninc)=.true.

          if(structure(ns)%e_from_hull.lt.comp_efromhull(ninc)) &
               comp_efromhull(ninc)=structure(ns)%e_from_hull

       end do
    end do

    ! * Output data on stable compounds

    ninc=0
    do nc=1,num_compositions
       if(.not.comp_include(nc)) cycle
       ninc=ninc+1
       indx=comp_index(nc)
       if(structure(indx)%num_copies==0) cycle

       ctemp2=''
       if(comp_onhull(ninc)) then
          ctemp2=' +'
       else
          ctemp2=' -'
       end if

       if(long) then
          fmt='(a40,f9.2,f10.3,f11.3,2f8.3,a2,i3,1x,a12,1x,a7,i5)'
          label=structure(indx)%structure_label
       else
          fmt='(a20,f9.2,f10.3,f11.3,2f8.3,a2,i3,1x,a12,1x,a7,i5)'
          label=structure(indx)%structure_label_short
       end if

       write (stdout,fmt) &
            label,&
            structure(indx)%pressure,&
            structure(indx)%volume/structure(comp_index(ninc))%num_form,&
            structure(indx)%enthalpy/structure(comp_index(ninc))%num_form,&
            comp_energy(ninc),&
            comp_efromhull(ninc),&
            trim(ctemp2),&
            structure(indx)%num_form,&
            structure(indx)%formula,&
            structure(indx)%symmetry,&
            structure(indx)%num_copies

    end do



    if(num_points.eq.2) then

       call write_hull()

       if(xmgrace) then
          ctemp="xmgrace -geometry 1180x920+0+0 hull.agr &"
          stat=0
          call system(ctemp)
          ctemp='Problem executing external command :: '//trim(ctemp)
          if (stat.ne.0) then
             write (stderr,'(a)') trim(ctemp)
             stop 
          end if
       endif

    else if(num_points.eq.3) then

       call write_ternary()

       if(xmgrace) then
          ctemp="xmgrace -geometry 1180x920+0+0 ternary.agr &"
          stat=0
          call system(ctemp)
          ctemp='Problem executing external command :: '//trim(ctemp)
          if (stat.ne.0) then
             write (stderr,'(a)') trim(ctemp)
             stop 
          end if
       endif

       call write_ternary_R()

    else if(num_points.eq.4) then

       call write_quaternary()

    end if

  end subroutine maxwell

  subroutine phull()

    write (stderr,'(a)') 'Constructing pressure hulls'

    do nc=1,num_compositions

       write (stderr,'(2a)') 'Composition: ',composition(nc)
       minvol=huge(1.0_dp)
       maxvol=-huge(1.0_dp)
       emin=huge(1.0_dp)
       emax=-huge(1.0_dp)
       num_struc=0
       do na=1,num_structures
          if(structure(na)%formula.ne.composition(nc)) cycle
          if(structure(na)%num_copies==0) cycle
          num_struc=num_struc+1
          if(structure(na)%volume/structure(na)%num_form.lt.minvol) minvol=structure(na)%volume/structure(na)%num_form
          if(structure(na)%volume/structure(na)%num_form.gt.maxvol) maxvol=structure(na)%volume/structure(na)%num_form
          if(structure(na)%enthalpy/structure(na)%num_form.lt.emin) emin=structure(na)%enthalpy/structure(na)%num_form
          if(structure(na)%enthalpy/structure(na)%num_form.gt.emax) emax=structure(na)%enthalpy/structure(na)%num_form
       end do

       open(unit=unit_qhull,form="formatted",file="qconvex.in")

       write (unit_qhull,*) 2,num_struc+2

       write (unit_qhull,'(2f17.8)') minvol,1000000.00_dp
       write (unit_qhull,'(2f17.8)') maxvol,1000000.00_dp

       do na=1,num_structures
          if(structure(na)%formula.ne.composition(nc)) cycle
          if(structure(na)%num_copies==0) cycle
          write (unit_qhull,'(2f17.8)')  structure(na)%volume/structure(na)%num_form,&
               structure(na)%enthalpy/structure(na)%num_form
       end do

       flush(unit_qhull)

       ! ** It is recommended the qconvex is called as an external program - use system
       ctemp="qconvex Fx n < qconvex.in > qconvex.out "
       stat=0
       call system(ctemp)
       ctemp='Problem executing external command :: '//trim(ctemp)
       if (stat.ne.0) then
          write (stderr,'(a)') trim(ctemp)
          stop 
       end if

       close(unit_qhull,status="delete")

       ! * Read the external data back in

       open(unit=unit_qhull,form="formatted",file="qconvex.out")
       read(unit_qhull,*) num_hull
       allocate(hull_index(num_hull))
       do n=1,num_hull
          read(unit_qhull,*) hull_index(n)
       end do
       hull_index=hull_index+1-num_points
       read(unit_qhull,*) hull_dim
       read(unit_qhull,*) num_normal
       allocate(hull_normal(hull_dim,num_normal))
       hull_normal=0.0_dp
       do n=1,num_normal
          read(unit_qhull,*) hull_normal(:,n)
       end do
       close(unit_qhull,status="delete")

       ! * Eliminate structures too far from the hull

       do na=1,num_structures
          if(structure(na)%formula.ne.composition(nc)) cycle
          if(structure(na)%num_copies==0) cycle

          ehull=-huge(1.0_dp)

          do n=2,num_normal

             if(abs(hull_normal(2,n)).gt.delta) then

                etemp = -hull_normal(hull_dim,n)/hull_normal(2,n)-&
                     hull_normal(1,n)*structure(na)%volume/structure(na)%num_form/hull_normal(2,n)
                if(etemp.gt.ehull) ehull=etemp

             end if

          end do

          structure(na)%e_from_hull=structure(na)%enthalpy/structure(na)%num_form-ehull

          if((structure(na)%e_from_hull).gt.(delta_e+1e-6_dp)) structure(na)%num_copies=0

          if((structure(na)%e_from_hull).lt.(1e-6_dp)) then
             structure(na)%on_hull=.true.
          else
             structure(na)%on_hull=.false.                 
          end if

       end do

       do na=1,num_structures
          if(structure(na)%formula.ne.composition(nc)) cycle
          if(structure(na)%num_copies==0) cycle


          ctemp2=''
          if(structure(na)%on_hull) then
             ctemp2=' +'
          else
             ctemp2=' -'
          end if

          if(long) then
             fmt='(a40,f9.2,f10.3,f11.3,f8.3,a2,i3,1x,a12,1x,a7,i5)'
             label=structure(na)%structure_label
          else
             fmt='(a20,f9.2,f10.3,f11.3,f8.3,a2,i3,1x,a12,1x,a7,i5)'
             label=structure(na)%structure_label_short
          end if

          write (stdout,fmt) &
               label,&
               structure(na)%pressure,&
               structure(na)%volume/structure(na)%num_form,&
               structure(na)%enthalpy/structure(na)%num_form-emin,&
               structure(na)%e_from_hull,&
               trim(ctemp2),&
               structure(na)%num_form,&
               structure(na)%formula,&
               structure(na)%symmetry,&
               structure(na)%num_copies


       end do

       call write_phull()

       deallocate(hull_index,hull_normal)

    end do

  end subroutine phull

  subroutine dos()

    allocate(egrid(ngrid),sdos(ngrid))

    sdos=0.0_dp

    emin=huge(1.0_dp)
    emax=-huge(1.0_dp)
    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       etemp=structure(ns)%enthalpy/structure(ns)%num_form
       if(etemp.lt.emin) emin=etemp
       if(etemp.gt.emax) emax=etemp
    end do

    emin=emin-10.0_dp*dos_smear
    emax=emax+10.0_dp*dos_smear

    do n=1,ngrid
       egrid(n)=emin+real(n-1,dp)/real(ngrid-1,dp)*(emax-emin)
    end do

    n=0
    do ns=1,num_structures
       if(structure(ns)%num_copies==0) cycle
       n=n+1
       etemp=structure(ns)%enthalpy/structure(ns)%num_form

       call  add_gaussian(etemp,egrid,sdos,dos_smear,1.0_dp)
    end do

    sdos=sdos/real(n,dp)

    call write_dos(egrid,sdos)

!!$    do n=1,ngrid
!!$       write (22,*) egrid(n),sdos(n)
!!$    end do

    deallocate(egrid,sdos)

  end subroutine dos

  subroutine write_enthalpy()

    integer :: unit_enthalpy=32
    real(kind=dp) :: xmin,xmax,ymin,ymax,xtick,ytick,low,high
    
    ! * Set the ranges etc

    xmin=minval(peinterp(:,1,:))
    xmax=maxval(peinterp(:,1,:))
    xtick=(xmax-xmin)/10.0_dp

    low=0.0_dp
    high=1.0_dp

    ymin=+huge(1.0_dp)
    ymax=-huge(1.0_dp)
    do n=1,num_interp
       if(minval(peinterp(n,2,:)-peinterp(n,2,nref(1))).lt.ymin) ymin=minval(peinterp(n,2,:)-peinterp(n,2,nref(1)))
       if(maxval(peinterp(n,2,:)-peinterp(n,2,nref(1))).gt.ymax) ymax=maxval(peinterp(n,2,:)-peinterp(n,2,nref(1)))
    end do
    ymin=min(ymin,-ymax*0.1_dp)

    ytick=(ymax-ymin)/10.0_dp

    open(unit=unit_enthalpy,form='formatted',file=trim(composition(nc))//'-enthalpy.agr')

    write(unit_enthalpy,'(a)') '@version 50109'
    write(unit_enthalpy,'(a)') '@default linewidth 2.0'
    write(unit_enthalpy,'(a)') '@default linestyle 1'
    write(unit_enthalpy,'(a)') '@g0 on'
    write(unit_enthalpy,'(a)') '@with g0'
    write(unit_enthalpy,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_enthalpy,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_enthalpy,'(a,f12.3)') '@    world xmin',xmin
    write(unit_enthalpy,'(a,f12.3)') '@    world xmax',xmax
    write(unit_enthalpy,'(a,f12.3)') '@    world ymin',ymin
    write(unit_enthalpy,'(a,f12.3)') '@    world ymax',ymax
    write(unit_enthalpy,'(a)') '@    view xmin 0.150000'
    write(unit_enthalpy,'(a)') '@    view xmax 0.85'
    write(unit_enthalpy,'(a)') '@    view ymin 0.200000'
    write(unit_enthalpy,'(a)') '@    view ymax 0.900000'
    write(unit_enthalpy,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_enthalpy,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_enthalpy,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_enthalpy,'(a,f10.2)') '@    xaxis  tick major',xtick
    write(unit_enthalpy,'(a)') '@    xaxis  label "Pressure (GPa)"'
    write(unit_enthalpy,'(a)') '@    xaxis  label font 4'
    write(unit_enthalpy,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_enthalpy,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_enthalpy,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_enthalpy,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_enthalpy,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_enthalpy,'(a)') '@    yaxis  ticklabel prec 2'
    write(unit_enthalpy,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_enthalpy,'(a)') '@    yaxis  label "Relative Enthalpy (eV/formula unit)"'
    write(unit_enthalpy,'(a)') '@    yaxis  label font 4'
    write(unit_enthalpy,'(a)') '@    yaxis  ticklabel font 4'

    write(unit_enthalpy,'(a)') '@    legend 0.9, 0.9'
    write(unit_enthalpy,'(a)') '@    legend box linewidth 2.0'
    write(unit_enthalpy,'(a)') '@    legend font 4'
    
    
    do na=1,num_pestruc

       write (ctemp,*) na-1
       write(unit_enthalpy,'(a,a,a,i5)') '@    s',adjustl(trim(ctemp)),' line color ',na
       
       do n=1,num_interp
          write (unit_enthalpy,*) peinterp(n,1,na),peinterp(n,2,na)-peinterp(n,2,nref(1))
       end do
       write(unit_enthalpy,'(a)') '&'

    end do
    
    do na=1,num_pestruc

       write (ctemp,*) na-1+num_pestruc
       write(unit_enthalpy,'(a,a,a,a,a)') '@    s',adjustl(trim(ctemp)),' legend  "',adjustl(trim(pestruct(na))),'"'
       write(unit_enthalpy,'(a,a,a)')     '@    s',adjustl(trim(ctemp)),' line type 0'
       write(unit_enthalpy,'(a,a,a)')     '@    s',adjustl(trim(ctemp)),' symbol 1'
       write(unit_enthalpy,'(a,a,a,i5)')  '@    s',adjustl(trim(ctemp)),' symbol color ',na
       write(unit_enthalpy,'(a,a,a)')     '@    s',adjustl(trim(ctemp)),' symbol size 0.5'
       write(unit_enthalpy,'(a,a,a)')     '@    s',adjustl(trim(ctemp)),' symbol fill color 0'
       write(unit_enthalpy,'(a,a,a)')     '@    s',adjustl(trim(ctemp)),' symbol fill pattern 1'

       do n=1,num_penth
          ns=sort_index(n) 
          if(struct_name(ns).eq.pestruct(na)) write (unit_enthalpy,*) pressure(ns),enthalpy(ns)-enref(ns,nref(1))
       end do
       write(unit_enthalpy,'(a)') '&'

    end do    

    close(unit=unit_enthalpy)

    if(xmgrace) then
       ctemp="xmgrace -geometry 1180x920+0+0 "//trim(composition(nc))//"-enthalpy.agr &"
       stat=0
       call system(ctemp)
       ctemp='Problem executing external command :: '//trim(ctemp)
       if (stat.ne.0) then
          write (stderr,'(a)') trim(ctemp)
          stop 
       end if
    endif

  end subroutine write_enthalpy

  subroutine write_dos(energy,sd)

    real(kind=dp), dimension(:), intent(in) :: energy
    real(kind=dp), dimension(:), intent(in) :: sd

    integer :: unit_dos=33,ngrd

    real(kind=dp) :: xmin,xmax,ymin,ymax,xtick,ytick

    ! * Set the ranges etc

    ngrd=size(energy)

    xmin=minval(energy)
    xmax=maxval(energy)
    xtick=(xmax-xmin)/5.0_dp

    ymin=0.0_dp
    ymax=maxval(sd)*1.1_dp
    ytick=(ymax-ymin)/10.0_dp

    open(unit=unit_dos,form='formatted',file='struc-dos.agr')

    write(unit_dos,'(a)') '@version 50109'
    write(unit_dos,'(a)') '@default linewidth 2.0'
    write(unit_dos,'(a)') '@default linestyle 1'
    write(unit_dos,'(a)') '@g0 on'
    write(unit_dos,'(a)') '@with g0'
    write(unit_dos,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_dos,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_dos,'(a)') '@map color 0 to (255, 255, 255), "white"'
    write(unit_dos,'(a)') '@map color 1 to (0, 0, 0), "black"'
    write(unit_dos,'(a)') '@map color 2 to (228, 26, 28), "red"'
    write(unit_dos,'(a)') '@map color 3 to (55, 126, 184), "blue"'
    write(unit_dos,'(a)') '@map color 4 to (77, 175, 74), "green"'
    write(unit_dos,'(a)') '@map color 5 to (152, 78, 163), "purple"'
    write(unit_dos,'(a)') '@map color 6 to (255, 127, 0), "orange"'
    write(unit_dos,'(a)') '@map color 7 to (255, 255, 51), "yellow"'
    write(unit_dos,'(a)') '@map color 8 to (166, 86, 40), "brown"'
    write(unit_dos,'(a)') '@map color 9 to (247, 129, 191), "pink"'
    write(unit_dos,'(a)') '@map color 10 to (153, 153, 153), "grey"'
    write(unit_dos,'(a)') '@map color 11 to (166, 206, 227), "lightblue"'
    write(unit_dos,'(a)') '@map color 12 to (178, 223, 138), "lightgreen"'
    write(unit_dos,'(a)') '@map color 13 to (251, 154, 153), "lightred"'
    write(unit_dos,'(a)') '@map color 14 to (253, 191, 111), "lightorange"'
    write(unit_dos,'(a)') '@map color 15 to (202, 178, 214), "lightpurple"'
    write(unit_dos,'(a,f12.3)') '@    world xmin',xmin
    write(unit_dos,'(a,f12.3)') '@    world xmax',xmax
    write(unit_dos,'(a,f12.3)') '@    world ymin',ymin
    write(unit_dos,'(a,f12.3)') '@    world ymax',ymax
    write(unit_dos,'(a)') '@    view xmin 0.150000'
    write(unit_dos,'(a)') '@    view xmax 0.85'
    write(unit_dos,'(a)') '@    view ymin 0.200000'
    write(unit_dos,'(a)') '@    view ymax 0.900000'
    write(unit_dos,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_dos,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_dos,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_dos,'(a,f10.2)') '@    xaxis  tick major',xtick
    write(unit_dos,'(a)') '@    xaxis  label "Energy (eV)"'
    write(unit_dos,'(a)') '@    xaxis  label font 4'
    write(unit_dos,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_dos,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_dos,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_dos,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_dos,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_dos,'(a)') '@    yaxis  ticklabel prec 2'
    write(unit_dos,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_dos,'(a)') '@    yaxis  label "Density of States"'
    write(unit_dos,'(a)') '@    yaxis  label font 4'
    write(unit_dos,'(a)') '@    yaxis  ticklabel font 4'
    write(unit_dos,'(a)') '@    s0 line color 3'
    write(unit_dos,'(a)') '@    s0 fill type 2'
    write(unit_dos,'(a)') '@    s0 fill rule 0'
    write(unit_dos,'(a)') '@    s0 fill color 11'
    write(unit_dos,'(a)') '@    s0 fill pattern 1'


    do n=1,ngrd
       write (unit_dos,*) energy(n),sd(n)
    end do
    write(unit_dos,'(a)') '&'

    !write(unit_dos,'(a)') '@default linestyle 0'

    close(unit=unit_dos)

    if(xmgrace) then
       ctemp="xmgrace -geometry 1180x920+0+0 struc-dos.agr &"
       stat=0
       call system(ctemp)
       ctemp='Problem executing external command :: '//trim(ctemp)
       if (stat.ne.0) then
          write (stderr,'(a)') trim(ctemp)
          stop 
       end if
    endif

  end subroutine write_dos

  subroutine write_hull()

    integer :: unit_hull=31
    real(kind=dp) :: xmin,xmax,ymin,ymax,xtick,ytick,ypos,low,high
    character(len=10) :: cform,cforma,cformb
    character(len=7)  :: csymm

    ! * Set the ranges etc

    xmin=0.0_dp
    xmax=1.0_dp
    xtick=0.1_dp


    low=minval(comp_energy)

    if(low.lt.0.0_dp) then
       low=min(-0.1_dp,low*1.1_dp)
       high=abs(low)/3.0_dp
    else
       low=0.0_dp
       high=1.0_dp
    end if

    ymin=low
    ymax=high
    ytick=0.1_dp*(high-low)

    open(unit=unit_hull,form='formatted',file='hull.agr')

    write(unit_hull,'(a)') '@version 50109'
    write(unit_hull,'(a)') '@default linewidth 1.5'
    write(unit_hull,'(a)') '@g0 on'
    write(unit_hull,'(a)') '@with g0'
    write(unit_hull,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_hull,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_hull,'(a,f10.3)') '@    world xmin',xmin
    write(unit_hull,'(a,f10.3)') '@    world xmax',xmax
    write(unit_hull,'(a,f10.3)') '@    world ymin',ymin
    write(unit_hull,'(a,f10.3)') '@    world ymax',ymax
    write(unit_hull,'(a)') '@    view xmin 0.150000'
    write(unit_hull,'(a)') '@    view xmax 0.85'
    write(unit_hull,'(a)') '@    view ymin 0.200000'
    write(unit_hull,'(a)') '@    view ymax 0.900000'
    write(unit_hull,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_hull,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_hull,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_hull,'(a,f10.2)') '@    xaxis  tick major',xtick
    cforma=adjustl(structure(comp_index_inc(hull_index(num_points+1)))%formula)
    cformb=adjustl(structure(comp_index_inc(hull_index(num_hull)))%formula)
    write(unit_hull,'(5a)')'@    xaxis  label "x in ',trim(cforma),'\s1-x\N',trim(cformb),'\sx\N"'
    write(unit_hull,'(a)') '@    xaxis  label font 4'
    write(unit_hull,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_hull,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_hull,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_hull,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_hull,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_hull,'(a)') '@    yaxis  ticklabel prec 2'
    write(unit_hull,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_hull,'(a)') '@    yaxis  label "Formation Enthalpy (eV/unit)"'
    write(unit_hull,'(a)') '@    yaxis  label font 4'
    write(unit_hull,'(a)') '@    yaxis  ticklabel font 4'
    write(unit_hull,'(a)') '@    s0 type xy'
    write(unit_hull,'(a)') '@    s0 symbol 1'
    write(unit_hull,'(a)') '@    s0 symbol size 0.500000'
    write(unit_hull,'(a)') '@    s0 symbol color 1'
    write(unit_hull,'(a)') '@    s0 symbol pattern 1'
    write(unit_hull,'(a)') '@    s0 symbol fill color 1'
    write(unit_hull,'(a)') '@    s0 symbol fill pattern 1'
    write(unit_hull,'(a)') '@    s0 line type 1'
    write(unit_hull,'(a)') '@    s0 line linestyle 1'
    write(unit_hull,'(a)') '@    s0 line linewidth 2.0'
    write(unit_hull,'(a)') '@    s0 line color 1'
    write(unit_hull,'(a)') '@    s1 type xy'
    write(unit_hull,'(a)') '@    s1 symbol 1'
    write(unit_hull,'(a)') '@    s1 symbol size 0.300000'
    write(unit_hull,'(a)') '@    s1 symbol color 2'
    write(unit_hull,'(a)') '@    s1 symbol pattern 1'
    write(unit_hull,'(a)') '@    s1 symbol fill color 2'
    write(unit_hull,'(a)') '@    s1 symbol fill pattern 1'
    write(unit_hull,'(a)') '@    s1 line type 0'
    write(unit_hull,'(a)') '@target G0.S0'
    write(unit_hull,'(a)') '@type xy'
    do n=num_points+1,num_hull
       write(unit_hull,'(2f10.5)') comp_frac(1,hull_index(n)),comp_energy(hull_index(n))
    end do
    write(unit_hull,'(a)') '&'
    write(unit_hull,'(a)') '@target G0.S1'
    write(unit_hull,'(a)') '@type xy'
    ninc=0
    do nc=1,num_compositions
       if(.not.comp_include(nc)) cycle
       ninc=ninc+1
       do ns=1,num_structures
          if(structure(ns)%num_copies==0) cycle
          if(structure(ns)%formula.ne.composition(nc)) cycle
          eref=0
          do n=1,num_points
             eref=eref+comp_frac(n,ninc)*point_energy(n)
          end do
          write(unit_hull,'(2f12.5)') comp_frac(1,ninc),&
               comp_factor(nc)*structure(ns)%enthalpy/real(structure(ns)%num_form*sum(comp_vec(:,nc)))-eref
       end do
    end do
    write(unit_hull,'(a)') '&'

    ypos=0.91_dp
    do n=num_points+1,num_hull
       write(unit_hull,'(a)') '@with string'
       write(unit_hull,'(a)') '@    string on'
       write(unit_hull,'(a)') '@    string loctype view'
       ypos=ypos-0.025
       write(unit_hull,'(a,f10.5)') '@    string 0.875,',ypos
       write(unit_hull,'(a)') '@    string color 1'
       write(unit_hull,'(a)') '@    string rot 0'
       write(unit_hull,'(a)') '@    string font 10'
       write(unit_hull,'(a)') '@    string just 0'
       write(unit_hull,'(a)') '@    string char size 0.7500000'
       cform=adjustl(structure(comp_index_inc(hull_index(n)))%formula)
       csymm=adjustl(structure(comp_index_inc(hull_index(n)))%symmetry)
       write(unit_hull,'(a,a10,i2,1x,a7,f6.3,a)') '@    string def "',&
            cform,&
            structure(comp_index_inc(hull_index(n)))%num_form,&
            csymm,&
            comp_frac(1,hull_index(n)),&
            '"'
    end do

    close(unit=unit_hull)

  end subroutine write_hull

  subroutine write_ternary()

    integer :: unit_ternary=32
    real(kind=dp) :: xmin,xmax,ymin,ymax,xtick,ytick,ypos,a,b,c

    ! * Set the ranges etc

    xmin=0.0_dp
    xmax=1.0_dp
    ymin=0
    ymax=1
    xtick=0.1_dp
    ytick=0.2_dp
    ypos=0.91_dp

    open(unit=unit_ternary,form='formatted',file='ternary.agr')

    write(unit_ternary,'(a)') '@version 50109'
    write(unit_ternary,'(a)') '@default linewidth 1.5'
    write(unit_ternary,'(a)') '@g0 on'
    write(unit_ternary,'(a)') '@with g0'
    write(unit_ternary,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_ternary,'(a,f10.3)') '@    world xmin',xmin
    write(unit_ternary,'(a,f10.3)') '@    world xmax',xmax
    write(unit_ternary,'(a,f10.3)') '@    world ymin',ymin
    write(unit_ternary,'(a,f10.3)') '@    world ymax',ymax
    write(unit_ternary,'(a)') '@    view xmin 0.150000'
    write(unit_ternary,'(a)') '@    view xmax 0.85'
    write(unit_ternary,'(a)') '@    view ymin 0.200000'
    write(unit_ternary,'(a)') '@    view ymax 0.900000'
    write(unit_ternary,'(a)') '@    xaxis  bar off'
    write(unit_ternary,'(a)') '@    xaxis  tick major off'
    write(unit_ternary,'(a)') '@    xaxis  tick minor off'
    write(unit_ternary,'(a,f10.2)') '@    xaxis  tick major',xtick
    write(unit_ternary,'(a)') '@    xaxis  label font 4'
    write(unit_ternary,'(a)') '@    xaxis  ticklabel off'
    write(unit_ternary,'(a)') '@    yaxis  bar off'
    write(unit_ternary,'(a)') '@    yaxis  tick major off'
    write(unit_ternary,'(a)') '@    yaxis  tick minor off'
    write(unit_ternary,'(a)') '@    yaxis  ticklabel off'
    write(unit_ternary,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_ternary,'(a)') '@    yaxis  label font 4'
    write(unit_ternary,'(a)') '@    yaxis  ticklabel font 4'
    write(unit_ternary,'(a)') '@    frame type 0'
    write(unit_ternary,'(a)') '@    frame linestyle 0'
    write(unit_ternary,'(a)') '@    s0 type xy'
    write(unit_ternary,'(a)') '@    s0 symbol 1'
    write(unit_ternary,'(a)') '@    s0 symbol size 0.7500000'
    write(unit_ternary,'(a)') '@    s0 symbol color 1'
    write(unit_ternary,'(a)') '@    s0 symbol pattern 1'
    write(unit_ternary,'(a)') '@    s0 symbol fill color 1'
    write(unit_ternary,'(a)') '@    s0 symbol fill pattern 1'
    write(unit_ternary,'(a)') '@    s0 line type 1'
    write(unit_ternary,'(a)') '@    s0 line linestyle 1'
    write(unit_ternary,'(a)') '@    s0 line linewidth 2.0'
    write(unit_ternary,'(a)') '@    s0 line color 1'
    write(unit_ternary,'(a)') '@    s1 type xy'
    write(unit_ternary,'(a)') '@    s1 symbol 1'
    write(unit_ternary,'(a)') '@    s1 symbol size 0.500000'
    write(unit_ternary,'(a)') '@    s1 symbol color 1'
    write(unit_ternary,'(a)') '@    s1 symbol pattern 1'
    write(unit_ternary,'(a)') '@    s1 symbol fill color 2'
    write(unit_ternary,'(a)') '@    s1 symbol fill pattern 1'
    write(unit_ternary,'(a)') '@    s1 line type 0'
    write(unit_ternary,'(a)') '@    s2 type xy'
    write(unit_ternary,'(a)') '@    s2 symbol 1'
    write(unit_ternary,'(a)') '@    s2 symbol size 0.300000'
    write(unit_ternary,'(a)') '@    s2 symbol color 2'
    write(unit_ternary,'(a)') '@    s2 symbol pattern 1'
    write(unit_ternary,'(a)') '@    s2 symbol fill color 2'
    write(unit_ternary,'(a)') '@    s2 symbol fill pattern 1'
    write(unit_ternary,'(a)') '@    s2 line type 0'
    write(unit_ternary,'(a)') '@target G0.S0'
    write(unit_ternary,'(a)') '@type xy'
    write(unit_ternary,'(a)') '0 0'
    write(unit_ternary,'(a)') '1 0'
    write(unit_ternary,'(a)') '0.5 0.8660254'
    write(unit_ternary,'(a)') '0 0'
    write(unit_ternary,'(a)') '&'
    write(unit_ternary,'(a)') '@target G0.S1'
    write(unit_ternary,'(a)') '@type xy'
    do n=num_points+1,num_hull
       a=comp_frac(1,hull_index(n)) ; b=comp_frac(2,hull_index(n)) ; c=comp_frac(3,hull_index(n))
       write(unit_ternary,'(2f10.5)') (2*b+c)/(a+b+c)/2.0_dp,sqrt(3.0_dp)/2.0_dp*c/(a+b+c)
    end do
    write(unit_ternary,'(a)') '&'
    write(unit_ternary,'(a)') '@target G0.S2'
    write(unit_ternary,'(a)') '@type xy'
    ninc=0
    do nc=1,num_compositions
       if(.not.comp_include(nc)) cycle
       ninc=ninc+1
       do ns=1,num_structures
          if(structure(ns)%num_copies==0) cycle
          if(structure(ns)%formula.ne.composition(nc)) cycle
          a=comp_frac(1,ninc) ; b=comp_frac(2,ninc) ; c=comp_frac(3,ninc)
          write(unit_ternary,'(2f12.5)') (2*b+c)/(a+b+c)/2.0_dp,sqrt(3.0_dp)/2.0_dp*c/(a+b+c)
       end do
    end do
    write(unit_ternary,'(a)') '&'
    do n=num_points+1,num_hull
       write(unit_ternary,'(a)') '@with string'
       write(unit_ternary,'(a)') '@    string on'
       write(unit_ternary,'(a)') '@    string loctype view'
       ypos=ypos-0.025
       write(unit_ternary,'(a,f10.5)') '@    string 0.875,',ypos
       write(unit_ternary,'(a)') '@    string color 1'
       write(unit_ternary,'(a)') '@    string rot 0'
       write(unit_ternary,'(a)') '@    string font 4'
       write(unit_ternary,'(a)') '@    string just 0'
       write(unit_ternary,'(a)') '@    string char size 0.7500000'
       write(unit_ternary,'(a,i4,1x,a,1x,a,a)') '@    string def "',structure(comp_index_inc(hull_index(n)))%num_form,&
            trim(adjustl(structure(comp_index_inc(hull_index(n)))%formula)),&
            trim(adjustl(structure(comp_index_inc(hull_index(n)))%symmetry)),&
            '"'
    end do

    write(unit_ternary,'(a)') '@with string'
    write(unit_ternary,'(a)') '@    string on'
    write(unit_ternary,'(a)') '@    string loctype view'
    write(unit_ternary,'(a,f10.5)') '@    string 0.12,0.15'
    write(unit_ternary,'(a)') '@    string color 1'
    write(unit_ternary,'(a)') '@    string rot 0'
    write(unit_ternary,'(a)') '@    string font 4'
    write(unit_ternary,'(a)') '@    string just 0'
    write(unit_ternary,'(a)') '@    string char size 1.0000000'
    write(unit_ternary,'(a,a,a)') '@    string def "',trim(comp_1),'"'

    write(unit_ternary,'(a)') '@with string'
    write(unit_ternary,'(a)') '@    string on'
    write(unit_ternary,'(a)') '@    string loctype view'
    write(unit_ternary,'(a,f10.5)') '@    string 0.875,0.15'
    write(unit_ternary,'(a)') '@    string color 1'
    write(unit_ternary,'(a)') '@    string rot 0'
    write(unit_ternary,'(a)') '@    string font 4'
    write(unit_ternary,'(a)') '@    string just 0'
    write(unit_ternary,'(a)') '@    string char size 1.0000000'
    write(unit_ternary,'(a,a,a)') '@    string def "',trim(comp_2),'"'

    write(unit_ternary,'(a)') '@with string'
    write(unit_ternary,'(a)') '@    string on'
    write(unit_ternary,'(a)') '@    string loctype view'
    write(unit_ternary,'(a,f10.5)') '@    string 0.49,0.85'
    write(unit_ternary,'(a)') '@    string color 1'
    write(unit_ternary,'(a)') '@    string rot 0'
    write(unit_ternary,'(a)') '@    string font 4'
    write(unit_ternary,'(a)') '@    string just 0'
    write(unit_ternary,'(a)') '@    string char size 1.0000000'
    write(unit_ternary,'(a,a,a)') '@    string def "',trim(comp_3),'"'

    close(unit=unit_ternary)

  end subroutine write_ternary

  subroutine write_ternary_R()

    integer :: unit_ternary=32,nhull
    real(kind=dp) :: x(num_compositions),y(num_compositions),z(num_compositions),energy(num_compositions)
    real(kind=dp) :: efromhull(num_compositions)
    character(len=20) :: form(num_compositions),onhull(num_compositions),spacegrp(num_compositions),nfu(num_compositions)
    logical :: anyplus=.false.


    nhull=0
    ninc=0
    do nc=1,num_compositions
       if(.not.comp_include(nc)) cycle
       ninc=ninc+1
!!$       if(comp_energy(ninc).gt.delta) 
       nhull=nhull+1
       x(nhull)=comp_frac(1,ninc)
       y(nhull)=comp_frac(2,ninc)
       z(nhull)=comp_frac(3,ninc)
       form(nhull)='"'//trim(adjustl(structure(comp_index_inc(ninc))%formula))//'"'
       write(nfu(nhull),*) structure(comp_index_inc(ninc))%num_form
       spacegrp(nhull)='"'//trim(adjustl(nfu(nhull)))//' '//trim(adjustl(structure(comp_index_inc(ninc))%symmetry))//'"'
       energy(nhull)=comp_energy(ninc)
       efromhull(nhull)=comp_efromhull(ninc)
       if(comp_onhull(ninc)) then
          onhull(nhull)='"On"'
       else
          onhull(nhull)='"Off"'
          if(comp_energy(ninc).gt.delta) then
             onhull(nhull)='""'
             energy(nhull)=0.0_dp
             anyplus=.true.
          end if
       end if
    end do

    open(unit=unit_ternary,form='formatted',file='ternary.r')

    write(unit_ternary,*) 'library(ggtern)'
    write(unit_ternary,*) 'maxwell <- data.frame('//trim(comp_1)//'=c('
    write(unit_ternary,'(f10.7,",")') x(1:nhull-1)
    write(unit_ternary,'(f10.7,"),")') x(nhull)
    write(unit_ternary,*) trim(comp_2)//'=c('
    write(unit_ternary,'(f10.7,",")') y(1:nhull-1)
    write(unit_ternary,'(f10.7,"),")') y(nhull)
    write(unit_ternary,*) trim(comp_3)//'=c('
    write(unit_ternary,'(f10.7,",")') z(1:nhull-1)
    write(unit_ternary,'(f10.7,"),")') z(nhull)
    write(unit_ternary,*) 'Energy=c('
    write(unit_ternary,'(f10.3,",")') energy(1:nhull-1)
    write(unit_ternary,'(f10.3,"),")') energy(nhull)
    write(unit_ternary,*) 'Efromhull=c('
    write(unit_ternary,'(f10.3,",")') efromhull(1:nhull-1)
    write(unit_ternary,'(f10.3,"),")') efromhull(nhull)
    write(unit_ternary,*) 'Hull=c('
    write(unit_ternary,'(1x,a20,",")') onhull(1:nhull-1)
    write(unit_ternary,'(1x,a20,"),")') onhull(nhull)
    write(unit_ternary,*) 'SpaceGrp=c('
    write(unit_ternary,'(1x,a20,",")') spacegrp(1:nhull-1)
    write(unit_ternary,'(1x,a20,"),")') spacegrp(nhull)
    write(unit_ternary,*) 'Formula=c('
    write(unit_ternary,'(1x,a20,",")') form(1:nhull-1)
    write(unit_ternary,'(1x,a20,"))")') form(nhull)
    write(unit_ternary,*) 'ggtern(data = maxwell, aes(x = '//trim(comp_1)//', y = '//trim(comp_2)//', z = '//trim(comp_3)//')) +'
    !write(unit_ternary,*) 'geom_density2d(aes(weight=Energy,color=abs(..level..))) +'
    do n=1,num_tri
       if(any(tri_index(:,n).lt.3)) cycle
       write(unit_ternary,*) 'geom_segment(aes(x=',comp_frac(1,tri_index(1,n)-2),',y=',&
            comp_frac(2,tri_index(1,n)-2),',z=',comp_frac(3,tri_index(1,n)-2),',xend=',&
            comp_frac(1,tri_index(2,n)-2),',yend=',comp_frac(2,tri_index(2,n)-2),&
            ',zend=',comp_frac(3,tri_index(2,n)-2),&
            '),size=0.2,linetype=1,colour="darkgreen",alpha=0.33) +'
       write(unit_ternary,*) 'geom_segment(aes(x=',comp_frac(1,tri_index(2,n)-2),',y=',&
            comp_frac(2,tri_index(2,n)-2),',z=',comp_frac(3,tri_index(2,n)-2),',xend=',&
            comp_frac(1,tri_index(3,n)-2),',yend=',comp_frac(2,tri_index(3,n)-2),&
            ',zend=',comp_frac(3,tri_index(3,n)-2),&
            '),size=0.2,linetype=1,colour="darkgreen",alpha=0.33) +'
       write(unit_ternary,*) 'geom_segment(aes(x=',comp_frac(1,tri_index(3,n)-2),&
            ',y=',comp_frac(2,tri_index(3,n)-2),',z=',comp_frac(3,tri_index(3,n)-2),&
            ',xend=',comp_frac(1,tri_index(1,n)-2),',yend=',comp_frac(2,tri_index(1,n)-2),&
            ',zend=',comp_frac(3,tri_index(1,n)-2),&
            '),size=0.2,linetype=1,colour="darkgreen",alpha=0.33) +'
    end do
    write(unit_ternary,*) 'geom_point(aes(fill=Energy, size=Efromhull, shape=Hull)) +'
    write(unit_ternary,*) 'geom_text(aes(label=Formula), vjust=0.5, size=0.5, color="gray") +'
    write(unit_ternary,*) 'geom_text(aes(label=SpaceGrp), vjust=-1, size=0.5, color="gray") +'
    if(anyplus) then
       write(unit_ternary,*) 'scale_shape_manual(values = c(46, 22, 21)) +'
    else
       if(any(onhull(1:nhull).eq.'"Off"')) then
          write(unit_ternary,*) 'scale_shape_manual(values = c(22, 21)) +'
       else
          write(unit_ternary,*) 'scale_shape_manual(values = c(21)) +'
       end if
    end if
    write(unit_ternary,*) 'scale_size_continuous(range = c(5, 0.25)) +'
    write(unit_ternary,*) 'scale_fill_gradient(low = "darkgreen", high = "white") +'
    write(unit_ternary,*) 'scale_color_gradient(low="darkgreen",high="white") +'
    write(unit_ternary,*) 'guides(color=FALSE) +'

    write(unit_ternary,*) 'theme_bw()'

    close(unit=unit_ternary)

  end subroutine write_ternary_R

  subroutine write_quaternary()

    integer :: unit_quaternary=33
    real(kind=dp) :: a,b,c,d

    open(unit=unit_quaternary,form='formatted',file='quaternary.plot')

    do n=num_points+1,num_hull
       a=comp_frac(1,hull_index(n)) ; b=comp_frac(2,hull_index(n)) ; c=comp_frac(3,hull_index(n)) ; d=comp_frac(4,hull_index(n))
       write(unit_quaternary,'(3f12.5)') (2.0_dp*b+c+d)/(a+b+c+d)/2.0_dp,&
            (sqrt(3.0_dp/4.0_dp)*c+d/sqrt(12.0_dp))/(a+b+c+d),&
            sqrt(2.0_dp/3.0_dp)*d/(a+b+c+d)
    end do

    write(unit_quaternary,*)
!!$
    ninc=0
    do nc=1,num_compositions
       if(.not.comp_include(nc)) cycle
       ninc=ninc+1
       do ns=1,num_structures
          if(structure(ns)%num_copies==0) cycle
          if(structure(ns)%formula.ne.composition(nc)) cycle
          a=comp_frac(1,ninc) ; b=comp_frac(2,ninc) ; c=comp_frac(3,ninc) ; d=comp_frac(4,ninc)
          write(unit_quaternary,'(3f12.5)') (2.0_dp*b+c+d)/(a+b+c+d)/2.0_dp,&
               (sqrt(3.0_dp/4.0_dp)*c+d/sqrt(12.0_dp))/(a+b+c+d),&
               sqrt(2.0_dp/3.0_dp)*d/(a+b+c+d)
       end do
    end do

    close(unit=unit_quaternary)

  end subroutine write_quaternary

  subroutine write_phull()

    integer :: unit_phull=34
    integer, allocatable, dimension(:) :: sind
    real(kind=dp) :: xmin,xmax,ymin,ymax,xtick,ytick,xpos,ypos,en0,en1,v0,v1
    real(kind=dp), allocatable, dimension(:) :: pest
    real(kind=dp), allocatable, dimension(:,:) :: plot_hull
!!$    character(len=10) :: cform
    character(len=7)  :: csymm

    ! * Set the ranges etc


    emin=huge(1.0_dp)
    emax=-huge(1.0_dp)
    minvol=huge(1.0_dp)
    maxvol=-huge(1.0_dp)
    do na=1,num_structures
       if(structure(na)%formula.ne.composition(nc)) cycle
       if(structure(na)%num_copies==0) cycle
       if(structure(na)%volume/structure(na)%num_form.lt.minvol) minvol=structure(na)%volume/structure(na)%num_form
       if(structure(na)%volume/structure(na)%num_form.gt.maxvol) maxvol=structure(na)%volume/structure(na)%num_form
       if(structure(na)%enthalpy/structure(na)%num_form.lt.emin) emin=structure(na)%enthalpy/structure(na)%num_form
       if(structure(na)%enthalpy/structure(na)%num_form.gt.emax) emax=structure(na)%enthalpy/structure(na)%num_form
    end do

    xmin=minvol-(maxvol-minvol)*0.1
    xmax=maxvol+(maxvol-minvol)*0.1
    xtick=(maxvol-minvol)/5.0_dp

    ymax=(emax-emin)*1.05_dp
    ymin=-(emax-emin)*0.05_dp
    ytick=(emax-emin)/10.0_dp

    open(unit=unit_phull,form='formatted',file=trim(composition(nc))//'-phull.agr')

    write(unit_phull,'(a)') '@version 50109'
    write(unit_phull,'(a)') '@default linewidth 1.5'
    write(unit_phull,'(a)') '@g0 on'
    write(unit_phull,'(a)') '@with g0'
    write(unit_phull,'(a)') '@map font 4 to "Helvetica", "Helvetica"'
    write(unit_phull,'(a)') '@map font 10 to "Courier-Bold", "Courier-Bold"'
    write(unit_phull,'(a)') '@map color 0 to (255, 255, 255), "white"'
    write(unit_phull,'(a)') '@map color 1 to (0, 0, 0), "black"'
    write(unit_phull,'(a)') '@map color 2 to (228, 26, 28), "red"'
    write(unit_phull,'(a)') '@map color 3 to (55, 126, 184), "blue"'
    write(unit_phull,'(a)') '@map color 4 to (77, 175, 74), "green"'
    write(unit_phull,'(a)') '@map color 5 to (152, 78, 163), "purple"'
    write(unit_phull,'(a)') '@map color 6 to (255, 127, 0), "orange"'
    write(unit_phull,'(a)') '@map color 7 to (255, 255, 51), "yellow"'
    write(unit_phull,'(a)') '@map color 8 to (166, 86, 40), "brown"'
    write(unit_phull,'(a)') '@map color 9 to (247, 129, 191), "pink"'
    write(unit_phull,'(a)') '@map color 10 to (153, 153, 153), "grey"'
    write(unit_phull,'(a)') '@map color 11 to (166, 206, 227), "lightblue"'
    write(unit_phull,'(a)') '@map color 12 to (178, 223, 138), "lightgreen"'
    write(unit_phull,'(a)') '@map color 13 to (251, 154, 153), "lightred"'
    write(unit_phull,'(a)') '@map color 14 to (253, 191, 111), "lightorange"'
    write(unit_phull,'(a)') '@map color 15 to (202, 178, 214), "lightpurple"'
    write(unit_phull,'(a,f10.3)') '@    world xmin',xmin
    write(unit_phull,'(a,f10.3)') '@    world xmax',xmax
    write(unit_phull,'(a,f10.3)') '@    world ymin',ymin
    write(unit_phull,'(a,f10.3)') '@    world ymax',ymax
    write(unit_phull,'(a)') '@    view xmin 0.150000'
    write(unit_phull,'(a)') '@    view xmax 0.85'
    write(unit_phull,'(a)') '@    view ymin 0.200000'
    write(unit_phull,'(a)') '@    view ymax 0.900000'
    write(unit_phull,'(a)') '@    xaxis  bar linewidth 1.5'
    write(unit_phull,'(a)') '@    xaxis  tick major linewidth 1.5'
    write(unit_phull,'(a)') '@    xaxis  tick minor linewidth 1.5'
    write(unit_phull,'(a,f10.2)') '@    xaxis  tick major',xtick
    write(unit_phull,'(a)') '@    xaxis  label "Volume (\cE\C\S3\N)"'
    write(unit_phull,'(a)') '@    xaxis  label font 4'
    write(unit_phull,'(a)') '@    xaxis  ticklabel font 4'
    write(unit_phull,'(a)') '@    yaxis  bar linewidth 1.5'
    write(unit_phull,'(a)') '@    yaxis  tick major linewidth 1.5'
    write(unit_phull,'(a)') '@    yaxis  tick minor linewidth 1.5'
    write(unit_phull,'(a)') '@    yaxis  ticklabel format decimal'
    write(unit_phull,'(a)') '@    yaxis  ticklabel prec 2'
    write(unit_phull,'(a,f10.2)') '@    yaxis  tick major',ytick
    write(unit_phull,'(a)') '@    yaxis  label "Enthalpy (eV/unit)"'
    write(unit_phull,'(a)') '@    yaxis  label font 4'
    write(unit_phull,'(a)') '@    yaxis  ticklabel font 4'
    write(unit_phull,'(a)') '@    s0 type xy'
    write(unit_phull,'(a)') '@    s0 symbol 1'
    write(unit_phull,'(a)') '@    s0 symbol size 0.500000'
    write(unit_phull,'(a)') '@    s0 symbol color 3'
    write(unit_phull,'(a)') '@    s0 symbol pattern 1'
    write(unit_phull,'(a)') '@    s0 symbol fill color 3'
    write(unit_phull,'(a)') '@    s0 symbol fill pattern 1'
    write(unit_phull,'(a)') '@    s0 line type 1'
    write(unit_phull,'(a)') '@    s0 line linestyle 1'
    write(unit_phull,'(a)') '@    s0 line linewidth 1.0'
    write(unit_phull,'(a)') '@    s0 line color 11'
    write(unit_phull,'(a)') '@    s1 type xy'
    write(unit_phull,'(a)') '@    s1 symbol 1'
    write(unit_phull,'(a)') '@    s1 symbol size 0.300000'
    write(unit_phull,'(a)') '@    s1 symbol color 11'
    write(unit_phull,'(a)') '@    s1 symbol pattern 1'
    write(unit_phull,'(a)') '@    s1 symbol fill color 11'
    write(unit_phull,'(a)') '@    s1 symbol fill pattern 1'
    write(unit_phull,'(a)') '@    s1 line type 0'
    write(unit_phull,'(a)') '@target G0.S0'
    write(unit_phull,'(a)') '@type xy'

    allocate(plot_hull(size(structure),2),sind(size(structure)),pest(size(structure)))
    n=0
    do ns=1,num_structures
       if(structure(ns)%formula.ne.composition(nc)) cycle
       if(structure(ns)%num_copies==0) cycle
       if(.not.structure(ns)%on_hull) cycle
       n=n+1
       plot_hull(n,1)=structure(ns)%volume/structure(ns)%num_form
       plot_hull(n,2)=structure(ns)%enthalpy/structure(ns)%num_form-emin
    end do

    do i=1,n
       sind(i)=i
    end do

    call heap_sort_index(n,plot_hull(1:n,1),sind(1:n))

    do i=1,n
       write(unit_phull,'(2f15.9)') plot_hull(sind(i),1),plot_hull(sind(i),2)
       en1=plot_hull(sind(i),2)
       v1=plot_hull(sind(i),1)
       if(i.gt.1) then
          pest(i)=-(en1-en0)/(v1-v0)*evbyang3
       end if
       en0=en1
       v0=v1
    end do


    write(unit_phull,'(a)') '&'
    write(unit_phull,'(a)') '@target G0.S1'
    write(unit_phull,'(a)') '@type xy'
    do ns=1,num_structures
       if(structure(ns)%formula.ne.composition(nc)) cycle
       if(structure(ns)%num_copies==0) cycle
       write(unit_phull,*) structure(ns)%volume/structure(ns)%num_form,structure(ns)%enthalpy/structure(ns)%num_form-emin
    end do
    write(unit_phull,'(a)') '&'



    ! ** Sort first

    allocate(ctext(num_structures),vol(num_structures))

    n=0
    do na=1,num_structures
       if(structure(na)%formula.ne.composition(nc)) cycle
       if(structure(na)%num_copies==0) cycle
       if(.not.structure(na)%on_hull) cycle
       n=n+1
       vol(n)=structure(na)%volume/structure(na)%num_form
       csymm=adjustl(structure(na)%symmetry)
       write(ctext(n),'(f7.3,i3,1x,a7,a)') &
            structure(na)%enthalpy/structure(na)%num_form-emin,&
            structure(na)%num_form,&
            csymm,&
            '"'
    end do

    ! ** Add pressures

    xpos=0.5
    ypos=0.5
    do i=2,n
       xpos=(plot_hull(sind(i),1)+plot_hull(sind(i-1),1)-2.0_dp*xmin)/2.0_dp/(xmax-xmin)*0.7_dp&
            +0.1_dp-abs(pest(i))/pest(i)*0.025_dp
       ypos=(plot_hull(sind(i),2)+plot_hull(sind(i-1),2)-2.0_dp*ymin)/2.0_dp/(ymax-ymin)*0.7_dp&
            +0.2_dp-0.01_dp
       write(unit_phull,'(a)') '@with string'
       write(unit_phull,'(a)') '@    string on'
       write(unit_phull,'(a)') '@    string loctype view'
       write(unit_phull,'(a,f10.5,a,f10.5)') '@    string ',xpos,',',ypos
       write(unit_phull,'(a)') '@    string color 1'
       write(unit_phull,'(a)') '@    string rot 0'
       write(unit_phull,'(a)') '@    string font 10'
       write(unit_phull,'(a)') '@    string just 0'
       write(unit_phull,'(a)') '@    string char size 0.400000'
       write(unit_phull,'(a,f7.1,a)') '@    string def "',pest(i),' GPa "'
    end do

    do i=1,n
       sind(i)=i
    end do

    call heap_sort_index(n,vol(1:n),sind(1:n))

    ! ** Now plot legend

    ypos=0.9175_dp

    do na=1,n
       write(unit_phull,'(a)') '@with string'
       write(unit_phull,'(a)') '@    string on'
       write(unit_phull,'(a)') '@    string loctype view'
       ypos=ypos-0.025
       write(unit_phull,'(a,f10.5)') '@    string 0.9,',ypos
       write(unit_phull,'(a)') '@    string color 1'
       write(unit_phull,'(a)') '@    string rot 0'
       write(unit_phull,'(a)') '@    string font 10'
       write(unit_phull,'(a)') '@    string just 0'
       write(unit_phull,'(a)') '@    string char size 0.7500000'
       write(unit_phull,'(a,f7.2,a)') '@    string def "',vol(sind(na)),trim(ctext(sind(na)))
    end do



    close(unit=unit_phull)

    if(xmgrace) then
       ctemp="xmgrace -geometry 1180x920+0+0 "//trim(composition(nc))//"-phull.agr &"
       stat=0
       call system(ctemp)
       ctemp='Problem executing external command :: '//trim(ctemp)
       if (stat.ne.0) then
          write (stderr,'(a)') trim(ctemp)
          stop 
       end if
    endif

    deallocate(plot_hull,sind,pest,ctext,vol)


  end subroutine write_phull

  subroutine get_arguments()

    integer :: iargc,num_args,na

    character(len=40), allocatable, dimension(:) :: argbuff

    logical :: usage=.false.

    num_args = iargc()

    allocate(argbuff(num_args+1))

    do na=1,num_args
       call getarg(na,argbuff(na))
    end do

    task=''
    num_tasks=0

    na = 0
    do while(na.lt.num_args)
       na=na+1
       select case (argbuff(na))
       case ('-g','--geometry')
          num_tasks=num_tasks+1
          task(num_tasks)='geometry'
          read(argbuff(na+1),*,end=110,err=110) thresh
          na=na+1
          cycle
110       thresh=0.1_dp
       case ('-r','--rank')
          num_tasks=num_tasks+1
          task(num_tasks)='rank'
       case ('-s','--summary')
          num_tasks=num_tasks+1
          task(num_tasks)='summary'    
       case ('-e','--enthalpy')
          num_tasks=num_tasks+1
          task(num_tasks)='enthalpy'
          read(argbuff(na+1),*,err=117,end=117) lscale
          na=na+1          
          cycle
117       lscale=1.0_dp
       case ('-f','--formula')
          num_tasks=num_tasks+1
          task(num_tasks)='formula'
          na=na+1
          formula=argbuff(na)
       case ('-fc','--formula_convert')
          num_tasks=num_tasks+1
          task(num_tasks)='formula_convert'
          na=na+1
          formula_convert=argbuff(na)
       case ('-c','--compare')
          num_tasks=num_tasks+1
          task(num_tasks)='compare'
          na=na+1
          read(argbuff(na),*,err=101) comp_threshold
          if(comp_threshold.lt.0.0_dp) then
             sloppy=.true.
             comp_threshold=abs(comp_threshold)
          end if
          na=na+1
          comp_structure=argbuff(na)
          cycle
101       usage=.true. 
       case ('-u','--unite')
          if(count(task=='eliminate')==0) then
             num_tasks=num_tasks+1
             task(num_tasks)='eliminate'
          end if
          na=na+1
          read(argbuff(na),*,err=104) elim_threshold
          if(elim_threshold.lt.0.0_dp) then
             sloppy=.true.
             elim_threshold=abs(elim_threshold)
          end if
          cycle
104       usage=.true.
       case ('-sc','--struct_comm')
          if(count(task=='struct_comm')==0) then
             num_tasks=num_tasks+1
             task(num_tasks)='struct_comm'
          end if
          na=na+1
          read(argbuff(na),*,err=114) sc_threshold        
          cycle
114       usage=.true.
       case ('-m','--maxwell')
          num_tasks=num_tasks+1
          task(num_tasks)='maxwell'
       case ('-ph','--pressure_hull')
          num_tasks=num_tasks+1
          task(num_tasks)='phull'
       case('-xg','--xmgrace')
          xmgrace=.true.
       case('-1')
          na=na+1
          read(argbuff(na),*,err=106) comp_1
          cycle
106       usage=.true.
       case('-2')
          na=na+1
          read(argbuff(na),*,err=107) comp_2 
          cycle
107       usage=.true.
       case('-3')
          na=na+1
          read(argbuff(na),*,err=113) comp_3 
          cycle
113       usage=.true.
       case('-4')
          na=na+1
          read(argbuff(na),*,err=115) comp_4
          cycle
115       usage=.true.
       case ('--delete') ! ** No shortcut for safety
          delete=.true.
       case ('-cl','--cluster')
          cluster=.true.
       case ('-ns','--notsymm')
          notsymm=.true.   
       case ('-x','--xyz')
          xyz=.true.
       case ('-o','--off')
          off=.true.
       case ('-cm','--community')
          comout=.true.
       case ('-wt','--weight')
          weight=.true.
       case ('-am','--adjacancymatrix')
          adj=.true.
       case ('-l','--long')
          long=.true.
       case ('-dr','--distance')
          na=na+1
          read(argbuff(na),*,end=108,err=108) rmax
          cycle
108       usage=.true.
       case ('-bl','--bondlength')
          na=na+1
          read(argbuff(na),*,end=109,err=109) bondlength
          if(bondlength.lt.0.0_dp) then
             bondlength=abs(bondlength)
             modular=.true.
          end if
          cycle
109       usage=.true.
       case ('-bs','--bondscale')
          na=na+1
          read(argbuff(na),*,end=116,err=116) bondscale
          if(bondscale.lt.0.0_dp) then
             bondscale=abs(bondscale)
             modular=.true.
          end if
          bondscale=bondscale*1.2_dp ! ** Adjust so that bs=1.0 connect C-C bonds
          cycle
116       usage=.true.
       case ('-dm','--deltamodularity')
          na=na+1
          read(argbuff(na),*,end=111,err=111) dmod
          cycle
111       usage=.true.
       case ('-al','--alpha')
          na=na+1
          read(argbuff(na),*,end=112,err=112) alpha
          cycle
112       usage=.true.
       case ('-t','--top')
          read(argbuff(na+1),*,end=103,err=103) num_top        
          na=na+1
          cycle
103       num_top=10
       case ('-de','--denergy')
          if(count(task=='eliminate')==0) then       
             num_tasks=num_tasks+1
             task(num_tasks)='eliminate'
          end if
          na=na+1
          read(argbuff(na),*,err=105) delta_e   
          cycle
105       usage=.true.
       case ('-n','--num_units')
          num_tasks=num_tasks+1
          task(num_tasks)='connect'
          na=na+1
          read(argbuff(na),*,end=118,err=118) num_units
          cycle
118       usage=.true.
       case ('-d','--dimensionality')
          num_tasks=num_tasks+1
          task(num_tasks)='connect'
          na=na+1
          read(argbuff(na),*,end=119,err=119) dimensionality
          cycle
119       usage=.true.
       case ('-sd','--struc_dos')
          num_tasks=num_tasks+1
          task(num_tasks)='dos'
          na=na+1
          read(argbuff(na),*,end=120,err=120) dos_smear
          cycle
120       usage=.true.
       case ('-p','--pressure')
          na=na+1
          read(argbuff(na),*,end=121,err=121) press_add
          cycle
121       usage=.true.
       case ('-fr','--formula_reference')
          na=na+1
          formula_ref=argbuff(na)
       case ('-h','--help','-?')
          usage=.true.
       case default
          usage=.true.
       end select
    end do

    if(num_tasks==0) usage=.true.

    if(usage) then

       print '(a)', '' 
       print '(a)', 'Usage: cryan [OPTIONS]' 
       print '(a)', '' 
       print '(a)', 'The structures are read from STDIN, for example:' 
       print '(a)', '' 
       print '(a)', '     cat *.res | cryan -s' 
       print '(a)', '     gunzip -c lots.res.gz | cryan -f H2O' 
       print '(a)', '     find . -name "*.res" | xargs cat | cryan -m' 
       print '(a)', '' 
       print '(a)', 'cryan options (note - the order matters):' 
       print '(a)', '' 
       print '(a)', ' -r,  --rank                          Rank all structures, of any composition' 
       print '(a)', ' -s,  --summary                       Summary, most stable from each composition' 
       print '(a)', ' -e,  --enthalpy <length_scale>       Plot enthalpy vs. pressure, interpolate with <length_scale>' 
       print '(a)', ' -f,  --formula <formula>             Select structures of a given composition' 
       print '(a)', ' -fc,  --formula_convert <formula>    Attempt to convert structure to this composition' 
       print '(a)', ' -t,  --top [num]                     Output top few results (default 10)' 
       print '(a)', ' -u,  --unite <thresh>                Unite similar structures' 
       print '(a)', ' -dr, --distance <rmax>               Distance threshold for structure comparison (default 20)' 
       print '(a)', ' -de, --delta_e <energy>              Ignore structures above energy (per atom)' 
       print '(a)', ' -sd, --struc_dos <smear>             Plot a structural density of states, smeared' 
       print '(a)', ' -p,  --pressure <pressure>           Additional pressure (default 0 GPa)' 
       print '(a)', ' -m,  --maxwell                       Extract the stable compositions' 
       print '(a)', ' -ph, --pressure_hull                 Extract the stable structures with pressure' 
       print '(a)', ' -<n>                                 Component <n>' 
       print '(a)', ' -xg, --xmgrace                       Plot output with xmgrace' 
       print '(a)', ' -c,  --compare <thresh > <structure> Compare structure to all others' 
       print '(a)', '      --delete                        Delete unwanted structures' 
       print '(a)', ' -g,  --geometry [thresh]             Calculate the atomic geometry for the structures (default 0.1)' 
       print '(a)', ' -n,  --num_units                     Only report structures with n separate units (default -1)' 
       print '(a)', ' -d,  --dimensionality                Only report structures with dimensionality of d (default -1.0)' 
       print '(a)', ' -cl, --cluster                       No periodic boundary conditions' 
       print '(a)', ' -bl, --bondlength                    Maximum bond length (default 0.0, negative for modularity)' 
       print '(a)', ' -bs, --bondscale                     Bond length scaling (default 1.0, negative for modularity)' 
       print '(a)', ' -dm, --deltamodularity               Modularity bias parameter' 
       print '(a)', ' -wt, --weight                        Weight the adjacancy matrix toward short contacts' 
       print '(a)', ' -ns, --notsymm                       Do not calculate point group of clusters' 
       print '(a)', ' -sc, --struct_comm <thresh>          Determine the community structure' 
       print '(a)', ' -cm, --community                     Output the community structure' 
       print '(a)', ' -am, --adjacancymatrix               Output the adjacancy matrix' 
       print '(a)', ' -x,  --xyz                           Output clusters in XYZ format' 
       print '(a)', ' -o,  --off                           Output polyhedra in OFF format' 
       print '(a)', ' -al, --alpha                         Construct alpha shapes' 
       print '(a)', ' -l,  --long                          Long names for structures' 
       print '(a)', ' -h,  --help, -?                      Print usage information and exit' 
       print '(a)', '' 

       stop

    end if

  end subroutine get_arguments

  subroutine calculate_distances(num,rmax,nosort,laplacian,superlap)

    integer,                                 intent(in)  :: num
    real(kind=dp),                           intent(in)  :: rmax
    logical,                       optional, intent(in)  :: nosort
    real(kind=dp),                 optional, intent(in)  :: laplacian
    real(kind=dp), dimension(:,:), optional, intent(out) :: superlap

    ! ------------- !

    integer :: n,nn,na,nb,nmax,n1max,n2max,n3max,n1min,n2min,n3min,n1,n2,n3,naa,nbb,na1,na2,na3,nb1,nb2,nb3
    real(kind=dp) :: veca(3),vecb(3),dist2,thresh2,lap2

    logical :: sort,lap

    if(present(nosort)) then
       sort = .not.nosort
    else
       sort = .true.
    end if

    if(present(laplacian)) then
       lap = .true.
       lap2 = laplacian**2
    else
       lap = .false.
       lap2 = 0.0_dp
    end if


    ! ------------- !      

    nmax = 1

102 n1max=0 ; n2max=0 ; n3max=0 ; n1min=0 ; n2min=0 ; n3min=0

    if(.not.cluster) then

       thresh2=rmax**2

       do n1=-nmax,nmax
          do n2=-nmax,nmax
             do n3=-nmax,nmax

                veca(:) = n1*structure(num)%lattice_car(:,1)+&
                     n2*structure(num)%lattice_car(:,2)+&
                     n3*structure(num)%lattice_car(:,3)

                dist2 = dot_product(veca,veca)

                if(dist2<thresh2) then
                   if(n1<n1min) n1min=n1
                   if(n2<n2min) n2min=n2
                   if(n3<n3min) n3min=n3
                   if(n1>n1max) n1max=n1
                   if(n2>n2max) n2max=n2
                   if(n3>n3max) n3max=n3
                end if

             end do
          end do
       end do

       if((n1max==nmax).or.(n2max==nmax).or.(n3max==nmax).or.(n1min==-nmax).or.(n2min==-nmax).or.(n3min==-nmax)) then
          nmax=nmax+1
          goto 102
       end if

    end if

    structure(num)%num_dist = structure(num)%num_ions*(n1max-n1min+1)*(n2max-n2min+1)*(n3max-n3min+1)

    if(.not.allocated(structure(num)%distall)) then
       allocate(structure(num)%distall(structure(num)%num_dist*structure(num)%num_ions))
       allocate(structure(num)%dist_names(2,structure(num)%num_dist*structure(num)%num_ions))
       allocate(structure(num)%dist_ions(2,structure(num)%num_dist*structure(num)%num_ions))
       allocate(structure(num)%dist_incell(structure(num)%num_dist*structure(num)%num_ions))
       allocate(structure(num)%dist_indx(structure(num)%num_dist*structure(num)%num_ions))
       if(lap) then
          allocate(structure(num)%laplacian(structure(num)%num_ions,structure(num)%num_ions))
          allocate(structure(num)%adjacancy(structure(num)%num_ions,structure(num)%num_ions))
          structure(num)%laplacian=0.0_dp
          structure(num)%adjacancy=0                  
       end if
    end if

    thresh2=(rmax/1.75_dp)**2

    nn=0
    do na=1,structure(num)%num_ions

       veca(:) = structure(num)%positions_frac(1,na)*structure(num)%lattice_car(:,1)+&
            structure(num)%positions_frac(2,na)*structure(num)%lattice_car(:,2)+&
            structure(num)%positions_frac(3,na)*structure(num)%lattice_car(:,3)

       structure(num)%positions_cart(:,na)=veca(:)

       ! ** Loop over an extended number of cells

       n=0
       do nb=1,structure(num)%num_ions

          if((bondscale.gt.0.0_dp).and.lap) then
             lap2=elements_radius(structure(num)%ion_names_indx(na))+&
                  elements_radius(structure(num)%ion_names_indx(nb))
             lap2=bondscale*lap2
             lap2=lap2**2

          end if

          do n1=n1min,n1max
             do n2=n2min,n2max
                do n3=n3min,n3max
                   n=n+1

                   vecb(:)=structure(num)%positions_frac(1,nb)*structure(num)%lattice_car(:,1)+&
                        structure(num)%positions_frac(2,nb)*structure(num)%lattice_car(:,2)+&
                        structure(num)%positions_frac(3,nb)*structure(num)%lattice_car(:,3)+&
                        n1*structure(num)%lattice_car(:,1)+n2*structure(num)%lattice_car(:,2)+&
                        n3*structure(num)%lattice_car(:,3)

                   dist2 = dot_product((veca-vecb),(veca-vecb))

                   if(dist2<thresh2) then
                      nn=nn+1
                      structure(num)%distall(nn) = sqrt(dist2)
                      structure(num)%dist_names(1,nn) = structure(num)%ion_names(na)                      
                      structure(num)%dist_names(2,nn) = structure(num)%ion_names(nb) 
                      structure(num)%dist_ions(1,nn) = na                     
                      structure(num)%dist_ions(2,nn) = nb    
                      if ((n1.eq.0).and.(n2.eq.0).and.(n3.eq.0)) then
                         structure(num)%dist_incell(nn) = .true.
                      else
                         structure(num)%dist_incell(nn) = .false.                         
                      end if
                   end if

                   if((lap.and.(dist2<lap2)).and.(na.ne.nb)) then
                      structure(num)%laplacian(na,nb)=-1.0_dp
                      if(weight) then
                         structure(num)%adjacancy(na,nb)=structure(num)%adjacancy(na,nb)+&
                              int(100.0_dp*(1.0_dp-sqrt(dist2/lap2)))
                      else
                         structure(num)%adjacancy(na,nb)=structure(num)%adjacancy(na,nb)+100
                      end if
                   end if

                end do
             end do
          end do

       end do

    end do

    if(lap) then

       do na=1,structure(num)%num_ions
          structure(num)%laplacian(na,na)=-sum(structure(num)%laplacian(na,:))
       end do

    end if

    if(sort) then
       call heap_sort(nn,structure(num)%distall(1:nn))
    else
       do n=1,nn
          structure(num)%dist_indx(n)=n
       end do
       call heap_sort_index(nn,structure(num)%distall(1:nn),structure(num)%dist_indx(1:nn))       
    end if

    structure(num)%num_dist = nn

    if(present(superlap)) then
       superlap=0.0_dp

       if(cluster) stop 'Dimensionality evaluation incompatible with cluster calculation'

       nmax = 1

103    n1max=0 ; n2max=0 ; n3max=0 ; n1min=0 ; n2min=0 ; n3min=0


       thresh2=rmax**2

       do n1=-nmax,nmax
          do n2=-nmax,nmax
             do n3=-nmax,nmax

                veca(:) = 2*n1*structure(num)%lattice_car(:,1)+&
                     2*n2*structure(num)%lattice_car(:,2)+&
                     2*n3*structure(num)%lattice_car(:,3)

                dist2 = dot_product(veca,veca)

                if(dist2<thresh2) then
                   if(n1<n1min) n1min=n1
                   if(n2<n2min) n2min=n2
                   if(n3<n3min) n3min=n3
                   if(n1>n1max) n1max=n1
                   if(n2>n2max) n2max=n2
                   if(n3>n3max) n3max=n3
                end if

             end do
          end do
       end do

       if((n1max==nmax).or.(n2max==nmax).or.(n3max==nmax).or.(n1min==-nmax).or.(n2min==-nmax).or.(n3min==-nmax)) then
          nmax=nmax+1
          goto 103
       end if

       naa=0
       do na=1,structure(num)%num_ions
          do na1=0,1
             do na2=0,1
                do na3=0,1
                   naa=naa+1

                   veca(:) = (structure(num)%positions_frac(1,na)+na1)*structure(num)%lattice_car(:,1)+&
                        (structure(num)%positions_frac(2,na)+na2)*structure(num)%lattice_car(:,2)+&
                        (structure(num)%positions_frac(3,na)+na3)*structure(num)%lattice_car(:,3)

                   ! ** Loop over an extended number of cells

                   nbb=0
                   do nb=1,structure(num)%num_ions
                      do nb1=0,1
                         do nb2=0,1
                            do nb3=0,1
                               nbb=nbb+1

                               if((bondscale.gt.0.0_dp).and.lap) then

                                  lap2=elements_radius(structure(num)%ion_names_indx(na))+&
                                       elements_radius(structure(num)%ion_names_indx(nb))
                                  lap2=bondscale*lap2
                                  lap2=lap2**2

                               end if

                               do n1=n1min,n1max
                                  do n2=n2min,n2max
                                     do n3=n3min,n3max
                                        n=n+1


                                        vecb(:)=(structure(num)%positions_frac(1,nb)+2*n1+nb1)*structure(num)%lattice_car(:,1)+&
                                             (structure(num)%positions_frac(2,nb)+2*n2+nb2)*structure(num)%lattice_car(:,2)+&
                                             (structure(num)%positions_frac(3,nb)+2*n3+nb3)*structure(num)%lattice_car(:,3)


                                        dist2 = dot_product((veca-vecb),(veca-vecb))

                                        if((lap.and.(dist2<lap2)).and.(naa.ne.nbb)) then
                                           superlap(naa,nbb)=-1.0_dp                         
                                        end if

                                     end do
                                  end do
                               end do

                            end do
                         end do
                      end do
                   end do

                end do
             end do
          end do
       end do


       do naa=1,size(superlap(:,1))
          superlap(naa,naa)=-sum(superlap(naa,:))
       end do

    end if

  end subroutine calculate_distances

  subroutine join_unit(pos,nat,indx)

    real(kind=dp), dimension(:,:), intent(inout) :: pos
    integer,                       intent(in)    :: nat
    integer,       dimension(:),   intent(in)    :: indx

    real(kind=dp) :: vtemp(3),stemp(3),bl2

    integer :: nna,nnb,nn,nn1,nn2,nn3
    logical :: moved(nat)

    ! *-

    if(nat.le.1) return

    moved=.true.
    nn=0
    do while(any(moved).and.(nn.lt.100))
       nn=nn+1
       moved=.false.

       do nna=1,nat

          do nnb=1,nat

             if(nna.eq.nnb) cycle

             vtemp(:)=pos(:,nnb)-pos(:,nna)

             if(bondscale.gt.0.0_dp) then
                bl2=bondscale*(elements_radius(indx(nna))+elements_radius(indx(nnb)))
                bl2=bl2**2
             else
                bl2=bondlength**2
             end if

             do nn1=-2,2
                do nn2=-2,2
                   do nn3=-2,2

                      if((nn1.eq.0).and.(nn2.eq.0).and.(nn3.eq.0)) cycle

                      stemp(:) = nn1*structure(ns)%lattice_car(:,1)+nn2*structure(ns)%lattice_car(:,2)+&
                           nn3*structure(ns)%lattice_car(:,3)

                      dist2=dot_product(vtemp+stemp,vtemp+stemp)

                      if((dist2.lt.bl2).and.(.not.moved(nnb))) then
                         pos(:,nnb) = pos(:,nnb)+stemp(:)
                         moved(nnb)=.true.
                      end if

                   end do
                end do
             end do

          end do

       end do

    end do

  end subroutine join_unit

  subroutine compact_unit(pos,nat,indx,nedge)

    real(kind=dp), dimension(:,:), intent(inout) :: pos
    integer,                       intent(in)    :: nat
    integer,       dimension(:),   intent(in)    :: indx    
    integer,                       intent(out)   :: nedge

    real(kind=dp), allocatable, dimension(:,:) :: poso
    real(kind=dp) :: pshift(3,27),bl2

    integer :: nshape(2),nna,nnb,nemax,ne,nshift,nn,nn1,nn2,nn3,nsmx,nemx,nnn

    logical :: changed

    ! *-

    if(nat.le.1) return

    nshape(:)=shape(pos)
    allocate(poso(nshape(1),nshape(2)))

    nn=0
    do nn1=-1,1
       do nn2=-1,1
          do nn3=-1,1
             nn=nn+1
             pshift(:,nn)=nn1*structure(ns)%lattice_car(:,1)+&
                  nn2*structure(ns)%lattice_car(:,2)+&
                  nn3*structure(ns)%lattice_car(:,3)
          end do
       end do
    end do

    poso=pos
    nemax=unit_edges(pos,nat,indx)
    nn=0
    do while(nn.lt.100)
       nn=nn+1

       ! * Greedy 

       nnn=0
       changed=.true.
       do while(changed.and.(nnn.lt.1000))
          changed=.false.
          nnn=nnn+1

          do nna=1,nat

             nemx=0
             do nshift=1,27

                ne=0
                do nnb=1,nat

                   if(bondscale.gt.0.0_dp) then
                      bl2=bondscale*(elements_radius(indx(nna))+elements_radius(indx(nnb)))
                      bl2=bl2**2
                   else
                      bl2=bondlength**2
                   end if

                   vtemp(:)=pos(:,nna)+pshift(:,nshift)-pos(:,nnb)

                   if(dot_product(vtemp,vtemp).lt.bl2) ne=ne+1

                end do
                if(ne.gt.nemx) then
                   nemx=ne
                   nsmx=nshift
                end if

             end do

             if(nsmx.ne.14) then
                pos(:,nna)=pos(:,nna)+pshift(:,nsmx)
                changed=.true.
             end if
          end do
       end do

       ! * Check

       ne=unit_edges(pos,nat,indx)
       if(ne.gt.nemax) then
          poso=pos
          nemax=ne
          nn=0
       else
          pos=poso
       end if

       ! * Shake

       do nna=1,nat
          if(random_single().gt.0.5_dp) cycle
          nshift=1+int(random_single()*27)
          pos(:,nna)=pos(:,nna)+pshift(:,nshift)
       end do

    end do

    pos=poso

    deallocate(poso)

    nedge=nemax

  end subroutine compact_unit

  function unit_edges(pos,nat,indx)

    real(kind=dp), dimension(:,:), intent(in) :: pos
    integer,                       intent(in) :: nat
    integer,       dimension(:),   intent(in) :: indx

    real(kind=dp) :: bl2

    integer :: unit_edges

    integer :: ii,jj

    unit_edges=0
    do ii=1,nat
       do jj=ii+1,nat
          if(bondscale.gt.0.0_dp) then
             bl2=bondscale*(elements_radius(indx(ii))+elements_radius(indx(jj)))
             bl2=bl2**2
          else
             bl2=bondlength**2
          end if
          if(dot_product(pos(:,ii)-pos(:,jj),pos(:,ii)-pos(:,jj)).lt.bl2) &
               unit_edges=unit_edges+1
       end do
    end do

  end function unit_edges

  function comp2vec(comp)

    character(len=*), intent(in) :: comp

    integer, dimension(num_species) :: comp2vec

    integer :: n,ns

    character(len=10), parameter :: numbers='1234567890'

    character(len=240) :: cwork,cwork2

    cwork=comp

    comp2vec=0

    do while (len_trim(cwork).gt.0)

       if(any(species_names.eq.cwork(1:2))) then
          do ns=1,num_species
             if(cwork(1:2).eq.species_names(ns)) exit
          end do
          cwork=cwork(3:)
          cwork2=cwork(:verify(cwork,numbers)-1)
          if(len_trim(cwork2).gt.0) then
             read(cwork2,*) n
          else
             n=1
          end if
          cwork=cwork(verify(cwork,numbers):)
          comp2vec(ns)=comp2vec(ns)+n
       else if(any(species_names.eq.cwork(1:1)))then
          do ns=1,num_species
             if(cwork(1:1).eq.species_names(ns)) exit
          end do
          cwork=cwork(2:)
          cwork2=cwork(:verify(cwork,numbers)-1)
          if(len_trim(cwork2).gt.0) then
             read(cwork2,*) n
          else
             n=1
          end if
          cwork=cwork(verify(cwork,numbers):)
          comp2vec(ns)=comp2vec(ns)+n
       else
          stop 'comp2vec : composition not understood'
       end if

    end do

  end function comp2vec

  function vec2decomp(M,V,factor)

    integer, dimension(:,:), intent(in)  :: M
    integer, dimension(:),   intent(in)  :: V
    integer,                 intent(out) :: factor

    integer, dimension(size(M(1,:))) :: vec2decomp

    logical :: success

    ! --

    integer :: n,V1(num_species),nf

    call linear_solve_int(M,V,vec2decomp,success)

    if(success) then

       nf = vec2decomp(1)
       do n = 2,size(M(1,:))
          if(vec2decomp(n).gt.0) nf = gcd_rec(nf,vec2decomp(n))
       end do

       vec2decomp=vec2decomp/nf

       V1=matmul(M,vec2decomp)

       factor=maxval(V1)/maxval(V)

    else

       vec2decomp=0
       factor=1

    end if

    return

  end function vec2decomp

  subroutine heap_sort(num_items,weight)

    !=========================================================================!
    ! This subroutine sorts the list of weights into descending order.        !
    !                                                                         !
    ! This is a heap sort                                                     !
    !-------------------------------------------------------------------------!
    ! Arguments:                                                              !
    !   num_items (input) :: The number of items to sort                      !
    !   weight (in/out)   :: The weights of each item. On exit these are      !
    !                        sorted into descending order.                    !
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
    ! Written by Chris Pickard 22nd May 2009                                  !
    !=========================================================================!

    implicit none

    ! Arguments

    integer, intent(in) :: num_items
    real(kind=dp), dimension(num_items), intent(inout) :: weight

    ! Local variables

    integer :: i,ir,j,l ! Loop counters
    real(kind=dp) :: wta

    if(num_items.lt.2) return

    l=num_items/2+1
    ir=num_items

    do
       if(l.gt.1) then
          l=l-1
          wta=weight(l)
       else
          wta=weight(ir)
          weight(ir)=weight(1)
          ir=ir-1
          if(ir.eq.1) then
             weight(1)=wta
             return
          end if
       end if
       i=l
       j=l+l
20     if(j.le.ir) then
          if(j.lt.ir) then
             if(weight(j).lt.weight(j+1)) j=j+1
          end if
          if(wta.lt.weight(j)) then
             weight(i)=weight(j)
             i=j
             j=j+j
          else
             j=ir+1
          end if
          goto 20
       end if
       weight(i)=wta
    end do

  end subroutine heap_sort

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
    ! Written by Chris Pickard 26th July 2012                                 !
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

  recursive function gcd_rec(u, v) result(gcd)
    integer             :: gcd
    integer, intent(in) :: u, v

    if (mod(u, v) /= 0) then
       gcd = gcd_rec(v, mod(u, v))
    else
       gcd = v
    end if
  end function gcd_rec

  subroutine permute(list)

    integer, dimension(:), intent(out) :: list

    integer :: i,j

    do i=1,size(list)
       j=1+int(random_single()*i)
       list(i)=list(j)
       list(j)=i
    end do

  end subroutine permute

  function inv(A) result(Ainv)

    real(dp), dimension(:,:), intent(in) :: A
    real(dp), dimension(size(A,1),size(A,2)) :: Ainv

    real(dp), dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
       stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
       stop 'Matrix inversion failed!'
    end if
  end function inv

  subroutine linear_solve_int(Aint,Bint,Dint,success)

    integer, dimension(:,:), intent(in)     :: Aint
    integer, dimension(:),   intent(in)     :: Bint
    integer, dimension(:),   intent(out)    :: Dint
    logical,                 intent(out)    :: success

    integer :: m,n
    integer :: lda,ldb
    integer :: lwork,info

    integer,parameter :: lwmax=1000    
    real(kind=dp)     :: work(lwmax)

    real(kind=dp), allocatable :: A(:,:),B(:,:)

    external DGELS

    m=size(Aint,1) ; n=size(Aint,2)
    lda=m ; ldb=size(Bint,1)

    allocate(A(lda,n),B(ldb,1))

    A=real(Aint,dp)
    B(:,1)=real(Bint,dp)

    lwork = -1
    call DGELS( 'N', m, n, 1, A, lda, B, ldb, work, lwork, info )
    lwork = min(lwmax,int(work(1)))
    call DGELS( 'N', m, n, 1, A, lda, B, ldb, work, lwork, info )

    ! ** Check we got a positive integer solution

    if(all(matmul(Aint,nint(abs(B(1:n,1)))).eq.Bint)) then
       Dint=nint(B(1:n,1))
       success=.true.
    else
       Dint=0
       success=.false.
    end if

  end subroutine linear_solve_int

  function least_square(n,x,y)

    integer, intent(in) :: n
    real(kind=dp), dimension(:), intent(in) :: x,y

    real(kind=dp), dimension(2) :: least_square

    integer :: i

    real(kind=dp) xm,ym,cov,var

    xm=sum(x(1:n))/real(n,dp)
    ym=sum(y(1:n))/real(n,dp)

    least_square=0.0_dp

    cov=0.0_dp
    var=0.0_dp
    do i=1,n
       cov=cov+(x(i)-xm)*(y(i)-ym)
       var=var+(x(i)-xm)**2
    end do

    least_square(1)=cov/var
    least_square(2)=ym-least_square(1)*xm

  end function least_square

  subroutine add_gaussian(x0,x,f,sig,wgt)

    real(kind=dp),               intent(in)    :: x0
    real(kind=dp), dimension(:), intent(in)    :: x
    real(kind=dp), dimension(:), intent(inout) :: f
    real(kind=dp),               intent(in)    :: sig
    real(kind=dp),               intent(in)    :: wgt

    integer :: i,n

    n=size(x)

    do i=1,n
       f(i)=f(i)+exp(-(x(i)-x0)**2/2.0_dp/sig**2)/sig*wgt/sqrt(tpi)
    end do

  end subroutine add_gaussian

end program cryan
