! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!==================================================================================!
!                                      Build                                       !
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
! This module stochastically builds the crystal structure                          !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module build

  use constants
  use cell
  use symmetry
  use rng

  implicit none

  private

  public :: gen_lattice
  public :: gen_position_check
  public :: check_separation
  public :: check_coordination
  public :: push_opt
  public :: assign_spin

  ! ** Parameters
  
  ! ** Module variables

  logical, public :: initialised=.false.
  
  ! ** Counters and temp variables

  integer :: i,j,n,m,n1,n2,n3,ni,nj,ns,jmin,mult,njoin,nm,nindx(3)

  integer, allocatable, dimension(:) :: ionlist,indx
  real(kind=dp), allocatable, dimension(:) :: distlist2
  real(kind=dp), allocatable, dimension(:,:) :: veclist

  
  real(kind=dp) :: vol,ang,angle,rn(3),an(3),V(3),VO(3),VN(3),min2,dist2,temp,mat(3,3),mirror(3,3)
  real(kind=dp) :: vol_tmp,lattice_tmp(6),s(3,3)

  real(kind=dp), dimension(:,:,:), allocatable :: thresh

  logical, dimension(:), allocatable :: taken
  
contains

  subroutine gen_lattice(stat)

    integer, intent(inout) :: stat
    
    if(stat/=0) return
    
    ! ** Reset the symmetry operators

    if(internalsymm) then
       call symm_by_number(symmetry_number)
       do ns=1,num_symm
          symm_ops(1:3,1:3,ns)=matmul(matmul(setting_perm,symm_ops(1:3,1:3,ns)),transpose(setting_perm))
          symm_ops(1:3,4,ns)=matmul(setting_perm,symm_ops(1:3,4,ns))
       end do
       cang=matmul(setting_perm,cang)
       clen=matmul(setting_perm,clen)
    end if
    
    if(cellamp.ge.0.0_dp) then

       ! ** We want the set centres to maintain their fractional position
       
       ! * Convert to fractional

       do n = 1,num_sets
          ion_set_centre(:,n) = matmul(lattice_rec_orig,ion_set_centre(:,n))
       end do

       ! ** Set the cell shake matrix

       rn=(random_triple()-0.5_dp)*2.0_dp*cellamp
       do i=1,3
          mat(i,i) = 1.0_dp+rn(i)
       end do
       rn=(random_triple()-0.5_dp)*2.0_dp*cellamp
       mat(1,2) = rn(1)/2.0_dp
       mat(2,1) = mat(1,2)
       mat(1,3) = rn(2)/2.0_dp
       mat(3,1) = mat(1,3)
       mat(2,3) = rn(3)/2.0_dp
       mat(3,2) = mat(2,3)
       
       ! ** Shake the matrix
    
       lattice_car = matmul(mat,lattice_car_orig)
       
       ! ** Convert to ABC

       lattice_abc(1) = sqrt(lattice_car(1,1)**2+lattice_car(2,1)**2+lattice_car(3,1)**2)
       lattice_abc(2) = sqrt(lattice_car(1,2)**2+lattice_car(2,2)**2+lattice_car(3,2)**2)
       lattice_abc(3) = sqrt(lattice_car(1,3)**2+lattice_car(2,3)**2+lattice_car(3,3)**2)
       lattice_abc(4) = acos(dot_product(lattice_car(:,2),lattice_car(:,3))/lattice_abc(2)/lattice_abc(3))/dgrd
       lattice_abc(5) = acos(dot_product(lattice_car(:,1),lattice_car(:,3))/lattice_abc(1)/lattice_abc(3))/dgrd
       lattice_abc(6) = acos(dot_product(lattice_car(:,1),lattice_car(:,2))/lattice_abc(1)/lattice_abc(2))/dgrd

       ! ** Put back into absolute coordinates
       do n = 1,num_sets
          ion_set_centre(:,n) = matmul(lattice_car,ion_set_centre(:,n))
       end do
      
    elseif(.not.fix_cell) then

       ! ** Make a random unit cell from scratch

       lattice_tmp=lattice_abc
       
222    ang = -1.0_dp
       do while(ang<0.0_dp)

          ! ** Modify cell constraints

          if(con_cell) call update_canglen(stat)
          
          ! ** Permute the lattice index
          
          do i=1,3
             nindx(i)=i
          end do
          call perm(nindx)

          ! ** Set the cell lengths

          rn=random_triple()
          
          lattice_abc(nindx(1)) = len_min_use
          lattice_abc(nindx(2)) = len_min_use+(len_max_use-len_min_use)*rn(1)
          lattice_abc(nindx(3)) = len_min_use+(lattice_abc(nindx(2))-len_min_use)*rn(2)

          do i=1,3
             if(clen(i)<0.0_dp) then 
                temp=lattice_abc(i)
                exit
             end if
          end do
          
          do i=1,3
             if(clen(i)<0.0_dp) lattice_abc(i)=temp
          end do
          
          ! ** Set the cell angles

          rn = random_triple()
          do i=1,3
             if(cang(i)<=0.0_dp) then
                lattice_abc(3+i) = ang_min+(ang_max-ang_min)*rn(i) 
             else if(cang(i)==999.0_dp) then
                temp = lattice_abc(i)*(lattice_abc(1)**2+lattice_abc(2)**2+lattice_abc(3)**2-2.0_dp*lattice_abc(i)**2)&
                     /2.0_dp/lattice_abc(1)/lattice_abc(2)/lattice_abc(3)
                if(abs(temp).gt.1.0_dp) goto 222
                angle = acos(temp)/dgrd ! ** This can't be constrained to max/min angle 
                lattice_abc(3+i) = angle
             else
                lattice_abc(3+i) = cang(i)
             end if
          end do
          
          if(any(cang==9999.0_dp)) then
             angle=huge(1.0_dp)
             do while((angle.lt.ang_min).or.(angle.gt.ang_max))
                temp=9999.0_dp
                do while(abs(temp).gt.1.0_dp)
                   rn = random_triple()
                   lattice_abc(3+1) = ang_min+(ang_max-ang_min)*rn(1)
                   lattice_abc(3+3) = ang_min+(ang_max-ang_min)*rn(3)
                   temp = -(cos(lattice_abc(3+1)*dgrd)+cos(lattice_abc(3+3)*dgrd)+1.0_dp)
                end do
                angle=acos(temp)/dgrd
             end do
             lattice_abc(3+2) = angle
          end if

          if(any(cang==99999.0_dp)) then

             angle=huge(1.0_dp)
             do while((angle.lt.ang_min).or.(angle.gt.ang_max))
                temp=99999.0_dp
                do while(abs(temp).gt.1.0_dp)
                   rn = random_triple()
                   lattice_abc(3+1) = ang_min+(ang_max-ang_min)*rn(1)
                   lattice_abc(3+2) = ang_min+(ang_max-ang_min)*rn(1)
                   temp = -(cos(lattice_abc(3+1)*dgrd)+cos(lattice_abc(3+2)*dgrd)+1.0_dp)
                end do
                angle=acos(temp)/dgrd 
             end do
             lattice_abc(3+3) = angle 

          end if

          do i=1,3
             if(cang(i)<0.0_dp) then
                temp=lattice_abc(3+i)
                exit
             end if
          end do

          do i=1,3
             if(cang(i)<0.0_dp) lattice_abc(3+i)=temp
          end do

          ang = 1.0_dp-cos(dgrd*lattice_abc(4))**2-cos(dgrd*lattice_abc(5))**2-cos(dgrd*lattice_abc(6))**2+&
               2.0_dp*cos(dgrd*lattice_abc(4))*cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6))

       end do
       
       ! ** Reject a cell that is too flat

       if(sqrt(abs(ang)).lt.acons) goto 222
       
       ! ** We should now have a valid cell - scale to the required volume
       
       vol = lattice_abc(1)*lattice_abc(2)*lattice_abc(3)*sqrt(abs(ang))
       rn = random_triple()
       if(fix_vol) then
          rn(1) = 1.0_dp
       else
          rn(1) = vol_min+(vol_max-vol_min)*rn(1)
       end if
       
       lattice_abc(1:3) = lattice_abc(1:3)*(num_symm*targvol*rn(1)/vol)**(1.0_dp/3.0_dp)

       if(fix_caxis) then
          vol = lattice_abc(1)*lattice_abc(2)*lattice_abc(3)*sqrt(abs(ang))
          lattice_abc(3)   = lattice_tmp(3)
          lattice_abc(4:5) = lattice_tmp(4:5)
          ang = 1.0_dp-cos(dgrd*lattice_abc(4))**2-cos(dgrd*lattice_abc(5))**2-cos(dgrd*lattice_abc(6))**2+&
               2.0_dp*cos(dgrd*lattice_abc(4))*cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6))
          vol_tmp = lattice_abc(1)*lattice_abc(2)*lattice_abc(3)*sqrt(abs(ang))
          lattice_abc(1:2) = sqrt(vol/vol_tmp)*lattice_abc(1:2)          
       end if

       if(fix_abaxis) then
          vol = lattice_abc(1)*lattice_abc(2)*lattice_abc(3)*sqrt(abs(ang))
          lattice_abc(1:2) = lattice_tmp(1:2)
          lattice_abc(6)   = lattice_tmp(6)
          ang = 1.0_dp-cos(dgrd*lattice_abc(4))**2-cos(dgrd*lattice_abc(5))**2-cos(dgrd*lattice_abc(6))**2+&
               2.0_dp*cos(dgrd*lattice_abc(4))*cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6))
          vol_tmp = lattice_abc(1)*lattice_abc(2)*lattice_abc(3)*sqrt(abs(ang))
          lattice_abc(3) = (vol/vol_tmp)*lattice_abc(3)          
       end if

       ! ** Convert to cartesian from ABC
       
       lattice_car(:,1) = (/lattice_abc(1),0.0_dp,0.0_dp/)
       lattice_car(:,2) = (/lattice_abc(2)*cos(dgrd*lattice_abc(6)),lattice_abc(2)*sin(dgrd*lattice_abc(6)),0.0_dp/)
       lattice_car(1,3) = lattice_abc(3)*cos(dgrd*lattice_abc(5))
       lattice_car(2,3) = lattice_abc(3)*(cos(dgrd*lattice_abc(4)) &
            -cos(dgrd*lattice_abc(5))*cos(dgrd*lattice_abc(6)))/sin(dgrd*lattice_abc(6))
       lattice_car(3,3) = sqrt(lattice_abc(3)**2-lattice_car(1,3)**2-lattice_car(2,3)**2)
       
       cell_consistent=.false.
       
    end if

    call compact_cell()

    call update_cell()

  end subroutine gen_lattice

  ! *******************************************************************************************

  subroutine gen_position_check(stat)

    integer, intent(inout) :: stat

    integer, allocatable, dimension(:)   :: p
    integer, allocatable, dimension(:,:) :: nimton

    real(kind=dp) :: max_ion_min2,occsum,temp

    real(kind=dp), allocatable, dimension(:)     :: distall
    real(kind=dp), allocatable, dimension(:,:,:) :: join_lookup

    integer :: npow = 2, npbc, num_image, nim, nim_min, nk,nc, nattempt
    
    if(stat/=0) return

    ni=0
    do n=1,num_sets
       do ns=1,num_symm
          do m=1,ninset(n)

             ni=ni+1

             ion_occ(ni)   = ion_set_occ(m,n)

          end do
       end do
    end do
    
    ! ** Initialise some arrays the first time through

    if(.not.initialised) then
       ni=0
       ion_min=0
       do n=1,num_sets
          do ns=1,num_symm
             do m=1,ninset(n)

                ni=ni+1

                ion_names(ni)   = ion_set_names(m,n)
                ion_coord(:,ni) = ion_set_coord(:,m,n)
                ion_nn(:,ni)    = ion_set_nn(:,m,n)
                ion_rad(ni)     = ion_set_rad(m,n)
                ion_min(ni)     = ion_set_min(m,n)
                ion_fix(ni)     = ion_set_fix(m,n)
                ion_bfix(ni)    = ion_set_bfix(m,n)
!!$                ion_occ(ni)     = ion_set_occ(m,n) !! ** See above - reset each time
                ion_mult(ni)    = ion_set_mult(m,n)

             end do
          end do
       end do

       if(allocated(thresh)) deallocate(thresh)
       allocate(thresh(3,num_ions*num_symm,num_ions*num_symm))
       
       max_ion_min2 = maxval(ion_min)**2

       do ni=1,num_ions*num_symm 
          do nj=1,num_ions*num_symm

             if((ion_rad(ni).gt.0.0_dp).and.(ion_rad(nj).gt.0.0_dp)) then
                thresh(1,ni,nj) = (ion_rad(ni)+ion_rad(nj))**2
             else
                thresh(1,ni,nj) = (min(ion_min(ni),ion_min(nj)))**2
             end if

             ! ** Defaults

             do m=1,num_pairs
                if(((strip(ion_names(ni)).eq.pair_names(1,m)).and.(strip(ion_names(nj)).eq.pair_names(2,m)))&
                     .or.((strip(ion_names(nj)).eq.pair_names(1,m)).and.(strip(ion_names(ni)).eq.pair_names(2,m)))) then
                   thresh(1,ni,nj)=pair_minsep(m)**2
                end if
             end do

             ! ** Modified

             do m=1,num_pairs
                if(((ion_names(ni).eq.pair_names(1,m)).and.(ion_names(nj).eq.pair_names(2,m)))&
                     .or.((ion_names(nj).eq.pair_names(1,m)).and.(ion_names(ni).eq.pair_names(2,m)))) then
                   thresh(1,ni,nj)=pair_minsep(m)**2
                end if
             end do

             ! ** Apply slack, and other radii

             temp = thresh(1,ni,nj)
             thresh(1,ni,nj) = (1.0_dp-slack)**2*temp
             thresh(2,ni,nj) = ((1.4_dp**2-1.25_dp**2)*slack**2+1.25_dp**2)*temp
             thresh(3,ni,nj) = 1.4_dp**2*temp

          end do
       end do
       

       ! ** Slacken the bond angles
       
       minbangle=(1.0_dp-slack)*minbangle
       maxbangle=maxbangle+(180.0_dp-maxbangle)*slack
       
       ! ** Work out the species

       species_names=''
       ion_species=0
       num_spec=0

       do ni=1,num_ions*num_symm

          if(.not.any(ion_names(ni).eq.species_names)) then
             num_spec=num_spec+1
             species_names(num_spec)=ion_names(ni)
          end if

          do m=1,num_spec
             if(ion_names(ni).eq.species_names(m)) ion_species(ni)=m
          end do

       end do

       ! ** Set the pairs

       if(num_pairs.eq.0) then
          do n=1,num_spec
             do m=n,num_spec
                num_pairs=num_pairs+1
                pair_names(1,num_pairs)=species_names(n)
                pair_names(2,num_pairs)=species_names(m)
                rn = random_triple()
                pair_minsep(num_pairs)=minmin+rn(1)*(minmax-minmin)
             end do
          end do
       end if

       ! ** Set flags

       initialised=.true.

    end if

    ! ** Generate the ionic positions

    ni=0
    do n=1,num_sets

       nattempt = 0

       ! ** Random shift vector

333    nattempt = nattempt + 1

       free=.false.
       if(posamp_set(n)<0.0_dp) then
          rn = random_triple()!-0.5_dp
          rn(:) = rn(1)*lattice_car(:,1) + rn(2)*lattice_car(:,2) + rn(3)*lattice_car(:,3)
          free=.true.
       else
          rn  = random_triple()
          an  = random_triple()
          mat = rotation_matrix(an)
          ! ** If npow=2 the distribution is uniform in 3D space
          rn(1) = minamp_set(n)**(npow+1)+rn(1)*(posamp_set(n)**(npow+1)-minamp_set(n)**(npow+1))
          rn(1) = rn(1)**(1.0_dp/(npow+1))
          rn(2) = 0.0_dp ; rn(3) = 0.0_dp
          rn = matmul(mat,rn)
          
          if(ellipsoid) then
             rn=rn*ellipsvec
          end if
          
       end if

       if((xamp_set(n).ge.0.0_dp).or.(yamp_set(n).ge.0.0_dp).or.(zamp_set(n).ge.0.0_dp)) then
          free=.false.
          rn = random_triple()-0.5_dp
          rn(:) = rn(1)*lattice_car(:,1) + rn(2)*lattice_car(:,2) + rn(3)*lattice_car(:,3)

          an = random_triple()-0.5_dp

          if(xamp_set(n).ge.0.0_dp) then
             rn(1) = xamp_set(n)*an(1)
          end if
          if(yamp_set(n).ge.0.0_dp) then
             rn(2) = yamp_set(n)*an(2)
          end if
          if(zamp_set(n).ge.0.0_dp) then
             rn(3) = zamp_set(n)*an(3)
          end if
       end if

       if (ninset(n).gt.1) then

          ! ** Random rotation matrix to rotate set

          an = random_triple()

          if(angamp_set(n).ge.0.0_dp) an = an*(angamp_set(n)/360.0_dp)

          mat = rotation_matrix(an)

          ! ** Random mirror image
          
          if(flip) then

             if(random_single().lt.0.5_dp) then
                m=1+int(random_single()*3.0_dp)
                mirror=ident
                mirror(m,m)=-1.0_dp
                mat=matmul(mirror,mat)
             end if
             
          end if

       end if

       ! ** Find the new positions (absolute)

       do ns=1,num_symm

          do m=1,ninset(n)
             ni=ni+1

             if(ninset(n).gt.1) then

                if(free) then
                   ion_new_positions(:,ni) = rn(:) + matmul(mat,ion_set_positions(:,m,n))
                else
                   ion_new_positions(:,ni) = ion_set_centre(:,n) + rn(:) + matmul(mat,ion_set_positions(:,m,n))
                end if

             else

                if((posamp_set(n)<0.0_dp).and.(xamp_set(n)<0.0_dp).and.(yamp_set(n)<0.0_dp).and.(zamp_set(n)<0.0_dp)) then
                   ion_new_positions(:,ni) = rn(:)
                else
                   ion_new_positions(:,ni) = ion_set_centre(:,n) + rn(:) + ion_set_positions(:,m,n)
                end if

             end if

             if(ns.gt.1) then

                V = 0.0_dp
                do n1=1,3
                   do n2=1,3
                      V(n1) = V(n1)+Sym(n1,n2,ns)*ion_new_positions(n2,ni)
                   end do
                end do

                ion_new_positions(:,ni) = V(:) + symm_ops(1,4,ns)*lattice_car(:,1) + &
                     symm_ops(2,4,ns)*lattice_car(:,2) + symm_ops(3,4,ns)*lattice_car(:,3) 

             end if
             
             ! ** Check, and reject if necessary

             if(.not.push) then

                do i=1,ni

                   VO(:) = ion_new_positions(:,i) - ion_new_positions(:,ni)      

                   dist2 = VO(1)**2+VO(2)**2+VO(3)**2

                   if ((ion_occ(i).gt.1.0_dp-delta).and.(ion_occ(ni).gt.1.0_dp-delta)) then

                      if(dist2<thresh(1,i,ni)) then
                         rejected(1)=rejected(1)+1
                         if(i<=ni-m) then
                            ni=ni-m-(ns-1)*ninset(n)
                            if(nattempt>100) then
                               stat=999
                               return
                            end if
                            goto 333
                         end if
                      end if

                   end if

                end do

             end if

          end do
       end do

    end do
    
    ! ** Combine fractional atoms in a symmetry star to make full atoms

    if(allocated(taken)) deallocate(taken)
    allocate(taken(num_symm))
    
    ni=0
    do n=1,num_ions

       if(ion_occ(ni+1).gt.1.0_dp-delta) then

          ! ** No reduction required

          ni=ni+num_symm

       elseif (ion_occ(ni+1).lt.-delta) then

          ! ** Symmetry applied approximately - general positions only

          ion_occ(ni+1)=1.0_dp
          ion_occ(ni+2:ni+num_symm)=0.0_dp

          ni=ni+num_symm

       elseif (num_symm.gt.1) then
          
          ! ** Generate lookup table of shifted vectors

          npbc=2
          if(cluster)  npbc=0
          num_image=(2*npbc+1)**3
          if(allocated(join_lookup)) deallocate(join_lookup,nimton,distall)
          allocate(join_lookup(3,num_image,num_symm),nimton(3,num_image),distall(num_image))

          do ns=1,num_symm
             nim=0
             do n1=-npbc,npbc
                do n2=-npbc,npbc
                   do n3=-npbc,npbc
                      nim=nim+1
                      nimton(1,nim)=n1 ; nimton(2,nim)=n2 ; nimton(3,nim)=n3
                      join_lookup(:,nim,ns) = ion_new_positions(:,ni+ns)+n1*lattice_car(:,1)&
                           +n2*lattice_car(:,2)+n3*lattice_car(:,3)  
                   end do
                end do
             end do

          end do

          ! ** Calculate the required multiplicities, and number of partial atoms to be joined 

          mult  = nint(ion_occ(ni+1)*num_symm)
          njoin = nint(1.0_dp/ion_occ(ni+1))
          
!!$          ! ** Set the new occupations 
!!$
!!$          ion_occ(ni+1:ni+mult)=1.0_dp
!!$          ion_occ(ni+mult+1:ni+num_symm)=0.0_dp

          ! ** Join up the atoms

          taken=.false.

          nim_min=0

          do nm=1,mult

             VO(:)=0.0_dp

             do nj=1,njoin

                min2=huge(1.0_dp)

                do ns=1,num_symm

                   if(.not.taken(ns)) then

                      do nim=1,num_image

                         distall(nim)=(VO(1)-join_lookup(1,nim,ns))*(VO(1)-join_lookup(1,nim,ns))&
                              +(VO(2)-join_lookup(2,nim,ns))*(VO(2)-join_lookup(2,nim,ns))&
                              +(VO(3)-join_lookup(3,nim,ns))*(VO(3)-join_lookup(3,nim,ns))

                         if(distall(nim).lt.min2) then
                            jmin=ns ; min2=distall(nim) ; nim_min=nim
                         end if

                      end do

                   end if

                end do

                VO(:) = (VO(:)*real(nj-1,dp)+join_lookup(:,nim_min,jmin))/real(nj,dp)

                taken(jmin)=.true.

             end do

             ! ** Set the joined atom

             ion_new_positions(:,ni+nm) = VO(:)

          end do

          ! ** Check that the joined up star obeys the symmetry operations
          
          if(.not.check_star(ion_new_positions(:,ni+1:ni+mult))) then
             stat=999
             rejected(2)=rejected(2)+1
             return
          end if

          ! ** Set the new occupations (ion_mult as well?)

          ion_occ(ni+1:ni+mult)=1.0_dp
          ion_occ(ni+mult+1:ni+num_symm)=0.0_dp
          
          ! ** Set the counter to the next star

          ni=ni+num_symm  

       else
          continue
       end if

    end do
    
    if(any(ion_occ.lt.1.0_dp-delta)) then
       
       ! ** Attempt to join atoms with fractional occupancy

       allocate(p(num_ions*num_symm))
       do ni=1,num_ions*num_symm
          p(ni)=ni
       end do

       call perm(p)

       do ni=1,num_ions*num_symm
          if(ion_occ(p(ni)).lt.delta) cycle
          occsum=0.0_dp
          do nj=ni,num_ions*num_symm
             if(ion_names(p(ni)).ne.ion_names(p(nj))) cycle
             occsum=occsum+ion_occ(p(nj))
             if(abs(nint(occsum)-occsum).lt.delta) exit
          end do

          if(abs(nint(occsum)-occsum).gt.delta) then
             stat=999
             return
          else
             nc=0
             do nk=ni,nj
                if(ion_names(p(ni)).ne.ion_names(p(nk))) cycle

                nc=nc+1

                if(nc.le.nint(occsum)) then
                   ion_occ(p(nk))=1.0_dp

                else
                   ion_occ(p(nk))=0.0_dp
                end if
             end do
          end if

       end do

       deallocate(p)

    end if

    ! ** Check that the structure has the correct symmetry

    call symm_equiv(stat)
    if(stat/=0) return   
    
    ! ** Apply a symmetry breaking random displacement

    if (breakamp.gt.0.0_dp) then
       do ni=1,num_ions*num_symm
          if(ion_occ(ni).lt.1.0_dp-delta) cycle
          rn  = random_triple()
          an  = random_triple()
          mat = rotation_matrix(an)
          ! ** If npow=2 the distribution is uniform in 3D space
          rn(1) = rn(1)*breakamp**(npow+1)
          rn(1) = rn(1)**(1.0_dp/(npow+1))
          rn(2) = 0.0_dp ; rn(3) = 0.0_dp
          rn = matmul(mat,rn)
          ion_new_positions(:,ni) = ion_new_positions(:,ni) + rn(:)
       end do
       ! ** PROBLEM - should change to P1?
    end if

    ! ** Check that there are no fractional occupancies

    do ni=1,num_ions*num_symm
       if((ion_occ(ni).gt.delta).and.(ion_occ(ni).lt.1.0_dp-delta)) then
          stat=999
          return
       end if
    end do
    
    ! ** Remove any atoms at the same site

    if (remove) then
       
       do ni=1,num_ions*num_symm
          if(ion_occ(ni).lt.delta) cycle
          do nj=1,num_ions*num_symm
             if((ni.eq.nj).or.(ion_occ(nj).lt.delta)) cycle
             VO=ion_new_positions(:,ni)-ion_new_positions(:,nj)
             VO = matmul(lattice_car,matmul(lattice_rec,VO)-floor(matmul(lattice_rec,VO)))
             do n=1,27
                if(dot_product(VO(:)-lookup(:,n),VO(:)-lookup(:,n)).lt.1e-6_dp) ion_occ(nj)=0.0_dp
             end do
             
          end do
       end do

    end if
       
    ! ** Count the actual atoms remaining

    num_ions_total=0
    do ni=1,num_ions*num_symm
       if(ion_occ(ni).lt.1.0_dp-delta) cycle
       num_ions_total=num_ions_total+1
    end do

!!$    write (stderr,*) num_ions,num_ions_total
    
  end subroutine gen_position_check

  ! *******************************************************************************************

  subroutine check_separation(stat)

    integer, intent(inout) :: stat

    integer :: nc3,m,nk

    real(kind=dp) :: vec(3),vep(3),v12(3),dminang2,dmaxang2,distpln,vn2
    real(kind=dp) :: cmax,cmin,origin(3),normal(3),distsph2,rad2,crad2,widthby2

    logical :: docoord

    if(stat/=0) return

    rad2=radius**2
    crad2=coreradius**2
    widthby2=width/2.0_dp

    normal=(/0.0_dp,0.0_dp,1.0_dp/)
    origin=(/0.0_dp,0.0_dp,shift/)

    cmax=cos(maxbangle*dgrd)
    cmin=cos(minbangle*dgrd)

    call translate_to_cell()

    ion_fails=0
    s=0.0_dp

    do ni=1,num_ions*num_symm ! ** No need to use symmetry
       if(ion_occ(ni).lt.1.0_dp-delta) cycle
       do nj=1,num_ions*num_symm ! ** No need to use symmetry
          if(ion_occ(ni).lt.1.0_dp-delta) cycle
          ion_ion_vec(:,ni,nj)=0.0_dp
       end do
    end do

    docoord=any(ion_coord.ge.0)
    
    do ni=1,num_ions*num_symm ! ** No need to use symmetry

       if(ion_occ(ni).lt.1.0_dp-delta) cycle

       nc3=0
       distlist2=huge(1.0_dp)
       do nj=1,num_ions*num_symm

          if(ion_occ(nj).lt.1.0_dp-delta) cycle

          if((.not.ion_moved(ni)).and.(.not.ion_moved(nj)).and.push.and.(.not.docoord)) cycle

          VO(:) = ion_new_positions(:,nj)-ion_new_positions(:,ni)

          do n1=1,27

             if((n1.ne.14).and.(cluster)) cycle

             if((n1.eq.14).and.(ion_ion_set(ni).eq.(ion_ion_set(nj)))) cycle

             vec(:)=lookup(:,n1)+VO(:)

             dist2 = dot_product(vec,vec)
             
             ! ** Reject or push if too close

             if(push) then
                ! ** Push
                if(dist2.lt.thresh(3,ni,nj)) then
                   if(dist2.lt.thresh(1,ni,nj)) then
                      stat=1
                      ion_fails(ni)=ion_fails(ni)+1
                      vn2=sqrt(dot_product(vec,vec))
                      if(vn2.lt.delta) then
                         vec=random_triple()*1e-8_dp
                         vn2=sqrt(dot_product(vec,vec))
                      end if
                      v12=vec*sqrt(thresh(1,ni,nj)-dist2)/2.0_dp/vn2
                      if(ion_bfix(ni)) v12=0.0_dp
                      if(ion_bfix(nj)) v12=v12*2.0_dp
                      ion_ion_vec(:,ni,nj)=ion_ion_vec(:,ni,nj)+v12(:)
                      call update_s(v12,vec(:))
                   end if
                   if((any(ion_coord(:,ni).ge.0))) then
                      nc3=nc3+1
                      distlist2(nc3)=dist2
                      ionlist(nc3)=nj
                      veclist(:,nc3)=vec(:)
                   end if
                end if
             else
                ! ** Just reject
                if(dist2.lt.thresh(1,ni,nj)) then
                   stat=1
                   rejected(3)=rejected(3)+1
                   return
                end if
             end if


          end do

       end do

       ! ** Add coordination push

       if(any(ion_coord(:,ni).ge.0).and..true.) then

          do n=1,nc3
             indx(n)=n
          end do

          call heap_sort_index(nc3,distlist2(1:nc3),indx(1:nc3))

          if(nc3.lt.ion_coord(1,ni)) ion_fails(ni)=ion_fails(ni)+1

          do n=1,nc3
             nj=ionlist(indx(n))
             if(n.le.ion_coord(1,ni)) then
                if (distlist2(indx(n)).ge.thresh(2,ni,nj)) then
                   ion_fails(ni)=ion_fails(ni)+1
                   vec=veclist(:,indx(n))
                   v12=-vec*sqrt(distlist2(indx(n))-thresh(2,ni,nj))/2.0_dp/sqrt(dot_product(vec,vec))
                   if(ion_bfix(ni)) v12=0.0_dp
                   if(ion_bfix(nj)) v12=v12*2.0_dp
                   ion_ion_vec(:,ni,nj)=ion_ion_vec(:,ni,nj)+v12(:)                      
                   call update_s(v12,vec)
                end if

                do m=1,nc3
                   if(m.eq.n) cycle
                   nk=ionlist(indx(m))
                   if(m.le.ion_coord(1,ni)) then

                      vep=-veclist(:,indx(m))+veclist(:,indx(n))
                      dist2=dot_product(vep,vep)

                      dminang2=distlist2(indx(n))+distlist2(indx(m))&
                           -2.0_dp*sqrt(distlist2(indx(n))*distlist2(indx(m)))*cmin
                      dmaxang2=distlist2(indx(n))+distlist2(indx(m))&
                           -2.0_dp*sqrt(distlist2(indx(n))*distlist2(indx(m)))*cmax

                      if(dist2.le.dminang2) then
                         ion_fails(ni)=ion_fails(ni)+1
                         v12=-vep(:)*sqrt(dminang2-dist2)/2.0_dp/sqrt(dot_product(vep,vep))
                         if(ion_bfix(nj)) v12=0.0_dp
                         if(ion_bfix(nk)) v12=v12*2.0_dp
                         ion_ion_vec(:,nj,nk)=ion_ion_vec(:,nj,nk)+v12(:)
                         call update_s(v12,vep)
                      end if

                      if(dist2.ge.dmaxang2) then
                         ion_fails(ni)=ion_fails(ni)+1
                         v12=vep(:)*sqrt(dist2-dminang2)/2.0_dp/sqrt(dot_product(vep,vep))
                         if(ion_bfix(nj)) v12=0.0_dp
                         if(ion_bfix(nk)) v12=v12*2.0_dp
                         ion_ion_vec(:,nj,nk)=ion_ion_vec(:,nj,nk)+v12(:)
                         call update_s(v12,vep)
                      end if

                   end if
                end do

             else
                if (distlist2(indx(n)).le.thresh(2,ni,nj)) then
                   ion_fails(ni)=ion_fails(ni)+1
                   vec=veclist(:,indx(n))
                   v12=vec*sqrt(thresh(3,ni,nj)-distlist2(indx(n)))/2.0_dp/sqrt(dot_product(vec,vec))
                   if(ion_bfix(ni)) v12=0.0_dp
                   if(ion_bfix(nj)) v12=v12*2.0_dp
                   ion_ion_vec(:,ni,nj)=ion_ion_vec(:,ni,nj)+v12(:)                      
                   call update_s(v12,vec)
                end if
             end if
          end do

       end if

       ! ** Add sphere push

       if(sphere) then

          ! * Shortest distance to sphere centre

          if(cluster) then
             origin=0.0_dp
          else
             origin=(lattice_car(1:3,1)+lattice_car(1:3,2)+lattice_car(1:3,3))/2.0_dp
          end if

          VO(:) = ion_new_positions(:,ni)-origin

          distsph2=huge(1.0_dp)
          do n1=1,27

             if((n1.ne.14).and.(cluster)) cycle

             vec(:)=VO(:)-lookup(:,n1)

             if(ellipsoid) vec(:)=vec(:)/ellipsvec(:)
             
             dist2 = dot_product(vec,vec)

             if(dist2.lt.distsph2) then
                vep=vec
                distsph2=dist2
             end if

          end do

          ! * Confining sphere

          if (distsph2.gt.rad2) then
             ion_fails(ni)=ion_fails(ni)+1
             v12=vep*sqrt(distsph2-rad2)/sqrt(dot_product(vep,vep))
             if(ion_bfix(ni)) v12=0.0_dp
             ion_ion_vec(:,ni,ni)=ion_ion_vec(:,ni,ni)+v12(:)                      
          end if

          ! * Repulsive core

          if ((distsph2.lt.crad2).and.core) then
             ion_fails(ni)=ion_fails(ni)+1
             v12=-vep*sqrt(crad2-distsph2)/sqrt(dot_product(vep,vep))
             if(ion_bfix(ni)) v12=0.0_dp
             ion_ion_vec(:,ni,ni)=ion_ion_vec(:,ni,ni)+v12(:)                      
          end if

       end if

       ! ** Add cylinder push

       if(cylinder) then
          
          if(cluster) then
             origin=0.0_dp
          else
             origin=(lattice_car(1:3,1)+lattice_car(1:3,2)+lattice_car(1:3,3))/2.0_dp
          end if

          VO(:) = ion_new_positions(:,ni)-origin

          distsph2=huge(1.0_dp)
          do n1=1,27

             if((n1.ne.14).and.(cluster)) cycle

             vec(:)=lookup(:,n1)+VO(:)

             dist2 = dot_product(vec(1:2),vec(1:2))

             if(dist2.lt.distsph2) then
                vep=vec
                distsph2=dist2
             end if

          end do

          ! * Confining cylinder - axis parallel to z

          if (distsph2.gt.rad2) then
             ion_fails(ni)=ion_fails(ni)+1
             v12(1:2)=vep(1:2)*sqrt(distsph2-rad2)/sqrt(dot_product(vep(1:2),vep(1:2)))
             if(ion_bfix(ni)) v12=0.0_dp
             ion_ion_vec(1:2,ni,ni)=ion_ion_vec(1:2,ni,ni)+v12(1:2)                      
          end if

          ! * Inner core cylinder - axis parallel to z

          if ((distsph2.lt.crad2).and.core) then
             ion_fails(ni)=ion_fails(ni)+1
             v12(1:2)=-vep(1:2)*sqrt(crad2-distsph2)/sqrt(dot_product(vep(1:2),vep(1:2)))
             if(ion_bfix(ni)) v12=0.0_dp
             ion_ion_vec(1:2,ni,ni)=ion_ion_vec(1:2,ni,ni)+v12(1:2)                      
          end if

       end if

       ! ** Add in the plane push

       if(width.ge.0.0_dp) then

          ! * Shortest distance to plane

          distpln=huge(1.0_dp)

          do n1=1,27
             
             if((n1.ne.14).and.(cluster)) cycle

             VO(:) = ion_new_positions(:,ni)-origin
             
             vec(:)=lookup(:,n1)+VO(:)
             dist2 = dot_product(vec,normal)
             if(dist2.lt.abs(distpln)) then
                vep=vec
                distpln=dist2
             end if

          end do
          
          if(abs(distpln).gt.widthby2) then
             ion_fails(ni)=ion_fails(ni)+1
             v12(3)=vep(3)*sqrt(distpln**2-widthby2**2)/abs(vep(3))
             if(ion_bfix(ni)) v12=0.0_dp
             ion_ion_vec(3,ni,ni)=ion_ion_vec(3,ni,ni)+v12(3)   
          end if

       end if

    end do

    if(sum(ion_fails).gt.0) stat=1

  end subroutine check_separation

  subroutine update_s(f12,r12)

    real(kind=dp), dimension(3), intent(in) :: f12,r12
    
    s(1,1)=s(1,1)+f12(1)*r12(1)
    s(2,2)=s(2,2)+f12(2)*r12(2)
    s(3,3)=s(3,3)+f12(3)*r12(3)
    
    s(2,3)=s(2,3)+(f12(2)*r12(3)+f12(3)*r12(2))/2.0_dp
    s(3,2)=s(3,2)+(f12(3)*r12(2)+f12(2)*r12(3))/2.0_dp
    
    s(3,1)=s(3,1)+(f12(3)*r12(1)+f12(1)*r12(3))/2.0_dp
    s(1,3)=s(1,3)+(f12(1)*r12(3)+f12(3)*r12(1))/2.0_dp
    
    s(1,2)=s(1,2)+(f12(1)*r12(2)+f12(2)*r12(1))/2.0_dp
    s(2,1)=s(2,1)+(f12(2)*r12(1)+f12(1)*r12(2))/2.0_dp
    
  end subroutine update_s
  
  subroutine check_coordination(stat)

    integer, intent(inout) :: stat

    integer :: nc

    real(kind=dp) :: arg

    stop 'check_coordination deactivated'
    
    if(stat/=0) return

    do ni=1,num_ions*num_symm ! ** Cannot just do one of each copy

       if(ion_occ(ni).lt.1.0_dp-delta) cycle

       if(any(ion_coord(:,ni).ge.0)) then

          nc=0

          do nj=1,num_ions*num_symm

             if(ion_occ(nj).lt.1.0_dp-delta) cycle

             if(any(ion_coord(:,nj).lt.0)) cycle

             VO(:) = ion_new_positions(:,nj)-ion_new_positions(:,ni)

             do n1=1,27

                if((n1.ne.14).and.(cluster)) cycle

                if((ni.eq.nj).and.(n1.eq.14)) cycle

                if(molecules.and.(n1.eq.14).and.(ion_ion_set(ni).eq.(ion_ion_set(nj)))) cycle

                dist2 = (lookup(1,n1)+VO(1))**2+(lookup(2,n1)+VO(2))**2+(lookup(3,n1)+VO(3))**2

                ! ** Reject if too close

                if(dist2.lt.thresh(1,ni,nj)) then
                   rejected(4)=rejected(4)+1
                   stat=1
                   return
                end if

                ! ** Clear a gap beyond the nearest neighbours

                if((dist2.ge.thresh(2,ni,nj)).and.(dist2.lt.thresh(3,ni,nj))) then
                   rejected(5)=rejected(5)+1
                   stat=1
                   return
                end if

                ! ** Count nearest neighbours

                if((dist2.ge.thresh(1,ni,nj)).and.(dist2.lt.thresh(2,ni,nj))) then
                   nc=nc+1
                   
                   ! * Reject if already too many
                   if((nc.gt.maxval(ion_coord(:,ni)))) then
                      rejected(6)=rejected(6)+1
                      stat=1
                      return
                   end if

                   ! * Reject if wrong type
                   if(((len_trim(ion_nn(2,ni)).gt.0).and.(ion_names(nj).ne.ion_nn(2,ni))).or.&
                        (ion_names(nj).eq.ion_nn(1,ni))) then
                      rejected(7)=rejected(7)+1
                      stat=1
                      return
                   end if

                   indx_nearestn(nc) = nj
                   nearestn(:,nc) = ion_new_positions(:,nj)+lookup(:,n1)

                end if

             end do

          end do

          ! ** Check we have the right coordination
          
          if((nc.lt.ion_coord(1,ni)).or.(nc.gt.ion_coord(2,ni))) then
             rejected(8)=rejected(8)+1
             stat=1
             return
          end if
          
          ! ** Check the bond angles

          do n1=1,nc

             V = nearestn(:,n1)-ion_new_positions(:,ni)
             do n2=n1+1,nc

                VN = nearestn(:,n2)-ion_new_positions(:,ni)

                arg=dot_product(V,VN)/sqrt(dot_product(V,V))/sqrt(dot_product(VN,VN))
                
                arg=max(-1.0_dp,arg)
                arg=min(1.0_dp,arg)
                
                angle=acos(arg)/dgrd

                if(angle.lt.minbangle) then                   
                   rejected(9)=rejected(9)+1
                   stat=1
                   return
                end if

                if(angle.gt.maxbangle) then                 
                   rejected(10)=rejected(10)+1
                   stat=1
                   return
                end if
                
             end do
          end do

       end if

    end do

  end subroutine check_coordination

  subroutine push_vec(eint,vec)

    integer,       intent(out) :: eint
    real(kind=dp), intent(out) :: vec(3,num_ions*num_symm)

    ! **

    real(kind=dp) :: hold_vec(3,num_ions*num_symm)

    eint=0
    do ni=1,num_ions*num_symm
       if(.not.ion_bfix(ni)) eint=eint+ion_fails(ni)
       if(ion_fails(ni).gt.0) then
          ion_moved(ni)=.true.
          vec(:,ni)=0.0_dp
          do nj=1,num_ions*num_symm
             vec(:,ni) = vec(:,ni)-ion_ion_vec(:,ni,nj)
          end do          
          vec(:,ni)=vec(:,ni)*(1.0_dp+delta)
       else
          ion_moved(ni)=.false.
          vec(:,ni)=0.0_dp
       end if

    end do

    ! ** Symmetrise
    
    hold_vec = vec
    vec = 0.0_dp

    do ni=1,num_ions*num_symm
       if(ion_moved(ni)) then
          do ns=1,num_symm
             vec(:,ni) = vec(:,ni) + matmul(Sym(1:3,1:3,ns),hold_vec(:,ion_equiv(ns,ni)))
          end do
          vec(:,ni) = vec(:,ni)/real(num_symm,dp)
       end if
    end do
    

    ! ** All in unit move together

    do ni=1,num_ions*num_symm

       do nj=ni+1,num_ions*num_symm
          
          if(ion_moved(ni).and.(ion_ion_set(ni).eq.(ion_ion_set(nj)))) ion_moved(nj)=.true.

       end do

    end do

    ! ** Apply brace
    
    call brace(vec)
    
  end subroutine push_vec
  
  subroutine push_opt(stat)

    integer, intent(inout)  :: stat

    ! **

    integer       :: e,e_min,steps,n
 
    real(kind=dp) :: g(3,num_ions*num_symm),step

    if(stat/=0) return
 
    n=num_ions*num_symm*27
    if(allocated(ionlist)) deallocate(ionlist,indx,distlist2,veclist)
    allocate(ionlist(n),indx(n),distlist2(n),veclist(3,n))

    ion_moved=.true.

    call translate_to_cell()
    
    call symm_equiv(stat)

    if(stat/=0) return

    call check_separation(stat)
    
    if(.not.push) return

    g     = 1.0_dp
    e     = huge(1)
    e_min = e
    steps = 0
    step  = pushstep
    
    do while((e.gt.0).and.(sum(abs(g)).gt.delta).and.(steps.lt.pushmax))
       steps=steps+1
       
       call push_vec(e,g)
       
       do ni=1,num_ions*num_symm
          ion_new_positions(:,ni)  = ion_new_positions(:,ni) + step*g(:,ni)
       end do
       
       stat=0
       
       call push_snapshot()
       
       if(e.lt.e_min) then
          call push_report(e)
          e_min=e
          steps=0
       end if
       
       call check_separation(stat)

    end do
    
    if((e.eq.0).and.(stat.eq.0)) then
       write(stderr,'(a)',advance='no') '~'
    else
       write(stderr,'(a)',advance='no') 'X'
    end if
    
  end subroutine push_opt

  subroutine push_report(e)

    integer, intent(in) :: e

    select case(e)
    case(1000:)
       write(stderr,'(a)',advance='no') '*'
    case(100:999)
       write(stderr,'(a)',advance='no') '|'
    case(10:99)
       write(stderr,'(a)',advance='no') ':'
    case(:9)
       write(stderr,'(a)',advance='no') '-'
    end select

  end subroutine push_report

  subroutine push_snapshot()

    integer :: ni,nat
    character(len=10) :: ctemp
    
    if(.true.) return
    
    nat=0
    do ni=1,num_ions*num_symm
       if(ion_occ(ni).lt.1.0_dp-delta) cycle
       nat=nat+1
    end do
    
    write (ctemp,'(i10)') nat
    write (stdout,'(a)') trim(adjustl(ctemp))
    write (stdout,'(a,9f11.6,a)') 'Lattice="',reshape(lattice_car,(/9/)),'" Properties=species:S:1:pos:R:3'
    do ni=1,num_ions*num_symm
       if(ion_occ(ni).lt.1.0_dp-delta) cycle
       write (stdout,'(a,3f18.13)') trim(adjustl(ion_names(ni))),ion_new_positions(:,ni)
    end do
    flush(stdout)
    
  end subroutine push_snapshot
  
  subroutine assign_spin(stat)
    
    integer, intent(inout)  :: stat

    !-------------------------------
    
    integer :: n,nn,count
    real(kind=dp) :: temp,s1,s2,grd,incpt,step
    real(kind=dp), allocatable, dimension(:) :: swork

    if(stat/=0) return
        
    if(.not.havespin) return
    
    stat=888 ! Stat return code for failure

    n=num_ions*num_symm
    

    if(spin.gt.spinmod) stop 'spin > spinmod'
    
    spinvec=0.0_dp
    afmvec=0.0_dp

    nn=0
    do i=1,n
       if(ion_occ(i).lt.1.0_dp-delta) cycle
       if((index(spinlist,trim(ion_names(i))//' ').gt.0).or.(len_trim(spinlist).eq.0)) then
          nn=nn+1
          rn = random_triple()       
          spinvec(i)=rn(1)
          afmvec(i)=rn(2)
       endif
    end do
    
    s1=sum(afmvec(1:n))/real(nn,dp)
    s2=spin/sum(spinvec(1:n))*real(nn,dp)
    
    do i=1,n
       if(ion_occ(i).lt.1.0_dp-delta) cycle
       if((index(spinlist,trim(ion_names(i))//' ').gt.0).or.(len_trim(spinlist).eq.0)) then
          afmvec(i) = afmvec(i)-s1
          spinvec(i) = spinvec(i)*s2
       end if
    end do
    
    ! * Adjust so that the integrated and total spins are as requested
    
    step=1e-1_dp
    temp=step
    count=0
    do count=1,10000
       temp=temp*(1.0_dp+step)
       if(sum(abs(spinvec(1:n)+temp*afmvec(1:n)))/nn.ge.(spinmod-epsilon(1.0_dp))) then
          grd=(sum(abs(spinvec(1:n)+temp*(1.0_dp+step)*afmvec(1:n)))/nn-&
               sum(abs(spinvec(1:n)+temp*afmvec(1:n)))/nn)/(temp*(1.0_dp+step)-temp)
          incpt=sum(abs(spinvec(1:n)+temp*afmvec(1:n)))/nn-grd*temp
          if(abs(grd).gt.0.0_dp) temp=(spinmod-incpt)/grd
          spinvec(1:n)=spinvec(1:n)+temp*afmvec(1:n)
          exit
       end if
    end do

    ! * Symmetrise spins

    if(num_symm.gt.1) then
       
       allocate(swork(size(spinvec)))
       
       swork=spinvec
       spinvec=0.0_dp
       
       do i=1,n
          if(ion_occ(i).lt.1.0_dp-delta) cycle
          if((index(spinlist,trim(ion_names(i))).gt.0).or.(len_trim(spinlist).eq.0)) then
             do ns=1,num_symm
                spinvec(i) = spinvec(i) +  swork(ion_equiv(ns,i))
             end do
             spinvec(i) = spinvec(i)/real(num_symm,dp)
          end if
       end do
       
       deallocate(swork)

    end if

    ! ** Final check

    if((abs(spin*nn-sum(spinvec)).lt.1e-6_dp).and.(abs(spinmod*nn-sum(abs(spinvec))).lt.1e-6_dp)) then
       stat=0
    end if

    return
    
  end subroutine assign_spin
    
  subroutine translate_to_cell()

    real(kind=dp) :: R0(3)
    integer :: nin
    
    if(.not.translate) return

    ! ** Translate to unit cell
    
    ni=0
    do n=1,num_sets
       do ns=1,num_symm
          
          R0=0.0_dp
          nin=0
          do m=1,ninset(n)
             ni=ni+1
             if(ion_occ(ni).lt.1.0_dp-delta) cycle
             nin=nin+1             
             R0=R0+ion_new_positions(:,ni)
          end do
          
          if(nin.gt.0) R0=R0/real(nin,dp)  
          
          R0(:) = floor(matmul(lattice_rec,R0))
          R0(:) = R0(1)*lattice_car(:,1) + R0(2)*lattice_car(:,2) + R0(3)*lattice_car(:,3)
          
          ni=ni-ninset(n)
          
          do m=1,ninset(n)
             ni=ni+1
             
             if(ion_occ(ni).lt.1.0_dp-delta) cycle
             
             ion_new_positions(:,ni)  = ion_new_positions(:,ni) - R0(:)
             
          end do
          
       end do
    end do
    
  end subroutine translate_to_cell

  subroutine brace(force)

    real(kind=dp), dimension(:,:), intent(inout) :: force

    real(kind=dp) :: ftemp(3),R0(3),T(3),R(3),F(3)

    integer :: ni,n,ns,m

    ! ** Calculate the torque on a unit

    ni=0
    do n=1,num_sets
       do ns=1,num_symm

          ! * Determine the reference position, R0

          R0=0.0_dp
          do m=1,ninset(n)
             ni=ni+1
             R0=R0+ion_new_positions(:,ni)
          end do
          R0=R0/real(ninset(n),dp)

          ! * Evaluate the torque with respect to R0

          ni=ni-ninset(n)
          T=0.0_dp
          do m=1,ninset(n)
             ni=ni+1
             R(:)=ion_new_positions(:,ni)-R0(:)
             F(:)=force(:,ni)
             T(1)=T(1)+R(2)*F(3)-R(3)*F(2)
             T(2)=T(2)+R(3)*F(1)-R(1)*F(3)
             T(3)=T(3)+R(1)*F(2)-R(2)*F(1)

          end do

       end do
    end do    

    ! ** Calculate the overall force on a unit

    ni=0
    do n=1,num_sets
       do ns=1,num_symm

          ftemp=0.0_dp
          do m=1,ninset(n)
             ni=ni+1
             ftemp=ftemp+force(:,ni)
          end do

          do m=ni-ninset(n)+1,ni
             force(:,m) = ftemp/real(ninset(n),dp)
          end do

       end do
    end do

    ! ** Fix atoms

    ni=0
    do n=1,num_sets
       do ns=1,num_symm

          do m=1,ninset(n)
             ni=ni+1
             if(ion_bfix(ni)) force(:,ni)=0.0_dp
          end do
          
       end do
    end do


  end subroutine brace

  subroutine symm_equiv(status)

    integer, intent(inout) :: status

    real(kind=dp) :: symm_thresh=1e-6_dp

    ion_equiv=0

    if(num_symm==1) then
       do ni=1,num_ions*num_symm
          ion_equiv(1,ni)=ni
       end do
       return 
    end if

    ! ** Which atoms are symmetry related to which?    

    do ni=1,num_ions*num_symm
       if(ion_occ(ni).lt.1.0_dp-delta) cycle

       do ns=1,num_symm

          do nj=1,num_ions*num_symm
             if(ion_occ(nj).lt.1.0_dp-delta) cycle
             if(ion_names(ni).ne.ion_names(nj)) cycle

             VO=ion_new_positions(:,ni)-matmul(Sym(1:3,1:3,ns),ion_new_positions(1:3,nj))-&
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

    ! ** Check that all symmetry operations map atoms to atoms

    do ni=1,num_ions*num_symm
       if(ion_occ(ni).lt.1.0_dp-delta) cycle
       if(any(ion_equiv(1:num_symm,ni).eq.0)) status=1
    end do
    
!!$    do ni=1,num_ions*num_symm
!!$       if(ion_occ(ni).lt.1.0_dp-delta) cycle
!!$       do ns=1,num_symm
!!$          write (stderr,'(a4,3i5)') ion_names(ni),ns,ni,ion_equiv(ns,ni)
!!$       end do
!!$    end do
!!$    write (stderr,*)
        
  end subroutine symm_equiv

  function rotation_matrix(rangle)
    
    real(kind=dp), intent(inout) :: rangle(3)

    real(kind=dp) :: rotation_matrix(3,3)

    !------

    real(kind=dp) :: R(3,3),VR(3),H(3,3)
    
    rangle(1) = tpi*rangle(1)
    rangle(2) = tpi*rangle(2)

    R = 0.0_dp
    
    R(1,1) =  cos(rangle(1)) ; R(2,1) = sin(rangle(1))
    R(1,2) = -sin(rangle(1)) ; R(2,2) = cos(rangle(1))
    R(3,3) =  1.0_dp
    
    VR(1) = cos(rangle(2))*sqrt(rangle(3))
    VR(2) = sin(rangle(2))*sqrt(rangle(3))
    VR(3) = sqrt(1.0_dp-rangle(3))
    
    H = 0.0_dp
    H(1,1) = 1.0_dp ; H(2,2) = 1.0_dp ; H(3,3) = 1.0_dp
    
    do i=1,3
       do j=1,3
          H(i,j) = H(i,j) - 2.0_dp*VR(i)*VR(j)
       end do
    end do
    
    ! ** Not totally sure about this, but need to make angamp=0 give me the identity
    
    H(3,:) = -1.0_dp*H(3,:)
    
    rotation_matrix = matmul(H,R)

  end function rotation_matrix

  function check_star(star)

    real(kind=dp), dimension(:,:), intent(in) :: star

    logical :: check_star

    real(kind=dp) :: V(3),V0(3),mind2

    integer :: nstar,ni,n1,n2,n3,npbc

    real(kind=dp) :: dthresh2=1e-10_dp

    if(cluster) then
       npbc=0
    else
       npbc=2
    end if
    
    check_star = .true.

    nstar = size(star(1,:))
    
    ni=1

    do ns=1,num_symm

       V0 = 0.0_dp
       do n1=1,3
          do n2=1,3
             V0(n1) = V0(n1)+Sym(n1,n2,ns)*star(n2,ni)
          end do
       end do
       
       V0(:) = V0(:) + symm_ops(1,4,ns)*lattice_car(:,1) + symm_ops(2,4,ns)*lattice_car(:,2) + &
            symm_ops(3,4,ns)*lattice_car(:,3) 

       mind2=huge(1.0_dp)
       do nj=1,nstar

          do n1=-npbc,npbc
             do n2=-npbc,npbc
                do n3=-npbc,npbc

                   V = V0-star(:,nj) &
                        - n1*lattice_car(:,1) - n2*lattice_car(:,2) - n3*lattice_car(:,3)  
                   if(dot_product(V,V)<mind2) mind2=dot_product(V,V)

                end do
             end do
          end do

       end do

       if(mind2.gt.dthresh2) then
          check_star=.false.
          return
       end if

    end do
    
  end function check_star

  subroutine perm(list)
    
    integer, dimension(:), intent(out) :: list

    integer :: i,j

    do i=1,size(list)
       j=1+int(random_single()*i)
       list(i)=list(j)
       list(j)=i
    end do
    
  end subroutine perm

  subroutine update_canglen(stat)

    integer, intent(inout) :: stat
    
    real(kind=dp) :: vec(3)
    
    integer :: n,m

    if(stat/=0) return

    if(any(cang.gt.180.0_dp)) then
       stat=1
       return
    end if
    
!!$    write (stderr,'(6f10.2)') clen,cang
!!$    write (stderr,'(6f10.2)') cellconvec
!!$    write (stderr,'(6f10.2)') 

    clen=min(clen,cellconvec(1:3))
    
    vec=cang
    
    do n=1,3
       if((vec(n).lt.delta).and.(cellconvec(3+n).lt.delta)) then
          cang(n)=min(vec(n),cellconvec(3+n))
       end if
       if((vec(n).gt.-delta).and.(cellconvec(3+n).gt.-delta)) then
          cang(n)=max(vec(n),cellconvec(3+n))
       end if
       if((vec(n).lt.-delta).and.(cellconvec(3+n).gt.delta)) then
          cang(n)=max(vec(n),cellconvec(3+n))
       end if
       if((vec(n).gt.delta).and.(cellconvec(3+n).gt.delta)) then
          if((abs(vec(n)-cang(n)).gt.delta).or.(abs(vec(n)-cellconvec(3+n)).gt.delta))then
             stat=1
             return
          end if
       end if
    end do
    
    do n=1,3
       if(cang(n).lt.-delta) then
          do m=1,n-1
             if(vec(m).lt.-delta) then
                cang(n)=cellconvec(3+m)
                exit
             end if
             ! * Remove this.Doesn't work for Rhom
!!$             if(cellconvec(3+m).lt.-delta) then
!!$                cang(n)=vec(m) !???
!!$                exit
!!$             end if
          end do
       end if
    end do

    do n=1,3
       if(cellconvec(3+n).lt.-delta) then
          do m=n+1,3
             if(cellconvec(3+m).lt.-delta) then
                if(abs(cang(n)-cang(m)).gt.delta) then
                   stat=1
                   return
                end if
             end if
          end do
       end if
       if(vec(n).lt.-delta) then
          do m=n+1,3
             if(vec(m).lt.-delta) then
                if(abs(cang(n)-cang(m)).gt.delta) then
                   stat=1
                   return
                end if
             end if
          end do
       end if
    end do
    
!!$    write (stderr,'(6f10.2)') clen,cang
!!$    write (stderr,*) '==='
!!$
!!$    stop
    
  end subroutine update_canglen
  
end module build
