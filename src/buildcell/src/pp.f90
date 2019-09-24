! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!==================================================================================!
!                                      pp                                          !
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
! This module reads, knows, and evaluates the pair potential                       !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module pp

  use constants
  use cell
  use symmetry

  implicit none

  private

  public :: init_pp
  public :: eval_pp

  ! ** Data defining the potential

  integer,                allocatable, dimension(:,:) :: PP_m, PP_n
  real(kind=dp)                                       :: range
  real(kind=dp), public,  allocatable, dimension(:,:) :: epsil,sigma,beta
  real(kind=dp),          allocatable, dimension(:,:) :: rcut, rcut_sq, E_shift, f_shift

  logical, public :: forceshift=.true.
  logical, public :: inipp=.false.

  !---------------------------------------------------!

  ! ** Private variables

  real(kind=dp) :: mod_r12, mod_r12_sq, inv_mod_r12_n, inv_mod_r12_m

contains

  subroutine init_pp

    integer       :: i,j,m
    real(kind=dp) :: mod_r12, inv_mod_r12_n, inv_mod_r12_m

    if(allocated(PP_m)) then
       deallocate(PP_m,PP_n,epsil,sigma,beta,E_shift,F_shift,rcut,rcut_sq)
    end if

    allocate(PP_m(num_spec,num_spec))
    allocate(PP_n(num_spec,num_spec))
    allocate(epsil(num_spec,num_spec))
    allocate(sigma(num_spec,num_spec))
    allocate(beta(num_spec,num_spec))
    allocate(E_shift(num_spec,num_spec))
    allocate(F_shift(num_spec,num_spec))
    allocate(rcut(num_spec,num_spec))
    allocate(rcut_sq(num_spec,num_spec))

    PP_m = 12
    PP_n = 6

    range = 1.0_dp       

    epsil = 1.0_dp
    sigma = minsep
    beta = 1.0_dp

    ! ** Set the the stripped species names

    do i=1,num_spec
       do j=1,num_spec
          do m=1,num_pairs
             if(((strip(species_names(i)).eq.pair_names(1,m)).and.(strip(species_names(j)).eq.pair_names(2,m)))&
                  .or.((strip(species_names(j)).eq.pair_names(1,m)).and.(strip(species_names(i)).eq.pair_names(2,m)))) then
                sigma(i,j)=pair_minsep(m)
             end if
          end do
       end do
    end do

    ! ** Set including the modified names

    do i=1,num_spec
       do j=1,num_spec
          do m=1,num_pairs
             if(((species_names(i).eq.pair_names(1,m)).and.(species_names(j).eq.pair_names(2,m)))&
                  .or.((species_names(j).eq.pair_names(1,m)).and.(species_names(i).eq.pair_names(2,m)))) then
                sigma(i,j)=pair_minsep(m)
             end if
          end do
       end do
    end do

    rcut = max(delta,range*sigma)

    do i=1,num_spec
       do j=1,num_spec
          if(sigma(i,j).lt.0.0_dp) then
             sigma(i,j) = abs(sigma(i,j))/(2.0_dp)**(1.0_dp/6.0_dp)
             rcut(i,j)  = 10.0_dp*sigma(i,j)
             epsil(i,j) = 10.0_dp
          end if
       end do
    end do

    rcut_sq=rcut**2

    ! ** Calculate shifts so that f=0 at |r1-r2|=rcut

    do i=1,num_spec
       do j=1,num_spec

          mod_r12=rcut(i,j)

          inv_mod_r12_n=(sigma(i,j)/mod_r12)**PP_n(i,j)*beta(i,j)

          inv_mod_r12_m=(sigma(i,j)/mod_r12)**PP_m(i,j)

          E_shift(i,j)=2.0_dp*epsil(i,j)*(inv_mod_r12_m-inv_mod_r12_n)

          if(forceshift) then
             f_shift(i,j)=2.0_dp*epsil(i,j)*(real(PP_m(i,j),dp)*inv_mod_r12_m/mod_r12-real(PP_n(i,j),dp)*inv_mod_r12_n/mod_r12)
          else
             f_shift=0.0_dp
          end if

       end do
    end do

    inipp=.true.

  end subroutine init_pp

  subroutine eval_pp(e,f,s)

    real(kind=dp),           intent(out) :: e
    real(kind=dp), optional, intent(out) :: f(3,num_ions*num_symm)
    real(kind=dp), optional, intent(out) :: s(3,3)

    integer       :: ion_i,ion_j,ion_k,ispeci,ispecj,ispeck,i
    integer       :: na,nna,nb,nnb,nc,nnc,ni,ns,nn,nnmx,nn0,n,m,mm
    real(kind=dp) :: radsph,R0(3),V0(3),dmin,ha,hb,hc,acb(3),acc(3),ccb(3),dec(num_ions*num_symm),ecoord,ethree,erep
    real(kind=dp) :: diag,de,rcsq,rcmod,mod_f12,work(3,num_ions*num_symm),fc(3,num_ions*num_symm,num_ions*num_symm)
    real(kind=dp) :: r12(3),r13(3),r23(3),rij(3),f12(3),rik(3),rjk(3),vec(3),vep(3)
    real(kind=dp) :: mod_r123,f_shift_123,inv_mod_r123_12,inv_mod_r123_6
    real(kind=dp) :: rijsq,rijmod,mod_f123,mod_r13,mod_r23,riksq,rikmod,rjksq,rjkmod
    real(kind=dp) :: normal(3),widthby2,origin(3),distpln,fplnshift,distsph,dtmp,ftmp

    real(kind=dp), allocatable :: lvec(:,:),ll(:),pos(:,:)

    logical       :: sameion,sameionik,sameionjk

    ! ** Translate atomic positions (units) to unit cell

    allocate(pos(3,num_ions*num_symm))

    if(.not.cluster) then

       pos(1:3,1:num_ions*num_symm)=ion_new_positions(1:3,1:num_ions*num_symm)

       ni=0
       do n=1,num_sets
          do ns=1,num_symm
             R0=0.0_dp
             do m=1,ninset(n)
                ni=ni+1
                R0=R0+ion_new_positions(:,ni)
             end do
             R0=R0/real(ninset(n),dp)
             V0=matmul(lattice_rec,R0)
             V0=V0-floor(V0)
             V0=matmul(lattice_car,V0)
             ni=ni-ninset(n)
             do m=1,ninset(n)
                ni=ni+1
                pos(:,ni)=ion_new_positions(:,ni)-R0(:)+V0(:)
             end do

          end do
       end do

    else

       ni=0
       do n=1,num_sets
          do ns=1,num_symm
             do m=1,ninset(n)
                ni=ni+1
                pos(:,ni)=ion_new_positions(:,ni)
             end do
          end do
       end do

    end if

    ! ** Initialise if needed

    if(.not.inipp) call init_pp()

    e=0.0_dp 
    if(present(f)) f=0.0_dp 
    if(present(s)) s=0.0_dp

    ! ** Determine the supercell required
    
    diag=sqrt(dot_product(sum(lattice_car,2),sum(lattice_car,2)))

    radsph=maxval(rcut)+maxval(ion_set_sphere)*2.0_dp+diag

    acb(1)=lattice_car(2,1)*lattice_car(3,2)-lattice_car(3,1)*lattice_car(2,2)
    acb(2)=lattice_car(3,1)*lattice_car(1,2)-lattice_car(1,1)*lattice_car(3,2)
    acb(3)=lattice_car(1,1)*lattice_car(2,2)-lattice_car(2,1)*lattice_car(1,2)

    hc=volume/sqrt(dot_product(acb,acb))

    acc(1)=lattice_car(2,1)*lattice_car(3,3)-lattice_car(3,1)*lattice_car(2,3)
    acc(2)=lattice_car(3,1)*lattice_car(1,3)-lattice_car(1,1)*lattice_car(3,3)
    acc(3)=lattice_car(1,1)*lattice_car(2,3)-lattice_car(2,1)*lattice_car(1,3)

    hb=volume/sqrt(dot_product(acc,acc))

    ccb(1)=lattice_car(2,3)*lattice_car(3,2)-lattice_car(3,3)*lattice_car(2,2)
    ccb(2)=lattice_car(3,3)*lattice_car(1,2)-lattice_car(1,3)*lattice_car(3,2)
    ccb(3)=lattice_car(1,3)*lattice_car(2,2)-lattice_car(2,3)*lattice_car(1,2)

    ha=volume/sqrt(dot_product(ccb,ccb))
   
    nna = int(radsph/ha)
    nnb = int(radsph/hb)
    nnc = int(radsph/hc)
    
    if(cluster) then
       nna=0 ; nnb=0 ; nnc=0
    end if

    ! ** Precalculate the lattice vectors and distances

    nnmx=(2*nna+1)*(2*nnb+1)*(2*nnc+1)
    allocate(lvec(3,nnmx),ll(nnmx))    

    nn=0
    do na=-nna,nna
       do nb=-nnb,nnb
          do nc=-nnc,nnc
             nn=nn+1
             lvec(1:3,nn)=na*lattice_car(1:3,1)+nb*lattice_car(1:3,2)+nc*lattice_car(1:3,3)
             ll(nn)=sqrt(lvec(1,nn)*lvec(1,nn)+lvec(2,nn)*lvec(2,nn)+lvec(3,nn)*lvec(3,nn))
             if(((na==0).and.(nb==0).and.(nc==0))) nn0=nn
          end do
       end do
    end do

    ! ** PP calculation of energy, forces and stress

    dec=0.0_dp
    fc=0.0_dp

    erep=0.0_dp
    ecoord=0.0_dp
    ethree=0.0_dp

    do ion_i=1,num_ions*num_symm
       if(ion_occ(ion_i).lt.1.0_dp-delta) cycle

       ispeci=ion_species(ion_i)


       do ion_j=1,num_ions*num_symm
          if(ion_occ(ion_j).lt.1.0_dp-delta) cycle

          ispecj=ion_species(ion_j)

          rij(:) = pos(:,ion_i)-pos(:,ion_j)

          rijsq = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
          rijmod=sqrt(rijsq)

          sameion=ion_ion_set(ion_i).eq.ion_ion_set(ion_j)

          rcsq=rcut_sq(ispeci,ispecj)
          rcmod=sqrt(rcsq)

          if(three.gt.0.0_dp) rcmod=max(rcmod,three)
          if(ion_coord(1,ion_i).ge.0) rcmod=max(rcmod,5.0_dp) ! ** HACK

          do nn=1,nnmx

             if((nn.eq.nn0).and.sameion) cycle

             if(abs(rijmod-ll(nn)).gt.rcmod) cycle

             r12(1:3) = rij(1:3)-lvec(1:3,nn)

             mod_r12_sq=r12(1)*r12(1)+r12(2)*r12(2)+r12(3)*r12(3)

             ! Check separation is not too big

             if (mod_r12_sq <= rcsq) then

                mod_r12=sqrt(mod_r12_sq)

                ! Calculate PP energy, forces and stress

                inv_mod_r12_n=(sigma(ispeci,ispecj)/mod_r12)**PP_n(ispeci,ispecj)*beta(ispeci,ispecj)
                inv_mod_r12_m=(sigma(ispeci,ispecj)/mod_r12)**PP_m(ispeci,ispecj)

                de=2.0_dp*epsil(ispeci,ispecj)*(inv_mod_r12_m-inv_mod_r12_n)
                de=de-E_shift(ispeci,ispecj)+(mod_r12-rcut(ispeci,ispecj))*f_shift(ispeci,ispecj)

                if(present(f)) then

                   mod_f12=2.0_dp*epsil(ispeci,ispecj)*&
                        (real(PP_m(ispeci,ispecj),dp)*inv_mod_r12_m/mod_r12-real(PP_n(ispeci,ispecj),dp)*inv_mod_r12_n/mod_r12)
                   mod_f12=mod_f12-f_shift(ispeci,ispecj) 
                   f12(:)=mod_f12*r12(:)/mod_r12

                end if

                ! Update energy in model

                e=e+de

                erep=erep+de
                
                ! Update forces in model

                if(present(f)) then
                   f(:,ion_i)=f(:,ion_i)+f12(:)
                end if

                ! Update stress

                if(present(s)) then
                   s(1,1)=s(1,1)+f12(1)*r12(1)
                   s(2,2)=s(2,2)+f12(2)*r12(2)
                   s(3,3)=s(3,3)+f12(3)*r12(3)

                   s(2,3)=s(2,3)+(f12(2)*r12(3)+f12(3)*r12(2))/2.0_dp
                   s(3,2)=s(3,2)+(f12(3)*r12(2)+f12(2)*r12(3))/2.0_dp

                   s(3,1)=s(3,1)+(f12(3)*r12(1)+f12(1)*r12(3))/2.0_dp
                   s(1,3)=s(1,3)+(f12(1)*r12(3)+f12(3)*r12(1))/2.0_dp

                   s(1,2)=s(1,2)+(f12(1)*r12(2)+f12(2)*r12(1))/2.0_dp
                   s(2,1)=s(2,1)+(f12(2)*r12(1)+f12(1)*r12(2))/2.0_dp
                end if

             end if

             ! ** Collect terms for coordination force

             if(ion_coord(1,ion_i).ge.0) then

                if(mod_r12_sq.lt.2.85_dp**2) then

                   mod_r12=sqrt(mod_r12_sq)

                   ! Calculate energy, forces and stress

                   dec(ion_i)=dec(ion_i)+softstep(mod_r12,sigma(ispeci,ispecj),sigma(ispeci,ispecj)*1.3_dp,mod_f12)

                   ! Update forces in model

                   if(present(f)) then
                      fc(:,ion_i,ion_j)=fc(:,ion_i,ion_j)+mod_f12*r12(:)/mod_r12
                      fc(:,ion_i,ion_i)=fc(:,ion_i,ion_i)+mod_f12*r12(:)/mod_r12
                   end if

                   ! Update stress - NOT DONE YET

                   if(present(s)) then
                      s=s+0.0_dp
                   end if

                end if

             end if

             ! ** Add the three body hard sphere force -- STRESS MISSING

             if (three.gt.0.0_dp) then

                do ion_k=1,num_ions*num_symm
                   if(ion_occ(ion_k).lt.1.0_dp-delta) cycle

                   ispeck=ion_species(ion_k)

                   sameionik=(ion_ion_set(ion_i).eq.ion_ion_set(ion_k))
                   sameionjk=(ion_ion_set(ion_j).eq.ion_ion_set(ion_k))

                   rij(1:3) = pos(:,ion_i)-pos(:,ion_j)
                   rik(1:3) = pos(:,ion_i)-pos(:,ion_k)
                   rjk(1:3) = pos(:,ion_j)-pos(:,ion_k)

                   rijsq = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
                   riksq = rik(1)*rik(1)+rik(2)*rik(2)+rik(3)*rik(3)
                   rjksq = rjk(1)*rjk(1)+rjk(2)*rjk(2)+rjk(3)*rjk(3)

                   rijmod=sqrt(rijsq)
                   rikmod=sqrt(riksq)
                   rjkmod=sqrt(rjksq)

                   do mm=1,nnmx

                      if((mm.eq.nn0).and.sameionik) cycle                      
                      if((mm.eq.nn).and.sameionjk) cycle

                      if(abs(rikmod-ll(mm)).gt.rcmod) cycle
                      if(abs(abs(rjkmod-ll(mm))-ll(nn)).gt.rcmod) cycle


                      r12(1:3) = pos(:,ion_i)-pos(:,ion_j)-lvec(1:3,nn)
                      r13(1:3) = pos(:,ion_i)-pos(:,ion_k)-lvec(1:3,mm)
                      r23(1:3) = pos(:,ion_j)-pos(:,ion_k)-lvec(1:3,mm)+lvec(1:3,nn)

                      mod_r12=sqrt(r12(1)*r12(1)+r12(2)*r12(2)+r12(3)*r12(3))
                      mod_r13=sqrt(r13(1)*r13(1)+r13(2)*r13(2)+r13(3)*r13(3))
                      mod_r23=sqrt(r23(1)*r23(1)+r23(2)*r23(2)+r23(3)*r23(3))

                      mod_r123=mod_r12+mod_r13+mod_r23

                      if(mod_r123.lt.three) then

                         de=hardrcore(mod_r123,three,mod_f123)                         

                         e=e+de

                         ethree=ethree+de

                         if(present(f)) then

                            f(:,ion_i)=f(:,ion_i)-mod_f123*(r12(:)/mod_r12+2*r13(:)/mod_r13)

                         end if

                      end if

                   end do

                end do


             end if



          end do !nn


       end do  !ion_j

       ! ** Add in the plane force

       if(width.ge.0.0_dp) then

          widthby2=width/2.0_dp          

          normal=(/0.0_dp,0.0_dp,1.0_dp/)

          origin=(/0.0_dp,0.0_dp,shift/)

          ! * Shortest distance to plane

          distpln=huge(1.0_dp)

          do nn=1,nnmx
             if(abs(dot_product(pos(:,ion_i)-lvec(:,nn)-origin,normal)).lt.abs(distpln)) &
                  distpln=dot_product(pos(:,ion_i)-lvec(:,nn)-origin,normal)
          end do

          if(abs(distpln).ge.widthby2) then
             f(3,ion_i)=f(3,ion_i)-distpln/abs(distpln)*(abs(distpln)-widthby2)
          end if

       end if

       ! ** Add in sphere force

       if(sphere) then

          ! * Shortest distance to sphere centre

          if(cluster) then
             origin=0.0_dp
          else
             origin=(lattice_car(1:3,1)+lattice_car(1:3,2)+lattice_car(1:3,3))/2.0_dp
          end if

          V0=pos(:,ion_i)-origin(:)

          distsph=huge(1.0_dp)
          do nn=1,nnmx
             
             vec=V0(:)-lvec(:,nn)
             
             if(ellipsoid) vec(:)=vec(:)/ellipsvec(:)

             dtmp=sqrt(dot_product(vec,vec))
             if(dtmp.lt.distsph) then
                distsph=dtmp
                vep=vec
                R0(:)=origin(:)+lvec(:,nn)
             end if
          end do

          if(radius.gt.0.0_dp) then

             if (distsph.gt.radius) then
                ftmp=-((12.0_dp*(distsph/radius)**12-6.0_dp*(distsph/radius)**6)-6.0_dp)/distsph**2
                f(:,ion_i)=f(:,ion_i)+ftmp*vep
             end if

          else

             ftmp=-4.0_dp*abs(radius)*distsph**2
             f(:,ion_i)=f(:,ion_i)+ftmp*vep

          end if


       end if

       if(cylinder) then

          ! * Shortest distance to line

          if(cluster) then
             origin=0.0_dp
          else
             origin=(lattice_car(1:3,1)+lattice_car(1:3,2)+lattice_car(1:3,3))/2.0_dp
          end if

          distsph=huge(1.0_dp)
          do nn=1,nnmx
             dtmp=sqrt(dot_product(pos(1:2,ion_i)-origin(1:2)-lvec(1:2,nn),pos(1:2,ion_i)-origin(1:2)-lvec(1:2,nn)))
             if(dtmp.lt.distsph) then
                distsph=dtmp
                R0(:)=origin(:)+lvec(:,nn)
             end if
          end do

          if(radius.gt.0.0_dp) then

             if (distsph.gt.radius) then
                ftmp=-((12.0_dp*(distsph/radius)**12-6.0_dp*(distsph/radius)**6)-6.0_dp)/distsph**2
                f(:,ion_i)=f(:,ion_i)+ftmp*(pos(:,ion_i)-origin(:))
             end if

          else

             ftmp=-4.0_dp*abs(radius)*distsph**2
             f(1:2,ion_i)=f(1:2,ion_i)+ftmp*(pos(1:2,ion_i)-origin(1:2))


          end if


       end if

    end do !ion_i

    ! ** Combine the coordination force and energy

    do ion_i=1,num_ions*num_symm
       if(ion_occ(ion_i).lt.1.0_dp-delta) cycle

       if(ion_coord(1,ion_i).ge.0) then

          ! Update force

          if(present(f)) then
             do ion_j=1,num_ions*num_symm
                if(ion_occ(ion_i).lt.1.0_dp-delta) cycle
                f(:,ion_i)=f(:,ion_i)+1.0_dp*(ion_coord(1,ion_j)-dec(ion_j))*fc(:,ion_i,ion_j)
             end do
          end if

          ! Update energy

          e=e+(ion_coord(1,ion_i)-dec(ion_i))**2
          ecoord=ecoord+(ion_coord(1,ion_i)-dec(ion_i))**2

       end if

       ! * No force on the fixed atoms
       
       if(ion_bfix(ion_i)) f(:,ion_i)=0.0_dp

    end do

    deallocate(lvec)

    ! ** Symmetrise forces

    if((num_symm.gt.1).and.present(f)) then

       work=f
       f=0.0_dp
       do ni=1,num_ions*num_symm
          if(ion_equiv(1,ni).eq.0) cycle
          do ns=1,num_symm
             f(:,ni) = f(:,ni) + matmul((Sym(1:3,1:3,ns)),work(:,ion_equiv(ns,ni)))
          end do
          f(:,ni) = f(:,ni)/real(num_symm,dp)
       end do

    end if

!!$    write (50,*) erep,ecoord,ethree
    
  end subroutine eval_pp

  function softstep(x,r,rc,ssp)

    real(kind=dp), intent(in)  :: x
    real(kind=dp), intent(in)  :: r
    real(kind=dp), intent(in)  :: rc
    real(kind=dp), intent(out) :: ssp

    real(kind=dp) :: softstep

    real(kind=dp), parameter :: a=20.0_dp

    real(kind=dp) :: xx

    ! ** Smoothstep function

!!$    if(x.lt.r) then
!!$       softstep=1.0_dp
!!$       ssp=0.0_dp
!!$    else if(x.gt.rc) then
!!$       softstep=0.0_dp
!!$       ssp=0.0_dp
!!$    else
!!$       xx=(r-x)/(r-rc)
!!$       softstep=1-(3*xx**2-2*xx**3)
!!$       ssp=6/(r-rc)*(xx-xx**2)
!!$    end if

    ! ** Smoother step function

    if(x.lt.r) then
       softstep=1.0_dp
       ssp=0.0_dp
    else if(x.gt.rc) then
       softstep=0.0_dp
       ssp=0.0_dp
    else
       xx=(r-x)/(r-rc)
       softstep=1-(6*xx**5-15*xx**4+10*xx**3)
       ssp=30/(r-rc)*(xx**4-2*xx**3+xx**2)
    end if


!!$    ! ** Smoothest step function
!!$
!!$    if(x.lt.r) then
!!$       softstep=1.0_dp
!!$       ssp=0.0_dp
!!$    else if(x.gt.rc) then
!!$       softstep=0.0_dp
!!$       ssp=0.0_dp
!!$    else
!!$       xx=(r-x)/(r-rc)
!!$       softstep=1-(-20*xx**7+70*xx**6-84*xx**5+35*xx**4)
!!$       ssp=140/(r-rc)*(-xx**6+3*xx**5-3*xx**4+xx**3)
!!$    end if


    ! ** Just a line

    if(x.lt.r) then
       softstep=1.0_dp
       ssp=0.0_dp
    else if(x.gt.rc) then
       softstep=0.0_dp
       ssp=0.0_dp
    else
       xx=(r-x)/(r-rc)
       softstep=1-xx
       ssp=1/(r-rc)
    end if
    
  end function softstep

  function softrcore(x,r,srcp)

    real(kind=dp), intent(in)  :: x
    real(kind=dp), intent(in)  :: r
    real(kind=dp), intent(out) :: srcp

    real(kind=dp) :: softrcore

    real(kind=dp) :: xx

    ! ** Smooth repulsive core function

    if(x.gt.r) then
       softrcore=0.0_dp
       srcp=0.0_dp
    else
       xx=x/r
       softrcore=(xx-1)**2
       srcp=2/r*(xx-1)
    end if

  end function softrcore

  function hardrcore(x,r,hrcp)

    real(kind=dp), intent(in)  :: x
    real(kind=dp), intent(in)  :: r
    real(kind=dp), intent(out) :: hrcp

    real(kind=dp) :: hardrcore

    real(kind=dp) :: xx

    ! ** Smooth repulsive core function

    if(x.gt.r) then
       hardrcore=0.0_dp
       hrcp=0.0_dp
    else
       xx=x/r
       hardrcore=1/xx**12-1/xx**6+6*(xx-1)
       hrcp=-12/xx**13+6/xx**7+6
       hrcp=hrcp/r
    end if

  end function hardrcore


end module pp
