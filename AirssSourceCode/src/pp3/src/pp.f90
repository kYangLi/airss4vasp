!==================================================================================!
!                                       pp                                         !
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

  implicit none

  private

  public :: read_pp
  public :: write_pp
  public :: init_pp
  public :: eval_pp

  ! ** Data defining the potential
  
  integer,                allocatable, dimension(:,:) :: PP_m, PP_n
  real(kind=dp)                                       :: range
  real(kind=dp), public,  allocatable, dimension(:,:) :: epsil,sigma,beta
  real(kind=dp),          allocatable, dimension(:,:) :: rcut, rcut_sq, E_shift, f_shift, mat

  character(len=6), allocatable, dimension(:) :: spec
  
  logical, public :: forceshift=.true.

  !---------------------------------------------------!

  ! ** Private variables

  integer, parameter :: unit_pp=22
  integer            :: nspec
  character(len=80)  :: ppfile
  
  real(kind=dp) :: mod_r12, mod_r12_sq, inv_mod_r12_n, inv_mod_r12_m, nelec

contains

  subroutine read_pp

    integer           :: ns,n,m,i,j
    character(len=80) :: ctemp

    nspec=num_spec

    ! * Allocate pp arrays

    allocate(PP_m(nspec,nspec))
    allocate(PP_n(nspec,nspec))
    allocate(epsil(nspec,nspec))
    allocate(sigma(nspec,nspec))
    allocate(beta(nspec,nspec))
    allocate(E_shift(nspec,nspec))
    allocate(F_shift(nspec,nspec))
    allocate(rcut(nspec,nspec))
    allocate(rcut_sq(nspec,nspec))

    ! ** Open the pp file

    write(ppfile,'(a)') trim(seedname)//".pp"

    open(unit=unit_pp,file=ppfile,status="old",err=999)

    ! ** Check the species count, and read exponents

    read(unit_pp,*,err=999) ns, PP_m(1,1), PP_n(1,1), range
    allocate(spec(ns),mat(ns,ns))
    read(unit_pp,*,err=999) spec


    PP_m = PP_m(1,1)
    PP_n = PP_n(1,1)

    if(ns.lt.nspec) stop 'PP is for the wrong number of species in read_pp'

    ! ** Epsil

    read(unit_pp,*,err=999) ctemp

    do n=1,ns
       read(unit_pp,*,err=999) (mat(n,m),m=n,ns)
       do m=n+1,ns
          mat(m,n)=mat(n,m)
       end do
    end do

    do n=1,nspec
       do i=1,ns
          if(spec(i).eq.ion_names_spec(n)) exit
       end do
       do m=1,nspec
          do j=1,ns
             if(spec(j).eq.ion_names_spec(m)) exit
          end do
          epsil(n,m)=mat(i,j)
       end do
    end do

    ! ** Sigma

    read(unit_pp,*,err=999) ctemp

    do n=1,ns
       read(unit_pp,*,err=999) (mat(n,m),m=n,ns)
       do m=n+1,ns
          mat(m,n)=mat(n,m)
       end do
    end do

    do n=1,nspec
       do i=1,ns
          if(spec(i).eq.ion_names_spec(n)) exit
       end do
       do m=1,nspec
          do j=1,ns
             if(spec(j).eq.ion_names_spec(m)) exit
          end do
          sigma(n,m)=mat(i,j)
       end do
    end do

    ! ** Beta

    read(unit_pp,'(a)',err=999,end=998) ctemp

998 continue

    if(index(ctemp,"Beta").gt.0) then

       do n=1,ns
          read(unit_pp,*,err=999) (mat(n,m),m=n,ns)
          do m=n+1,ns
             mat(m,n)=mat(n,m)
          end do
       end do

       do n=1,nspec
          do i=1,ns
             if(spec(i).eq.ion_names_spec(n)) exit
          end do
          do m=1,nspec
             do j=1,ns
                if(spec(j).eq.ion_names_spec(m)) exit
             end do
             beta(n,m)=mat(i,j)
          end do
       end do
       
    else
       beta=1.0_dp
    end if

    ! * Set to scaled means if negative epsil or sigma provided

    do n=2,nspec
       if(epsil(n,n).lt.0.0_dp) epsil(n,n)=epsil(1,1)*abs(epsil(n,n))
       if(sigma(n,n).lt.0.0_dp) sigma(n,n)=sigma(1,1)*abs(sigma(n,n))
    end do

    do n=1,nspec
       do m=1,nspec
          if(n.ne.m) then
             if(epsil(m,n).lt.0.0_dp) epsil(m,n)=abs(epsil(m,n))*sqrt(abs(epsil(n,n)*epsil(m,m)))
             if(sigma(m,n).lt.0.0_dp) sigma(m,n)=abs(sigma(m,n))*(abs(sigma(n,n))+abs(sigma(m,m)))/2.0_dp
          end if
       end do
    end do

    close(unit_pp)
    
    return

999 stop 'There is a problem reading the pp information. Stopping.'

  end subroutine read_pp

  subroutine write_pp

    integer :: m,n

    ! ** Open the new pp file

    write(ppfile,'(a)') trim(seedname)//".new.pp"

    open(unit=unit_pp,file=ppfile,status="unknown",err=998)

    write (unit_pp,'(3i3,f10.5)')  nspec, PP_m(1,1), PP_n(1,1), range
    write (unit_pp,'(a)') "# Epsilon"
    do n=1,nspec
       write(unit_pp,*) (epsil(n,m),m=n,nspec)
    end do
    write (unit_pp,'(a)') "# Sigma"
    do n=1,nspec
       write(unit_pp,*) (sigma(n,m),m=n,nspec)
    end do
    write (unit_pp,'(a)') "# Beta"
    do n=1,nspec
       write(unit_pp,'(10i3)') (int(beta(n,m)),m=n,nspec)
    end do

    close(unit_pp)

    return

998 stop 'There is a problem opening the new pp file. Stopping.'

  end subroutine write_pp

  subroutine init_pp

    integer       :: i,j,n
    real(kind=dp) :: mod_r12, inv_mod_r12_n, inv_mod_r12_m

    rcut = range*sigma
    
    rcut_sq=rcut**2

    ! ** Calculate shifts so that f=0 at |r1-r2|=rcut

    do i=1,nspec
       do j=1,nspec

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

  end subroutine init_pp

  subroutine eval_pp(e,f,s)

    real(kind=dp),           intent(out) :: e
    real(kind=dp), optional, intent(out) :: f(3,num_ions)
    real(kind=dp), optional, intent(out) :: s(3,3)

    integer       :: ion_i,ion_j,ispeci,ispecj
    integer       :: na,nna,nb,nnb,nc,nnc,ni,ns,nn,nnmx,nn0
    
    real(kind=dp) :: r12(3),rij(3),f12(3)
    real(kind=dp) :: de,rcsq,rcmod,rijsq,rijmod,mod_f12,work(3,num_ions*num_symm),swork(3,3)
    real(kind=dp), allocatable :: lvec(:,:),ll(:)

    logical       :: sameion
    
    e=0.0_dp 
    if(present(f)) f=0.0_dp 
    if(present(s)) s=0.0_dp

    ! Determine the supercell required

    nna = (int(2*maxval(rcut)/lattice_abc(1))+1)/2+1
    nnb = (int(2*maxval(rcut)/lattice_abc(2))+1)/2+1
    nnc = (int(2*maxval(rcut)/lattice_abc(3))+1)/2+1
    
    if(cluster) then
       nna=0 ; nnb=0 ; nnc=0
    end if

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

    do ion_i=1,num_ions

       ispeci=ion_species(ion_i)

       do ion_j=ion_i,num_ions

          ispecj=ion_species(ion_j)

          rij(:) = ion_positions(:,ion_i)-ion_positions(:,ion_j)

          rijsq = rij(1)*rij(1)+rij(2)*rij(2)+rij(3)*rij(3)
          rijmod=sqrt(rijsq)

          sameion=ion_i.eq.ion_j

          rcsq=rcut_sq(ispeci,ispecj)
          rcmod=sqrt(rcsq)

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
                        (real(PP_m(ispeci,ispecj),dp)*inv_mod_r12_m/mod_r12-&
                        real(PP_n(ispeci,ispecj),dp)*inv_mod_r12_n/mod_r12)
                   mod_f12=mod_f12-f_shift(ispeci,ispecj) 
                   f12(:)=mod_f12*r12(:)/mod_r12

                end if

                if(ion_j.ne.ion_i) then
                   if(present(f)) f12=f12*2.0_dp
                   de=de*2.0_dp
                end if

                ! Update energy in model

                e=e+de

                ! Update forces in model

                if(present(f)) then
                   f(:,ion_i)=f(:,ion_i)+f12(:)
                   if(ion_j.ne.ion_i) then                         
                      f(:,ion_j)=f(:,ion_j)-f12(:)
                   end if
                end if

                ! Update stress

                if(present(s)) then
                   s(1,1)=s(1,1)+f12(1)*r12(1)
                   s(2,2)=s(2,2)+f12(2)*r12(2)
                   s(3,3)=s(3,3)+f12(3)*r12(3)
                   s(2,3)=s(2,3)+(f12(2)*r12(3)+f12(3)*r12(2))/2.0_dp
                   s(3,1)=s(3,1)+(f12(3)*r12(1)+f12(1)*r12(3))/2.0_dp
                   s(1,2)=s(1,2)+(f12(1)*r12(2)+f12(2)*r12(1))/2.0_dp
                end if

             end if


          end do


       end do  !ion_j
    end do !ion_i

    deallocate(lvec)
    
    if(present(s)) then
       s(3,2) = s(2,3)
       s(1,3) = s(3,1)
       s(2,1) = s(1,2)
       s=s/volume       ! * Convert stress to absolute form
    end if

    ! * Balance forces

    if(present(f)) then
       f12(1) = sum(f(1,:))/real(num_ions,dp)
       f12(2) = sum(f(2,:))/real(num_ions,dp)
       f12(3) = sum(f(3,:))/real(num_ions,dp)
       
       f(1,:) = f(1,:)-f12(1)
       f(2,:) = f(2,:)-f12(2)
       f(3,:) = f(3,:)-f12(3)
    end if

    if(num_symm.gt.1) then
       
       if(present(f)) then
          
          ! * Symmetrise forces
          
          work=f
          f=0.0_dp
          
          do ni=1,num_ions
             do ns=1,num_symm
                f(:,ni) = f(:,ni) +  matmul(symm_ops(1:3,1:3,ns),work(:,ion_equiv(ns,ni))) 
             end do
             f(:,ni) = f(:,ni)/real(num_symm,dp)
          end do
          
       end if
       
       if(present(s)) then

          ! * Symmetrise stresses

          swork = s/real(num_symm,dp)
          s=0.0_dp
          
          do ns=1,num_symm
             s = s + matmul(symm_ops(1:3,1:3,ns),matmul(swork,transpose(symm_ops(1:3,1:3,ns))))
          end do

       end if

    end if

    ! ** Constrain ions

    do ni=1,num_ions
       if(ion_cons(ni)) f(:,ni)=0.0_dp
    end do
    

  end subroutine eval_pp

end module pp
