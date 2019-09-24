! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
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
  use symmetry
  use rng

  implicit none

  private

  public :: rash
  public :: opt_tpsd

  integer, public :: maxsteps=1000

  !---------------------------------------------------!

  ! ** Private data

  integer :: nvec

  real(kind=dp) :: mu              ! ** Parameter for augmented Lagrangian

  real(kind=dp) :: thresh=1e-4_dp ! ** Convergence parameter
  
  real(kind=dp) :: alpha           ! ** The Lagrange multiplier for volume conservation
  real(kind=dp) :: volume0         ! ** The initial volume
  real(kind=dp) :: orig_car(3,3)   ! ** The original lattice vectors
  real(kind=dp) :: orig_rec(3,3)   ! ** The original reciprocal lattice vectors

  
  real(kind=dp), allocatable, dimension(:,:,:) :: unit_positions

  logical :: keep_going
  
contains

  subroutine rash(converged,over)
    
    logical,        intent(out) :: converged
    real(kind=dp), intent(out) :: over
    
    integer, parameter :: count_max=10
    
    integer :: nn,nc,n,m,ni,npow
    real(kind=dp), allocatable :: pos0(:,:),set_positions(:,:)
    real(kind=dp) :: o0,rn(3),minamp,an(3),mat(3,3),centre(3),shift(3),mrot(3,3),lat0(3,3)
    
    nn=max(num_ions*num_symm*num_cells,num_ions*num_form*num_symm)
    allocate(pos0(3,nn),set_positions(3,nn))

    o0=huge(1.0_dp)
    pos0=ion_new_positions
    lat0=lattice_car
    nc=0
    do
       nc=nc+1
       
       if(nc.gt.count_max) exit       
       
       ! ** Relax
       
       call opt_tpsd(converged,over)
       
       ! ** Store
       
       if((over.le.max(o0-thresh*10.0_dp,0.0_dp)).and.converged) then
          o0=over
          pos0=ion_new_positions
          lat0=lattice_car
          nc=0
       end if

       write (stderr,'(i4,2f10.5)') nc,over/real(num_ions_total),o0/real(num_ions_total)
       
       ! ** Check if we have already reached the required overlap

       if(converged.and.((over/real(num_ions_total).lt.overlap+delta))) exit

       ! ** Shake

       ion_new_positions=pos0
       lattice_car=lat0

       !call compact_cell()
       call update_cell()

       ni=0
       do n=1,num_sets

          ! ** Random shift
          
          rn  = random_triple()
          an  = random_triple()
          mat = rotation_matrix(an)
          
          npow=2 ! * If npow=2 the distribution is uniform in 3D space
          minamp=0.0_dp
          rn(1) = minamp**(npow+1)+rn(1)*(rash_posamp**(npow+1)-minamp**(npow+1))
          rn(1) = rn(1)**(1.0_dp/(npow+1))
          rn(2) = 0.0_dp ; rn(3) = 0.0_dp
          shift = matmul(mat,rn)

          ! ** Random rotation

          if (ninset(n).gt.1) then

             an = random_triple()

             an = an*(rash_angamp/360.0_dp)
             
             mrot = rotation_matrix(an)

             centre=0.0_dp
             
             do m=1,ninset(n)

                set_positions(:,m)=ion_new_positions(:,ni+m)
                
                centre(:)=centre(:)+set_positions(:,m)

             end do

             centre=centre/real(ninset(n),dp)

             do m=1,ninset(n)

                set_positions(:,m)=set_positions(:,m)-centre(:)
                
             end do
             
          end if
             
             
          do m=1,ninset(n)
             
             ni=ni+1
             if (ninset(n).gt.1) then
                
                ion_new_positions(:,ni) = centre(:) + shift(:) + matmul(mrot,set_positions(:,m))
                
             else
                
                ion_new_positions(:,ni)=ion_new_positions(:,ni)+shift(:)
                
             end if
             
             
          end do

!!$          nii=ni-ninset(n)
!!$          
!!$          do ns=2,num_symm
!!$
!!$             do m=1,ninset(n)
!!$             
!!$                ni=ni+1
!!$
!!$                V = 0.0_dp
!!$                do n1=1,3
!!$                   do n2=1,3
!!$                      V(n1) = V(n1)+Sym(n1,n2,ns)*ion_new_positions(n2,nii+m)
!!$                   end do
!!$                end do
!!$                
!!$                ion_new_positions(:,ni) = V(:) + symm_ops(1,4,ns)*lattice_car(:,1) + &
!!$                     symm_ops(2,4,ns)*lattice_car(:,2) + symm_ops(3,4,ns)*lattice_car(:,3) 
!!$                
!!$             end do
!!$             
!!$          end do
          
       end do
          
       
    end do

    ! ** Return to best
    
    ion_new_positions=pos0

    lattice_car=lat0

    !call compact_cell()
    call update_cell()
    
    call opt_tpsd(converged,over)
    
    deallocate(pos0,set_positions)

    
  end subroutine rash
  
  subroutine opt_tpsd(converged,over)

    !-------------------------------------------------------------!
    ! Two-point Step Size Gradient Methods - Barzilai and Borwein !
    ! IMA Journal of Numerical Analysis (1988) 8, 141-148         !
    !-------------------------------------------------------------!

    logical,                 intent(out) :: converged
    real(kind=dp), optional, intent(out) :: over

    real(kind=dp) :: step,o,o0,xg,gg,gdg

    real(kind=dp), allocatable, dimension(:) :: x,x0,g,g0,s,y
    real(kind=dp), allocatable, dimension(:) :: xxx,xxx0,ggg,ggg0,gggd

    integer :: n,nsteps

    logical :: first

    keep_going=.true.
    
    ! ** Force initialisation of pair potentials

    inipp=.false.

    ! ** Ensure everything is consistent for this cell

    !call compact_cell() 
    call update_cell()

    ! ** Store the original set of lattice vectors

    orig_car=lattice_car
    orig_rec=lattice_rec

    ! ** Store the volume, which will be maintained during optimisation

    volume0=volume

    ! ** Set up the vectors and gradients used in optimisation

    nvec=0
    do n=1,num_sets
       if(ninset(n).gt.1) nvec=nvec+1
    end do
    nvec=(num_sets*3+nvec*4)*num_symm

    if(celladapt) nvec=nvec+10 ! * 9 for cell, 1 for lagrange multiplier (alpha)

    allocate(x(nvec),x0(nvec),g(nvec),g0(nvec))
    allocate(xxx(nvec),xxx0(nvec),ggg(nvec),ggg0(nvec),gggd(nvec))
    allocate(s(nvec),y(nvec))

    call pos2vec(x0)

    if(celladapt) x0(nvec)=5.0_dp

!!$    call test_derivatives(x0,nvec)
!!$    stop
    
    mu=0.1_dp
    first  = .true.
    o0=huge(1.0_dp)
    gdg=huge(1.0_dp)
    converged = .true.
    step=0.0_dp
    g0=0.0_dp

    nsteps = 0
    do while((gdg.gt.thresh).and.keep_going)
       nsteps=nsteps+1
       if(nsteps.gt.maxsteps) exit

       x=x0+step*g0
       
       call force2vec(x,g,o)
       
       gdg = sqrt(dot_product(g,g)/real(nvec,dp))

!!$       write (40,*) nsteps,log10(max(gdg,epsilon(1.0_dp))),o
!!$       flush(40)

       if(.not.first) then

          s=x-x0
          y=g-g0

          xg=vecdotvec(s,y)
          gg=vecdotvec(y,y)

          if(gg.gt.0.0_dp) then
             step=abs(xg/gg)
          else
             step=0.0_dp
          end if
       
          step=max(step,epsilon(1.0_dp))

       else

          step=delta          
          first=.false.

       end if

       g0 = g
       x0 = x
       o0 = o

    end do

!!$    write (40,*)

    converged=(gdg.lt.thresh)

    if(present(over)) over=o-alpha*(volume-volume0)-mu*(volume-volume0)**2/2.0_dp

    
  end subroutine opt_tpsd

  subroutine pos2vec(vec)

    real(kind=dp), dimension(:), intent(out) :: vec
    
    real(kind=dp) :: R0(3)
    integer       :: ni,nv,n,i,m,ns,nu

    ! ** Allocate internal arrays
    
    if(allocated(unit_positions)) deallocate(unit_positions)
    allocate(unit_positions(3,maxval(ninset),num_sets*num_symm))
    
    vec=0.0_dp
    unit_positions=0.0_dp

    ni=0
    nv=0
    nu=0
    do n=1,num_sets
       do ns=1,num_symm
          
          nu=nu+1

          R0=0.0_dp
          do m=1,ninset(n)
             ni=ni+1
             R0=R0+ion_new_positions(:,ni)
          end do
          R0=R0/real(ninset(n),dp)
          
          ! * Orientation quaternion of unit

          if(ninset(n).gt.1) then

             do i=1,1
                nv=nv+1
                vec(nv)=1.0_dp
             end do

             do i=2,4
                nv=nv+1
                vec(nv)=0.0_dp
             end do

          end if

          ! * Centre of unit, convert to fractional

          vec(nv+1:nv+3)=matmul(lattice_rec,R0) 
          nv=nv+3

          ! * Store the units, in absolute coordinates

          ni=ni-ninset(n)
          do m=1,ninset(n)
             ni=ni+1
             unit_positions(:,m,nu)=ion_new_positions(:,ni)-R0(:)
          end do

       end do
    end do

    if(celladapt) then

       !call compact_cell()  
       call update_cell()
       
       vec(nv+1)=lattice_car(1,1)
       vec(nv+2)=lattice_car(2,1)
       vec(nv+3)=lattice_car(3,1)
       vec(nv+4)=lattice_car(1,2)
       vec(nv+5)=lattice_car(2,2)
       vec(nv+6)=lattice_car(3,2)
       vec(nv+7)=lattice_car(1,3)
       vec(nv+8)=lattice_car(2,3)
       vec(nv+9)=lattice_car(3,3)
       
       vec(nvec)=alpha
       
    end if
    
  end subroutine pos2vec

  subroutine vec2pos(vec)

    real(kind=dp), dimension(:),   intent(inout)  :: vec

    real(kind=dp) :: v(3),rot(3,3),symmerr,vol,vol0,lc0(3,3)
    integer :: ni,nv,n,m,ns,nu

    ! ** New lattice vectors

    if(celladapt) then

       nv=nvec-10

       lc0=lattice_car
       
       vol0=volume_cell(lc0)
       
       lattice_car(1,1)=vec(nv+1)
       lattice_car(2,1)=vec(nv+2)
       lattice_car(3,1)=vec(nv+3)
       lattice_car(1,2)=vec(nv+4)
       lattice_car(2,2)=vec(nv+5)
       lattice_car(3,2)=vec(nv+6)
       lattice_car(1,3)=vec(nv+7)
       lattice_car(2,3)=vec(nv+8)
       lattice_car(3,3)=vec(nv+9)

       vol=volume_cell(lattice_car)
       
       if(abs((vol-vol0)/vol0).gt.(1.0_dp)) then
          write (stderr,'(a)') 'Volume change too large'
          keep_going=.false.
       end if
       
       alpha=vec(nvec)

       cell_consistent=.false.

       !call symm_cell()

       call update_cell()

    end if

    ! ** New ionic positions

    ion_new_positions=0.0_dp

    ni=0 ; nv=0 ; nu=0
    do n=1,num_sets
       do ns=1,num_symm
          nu=nu+1

          if(ninset(n).gt.1) then
             rot(:,:) = q2rot(vec(nv+1:nv+4))
             v(:)=matmul(lattice_car,vec(nv+5:nv+7))
             nv=nv+7
          else

             v(:)=matmul(lattice_car,vec(nv+1:nv+3))
             nv=nv+3

          end if

          if(ninset(n).gt.1) then

             do m=1,ninset(n)
                ni=ni+1
                ion_new_positions(1:3,ni)= v(:) + matmul(rot,unit_positions(1:3,m,nu))
             end do

          else

             ni=ni+1
             ion_new_positions(1:3,ni) = v(:)

          end if

       end do
    end do

    ! ** Enforce symmetry

    if(num_symm.gt.1) then

       ion_temp_positions=0.0_dp
       symmerr=0.0_dp
       do ni=1,num_ions*num_symm
          if(ion_occ(ni).lt.1.0_dp-delta) cycle
          do ns=1,num_symm

             v=ion_new_positions(1:3,ni)-matmul(Sym(1:3,1:3,ns),ion_new_positions(1:3,ion_equiv(ns,ni)))-&
                  symm_ops(1,4,ns)*lattice_car(1:3,1) - &
                  symm_ops(2,4,ns)*lattice_car(1:3,2) - &
                  symm_ops(3,4,ns)*lattice_car(1:3,3)

             v = matmul(lattice_rec,v)
             v = v - nint(v)
             v = matmul(lattice_car,v)

             ion_temp_positions(1:3,ni)=ion_temp_positions(1:3,ni)+v/real(num_symm,dp)

          end do
          symmerr=symmerr+dot_product(ion_temp_positions(1:3,ni),ion_temp_positions(1:3,ni))
       end do

       symmerr=sqrt(symmerr)

       ion_new_positions=ion_new_positions-ion_temp_positions

       if(symmerr.gt.1e-4_dp) then
          write (stderr,'(a,f20.15)') 'Lost symmetry',symmerr
          keep_going=.false.
       end if

       ! ** Make the vector consistent if deviation significant

       if(symmerr.gt.1e3_dp*epsilon(symmerr)) call pos2vec(vec)

    end if

  end subroutine vec2pos

  subroutine force2vec(vec,fvec,over)

    real(kind=dp), dimension(:),   intent(inout)  :: vec
    real(kind=dp), dimension(:),   intent(out) :: fvec
    real(kind=dp),                 intent(out) :: over
    
    real(kind=dp) :: F0(3),T(3),Tq(4),R(3)
    real(kind=dp) :: f(3,num_ions*num_symm),s(3,3),sv(3,3),ss(3,3)

    integer :: ni,nv,n,m,ns,nu,i

    ! ** Convert the vector to ionic positions and current cell
    
    call vec2pos(vec)
    
    ! ** Calculate the forces and stresses
    
    call eval_pp(over,f,s)

    if(celladapt) over=over+alpha*(volume-volume0)+mu*(volume-volume0)**2/2.0_dp

    ! ** Construct the vector
    
    fvec=0.0_dp

    ni=0
    nv=0
    nu=0

    do n=1,num_sets
       do ns=1,num_symm
          nu=nu+1
          
          T=0.0_dp
          F0=0.0_dp
          do m=1,ninset(n)
             ni=ni+1
             F0=F0+f(:,ni)
             if(ninset(n).gt.1) then
               
                R(:)=matmul(q2rot(vec(nv+1:nv+4)),unit_positions(:,m,nu))
            
                T(1)=T(1)+R(2)*f(3,ni)-R(3)*f(2,ni)
                T(2)=T(2)+R(3)*f(1,ni)-R(1)*f(3,ni)
                T(3)=T(3)+R(1)*f(2,ni)-R(2)*f(1,ni)
                
             end if
          end do
          
          F0=F0/real(ninset(n),dp)

          ! * Torque quaternion

          if(ninset(n).gt.1) then
             
             T=T/real(ninset(n),dp)
             
             Tq(1)=0.0_dp
             Tq(2:4)=T(1:3)

             fvec(nv+1:nv+4)=0.5_dp*qmult(Tq,vec(nv+1:nv+4))
             nv=nv+4

          end if

          ! * Force, convert to fractional

          !fvec(nv+1:nv+3)=matmul(transpose(lattice_car),F0)*2.0_dp
          fvec(nv+1:nv+3)=matmul(lattice_rec,F0)*2.0_dp ! ** This factor of 2 lives somewhere else
          nv=nv+3
            
       end do
    end do
   
    if(celladapt) then
       
       ss=transpose(matmul(lattice_rec,s))

       sv(1,1)=(lattice_car(2,2)*lattice_car(3,3)-lattice_car(2,3)*lattice_car(3,2))
       sv(1,2)=(lattice_car(2,3)*lattice_car(3,1)-lattice_car(2,1)*lattice_car(3,3))
       sv(1,3)=(lattice_car(2,1)*lattice_car(3,2)-lattice_car(3,1)*lattice_car(2,2))
       sv(2,1)=(lattice_car(1,3)*lattice_car(3,2)-lattice_car(1,2)*lattice_car(3,3))
       sv(2,2)=(lattice_car(1,1)*lattice_car(3,3)-lattice_car(1,3)*lattice_car(3,1))
       sv(2,3)=(lattice_car(1,2)*lattice_car(3,1)-lattice_car(1,1)*lattice_car(3,2))
       sv(3,1)=(lattice_car(1,2)*lattice_car(2,3)-lattice_car(1,3)*lattice_car(2,2))
       sv(3,2)=(lattice_car(1,3)*lattice_car(2,1)-lattice_car(1,1)*lattice_car(2,3))
       sv(3,3)=(lattice_car(1,1)*lattice_car(2,2)-lattice_car(1,2)*lattice_car(2,1))

       ss=ss-(alpha+mu*(volume-volume0))*sv

       fvec(nv+1)=ss(1,1)
       fvec(nv+2)=ss(2,1)
       fvec(nv+3)=ss(3,1)
       fvec(nv+4)=ss(1,2)
       fvec(nv+5)=ss(2,2)
       fvec(nv+6)=ss(3,2)
       fvec(nv+7)=ss(1,3)
       fvec(nv+8)=ss(2,3)
       fvec(nv+9)=ss(3,3)
       
       fvec(nvec)=(volume-volume0)

    end if
    
  end subroutine force2vec

  function vecdotvec(a,b)

    real(kind=dp), dimension(:), intent(in) :: a,b
    
    real(kind=dp) :: vecdotvec

    real(kind=dp) :: cartcar(3,3)

    integer :: n,ns,nv
    
    cartcar=matmul(transpose(lattice_car),lattice_car)
    
    vecdotvec=0.0_dp

    nv=0
    do n=1,num_sets
       do ns=1,num_symm
                   
          if(ninset(n).gt.1) then
             vecdotvec=vecdotvec+dot_product(a(nv+1:nv+4),b(nv+1:nv+4))
             nv=nv+4
          end if
          
          vecdotvec=vecdotvec+dot_product(a(nv+1:nv+3),matmul(cartcar,b(nv+1:nv+3)))
          nv=nv+3
        
       end do
    end do
       
    if(celladapt) vecdotvec=vecdotvec+dot_product(a(nv+1:nv+10),b(nv+1:nv+10))
    
  end function vecdotvec

  subroutine test_derivatives(vec0,N)

    real(kind=dp), dimension(:), intent(inout) :: vec0
    
    integer,intent(in) :: N

    real(kind=dp) :: vec(N),grad(N),gradn(N),grad0(N),func0,func

    real(kind=dp) :: dg=1e-11_dp

    integer :: i
    
    call force2vec(vec0,grad0,func0)

    grad0=grad0*5
    
    write (stderr,*) func0
    
    do i=1,N
       vec=vec0
       vec(i)=vec(i)+dg
       call force2vec(vec,grad,func)
       gradn(i)=-(func-func0)/dg/5
       write (stderr,'(i4,3f15.6)') i,grad0(i),gradn(i),grad0(i)-gradn(i)
    end do
    
  end subroutine test_derivatives
  
  function qmult(q1,q2)
    
    real(kind=dp), intent(in) :: q1(4),q2(4)

    real(kind=dp) :: qmult(4)

    qmult(1) = q1(1)*q2(1)-q1(2)*q2(2)-q1(3)*q2(3)-q1(4)*q2(4)
    qmult(2) = q1(1)*q2(2)+q1(2)*q2(1)+q1(3)*q2(4)-q1(4)*q2(3)
    qmult(3) = q1(1)*q2(3)-q1(2)*q2(4)+q1(3)*q2(1)+q1(4)*q2(2)
    qmult(4) = q1(1)*q2(4)+q1(2)*q2(3)-q1(3)*q2(2)+q1(4)*q2(1)

  end function qmult

  function q2rot(q)
    
    real(kind=dp), intent(in) :: q(4)

    real(kind=dp) :: q2rot(3,3)

    ! ---
    
    real(kind=dp) :: nq

    nq=q(1)**2+q(2)**2+q(3)**2+q(4)**2

    
    if(nq.gt.0.0_dp) then
       
       q2rot(1,1) = 1.0_dp-(2.0_dp*q(3)**2+2.0_dp*q(4)**2)/nq
       q2rot(2,2) = 1.0_dp-(2.0_dp*q(2)**2+2.0_dp*q(4)**2)/nq
       q2rot(3,3) = 1.0_dp-(2.0_dp*q(2)**2+2.0_dp*q(3)**2)/nq
       
       q2rot(1,2) = (2.0_dp*q(2)*q(3)-2.0_dp*q(4)*q(1))/nq
       q2rot(2,1) = (2.0_dp*q(2)*q(3)+2.0_dp*q(4)*q(1))/nq
       
       q2rot(1,3) = (2.0_dp*q(2)*q(4)+2.0_dp*q(3)*q(1))/nq
       q2rot(3,1) = (2.0_dp*q(2)*q(4)-2.0_dp*q(3)*q(1))/nq
       
       q2rot(2,3) = (2.0_dp*q(3)*q(4)-2.0_dp*q(2)*q(1))/nq
       q2rot(3,2) = (2.0_dp*q(3)*q(4)+2.0_dp*q(2)*q(1))/nq

    else
       
       q2rot = 0.0_dp
       q2rot(1,1) = 1.0_dp
       q2rot(2,2) = 1.0_dp
       q2rot(3,3) = 1.0_dp

    end if

  end function q2rot

  function rotation_matrix(rangle)

    real(kind=dp), intent(inout) :: rangle(3)

    real(kind=dp) :: rotation_matrix(3,3)

    !------

    integer :: i,j
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
  
end module opt
