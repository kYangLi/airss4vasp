!==================================================================================!
!                                   modularity                                     !
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
! This module optimises the modularity of a network                                !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module modularity

  use constants
  use rng

  implicit none

  private

  integer                                                      :: nvert
  integer                                                      :: nedge
  integer(kind=ki)                               :: m2
  integer(kind=ki), allocatable, dimension(:,:)  :: adjacancy
  integer(kind=ki), allocatable, dimension(:,:)  :: ddia,ddiax
  integer(kind=ki), allocatable, dimension(:,:)  :: eab
  integer, allocatable, dimension(:,:)                         :: neighbours
  integer, allocatable, dimension(:,:)                         :: incom
  integer(kind=ki), allocatable, dimension(:)    :: degree
  integer(kind=ki), allocatable, dimension(:)    :: degree2
  integer, allocatable, dimension(:)                           :: nincom
  integer(kind=ki), allocatable, dimension(:)    :: dcom,dcomx
  integer, allocatable, dimension(:)                           :: numneigh
  integer, allocatable, dimension(:)                           :: edgwgt
  integer, allocatable, dimension(:,:)                         :: edges
 
  integer, public                               :: nunit
  integer, public                               :: rhigh
  integer(kind=ki), public        :: rasym

  integer, public, allocatable, dimension(:)    :: com
  real(kind=dp), public                         :: q

  public :: init_modularity
  public :: optimise_modularity
  public :: optimise_rash
  public :: modules_evaluate

contains

  subroutine init_modularity(adj,r)

    integer, intent(in), dimension(:,:) :: adj
    real(kind=dp), intent(in)           :: r

    ! *-

    integer :: i,j
    logical :: finish

    ! ** Initialise pseudo random number generator

    call init_pseudorandom()

    ! ** Allocate network data

    nvert=size(adj,1)

    if(allocated(adjacancy)) &
         deallocate(adjacancy,neighbours,degree,degree2,com,numneigh,incom,nincom)

    allocate(adjacancy(nvert,nvert),neighbours(nvert,nvert))
    allocate(degree(nvert),degree2(nvert),com(nvert),numneigh(nvert))
    allocate(incom(nvert,nvert),nincom(nvert))

    ! ** Set the internal adjacancy matrix to the one supplied

    adjacancy=adj

    ! ** Find rhigh and rasym

    rasym=0
    do i=1,nvert
       degree(i)=0
       do j=1,nvert
          if(i.ne.j) degree(i)=degree(i)+adjacancy(j,i)
       end do
       rasym=rasym+degree(i)
    end do
    m2=rasym

    rasym=-rasym/nvert

    rhigh=0
    finish=.false.
    do while(.not.finish)
       rhigh=rhigh+1
       finish=.true.
       do i=1,nvert
          do j=i+1,nvert         
             if(adjacancy(i,j)*(m2+rhigh*nvert).gt.(degree(i)+rhigh)*(degree(j)+rhigh)) then
                finish=.false.
                goto 313
             end if
          end do
       end do
313    continue
    end do

    ! ** Set the diagonal element (self loops - negative allowed)

    if(r.ge.0.0_dp) then
       do i=1,nvert
          adjacancy(i,i)=int((1.0_dp-r)*rasym)
       end do
    else
       do i=1,nvert
          adjacancy(i,i)=0
       end do
    end if
    
    ! ** Check for floating nodes - loop to self

    do i=1,nvert
       if((sum(adjacancy(:,i))-adjacancy(i,i)).eq.0) adjacancy(i,i)=1
    end do

    ! ** Count the degrees of each node, and total

    m2=0
    do i=1,nvert
       degree(i)=sum(adjacancy(:,i))
       degree2(i)=degree(i)*degree(i)
       m2=m2+degree(i)
    end do

    if(m2.le.0) then
       write (stderr,*) 'WARNING - modularity: increase r'
    end if

    ! ** Set the neighbourlist

    numneigh=0
    neighbours=0
    do i=1,nvert
       do j=1,nvert
          if((adjacancy(i,j).ne.0)) then
             numneigh(i)=numneigh(i)+1
             neighbours(numneigh(i),i)=j
          end if
       end do
    end do

    ! ** Identify edges and weights
    
    nedge=0
    do i=1,nvert
       do j=i,nvert
          if(adjacancy(i,j).ne.0) nedge=nedge+1
       end do
    end do

    if(allocated(edges)) deallocate(edges,edgwgt)
    allocate(edges(2,nedge),edgwgt(nedge))

    nedge=0
    edgwgt=0
    do i=1,nvert
       do j=i,nvert
          if(adjacancy(i,j).ne.0) then
             nedge=nedge+1
             edgwgt(nedge)=adjacancy(i,j)
             edges(1,nedge)=i
             edges(2,nedge)=j
             if(i.ne.j) edgwgt(nedge)=edgwgt(nedge)*2
          end if
       end do
    end do

  end subroutine init_modularity

  subroutine init_communities(com_init,com_max)

    character(len=*),  intent(in)       :: com_init
    integer, optional, intent(in)       :: com_max

    integer :: i

    ! ** Initialise the communities

    do i=1,nvert
       select case(com_init)
       case('RAND');
          if(present(com_max)) then
             com(i)=1+int(random_single()*com_max)
          else
             com(i)=1+int(random_single()*nvert)
          end if
       case('ALL');
          com(i)=i
       case('ONE');
          com(i)=1
       case default;
          stop 'Community initialisation not recognised'
       end select
    end do

  end subroutine init_communities

  subroutine modules_evaluate(n,mods,mods_q)

    integer, intent(in) :: n
    integer, intent(in), dimension(:,:) :: mods
    real(kind=dp), intent(out), dimension(:) :: mods_q

    integer :: i

    do i=1,n
       com=mods(:,i)
       call q_init(calc_eab=.false.,inunit=.true.) 
    
       mods_q(i)=q_modularity()
    end do

  end subroutine modules_evaluate

  subroutine optimise_rash()

    integer, parameter :: nmax=200
    real(kind=dp),parameter :: acc=0.2_dp
    integer :: n,como(nvert),i,nc
    real(kind=dp) :: qo

    qo=-huge(1.0_dp)
    do n=1,10
!!$       call init_communities("RAND")
       call init_communities("ALL")

       call optimise_modularity()

       if(q.gt.qo) then
          como=com
          qo=q
       end if

    end do

    com=como
    n=0
    write (stderr,'(a,f15.10,i5)') ' }',qo,nunit;flush(stdout)

    do while(n.lt.nmax)
       n=n+1

       ! * Shake

       do i=1,nvert
          if(random_single().lt.acc) then
             nc=1+int(random_single()*nunit*2)
             com(i)=nc
          end if
       end do

       call optimise_modularity()

       if(q.gt.qo) then
          como=com
          qo=q
          write (stderr,'(a,f15.10,2i5)') ' :',qo,nunit,n;flush(stdout)
          n=0
       else
          com=como
       end if
       
    end do

    q=qo
    com=como

    call q_init(calc_eab=.false.,inunit=.true.) 

  end subroutine optimise_rash

  subroutine optimise_modularity()

    integer,       parameter :: niter=100
    real(kind=dp), parameter :: pacc=0.1_dp

    real(kind=dp) :: qo
    integer       :: i,n,nc,como(nvert)
    logical       :: mgd,changed

    call q_init(calc_eab=.false.,inunit=.true.) 

    changed=.true.
    do while(changed)
    
       qo=q_modularity()
       
       como=com
       dcomx=dcom
       ddiax=ddia

       n=0
       do while(n.lt.niter)                 
          n=n+1

          ! * Relax

          call optimise_local()

          ! * Update

          q=q_modularity()

          if(q.gt.qo) then
             qo=q
             call q_init(calc_eab=.false.,inunit=.true.)
             !write (stderr,'(i5,a1)',advance='no') n,">";flush(stdout)
             como=com
             dcomx=dcom
             ddiax=ddia
             n=0
          else
             com=como
             dcom=dcomx
             ddia=ddiax
          end if

          ! * Shake

          do i=1,nvert
             if(random_single().lt.pacc) then
                if(numneigh(i).gt.0) then
                   nc=1+int(random_single()*numneigh(i))
                   call q_update(i,com(neighbours(nc,i)))
                end if
             end if
          end do         

       end do
       !write (stderr,*) '+';flush(stdout)
       com=como
       
       changed=.false.
       do

          call optimise_merge(mgd)
          if(mgd) then
             changed=.true.
          else
             exit
          end if

       end do

       call q_init(calc_eab=.false.,inunit=.true.) 

    end do
    
    q=q_modularity()

  end subroutine optimise_modularity

  subroutine optimise_local()

    integer :: i,nc,pi,perm(nvert)
    integer(kind=ki) :: vtemp(nunit),temp
    logical :: changed

    changed=.true.
    do while(changed)
       changed=.false.

       call permute(perm)
       do i=1,nvert
          pi=perm(i)

          
          vtemp=ddia(:,pi)-degree(pi)*dcom(:)

          temp=vtemp(com(pi))+degree2(pi)

          if(any(vtemp.gt.temp)) then
             nc=maxloc(vtemp,1)
             call q_update(pi,nc)
             changed=.true.
          end if

       end do

    end do

  end subroutine optimise_local

  subroutine optimise_merge(merged)

    logical, intent(out) :: merged

    integer                        :: i,na,nb,indx(2)
    integer(kind=ki) :: mdel(nunit,nunit)

    logical :: changed,avail(nunit,nunit)

    call q_init(calc_eab=.true.,inunit=.true.) 
    merged=.false.
    changed=.true.
    do while(changed)
       changed=.false.

       avail=.false.

       do na=1,nunit
          do nb=1,nunit
             mdel(na,nb)=q_merge(na,nb)
             if(mdel(na,nb).gt.0) avail(na,nb)=.true.
          end do
       end do

       do while(any(avail))
          changed=.true.
          merged=.true.
          indx=maxloc(mdel,avail)
          
          do i=1,nvert
             if(com(i).eq.indx(2)) com(i)=indx(1)
          end do

          avail(indx(1),:)=.false.
          avail(:,indx(1))=.false.
          avail(indx(2),:)=.false.
          avail(:,indx(2))=.false.

       end do

       if(changed) call q_init(calc_eab=.true.,inunit=.true.) 

    end do    

  end subroutine optimise_merge

  subroutine q_init(calc_eab,inunit)

    logical, intent(in) :: calc_eab
    logical, intent(in) :: inunit

    integer :: n1,n2,n3

    if(inunit) call q_inunit_pack()

    if(allocated(dcom)) deallocate(dcom,dcomx,ddia,ddiax)
    allocate(dcom(nunit),dcomx(nunit),ddia(nunit,nvert),ddiax(nunit,nvert))

    dcom=0
    do n1=1,nvert
       dcom(com(n1))=dcom(com(n1))+degree(n1)
    end do

    ddia=0
    do n1=1,nvert
       do n2=1,numneigh(n1)
          if(neighbours(n2,n1).eq.n1) cycle
          ddia(com(neighbours(n2,n1)),n1)=ddia(com(neighbours(n2,n1)),n1)+adjacancy(n1,neighbours(n2,n1))*m2
       end do
    end do

    if(calc_eab) then

       if(allocated(eab)) deallocate(eab)
       allocate(eab(nunit,nunit))
       
       eab=0
       do n1=1,nvert
          n2=com(n1)
          do n3=1,nvert
             if(adjacancy(n3,n1).ne.0) then
                eab(n2,com(n3))=eab(n2,com(n3))+adjacancy(n1,n3)
             end if
          end do
       end do
       eab=eab*m2
       
    end if

  end subroutine q_init

  subroutine q_update(vert,newcom)

    integer, intent(in) :: vert
    integer, intent(in) :: newcom

    integer :: n1,n2,oldcom

    oldcom=com(vert)
    com(vert)=newcom

    dcom(oldcom)=dcom(oldcom)-degree(vert)
    dcom(newcom)=dcom(newcom)+degree(vert)

    do n1=1,numneigh(vert)
       n2=neighbours(n1,vert)
       if(n2.eq.vert) cycle
       ddia(oldcom,n2)=ddia(oldcom,n2)-adjacancy(n2,vert)*m2
       ddia(newcom,n2)=ddia(newcom,n2)+adjacancy(n2,vert)*m2
    end do

  end subroutine q_update

  subroutine q_inunit_pack()

    integer :: n1,n2,iwork(nvert),iwork2(nvert)

    nunit=0
    iwork=0
    do n1=1,nvert
       if(.not.any((com(n1).eq.(com(1:n1-1))))) then
          nunit=nunit+1
          iwork(nunit)=com(n1)
       end if
    end do

    iwork2=0
    do n1=1,nvert
       do n2=1,nunit
          if(iwork(n2).eq.com(n1)) iwork2(n1)=n2
       end do
    end do
    com=iwork2

  end subroutine q_inunit_pack

  function q_merge(na,nb)
    
    integer, intent(in) :: na,nb
    integer(kind=ki) :: q_merge

    if(na.eq.nb) then
       q_merge=0
    else
       q_merge=eab(na,nb)-dcom(na)*dcom(nb)
    end if

  end function q_merge

  function q_modularity()

    real(kind=dp) :: q_modularity

    integer :: n1,n2
    integer(kind=ki) :: imod,ipart
 
    nincom=0
    do n1=1,nvert
       nincom(com(n1))=nincom(com(n1))+1
       incom(nincom(com(n1)),com(n1))=n1
    end do

    imod=0
    do n1=1,nedge
       if(com(edges(1,n1)).eq.com(edges(2,n1))) imod=imod+edgwgt(n1)*m2
    end do

    do n1=1,nunit

       ipart=0
       do n2=1,nincom(n1)
          ipart=ipart+degree(incom(n2,n1))
       end do
       imod=imod-ipart**2

    end do

    q_modularity=real(imod,dp)/real(m2,dp)**2

  end function q_modularity

  subroutine permute(list)
    
    integer, dimension(:), intent(out) :: list

    integer :: i,j

    do i=1,size(list)
       j=1+int(random_single()*i)
       list(i)=list(j)
       list(j)=i
    end do
    
  end subroutine permute

end module modularity
