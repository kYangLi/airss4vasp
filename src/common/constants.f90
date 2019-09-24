!==================================================================================!
!                                    constants                                     !
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
! This module contains useful constants                                            !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module constants

  implicit none

  private

  integer, public, parameter :: dp=8 ! Single = 4 Double = 8 Quad = 16
  integer, public, parameter :: ki=8
  
  integer, public, parameter :: stdin=5
  integer, public, parameter :: stdout=6
  integer, public, parameter :: stderr=0

  real(kind=dp), public, parameter :: pi=3.141592653589793238462643383279502884197_dp
  real(kind=dp), public, parameter :: tpi=2.0_dp*pi
  real(kind=dp), public, parameter :: gr=(sqrt(5.0_dp)+1.0_dp)/2.0_dp
  real(kind=dp), public, parameter :: dgrd = pi/180.0_dp
  real(kind=dp), public, parameter :: evbyang3=160.2176487_dp
  real(kind=dp), public, parameter :: delta = 1e-13_dp

  real(kind=dp), public, parameter, dimension(3,3) :: &
       ident=reshape((/1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp/),(/3,3/))
  

  character(len=3), public, parameter, dimension(118) :: elements_alpha=(/&
       & 'Ac ','Ag ','Al ','Am ','Ar ','As ','At ','Au ','B  ','Ba '&
       &,'Be ','Bh ','Bi ','Bk ','Br ','C  ','Ca ','Cd ','Ce ','Cf '&
       &,'Cl ','Cm ','Co ','Cn ','Cr ','Cs ','Cu ','Db ','Ds ','Er '&
       &,'Es ','Eu ','F  ','Fe ','Fm ','Fr ','Ga ','Gd ','Ge ','H  '&
       &,'He ','Hf ','Hg ','Ho ','Hs ','I  ','In ','Ir ','K  ','Kr '&
       &,'La ','Li ','Lr ','Lu ','Md ','Mg ','Mn ','Mo ','Mt ','N  '&
       &,'Na ','Nb ','Nd ','Ne ','Ni ','No ','Np ','O  ','Os ','P  '&
       &,'Pa ','Pb ','Pd ','Pm ','Po ','Pr ','Pt ','Pu ','Ra ','Rb '&
       &,'Re ','Rf ','Rg ','Rh ','Rn ','Ru ','S  ','Sb ','Sc ','Se '&
       &,'Sg ','Si ','Sm ','Sn ','Sr ','Ta ','Tb ','Tc ','Te ','Th '&
       &,'Ti ','Tl ','Tm ','U  ','Uuh','Uun','UUo','Uup','Uuq','Uut'&
       &,'Uuu','V  ','W  ','Xe ','Y  ','Yb ','Zn ','Zr '/)

  integer, public, parameter, dimension(118) :: elements_valence=(/&
       & 0, 0, 3, 0, 8, 5, 0, 0, 3, 2 &
       &,2, 0, 0, 0, 0, 4, 2, 0, 0, 0 &
       &,7, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &,0, 0, 7, 0, 0, 0, 3, 0, 4, 1 &
       &,2, 0, 0, 0, 0, 7, 0, 0, 1, 8 &
       &,0, 1, 0, 0, 0, 2, 0, 0, 0, 5 &
       &,1, 0, 0, 8, 0, 0, 0, 6, 0, 5 &
       &,0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &,0, 0, 0, 0, 0, 0, 6, 0, 0, 6 &
       &,0, 4, 0, 4, 0, 0, 0, 0, 0, 0 &
       &,0, 0, 0, 0, 0, 0, 0, 0, 0, 0 &
       &,0, 0, 0, 8, 0, 0, 0 ,0/)
  
  
  ! ** Radii taken from Fig 3 Clementi, E.; Raimond, D. L.; Reinhardt, W. P. (1967).
  !    "Atomic Screening Constants from SCF Functions. II. Atoms with 37 to 86 Electrons".
  !    Journal of Chemical Physics. 47 (4): 1300–1307.
  !    H=0.39 Missing elements set to 2.00

  real(kind=dp), public, parameter, dimension(118) :: elements_radius=(/&
       &  2.00_dp,1.65_dp,1.18_dp,2.00_dp,0.71_dp,1.14_dp,1.27_dp,1.74_dp,0.87_dp,2.53_dp&
       &, 1.12_dp,2.00_dp,1.43_dp,2.00_dp,0.94_dp,0.67_dp,1.94_dp,1.61_dp,2.00_dp,2.00_dp&
       &, 0.79_dp,2.00_dp,1.52_dp,2.00_dp,1.66_dp,2.98_dp,1.45_dp,2.00_dp,2.00_dp,2.26_dp&
       &, 2.00_dp,2.31_dp,0.42_dp,1.56_dp,2.00_dp,2.00_dp,1.36_dp,2.33_dp,1.25_dp,0.39_dp&
       &, 0.31_dp,2.08_dp,1.71_dp,2.26_dp,2.00_dp,1.15_dp,1.56_dp,1.80_dp,2.43_dp,0.88_dp&
       &, 2.00_dp,1.67_dp,2.00_dp,2.17_dp,2.00_dp,1.45_dp,1.61_dp,1.90_dp,2.00_dp,0.56_dp&
       &, 1.90_dp,1.98_dp,2.06_dp,0.38_dp,1.49_dp,2.00_dp,2.00_dp,0.48_dp,1.85_dp,0.98_dp&
       &, 2.00_dp,1.54_dp,1.69_dp,2.05_dp,1.35_dp,2.47_dp,1.77_dp,2.00_dp,2.00_dp,2.00_dp&
       &, 1.88_dp,2.00_dp,2.00_dp,1.73_dp,1.20_dp,1.78_dp,0.88_dp,1.33_dp,1.84_dp,1.03_dp&
       &, 2.00_dp,1.11_dp,1.38_dp,1.45_dp,2.19_dp,2.00_dp,2.25_dp,1.83_dp,1.23_dp,2.00_dp&
       &, 1.76_dp,1.56_dp,2.22_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp&
       &, 2.00_dp,1.71_dp,1.93_dp,1.08_dp,2.12_dp,2.22_dp,1.42_dp,2.06_dp/)

  ! ** van der Waals raii taken from Table 12 of Mantina et al,
  !    "Consistent van der Waals Radii for the Whole Main Group",
  !    J Phys Chem A, 113(19), 5806-5812, 2009
  !    Missing elements set to 2.00

  real(kind=dp), public, parameter, dimension(118) :: elements_vdw=(/&
       &  2.00_dp,2.00_dp,1.84_dp,2.00_dp,1.88_dp,1.85_dp,2.02_dp,2.00_dp,1.92_dp,2.68_dp&
       &, 1.53_dp,2.00_dp,2.07_dp,2.00_dp,1.83_dp,1.70_dp,2.31_dp,2.00_dp,2.00_dp,2.00_dp&
       &, 1.75_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,3.43_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp&
       &, 2.00_dp,2.00_dp,1.47_dp,2.00_dp,2.00_dp,3.48_dp,1.87_dp,2.00_dp,2.11_dp,1.10_dp&
       &, 1.40_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,1.98_dp,1.93_dp,2.00_dp,2.75_dp,2.02_dp&
       &, 2.00_dp,1.81_dp,2.00_dp,2.00_dp,2.00_dp,1.73_dp,2.00_dp,2.00_dp,2.00_dp,1.55_dp&
       &, 2.27_dp,2.00_dp,2.00_dp,1.54_dp,2.00_dp,2.00_dp,2.00_dp,1.52_dp,2.00_dp,1.80_dp&
       &, 2.00_dp,2.02_dp,2.00_dp,2.00_dp,1.97_dp,2.00_dp,2.00_dp,2.00_dp,2.83_dp,3.03_dp&
       &, 2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.20_dp,2.00_dp,1.80_dp,2.06_dp,2.00_dp,1.90_dp&
       &, 2.00_dp,2.10_dp,2.00_dp,2.17_dp,2.49_dp,2.00_dp,2.00_dp,2.00_dp,2.06_dp,2.00_dp&
       &, 2.00_dp,1.96_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp& 
       &, 2.00_dp,2.00_dp,2.00_dp,2.16_dp,2.00_dp,2.00_dp,2.00_dp,2.00_dp/)
  
end module constants
