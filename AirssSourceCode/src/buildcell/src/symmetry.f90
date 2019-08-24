! -*- mode: F90 ; mode: font-lock ; column-number-mode: true ; vc-back-end: RCS -*-
!==================================================================================!
!                                    symmetry                                      !
!==================================================================================!
!
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
! This module stores the symmetry information                                      !
!----------------------------------------------------------------------------------!
! Written by Chris Pickard, Copyright (c) 2005-2018                                !
!----------------------------------------------------------------------------------!
!                                                                                  !
!==================================================================================!

module symmetry

  use constants

  implicit none

  private

  public :: symm_by_name
  public :: symm_by_number

  integer,           public :: num_symm,symmetry_number
  real(kind=dp),     public :: clen(3), cang(3), symm_ops(3,4,48), Sym(3,4,48)
  character(len=10), public :: symmetry_name

  ! ** Space group ranking data for CCDC database of molecular crystals (2015)
  
  integer, dimension(230), public :: sgrank=(/14,2,15,19,4,61,33,62,9,1,60,5,29,13,148,12,11,7,18,88,56,43,92,&
       20,96,167,36,82,86,146,114,64,147,176,41,70,52,57,63,161,76,205,152,85,78,144,58,145,165,87,154,173,122,&
       31,169,45,225,198,170,142,54,68,155,130,166,72,163,8,110,217,159,34,59,220,160,73,55,178,79,113,80,143,23,&
       81,129,194,179,126,74,118,32,121,141,139,94,3,40,197,136,26,190,229,123,227,218,46,221,228,186,30,10,53,65,&
       222,37,182,137,192,204,84,104,138,206,66,71,164,98,203,77,219,106,128,17,135,50,120,150,230,124,158,196,44,&
       42,140,127,91,171,21,24,48,67,69,90,95,117,180,51,199,193,215,216,39,97,223,172,226,134,191,75,213,83,202,210,&
       16,109,181,131,133,47,175,209,212,108,112,116,119,201,22,185,102,151,189,214,224,38,168,174,211,6,103,125,49,&
       162,27,153,187,195,200,28,107,132,157,25,188,149,184,93,156,89,177,35,101,111,207,208,115,183,99,100,105/)
  
  integer       :: n,m,nn
  real(kind=dp) :: ang

contains

  subroutine symm_by_name(symm_name)

    character(len=10), intent(in) :: symm_name

    select case(symm_name)

    case('P1')
       call symm_by_number(1)
    case('P-1')
       call symm_by_number(2)
    case('P2')
       call symm_by_number(3)
    case('P21')
       call symm_by_number(4)
    case('C2')
       call symm_by_number(5)
    case('Pm')
       call symm_by_number(6)
    case('Pc')
       call symm_by_number(7)
    case('Cm')
       call symm_by_number(8)
    case('Cc')
       call symm_by_number(9)
    case('P2m')
       call symm_by_number(10)
    case('P21m')
       call symm_by_number(11)
    case('C2m')
       call symm_by_number(12)
    case('P2c')
       call symm_by_number(13)
    case('P21c')
       call symm_by_number(14)
    case('C2c')
       call symm_by_number(15)
    case('P222')
       call symm_by_number(16)
    case('P2221')
       call symm_by_number(17)
    case('P21212')
       call symm_by_number(18)
    case('P212121')
       call symm_by_number(19)
    case('C2221')
       call symm_by_number(20)
    case('C222')
       call symm_by_number(21)
    case('F222')
       call symm_by_number(22)
    case('I222')
       call symm_by_number(23)
    case('I212121')
       call symm_by_number(24)
    case('Pmm2')
       call symm_by_number(25)
    case('Pmc21')
       call symm_by_number(26)
    case('Pcc2')
       call symm_by_number(27)
    case('Pma2')
       call symm_by_number(28)
    case('Pca21')
       call symm_by_number(29)
    case('Pnc2')
       call symm_by_number(30)
    case('Pmn21')
       call symm_by_number(31)
    case('Pba2')
       call symm_by_number(32)
    case('Pna21')
       call symm_by_number(33)
    case('Pnn2')
       call symm_by_number(34)
    case('Cmm2')
       call symm_by_number(35)
    case('Cmc21')
       call symm_by_number(36)
    case('Ccc2')
       call symm_by_number(37)
    case('Amm2')
       call symm_by_number(38)
    case('Abm2')
       call symm_by_number(39)
    case('Ama2')
       call symm_by_number(40)
    case('Aba2')
       call symm_by_number(41)
    case('Fmm2')
       call symm_by_number(42)
    case('Fdd2')
       call symm_by_number(43)
    case('Imm2')
       call symm_by_number(44)
    case('Iba2')
       call symm_by_number(45)
    case('Ima2')
       call symm_by_number(46)
    case('Pmmm')
       call symm_by_number(47)
    case('Pnnn')
       call symm_by_number(48)
    case('Pccm')
       call symm_by_number(49)
    case('Pban')
       call symm_by_number(50)
    case('Pmma')
       call symm_by_number(51)
    case('Pnna')
       call symm_by_number(52)
    case('Pmna')
       call symm_by_number(53)
    case('Pcca')
       call symm_by_number(54)
    case('Pbam')
       call symm_by_number(55)
    case('Pccn')
       call symm_by_number(56)
    case('Pbcm')
       call symm_by_number(57)
    case('Pnnm')
       call symm_by_number(58)
    case('Pmmn')
       call symm_by_number(59)
    case('Pbcn')
       call symm_by_number(60)
    case('Pbca')
       call symm_by_number(61)
    case('Pnma')
       call symm_by_number(62)
    case('Cmcm')
       call symm_by_number(63)
    case('Cmca')
       call symm_by_number(64)
    case('Cmmm')
       call symm_by_number(65)
    case('Cccm')
       call symm_by_number(66)
    case('Cmma')
       call symm_by_number(67)
    case('Ccca')
       call symm_by_number(68)
    case('Fmmm')
       call symm_by_number(69)
    case('Fddd')
       call symm_by_number(70)
    case('Immm')
       call symm_by_number(71)
    case('Ibam')
       call symm_by_number(72)
    case('Ibca')
       call symm_by_number(73)
    case('Imma')
       call symm_by_number(74)
    case('P4')
       call symm_by_number(75)
    case('P41')
       call symm_by_number(76)
    case('P42')
       call symm_by_number(77)
    case('P43')
       call symm_by_number(78)
    case('I4')
       call symm_by_number(79)
    case('I41')
       call symm_by_number(80)
    case('P-4')
       call symm_by_number(81)
    case('I-4')
       call symm_by_number(82)
    case('P4m')
       call symm_by_number(83)
    case('P42m')
       call symm_by_number(84)
    case('P4n')
       call symm_by_number(85)
    case('P42n')
       call symm_by_number(86)
    case('I4m')
       call symm_by_number(87)
    case('I41a')
       call symm_by_number(88)
    case('P422')
       call symm_by_number(89)
    case('P4212')
       call symm_by_number(90)
    case('P4122')
       call symm_by_number(91)
    case('P41212')
       call symm_by_number(92)
    case('P4222')
       call symm_by_number(93)
    case('P42212')
       call symm_by_number(94)
    case('P4322')
       call symm_by_number(95)
    case('P43212')
       call symm_by_number(96)
    case('I422')
       call symm_by_number(97)
    case('I4122')
       call symm_by_number(98)
    case('P4mm')
       call symm_by_number(99)
    case('P4bm')
       call symm_by_number(100)
    case('P42cm')
       call symm_by_number(101)
    case('P42nm')
       call symm_by_number(102)
    case('P4cc')
       call symm_by_number(103)
    case('P4nc')
       call symm_by_number(104)
    case('P42mc')
       call symm_by_number(105)
    case('P42bc')
       call symm_by_number(106)
    case('14mm')
       call symm_by_number(107)
    case('I4cm')
       call symm_by_number(108)
    case('I41md')
       call symm_by_number(109)
    case('I41cd')
       call symm_by_number(110)
    case('P-42m')
       call symm_by_number(111)
    case('P-42c')
       call symm_by_number(112)
    case('P-421m')
       call symm_by_number(113)
    case('P-421c')
       call symm_by_number(114)
    case('P-4m2')
       call symm_by_number(115)
    case('P-4c2')
       call symm_by_number(116)
    case('P-4b2')
       call symm_by_number(117)
    case('P-4n2')
       call symm_by_number(118)
    case('I-4m2')
       call symm_by_number(119)
    case('I-4c2')
       call symm_by_number(120)
    case('I-42m')
       call symm_by_number(121)
    case('I-42d')
       call symm_by_number(122)
    case('P4mmm')
       call symm_by_number(123)
    case('P4mcc')
       call symm_by_number(124)
    case('P4nbm')
       call symm_by_number(125)
    case('P4nnc')
       call symm_by_number(126)
    case('P4mbm')
       call symm_by_number(127)
    case('P4mnc')
       call symm_by_number(128)
    case('P4nmm')
       call symm_by_number(129)
    case('P4ncc')
       call symm_by_number(130)
    case('P42mmc')
       call symm_by_number(131)
    case('P42mcm')
       call symm_by_number(132)
    case('P42nbc')
       call symm_by_number(133)
    case('P42nnm')
       call symm_by_number(134)
    case('P42mbc')
       call symm_by_number(135)
    case('P42mnm')
       call symm_by_number(136)
    case('P42nmc')
       call symm_by_number(137)
    case('P42ncm')
       call symm_by_number(138)
    case('I4mmm')
       call symm_by_number(139)
    case('I4mcm')
       call symm_by_number(140)
    case('I41amd')
       call symm_by_number(141)
    case('I41acd')
       call symm_by_number(142)
    case('P3')
       call symm_by_number(143)
    case('P31')
       call symm_by_number(144)
    case('P32')
       call symm_by_number(145)
    case('R3')
       call symm_by_number(146)
    case('P-3')
       call symm_by_number(147)
    case('R-3')
       call symm_by_number(148)
    case('P312')
       call symm_by_number(149)
    case('P321')
       call symm_by_number(150)
    case('P3112')
       call symm_by_number(151)
    case('P3121')
       call symm_by_number(152)
    case('P3212')
       call symm_by_number(153)
    case('P3221')
       call symm_by_number(154)
    case('R32')
       call symm_by_number(155)
    case('P3m1')
       call symm_by_number(156)
    case('P31m')
       call symm_by_number(157)
    case('P3c1')
       call symm_by_number(158)
    case('P31c')
       call symm_by_number(159)
    case('R3m')
       call symm_by_number(160)
    case('R3c')
       call symm_by_number(161)
    case('P-31m')
       call symm_by_number(162)
    case('P-31c')
       call symm_by_number(163)
    case('P-3m1')
       call symm_by_number(164)
    case('P-3c1')
       call symm_by_number(165)
    case('R-3m')
       call symm_by_number(166)
    case('R-3c')
       call symm_by_number(167)
    case('P6')
       call symm_by_number(168)
    case('P61')
       call symm_by_number(169)
    case('P65')
       call symm_by_number(170)
    case('P62')
       call symm_by_number(171)
    case('P64')
       call symm_by_number(172)
    case('P63')
       call symm_by_number(173)
    case('P-6')
       call symm_by_number(174)
    case('P6m')
       call symm_by_number(175)
    case('P63m')
       call symm_by_number(176)
    case('P622')
       call symm_by_number(177)
    case('P6122')
       call symm_by_number(178)
    case('P6522')
       call symm_by_number(179)
    case('P6222')
       call symm_by_number(180)
    case('P6422')
       call symm_by_number(181)
    case('P6322')
       call symm_by_number(182)
    case('P6mm')
       call symm_by_number(183)
    case('P6cc')
       call symm_by_number(184)
    case('P63cm')
       call symm_by_number(185)
    case('P63mc')
       call symm_by_number(186)
    case('P-6m2')
       call symm_by_number(187)
    case('P-6c2')
       call symm_by_number(188)
    case('P-62m')
       call symm_by_number(189)
    case('P-62c')
       call symm_by_number(190)
    case('P6mmm')
       call symm_by_number(191)
    case('P6mcc')
       call symm_by_number(192)
    case('P63mcm')
       call symm_by_number(193)
    case('P63mmc')
       call symm_by_number(194)
    case('P23')
       call symm_by_number(195)
    case('F23')
       call symm_by_number(196)
    case('I23')
       call symm_by_number(197)
    case('P213')
       call symm_by_number(198)
    case('I213')
       call symm_by_number(199)
    case('Pm-3')
       call symm_by_number(200)
    case('Pn-3')
       call symm_by_number(201)
    case('Fm-3')
       call symm_by_number(202)
    case('Fd-3')
       call symm_by_number(203)
    case('Im-3')
       call symm_by_number(204)
    case('Pa-3')
       call symm_by_number(205)
    case('Ia-3')
       call symm_by_number(206)
    case('P432')
       call symm_by_number(207)
    case('P4232')
       call symm_by_number(208)
    case('F432')
       call symm_by_number(209)
    case('F4132')
       call symm_by_number(210)
    case('I432')
       call symm_by_number(211)
    case('P4332')
       call symm_by_number(212)
    case('P4132')
       call symm_by_number(213)
    case('I4132')
       call symm_by_number(214)
    case('P-43m')
       call symm_by_number(215)
    case('F-43m')
       call symm_by_number(216)
    case('I-43m')
       call symm_by_number(217)
    case('P-43n')
       call symm_by_number(218)
    case('F-43c')
       call symm_by_number(219)
    case('I-43d')
       call symm_by_number(220)
    case('Pm-3m')
       call symm_by_number(221)
    case('Pn-3n')
       call symm_by_number(222)
    case('Pm-3n')
       call symm_by_number(223)
    case('Pn-3m')
       call symm_by_number(224)
    case('Fm-3m')
       call symm_by_number(225)
    case('Fm-3c')
       call symm_by_number(226)
    case('Fd-3m')
       call symm_by_number(227)
    case('Fd-3c')
       call symm_by_number(228)
    case('Im-3m')
       call symm_by_number(229)
    case('Ia-3d')
       call symm_by_number(230)
    Case default
       write (stderr,'(a,a)') 'Space group symmetry not recognised: ',trim(symm_name)
       stop
    end select
  end subroutine symm_by_name

  subroutine symm_by_number(symm_number)

    integer, intent(in) :: symm_number

    symmetry_number=symm_number
    
    select case(symm_number)

    case(1)

       symmetry_name="P1"

       num_symm = 1

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(2)

       symmetry_name="P-1"

       num_symm = 2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(3)

       symmetry_name="P2"

       num_symm = 2

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(4)

       symmetry_name="P21"

       num_symm = 2

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.5_dp,0.0_dp/)

    case(5)

       symmetry_name="C2"

       num_symm = 2

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(6)

       symmetry_name="Pm"

       num_symm = 2

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(7)

       symmetry_name="Pc"

       num_symm = 2

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(8)

       symmetry_name="Cm"

       num_symm = 2

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(9)

       symmetry_name="Cc"

       num_symm = 2

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(10)

       symmetry_name="P2m"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(11)

       symmetry_name="P21m"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.5_dp,0.0_dp/)

    case(12)

       symmetry_name="C2m"

       num_symm = 4

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(13)

       symmetry_name="P2c"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(14)

       symmetry_name="P21c"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 0.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.5_dp,0.5_dp/)

    case(15)

       symmetry_name="C2c"

       num_symm = 4

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(16)

       symmetry_name="P222"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(17)

       symmetry_name="P2221"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(18)

       symmetry_name="P21212"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(19)

       symmetry_name="P212121"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(20)

       symmetry_name="C2221"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(21)

       symmetry_name="C222"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(22)

       symmetry_name="F222"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 999.0_dp
       cang(2) = 999.0_dp
       cang(3) = 999.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(23)

       symmetry_name="I222"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(24)

       symmetry_name="I212121"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

    case(25)

       symmetry_name="Pmm2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(26)

       symmetry_name="Pmc21"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(27)

       symmetry_name="Pcc2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(28)

       symmetry_name="Pma2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.0_dp,0.0_dp/)

    case(29)

       symmetry_name="Pca21"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.0_dp,0.5_dp/)

    case(30)

       symmetry_name="Pnc2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.5_dp,0.5_dp/)

    case(31)

       symmetry_name="Pmn21"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(32)

       symmetry_name="Pba2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(33)

       symmetry_name="Pna21"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(34)

       symmetry_name="Pnn2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(35)

       symmetry_name="Cmm2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(36)

       symmetry_name="Cmc21"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(37)

       symmetry_name="Ccc2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(38)

       symmetry_name="Amm2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(39)

       symmetry_name="Abm2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,-0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,-0.5_dp,0.0_dp/)

    case(40)

       symmetry_name="Ama2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(41)

       symmetry_name="Aba2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,-0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,-0.5_dp,0.5_dp/)

    case(42)

       symmetry_name="Fmm2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 999.0_dp
       cang(2) = 999.0_dp
       cang(3) = 999.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(43)

       symmetry_name="Fdd2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 999.0_dp
       cang(2) = 999.0_dp
       cang(3) = 999.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,4) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.25_dp,0.25_dp,0.25_dp/)

    case(44)

       symmetry_name="Imm2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(45)

       symmetry_name="Iba2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

    case(46)

       symmetry_name="Ima2"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.5_dp,0.5_dp/)

    case(47)

       symmetry_name="Pmmm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(48)

       symmetry_name="Pnnn"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(49)

       symmetry_name="Pccm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(50)

       symmetry_name="Pban"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(51)

       symmetry_name="Pmma"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.0_dp,0.0_dp/)

    case(52)

       symmetry_name="Pnna"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.5_dp,0.5_dp/)

    case(53)

       symmetry_name="Pmna"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(54)

       symmetry_name="Pcca"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.0_dp,0.5_dp/)

    case(55)

       symmetry_name="Pbam"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(56)

       symmetry_name="Pccn"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.0_dp,0.5_dp/)

    case(57)

       symmetry_name="Pbcm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.5_dp,0.0_dp/)

    case(58)

       symmetry_name="Pnnm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(59)

       symmetry_name="Pmmn"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(60)

       symmetry_name="Pbcn"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(61)

       symmetry_name="Pbca"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(62)

       symmetry_name="Pnma"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(63)

       symmetry_name="Cmcm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(64)

       symmetry_name="Cmca"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(65)

       symmetry_name="Cmmm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(66)

       symmetry_name="Cccm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(67)

       symmetry_name="Cmma"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(68)

       symmetry_name="Ccca"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 0.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(69)

       symmetry_name="Fmmm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 999.0_dp
       cang(2) = 999.0_dp
       cang(3) = 999.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(70)

       symmetry_name="Fddd"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 999.0_dp
       cang(2) = 999.0_dp
       cang(3) = 999.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,8) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.25_dp,0.25_dp,0.25_dp/)

    case(71)

       symmetry_name="Immm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(72)

       symmetry_name="Ibam"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,1.0_dp/)

    case(73)

       symmetry_name="Ibca"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,1.0_dp/)

    case(74)

       symmetry_name="Imma"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 9999.0_dp
       cang(2) = 9999.0_dp
       cang(3) = 9999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(75)

       symmetry_name="P4"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(76)

       symmetry_name="P41"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.25_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.75_dp/)

    case(77) 

       symmetry_name="P42"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(78) 

       symmetry_name="P43"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.75_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.25_dp/)

    case(79)

       symmetry_name="I4"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(80)

       symmetry_name="I41"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.75_dp,1.25_dp,0.5_dp/)

    case(81)

       symmetry_name="P-4"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(82)

       symmetry_name="I-4"

       num_symm = 4

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(83)

       symmetry_name="P4m"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(84)

       symmetry_name="P42m"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(85)

       symmetry_name="P4n"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(86)

       symmetry_name="P42n"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(87)

       symmetry_name="I4m"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(88)

       symmetry_name="I41a"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/1.0_dp,1.0_dp,1.0_dp/)

    case(89)

       symmetry_name="P422"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(90)

       symmetry_name="P4212"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(91)

       symmetry_name="P4122"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.25_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.75_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.75_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.25_dp/)

    case(92)

       symmetry_name="P41212"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.25_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.75_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.25_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.75_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(93)

       symmetry_name="P4222"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(94)

       symmetry_name="P42212"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(95)

       symmetry_name="P4322"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.75_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.25_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.25_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.75_dp/)

    case(96)

       symmetry_name="P43212"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.75_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.25_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.75_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.25_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(97)

       symmetry_name="I422"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(98)

       symmetry_name="I4122"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(99)

       symmetry_name="P4mm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(100)

       symmetry_name="P4bm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(101)

       symmetry_name="P42cm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(102)

       symmetry_name="P42nm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(103)

       symmetry_name="P4cc"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(104)

       symmetry_name="P4nc"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(105)

       symmetry_name="P42mc"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(106)

       symmetry_name="P42bc"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(107)

       symmetry_name="I4mm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(108)

       symmetry_name="I4cm"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(109)

       symmetry_name="I41md"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.75_dp,1.25_dp,0.5_dp/)

    case(110)

       symmetry_name="I41cd"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/1.25_dp,0.75_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.25_dp,0.75_dp,0.5_dp/)

    case(111)

       symmetry_name="P-42m"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(112)

       symmetry_name="P-42c"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(113)

       symmetry_name="P-421m"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(114)

       symmetry_name="P-421c"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(115)

       symmetry_name="P-4m2"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(116)

       symmetry_name="P-4c2"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(117)

       symmetry_name="P-4b2"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(118)

       symmetry_name="P-4n2"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(119)

       symmetry_name="I-4m2"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(120)

       symmetry_name="I-4c2"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(121)

       symmetry_name="I-42m"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(122)

       symmetry_name="I-42d"

       num_symm = 8

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.75_dp,1.25_dp,0.5_dp/)

    case(123)

       symmetry_name="P4mmm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(124)

       symmetry_name="P4mcc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(125)

       symmetry_name="P4nbm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(126)

       symmetry_name="P4nnc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(127)

       symmetry_name="P4mbm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(128)

       symmetry_name="P4mnc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(129)

       symmetry_name="P4nmm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(130)

       symmetry_name="P4ncc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(131)

       symmetry_name="P42mmc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(132)

       symmetry_name="P42mcm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(133)

       symmetry_name="P42nbc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(134)

       symmetry_name="P42nnm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(135)

       symmetry_name="P42mbc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(136)

       symmetry_name="P42mnm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(137)

       symmetry_name="P42nmc"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(138)

       symmetry_name="P42ncm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(139)

       symmetry_name="I4mmm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,14) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(140)

       symmetry_name="I4mcm"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,14) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.0_dp/)

    case(141)

       symmetry_name="I41amd"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,13) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,14) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.75_dp,0.25_dp,0.5_dp/)

    case(142)

       symmetry_name="I41acd"

       num_symm = 16

       ! * Cell constraints

       cang(1) = 99999.0_dp
       cang(2) = 99999.0_dp
       cang(3) = 99999.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.25_dp,0.75_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/1.25_dp,0.75_dp,0.5_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.75_dp,0.25_dp,0.5_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.75_dp,1.25_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/1.0_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,14) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.25_dp,0.75_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/1.25_dp,0.75_dp,0.5_dp/)

    case(143)

       symmetry_name="P3"

       num_symm = 3

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(144)

       symmetry_name="P31"

       num_symm = 3

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.3333333333333_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.6666666666667_dp/)

    case(145)

       symmetry_name="P32"

       num_symm = 3

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.6666666666667_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.3333333333333_dp/)

    case(146)

       symmetry_name="R3"

       num_symm = 3

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = -1.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(147)

       symmetry_name="P-3"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(148)

       symmetry_name="R-3"

       num_symm = 6

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = -1.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(149)

       symmetry_name="P312"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(150)

       symmetry_name="P321"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(151)

       symmetry_name="P3112"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(152)

       symmetry_name="P3121"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

    case(153)

       symmetry_name="P3212"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(154)

       symmetry_name="P3221"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.333333333333333333333333_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.666666666666666666666667_dp/)

    case(155)

       symmetry_name="R32"

       num_symm = 6

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = -1.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(156)

       symmetry_name="P3m1"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(157)

       symmetry_name="P31m"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(158)

       symmetry_name="P3c1"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(159)

       symmetry_name="P31c"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(160)

       symmetry_name="R3m"

       num_symm = 6

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = -1.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(161)

       symmetry_name="R3c"

       num_symm = 6

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = -1.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(162)

       symmetry_name="P-31m"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(163)

       symmetry_name="P-31c"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(164)

       symmetry_name="P-3m1"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(165)

       symmetry_name="P-3c1"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(166)

       symmetry_name="R-3m"

       num_symm = 12

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = -1.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(167)

       symmetry_name="R-3c"

       num_symm = 12

       ! * Cell constraints

       cang(1) = -1.0_dp
       cang(2) = -1.0_dp
       cang(3) = -1.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(168)

       symmetry_name="P6"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(169)

       symmetry_name="P61"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.33333333333333333333333_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.66666666666666666666667_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.83333333333333333333333_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.16666666666666666666667_dp/)

    case(170)

       symmetry_name="P65"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.66666666666666666666667_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.33333333333333333333333_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.16666666666666666666667_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.83333333333333333333333_dp/)

    case(171)

       symmetry_name="P62"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.66666666666666666666667_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.33333333333333333333333_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.66666666666666666666667_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.33333333333333333333333_dp/)

    case(172)

       symmetry_name="P64"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.33333333333333333333333_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.66666666666666666666667_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.33333333333333333333333_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.66666666666666666666667_dp/)

    case(173)

       symmetry_name="P63"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(174)

       symmetry_name="P-6"

       num_symm = 6

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(175)

       symmetry_name="P6m"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(176)

       symmetry_name="P63m"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(177)

       symmetry_name="P622"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(178)

       symmetry_name="P6122"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.833333333333333333333333333_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.166666666666666666666666667_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.83333333333333333333333333_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.16666666666666666666666667_dp/)

    case(179)

       symmetry_name="P6522"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.166666666666666666666666667_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.833333333333333333333333333_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.16666666666666666666666667_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.83333333333333333333333333_dp/)

    case(180)

       symmetry_name="P6222"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.66666666666666666666666667_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.33333333333333333333333333_dp/)

    case(181)

       symmetry_name="P6422"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.333333333333333333333333333_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.666666666666666666666666667_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.33333333333333333333333333_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.66666666666666666666666667_dp/)

    case(182)

       symmetry_name="P6322"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(183)

       symmetry_name="P6mm"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(184)

       symmetry_name="P6cc"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(185)

       symmetry_name="P63cm"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(186)

       symmetry_name="P63mc"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(187)

       symmetry_name="P-6m2"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(188)

       symmetry_name="P-6c2"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(189)

       symmetry_name="P-62m"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(190)

       symmetry_name="P-62c"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(191)

       symmetry_name="P6mmm"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(192)

       symmetry_name="P6mcc"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,23) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(193)

       symmetry_name="P63mcm"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(194)

       symmetry_name="P63mmc"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 120.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,6) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,9) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,23) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.5_dp/)

    case(195)

       symmetry_name="P23"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(196)

       symmetry_name="F23"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(197)

       symmetry_name="I23"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(198)

       symmetry_name="P213"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.0_dp,0.5_dp/)

    case(199)

       symmetry_name="I213"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,1.0_dp,0.5_dp/)

    case(200)

       symmetry_name="Pm-3"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,19) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(201)

       symmetry_name="Pn-3"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,19) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,24) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(202)

       symmetry_name="Fm-3"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(203)

       symmetry_name="Fd-3"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,19) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,20) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,23) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,24) = (/0.25_dp,0.25_dp,0.25_dp/)

    case(204)

       symmetry_name="Im-3"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(205)

       symmetry_name="Pa-3"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,19) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.0_dp,0.5_dp/)

    case(206)

       symmetry_name="Ia-3"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,20) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,20) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,22) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,24) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,1.0_dp,0.5_dp/)

    case(207)

       symmetry_name="P432"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(208)

       symmetry_name="P4232"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(209)

       symmetry_name="F432"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(210)

       symmetry_name="F4132"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,17) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.25_dp,0.25_dp,0.25_dp/)

    case(211)

       symmetry_name="I432"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(212)

       symmetry_name="P4332"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.25_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.75_dp,0.75_dp,0.25_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.75_dp,0.25_dp,0.75_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.25_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.75_dp,0.25_dp,0.75_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.75_dp,0.75_dp,0.25_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.25_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.75_dp,0.75_dp,0.25_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.75_dp,0.25_dp,0.75_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.25_dp,0.25_dp,0.25_dp/)

    case(213)

       symmetry_name="P4132"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,0.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.75_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.75_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.25_dp,0.25_dp,0.75_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.25_dp,0.75_dp,0.25_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.75_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.25_dp,0.75_dp,0.25_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.75_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.25_dp,0.25_dp,0.75_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.75_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.25_dp,0.25_dp,0.75_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.25_dp,0.75_dp,0.25_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.75_dp,0.75_dp,0.75_dp/)

    case(214)

       symmetry_name="I4132"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/1.5_dp,1.5_dp,1.5_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/1.0_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/1.0_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/1.0_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/1.5_dp,1.5_dp,1.5_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/1.0_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/1.0_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/1.0_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/1.5_dp,1.5_dp,1.5_dp/)

    case(215)

       symmetry_name="P-43m"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(216)

       symmetry_name="F-43m"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,19) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(217)

       symmetry_name="I-43m"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,19) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(218)

       symmetry_name="P-43n"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(219)

       symmetry_name="F-43c"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,19) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,19) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,24) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(220)

       symmetry_name="I-43d"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,14) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/1.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/1.0_dp,1.5_dp,1.0_dp/)

       symm_ops(:,1,16) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/1.0_dp,1.0_dp,1.5_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/1.0_dp,1.0_dp,1.5_dp/)

       symm_ops(:,1,19) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,19) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/1.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/1.0_dp,1.5_dp,1.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/1.0_dp,1.5_dp,1.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/1.0_dp,1.0_dp,1.5_dp/)

       symm_ops(:,1,24) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,24) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/1.5_dp,1.0_dp,1.0_dp/)
       
    case(221)

       symmetry_name="Pm-3m"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,26) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,26) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,27) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,27) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,27) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,27) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,28) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,28) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,30) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,30) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,31) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,31) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,31) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,31) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,32) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,32) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,32) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,34) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,34) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,35) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,35) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,35) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,35) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,36) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,36) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,36) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,37) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,37) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,39) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,39) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,39) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,40) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,40) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,40) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,40) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,41) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,41) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,41) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,41) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,42) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,42) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,44) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,44) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,45) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,45) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,45) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,46) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,46) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,46) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,46) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,47) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,47) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(222)

       symmetry_name="Pn-3n"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,26) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,26) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,27) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,27) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,27) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,27) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,28) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,28) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,30) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,30) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,31) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,31) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,31) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,31) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,32) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,32) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,32) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,34) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,34) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,35) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,35) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,35) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,35) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,36) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,36) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,36) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,37) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,37) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,39) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,39) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,39) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,40) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,40) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,40) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,40) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,41) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,41) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,41) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,41) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,42) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,42) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,44) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,44) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,45) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,45) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,45) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,46) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,46) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,46) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,46) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,47) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,47) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(223)

       symmetry_name="Pm-3n"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,26) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,26) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,27) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,27) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,27) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,27) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,28) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,28) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,30) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,30) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,31) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,31) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,31) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,31) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,32) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,32) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,32) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,34) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,34) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,35) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,35) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,35) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,35) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,36) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,36) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,36) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,37) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,37) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,39) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,39) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,39) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,40) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,40) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,40) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,40) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,41) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,41) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,41) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,41) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,42) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,42) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,44) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,44) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,45) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,45) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,45) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,46) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,46) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,46) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,46) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,47) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,47) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(224)

       symmetry_name="Pn-3m"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 90.0_dp
       cang(2) = 90.0_dp
       cang(3) = 90.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,6) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,7) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,8) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,10) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,11) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,12) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,15) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,16) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,17) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,17) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,20) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,20) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,21) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,22) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,23) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,26) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,26) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,27) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,27) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,27) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,27) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,28) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,28) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,30) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,30) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,31) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,31) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,31) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,31) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,32) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,32) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,32) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,34) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,34) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,35) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,35) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,35) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,35) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,36) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,36) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,36) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,37) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,37) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,39) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,39) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,39) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,40) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,40) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,40) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,40) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,41) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,41) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,41) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,41) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,42) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,42) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,44) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,44) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,45) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,45) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,45) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,46) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,46) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,46) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,46) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,47) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,47) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(225)

       symmetry_name="Fm-3m"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,26) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,26) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,27) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,27) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,27) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,27) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,28) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,28) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,30) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,30) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,31) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,31) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,31) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,31) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,32) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,32) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,32) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,34) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,34) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,35) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,35) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,35) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,35) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,36) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,36) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,36) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,37) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,37) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,39) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,39) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,39) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,40) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,40) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,40) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,40) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,41) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,41) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,41) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,41) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,42) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,42) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,44) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,44) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,45) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,45) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,45) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,46) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,46) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,46) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,46) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,47) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,47) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(226)

       symmetry_name="Fm-3c"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,26) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,26) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,27) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,27) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,27) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,27) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,28) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,28) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,30) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,30) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,31) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,31) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,31) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,31) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,32) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,32) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,32) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,34) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,34) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,35) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,35) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,35) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,35) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,36) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,36) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,36) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,37) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,37) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,39) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,39) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,39) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,40) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,40) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,40) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,40) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,41) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,41) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,41) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,41) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,42) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,42) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,44) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,44) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,45) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,45) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,45) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,46) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,46) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,46) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,46) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,47) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,47) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.5_dp,0.5_dp,0.5_dp/)

    case(227)

       symmetry_name="Fd-3m"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,17) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,26) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,26) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,27) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,27) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,27) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,27) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,28) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,28) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,30) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,30) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,31) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,31) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,31) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,31) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,32) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,32) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,32) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,34) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,34) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,35) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,35) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,35) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,35) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,36) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,36) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,36) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,37) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,37) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,39) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,39) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,39) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,40) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,40) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,40) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,40) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,41) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,41) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,41) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,41) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,42) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,42) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,44) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,44) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,45) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,45) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,45) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,46) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,46) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,46) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,46) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,47) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,47) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(228)

       symmetry_name="Fd-3c"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 60.0_dp
       cang(2) = 60.0_dp
       cang(3) = 60.0_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,7) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,8) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,10) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,1.0_dp/)

       symm_ops(:,1,11) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,1.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,12) = (/1.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,15) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,15) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,16) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,17) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,17) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,18) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,18) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,20) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,20) = (/1.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,21) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,21) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.25_dp,1.25_dp,0.25_dp/)

       symm_ops(:,1,23) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,23) = (/0.25_dp,0.25_dp,1.25_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.25_dp,0.25_dp,0.25_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.75_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,26) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,26) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,26) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,26) = (/-0.25_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,27) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,27) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,27) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,27) = (/0.75_dp,0.75_dp,-0.25_dp/)

       symm_ops(:,1,28) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,28) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,28) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,28) = (/0.75_dp,-0.25_dp,0.75_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.75_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,30) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,30) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,30) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,30) = (/0.75_dp,-0.25_dp,0.75_dp/)

       symm_ops(:,1,31) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,31) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,31) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,31) = (/-0.25_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,32) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,32) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,32) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,32) = (/0.75_dp,0.75_dp,-0.25_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.75_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,34) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,34) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,3,34) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,34) = (/0.75_dp,0.75_dp,-0.25_dp/)

       symm_ops(:,1,35) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,2,35) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,35) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,35) = (/0.75_dp,-0.25_dp,0.75_dp/)

       symm_ops(:,1,36) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,36) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,36) = (/1.0_dp,1.0_dp,1.0_dp/)
       symm_ops(:,4,36) = (/-0.25_dp,0.75_dp,0.75_dp/)

       symm_ops(:,1,37) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,37) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,37) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,37) = (/0.5_dp,-0.5_dp,0.5_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,39) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,39) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,39) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,39) = (/-0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,40) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,40) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,40) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,40) = (/0.5_dp,0.5_dp,-0.5_dp/)

       symm_ops(:,1,41) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,41) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,41) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,41) = (/0.5_dp,-0.5_dp,0.5_dp/)

       symm_ops(:,1,42) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,42) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,42) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,42) = (/0.5_dp,0.5_dp,-0.5_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,44) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,44) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,44) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,44) = (/-0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,45) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,45) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,45) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,4,45) = (/0.5_dp,-0.5_dp,0.5_dp/)

       symm_ops(:,1,46) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,46) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,3,46) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,46) = (/-0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,47) = (/-1.0_dp,-1.0_dp,-1.0_dp/)
       symm_ops(:,2,47) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,47) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,47) = (/0.5_dp,0.5_dp,-0.5_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.5_dp,0.5_dp,0.5_dp/)
 
    case(229)

       symmetry_name="Im-3m"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,26) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,26) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,26) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,26) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,27) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,27) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,27) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,27) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,28) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,28) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,28) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,28) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,30) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,30) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,30) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,30) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,31) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,31) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,31) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,31) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,32) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,32) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,32) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,32) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,34) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,34) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,34) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,34) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,35) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,35) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,35) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,35) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,36) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,36) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,36) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,36) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,37) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,37) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,37) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,37) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,39) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,39) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,39) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,39) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,40) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,40) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,40) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,40) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,41) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,41) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,41) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,41) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,42) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,42) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,42) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,42) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,44) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,44) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,44) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,44) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,45) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,45) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,45) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,45) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,46) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,46) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,46) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,46) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,47) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,47) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,47) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,47) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(230)

       symmetry_name="Ia-3d"

       num_symm = 48

       ! * Cell constraints

       cang(1) = 109.4712206344907_dp
       cang(2) = 109.4712206344907_dp
       cang(3) = 109.4712206344907_dp

       clen(1) = -1.0_dp
       clen(2) = -1.0_dp
       clen(3) = -1.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,3) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,4) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,4) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,5) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,5) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,6) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,7) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,7) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,7) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,7) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,8) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,8) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,8) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,8) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,9) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,9) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,9) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,9) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,10) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,10) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,10) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,10) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,11) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,11) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,11) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,11) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,12) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,12) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,12) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,12) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,13) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,14) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,14) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,14) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,14) = (/1.5_dp,1.5_dp,1.5_dp/)

       symm_ops(:,1,15) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,15) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,15) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,15) = (/1.0_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,16) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,16) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,16) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,16) = (/1.0_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,17) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,17) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,17) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,17) = (/0.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,18) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,18) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,18) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,18) = (/1.0_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,19) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,19) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,19) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,19) = (/1.5_dp,1.5_dp,1.5_dp/)

       symm_ops(:,1,20) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,20) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,20) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,20) = (/1.0_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,21) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,21) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,21) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,21) = (/0.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,22) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,22) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,22) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,22) = (/1.0_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,23) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,23) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,23) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,23) = (/1.0_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,24) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,24) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,24) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,24) = (/1.5_dp,1.5_dp,1.5_dp/)

       symm_ops(:,1,25) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,25) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,25) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,25) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,26) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,26) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,26) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,26) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,27) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,27) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,27) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,27) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,28) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,28) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,28) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,28) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,29) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,29) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,29) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,29) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,30) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,30) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,30) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,30) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,31) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,31) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,31) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,31) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,32) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,32) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,32) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,32) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,33) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,33) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,33) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,33) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,34) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,34) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,34) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,34) = (/1.0_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,35) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,35) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,35) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,35) = (/0.5_dp,0.5_dp,1.0_dp/)

       symm_ops(:,1,36) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,36) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,36) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,36) = (/0.5_dp,1.0_dp,0.5_dp/)

       symm_ops(:,1,37) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,37) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,3,37) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,37) = (/1.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,38) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,38) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,38) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,38) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,39) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,39) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,39) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,39) = (/1.0_dp,1.0_dp,1.5_dp/)

       symm_ops(:,1,40) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,40) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,3,40) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,40) = (/1.0_dp,1.5_dp,1.0_dp/)

       symm_ops(:,1,41) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,2,41) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,41) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,4,41) = (/1.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,42) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,42) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,42) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,42) = (/1.0_dp,1.5_dp,1.0_dp/)

       symm_ops(:,1,43) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,43) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,43) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,43) = (/0.5_dp,0.5_dp,0.5_dp/)

       symm_ops(:,1,44) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,2,44) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,44) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,44) = (/1.0_dp,1.0_dp,1.5_dp/)

       symm_ops(:,1,45) = (/-1.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,45) = (/-1.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,45) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,45) = (/1.5_dp,1.0_dp,1.0_dp/)

       symm_ops(:,1,46) = (/1.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,2,46) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,46) = (/0.0_dp,-1.0_dp,1.0_dp/)
       symm_ops(:,4,46) = (/1.0_dp,1.0_dp,1.5_dp/)

       symm_ops(:,1,47) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,2,47) = (/1.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,3,47) = (/0.0_dp,1.0_dp,-1.0_dp/)
       symm_ops(:,4,47) = (/1.0_dp,1.5_dp,1.0_dp/)

       symm_ops(:,1,48) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,48) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,48) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,48) = (/0.5_dp,0.5_dp,0.5_dp/)

       ! ** Point groups

    case(-1)

       ! * C1

       symmetry_name="C1"

       num_symm = 1

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(-2)

       ! * Cs

       symmetry_name="Cs"

       num_symm = 2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(-3)

       ! * Ci

       symmetry_name="Ci"

       num_symm = 2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,2) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

    case(-14:-4)
       
       ! * Cn

       symmetry_name="Cn"

       num_symm = -symm_number-2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       do n=1,num_symm

          ang=tpi/real(num_symm,dp)*real(n-1,dp)

          symm_ops(:,1,n) = (/cos(ang),-sin(ang),0.0_dp/)
          symm_ops(:,2,n) = (/sin(ang),cos(ang),0.0_dp/)
          symm_ops(:,3,n) = (/0.0_dp,0.0_dp,1.0_dp/)
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)

       end do

    case(-19:-15)
       
       ! * Dn

       symmetry_name="Dn"

       num_symm = (-symm_number-13)*2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       do n=1,num_symm/2

          ang=tpi/real(num_symm/2,dp)*real(n-1,dp)

          symm_ops(:,1,n) = (/cos(ang),-sin(ang),0.0_dp/)
          symm_ops(:,2,n) = (/sin(ang),cos(ang),0.0_dp/)
          symm_ops(:,3,n) = (/0.0_dp,0.0_dp,1.0_dp/)
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)

       end do

       symm_ops(:,1,num_symm/2+1) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,num_symm/2+1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,num_symm/2+1) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,num_symm/2+1) = (/0.0_dp,0.0_dp,0.0_dp/)

       do n=num_symm/2+2,num_symm
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-num_symm/2),symm_ops(1:3,1:3,num_symm/2+1))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

    case(-24:-20)
       
       ! * Cnv

       symmetry_name="Cnv"

       num_symm = (-symm_number-18)*2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       do n=1,num_symm/2

          ang=tpi/real(num_symm/2,dp)*real(n-1,dp)

          symm_ops(:,1,n) = (/cos(ang),-sin(ang),0.0_dp/)
          symm_ops(:,2,n) = (/sin(ang),cos(ang),0.0_dp/)
          symm_ops(:,3,n) = (/0.0_dp,0.0_dp,1.0_dp/)
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)

       end do

       symm_ops(:,1,num_symm/2+1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,num_symm/2+1) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,num_symm/2+1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,num_symm/2+1) = (/0.0_dp,0.0_dp,0.0_dp/)

       do n=num_symm/2+2,num_symm
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-num_symm/2),symm_ops(1:3,1:3,num_symm/2+1))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

    case(-29:-25)
       
       ! * Cnh

       symmetry_name="Cnh"

       num_symm = (-symm_number-23)*2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       do n=1,num_symm/2

          ang=tpi/real(num_symm/2,dp)*real(n-1,dp)

          symm_ops(:,1,n) = (/cos(ang),-sin(ang),0.0_dp/)
          symm_ops(:,2,n) = (/sin(ang),cos(ang),0.0_dp/)
          symm_ops(:,3,n) = (/0.0_dp,0.0_dp,1.0_dp/)
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)

       end do

       symm_ops(:,1,num_symm/2+1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,num_symm/2+1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,num_symm/2+1) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,num_symm/2+1) = (/0.0_dp,0.0_dp,0.0_dp/)

       do n=num_symm/2+2,num_symm
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-num_symm/2),symm_ops(1:3,1:3,num_symm/2+1))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

    case(-31:-30)
       
       ! * Dnh

       symmetry_name="Dnh"

       num_symm = (-symm_number-28)*4

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       do n=1,num_symm/4

          ang=tpi/real(num_symm/4,dp)*real(n-1,dp)

          symm_ops(:,1,n) = (/cos(ang),-sin(ang),0.0_dp/)
          symm_ops(:,2,n) = (/sin(ang),cos(ang),0.0_dp/)
          symm_ops(:,3,n) = (/0.0_dp,0.0_dp,1.0_dp/)
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)

       end do

       symm_ops(:,1,num_symm/4+1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,num_symm/4+1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,num_symm/4+1) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,num_symm/4+1) = (/0.0_dp,0.0_dp,0.0_dp/)

       do n=num_symm/4+2,num_symm/2
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-num_symm/4),symm_ops(1:3,1:3,num_symm/4+1))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

       symm_ops(:,1,num_symm/2+1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,num_symm/2+1) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,num_symm/2+1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,num_symm/2+1) = (/0.0_dp,0.0_dp,0.0_dp/)
      
       do n=num_symm/2+2,num_symm
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-num_symm/2),symm_ops(1:3,1:3,num_symm/2+1))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

    case(-33:-32)
       
       ! * Dnd

       symmetry_name="Dnd"

       num_symm = (-symm_number-30)*4

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       do n=1,num_symm/4

          ang=tpi/real(num_symm/4,dp)*real(n-1,dp)

          symm_ops(:,1,n) = (/cos(ang),-sin(ang),0.0_dp/)
          symm_ops(:,2,n) = (/sin(ang),cos(ang),0.0_dp/)
          symm_ops(:,3,n) = (/0.0_dp,0.0_dp,1.0_dp/)
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)

       end do

       symm_ops(:,1,num_symm/4+1) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,num_symm/4+1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,num_symm/4+1) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,num_symm/4+1) = (/0.0_dp,0.0_dp,0.0_dp/)

       do n=num_symm/4+2,num_symm/2
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-num_symm/4),symm_ops(1:3,1:3,num_symm/4+1))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

       if(num_symm==8) then

          symm_ops(:,1,num_symm/2+1) = (/0.0_dp,1.0_dp,0.0_dp/)
          symm_ops(:,2,num_symm/2+1) = (/1.0_dp,0.0_dp,0.0_dp/)
          symm_ops(:,3,num_symm/2+1) = (/0.0_dp,0.0_dp,1.0_dp/)
          symm_ops(:,4,num_symm/2+1) = (/0.0_dp,0.0_dp,0.0_dp/)
      
       else if(num_symm==12) then
          
          symm_ops(:,1,num_symm/2+1) = (/-1.0_dp,0.0_dp,0.0_dp/)
          symm_ops(:,2,num_symm/2+1) = (/0.0_dp,-1.0_dp,0.0_dp/)
          symm_ops(:,3,num_symm/2+1) = (/0.0_dp,0.0_dp,-1.0_dp/)
          symm_ops(:,4,num_symm/2+1) = (/0.0_dp,0.0_dp,0.0_dp/)

       else

          stop 'Dnd - symmetry not coded'
          
       end if

       do n=num_symm/2+2,num_symm
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-num_symm/2),symm_ops(1:3,1:3,num_symm/2+1))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

   case(-38:-34)
       
       ! * Sn, n even

      symmetry_name="Sn"

       num_symm = (-symm_number-32)*2

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators

       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       ang=tpi/real(num_symm,dp)
       
       symm_ops(:,1,2) = (/cos(ang),-sin(ang),0.0_dp/)
       symm_ops(:,2,2) = (/sin(ang),cos(ang),0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       do n=3,num_symm
          symm_ops(1:3,1:3,n) = matmul(symm_ops(1:3,1:3,n-1),symm_ops(1:3,1:3,2))
          symm_ops(:,4,n) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

    case(-39)
       
       ! * T

       symmetry_name="T"

       num_symm = 12

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators
       
       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)


       nn=6
       do n=2,3
          do m=4,6
             nn=nn+1
             symm_ops(1:3,1:3,nn)=matmul(symm_ops(1:3,1:3,n),symm_ops(1:3,1:3,m))
             symm_ops(:,4,nn) = (/0.0_dp,0.0_dp,0.0_dp/)
          end do

       end do

    case(-40)
       
       ! * Td

       symmetry_name="Td"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators
       
       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)


       nn=6
       do n=2,3
          do m=4,6
             nn=nn+1
             symm_ops(1:3,1:3,nn)=matmul(symm_ops(1:3,1:3,n),symm_ops(1:3,1:3,m))
             symm_ops(:,4,nn) = (/0.0_dp,0.0_dp,0.0_dp/)
          end do

       end do

       symm_ops(:,1,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       nn=13
       do n=2,12
          nn=nn+1
          symm_ops(1:3,1:3,nn)=matmul(symm_ops(1:3,1:3,n),symm_ops(1:3,1:3,13))
          symm_ops(:,4,nn) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do

    case(-41)
       
       ! * Th

       symmetry_name="Th"

       num_symm = 24

       ! * Cell constraints

       cang(1) = 0.0_dp
       cang(2) = 0.0_dp
       cang(3) = 0.0_dp

       clen(1) = 0.0_dp
       clen(2) = 0.0_dp
       clen(3) = 0.0_dp

       ! * Symmetry operators
       
       symm_ops(:,1,1) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,1) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,1) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,1) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,2) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,2,2) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,2) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,4,2) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,3) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,3) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,3,3) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,4,3) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,4) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,4) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,4) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,4) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,5) = (/-1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,5) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,3,5) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,5) = (/0.0_dp,0.0_dp,0.0_dp/)

       symm_ops(:,1,6) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,2,6) = (/0.0_dp,-1.0_dp,0.0_dp/)
       symm_ops(:,3,6) = (/0.0_dp,0.0_dp,-1.0_dp/)
       symm_ops(:,4,6) = (/0.0_dp,0.0_dp,0.0_dp/)


       nn=6
       do n=2,3
          do m=4,6
             nn=nn+1
             symm_ops(1:3,1:3,nn)=matmul(symm_ops(1:3,1:3,n),symm_ops(1:3,1:3,m))
             symm_ops(:,4,nn) = (/0.0_dp,0.0_dp,0.0_dp/)
          end do

       end do

       symm_ops(:,1,13) = (/0.0_dp,1.0_dp,0.0_dp/)
       symm_ops(:,2,13) = (/1.0_dp,0.0_dp,0.0_dp/)
       symm_ops(:,3,13) = (/0.0_dp,0.0_dp,1.0_dp/)
       symm_ops(:,4,13) = (/0.0_dp,0.0_dp,0.0_dp/)

       nn=13
       do n=2,12
          nn=nn+1
          symm_ops(1:3,1:3,nn)=matmul(symm_ops(1:3,1:3,n),symm_ops(1:3,1:3,13))
          symm_ops(:,4,nn) = (/0.0_dp,0.0_dp,0.0_dp/)
       end do
 
    case default
       write (stderr,*) 'Point group symmetry not recognised'
       stop
    end select

    do n=1,num_symm
       symm_ops(1:3,1:3,n)=transpose(symm_ops(1:3,1:3,n))
    end do
    
  end subroutine symm_by_number

end module symmetry
