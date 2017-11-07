!============================== License GPLv3 ===================================
!	 This file is derived from the project:

!	 VTKFortran A KISS pure Fortran Library to parse and emit files conforming VTK (XML) standard
!    Copyright (C) 2012 Stefano Zaghi GNU Public License version 3
!	 https://github.com/szaghi/VTKFortran

!    and is published as a part of:

!    ListLBM is a lattice Boltzmann solver for multiphase flow in porous media
!    Copyright (C) 2017  Sorush Khajepor sk451@hw.ac.uk
!    Google Scholar:  https://scholar.google.co.uk/citations?user=16Y7-NsAAAAJ&hl
!	 Youtube channel: https://www.youtube.com/channel/UC5JVOcl3BjA1EkPMuKDZtGQ

!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.

!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.

!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <https://www.gnu.org/licenses/>.
!================================================================================

! To define a variable:
!        Real( Kind = WP ) :: var
! To define a constant:
!        const = 3.14_wp
! Type conversion:
!   a = Real( <integer value>, WP)
!   i = Int ( <real    value>, I2P)

Module kinds_f
Implicit None
! 	The following are the portable kind parameters available.
! 	Real precision definitions:
	integer, parameter:: R16P = selected_real_kind(33,4931) !< 33  digits, range \f$[10^{-4931}, 10^{+4931} - 1]\f$; 128 bits.
	integer, parameter:: R8P  = selected_real_kind(15,307)  !< 15  digits, range \f$[10^{-307} , 10^{+307}  - 1]\f$; 64 bits.
	integer, parameter:: R4P  = selected_real_kind(6,37)    !< 6   digits, range \f$[10^{-37}  , 10^{+37}   - 1]\f$; 32 bits.
	integer, parameter:: WP  = R8P                         !< Default real precision.

! 	Integer precision definitions:
	integer, parameter:: I8P  = selected_int_kind(18) !< Range \f$[-2^{63},+2^{63} - 1]\f$, 19 digits plus sign; 64 bits.
	integer, parameter:: I4P  = selected_int_kind(9)  !< Range \f$[-2^{31},+2^{31} - 1]\f$, 10 digits plus sign; 32 bits.
	integer, parameter:: I2P  = selected_int_kind(4)  !< Range \f$[-2^{15},+2^{15} - 1]\f$, 5  digits plus sign; 16 bits.
	integer, parameter:: I1P  = selected_int_kind(2)  !< Range \f$[-2^{7} ,+2^{7}  - 1]\f$, 3  digits plus sign; 8  bits.
	integer, parameter:: IP  = I4P                   !< Default integer precision.

! 	Format parameters useful for writing in a well-ascii-format numeric variables.
! 	Real output formats:
	character(10), parameter:: FR16P = '(E42.33E4)' !< Output format for kind=R16P variable.
	character(10), parameter:: FR8P  = '(E23.15E3)' !< Output format for kind=R8P variable.
	character(9),  parameter:: FR4P  = '(E13.6E2)'  !< Output format for kind=R4P variable.
	character(10), parameter:: FWP  = FR8P         !< Output format for kind=R_P variable.

! 	Real number of digits of output formats:
	integer, parameter:: DR16P = 42   !< Number of digits of output format FR16P.
	integer, parameter:: DR8P  = 23   !< Number of digits of output format FR8P.
	integer, parameter:: DR4P  = 13   !< Number of digits of output format FR4P.
	integer, parameter:: DR_P  = DR8P !< Number of digits of output format FR_P.

! 	Integer output formats
	character(5), parameter:: FI8P   = '(I20)'    !< Output format                     for kind=I8P variable.
	character(8), parameter:: FI8PZP = '(I20.19)' !< Output format with zero prefixing for kind=I8P variable.
	character(5), parameter:: FI4P   = '(I11)'    !< Output format                     for kind=I4P variable.
	character(8), parameter:: FI4PZP = '(I11.10)' !< Output format with zero prefixing for kind=I4P variable.
	character(4), parameter:: FI2P   = '(I6)'     !< Output format                     for kind=I2P variable.
	character(6), parameter:: FI2PZP = '(I6.5)'   !< Output format with zero prefixing for kind=I2P variable.
	character(4), parameter:: FI1P   = '(I4)'     !< Output format                     for kind=I1P variable.
	character(6), parameter:: FI1PZP = '(I4.3)'   !< Output format with zero prefixing for kind=I1P variable.
	character(5), parameter:: FI_P   = FI4P       !< Output format                     for kind=I_P variable.
	character(8), parameter:: FI_PZP = FI4PZP     !< Output format with zero prefixing for kind=I_P variable.

! 	Integer number of digits of output formats:
	integer, parameter:: DI8P = 20   !< Number of digits of output format I8P.
	integer, parameter:: DI4P = 11   !< Number of digits of output format I4P.
	integer, parameter:: DI2P = 6    !< Number of digits of output format I2P.
	integer, parameter:: DI1P = 4    !< Number of digits of output format I1P.
	integer, parameter:: DI_P = DI4P !< Number of digits of output format I_P.

! Epsilon is smallest number of its kind
	real(kind=wp), parameter :: eps_wp = epsilon(1._wp)

End Module kinds_f

