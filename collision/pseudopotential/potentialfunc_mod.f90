!============================== License GPLv3 ===================================
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

module potentialfunc_mod
    use kinds_f
    implicit none

!	A type for storing pseudopotential parameters and potential function
    type :: potparam_t
!		Pseudopotential parameters
       	real(kind=wp)                        :: G,rho_0,a,b,Tc,T,Tr
!		A pointer for pseudopotential function
		procedure(getpottemplate),   pointer :: getpot
    end type

contains

!	Template of a pseudopotential function
    function getpottemplate (this,rho)
        class(potparam_t) :: this
        real(kind = wp), intent(IN):: rho
        real(kind = wp)::  getpottemplate
        getpottemplate = 0._wp
    end function

!	The potential for Carnahan-Starling equation of state
    function pot_CS (this,rho)
        class(potparam_t) :: this
        real(kind = wp), intent(in):: rho
        real(kind = wp):: pot_CS,tmp_p, tmp,tmp_sq
        associate (G=>this%G,a=>this%a,b=>this%b,T=>this%T)
        tmp = b*rho/4.0_wp
        tmp_p = rho*T * (1.0_wp + tmp + tmp*tmp - tmp*tmp*tmp) / ((1.0_wp-tmp)**3)  - a*rho*rho
        tmp_sq = 2.0_wp*(tmp_p - rho/3.0_wp)/G
        pot_CS = sqrt( abs(tmp_sq) )
        end associate
    end function pot_CS

!	The potential for van der Waals equation of state
    function pot_VW (this,rho)
        class(potparam_t) :: this
        real(kind = wp), intent(in) :: rho
        Real(kind = wp):: pot_VW, tmp_p,tmp
        associate (G=>this%G,a=>this%a,b=>this%b,T=>this%T)
        tmp_p = rho*T / (1.0_wp -b*rho) - a*(rho*rho)
        tmp = 2.0_wp * (tmp_p-rho/3.0_wp)/G
        pot_VW = sqrt(tmp)
        end associate
    end function

!	The Shan-Chen potential
    function pot_SC (this,rho)
        class(potparam_t) :: this
        real(kind = wp), intent(in):: rho
        real(kind = wp):: pot_SC
        pot_SC = this%rho_0 *(1-exp(-rho/this%rho_0))
    end function

!	The potential = rho
    function pot_givesrho (this,rho)
        class(potparam_t) :: this
        real(kind = wp), intent(IN):: rho
        real(kind = wp)::  pot_givesrho
        pot_givesrho = rho
    end function


end module
