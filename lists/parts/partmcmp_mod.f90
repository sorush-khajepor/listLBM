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
module partmcmp_mod
    use partppot_mod
    implicit none

    type,   extends(partppot2d_t) :: partmcmp2d_t
    contains
        procedure :: init01     => init01_partmcmp2d
        procedure :: allocmcpot => allocmcpot_partmcmp2d
    end type

    type,   extends(partmcmp2d_t) :: partmcmp3d_t
    contains
        procedure :: allocmainlist  => allocmainlist_partmcmp3d
    end type


    type :: partmcmppt_t
        class(partmcmp2d_t),pointer :: pt
    end type


contains

subroutine init01_partmcmp2d(this,lislim,dims)
    class(partmcmp2d_t)                 :: this
    class(listlimits2d_t)               :: lislim
    integer, intent(in)                 :: dims(1:numd)

    call this%setlists(lislim,dims)

end subroutine

subroutine allocmcpot_partmcmp2d(this,ncomps)
    class(partmcmp2d_t)     :: this
    integer, intent(in)     :: ncomps

    call this%main%allocmcpot(ncomps)

end subroutine


subroutine allocmainlist_partmcmp3d(this)
    class(partmcmp3d_t)   :: this
    allocate(mainlist3d_t::this%main)
end subroutine

end module
