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
module list_mod
    use kinds_f
    use lattice_mod
    implicit none
    type, abstract :: list_t
        integer :: ist = 0
        integer :: ind = 0
    contains
        procedure(getijk_list), deferred :: getijk
        procedure :: getnumocells => getnumocells_list
    end type

interface getijk_fun
function getijk_list(this,icell) result(ijk)
    import :: list_t,numd
    class (list_t) :: this
    integer :: ijk(1:numd)
    integer,intent(in) :: icell
end function
end interface

contains

elemental integer function getnumocells_list(this) result(num)
    class(list_t),intent(in) :: this
    num = this%ind-this%ist+1
end function


end module
