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
module listlimits_mod
    use lattice_mod
    implicit none

    type listlimits2d_t
        integer :: blackwall,blackghost,fluid,ghostwall,ghostfluid,&
                   ghost,wall,influid,total,indomain
        integer, allocatable ::    bcfluid(:)
    contains
        procedure :: init => init_listlimits2d
    end type

    type, extends(listlimits2d_t) :: listlimits3d_t
    contains
        procedure :: init => init_listlimits3d
    end type

contains

subroutine init_listlimits2d(this,dims)
    class(listlimits2d_t) :: this
    integer :: dims(:)
    allocate(this%bcfluid(1:nbcsit))
    this%indomain = dims(1)*dims(2)
    this%total = (dims(1)+2)*(dims(2)+2)
end subroutine

subroutine init_listlimits3d(this,dims)
    class(listlimits3d_t) :: this
    integer :: dims(:)
    allocate(this%bcfluid(1:nbcsit))
    this%indomain = dims(1)*dims(2)*dims(3)
    this%total = (dims(1)+2)*(dims(2)+2)*(dims(3)+2)
end subroutine

end module
