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
module partppot_mod
    use part_mod
    implicit none

    type,   extends(part2d_t) :: partppot2d_t
    contains
        procedure :: init01     => init01_partppot2d
        procedure :: collide    => collide_partppot2d
        procedure :: updatepot  => updatepot_partppot2d

    end type

    type partppotpt_t
        type(partppot2d_t), pointer :: pt
    end type

    type,   extends(partppot2d_t) :: partppot3d_t
    contains
        procedure :: allocmainlist  => allocmainlist_partppot3d
    end type

contains


subroutine init01_partppot2d(this,lislim,dims)
    class(partppot2d_t)                 :: this
    class(listlimits2d_t)               :: lislim
    integer, intent(in)                 :: dims(1:numd)

    call this%setlists(lislim,dims)
    if (.not.associated(this%main%pot)) &
            call this%main%allocpot()

end subroutine


subroutine updatepot_partppot2d(this)
    class(partppot2d_t) :: this
    integer :: n

    ! First step, all potentials in this part are updated.
    do n = 1, this%nslists
       call this%slists(n)%pt%updatepot()
    end do
end subroutine

subroutine collide_partppot2d(this)
    class(partppot2d_t) :: this
    integer :: n

    ! First step, all potentials in this part are updated.
    call this%updatepot()
    ! Then the main collision is happening.
    call this%part2d_t%collide()

end subroutine


subroutine allocmainlist_partppot3d(this)
    class(partppot3d_t)   :: this
    allocate(mainlist3d_t::this%main)
end subroutine

end module
