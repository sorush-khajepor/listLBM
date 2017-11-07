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
module multipart_mod
    use partmcmp_mod
    implicit none

    type, abstract :: multipart_t
        class(partmcmppt_t), pointer :: parts(:)
        integer :: nparts
    contains
        procedure(stream_multipart),    deferred    :: stream
        procedure(collide_multipart),   deferred    :: collide
    end type



    type, extends(multipart_t):: multipart2d_t
    contains
        procedure :: init01_multipart2d => init01_multipart2d
        generic   :: init               => init01_multipart2d
        procedure :: stream => stream_multipart2d
        procedure :: collide => collide_multipart2d
    end type



interface collide_multipartsub
subroutine collide_multipart(this)
    import :: multipart_t
    class(multipart_t) :: this
end subroutine
end interface

interface stream_multipartsub
subroutine stream_multipart(this)
    import  :: multipart_t
    class(multipart_t) :: this
end subroutine
end interface



contains

subroutine init01_multipart2d(this,parts)
    class(multipart2d_t)            :: this
    class(partmcmppt_t),        target  :: parts(:)
    integer                         :: n

    this%parts => parts
    this%nparts = ubound(parts,dim=1)-lbound(parts,dim=1)+1
    do n = 1, this%nparts
        this%parts(n)%pt%partid = n
        this%parts(n)%pt%main%partid = n
        call this%parts(n)%pt%allocmcpot(this%nparts)
    end do

end subroutine

subroutine stream_multipart2d(this)
    class(multipart2d_t) :: this
    integer :: n

    do n = 1, this%nparts
        call this%parts(n)%pt%stream()
    end do

end subroutine


subroutine collide_multipart2d(this)
    class(multipart2d_t) :: this
    integer :: n

    do n = 1, this%nparts
        call this%parts(n)%pt%updatepot()
    end do

   ! print*,this%parts(1)%pt%main%mcpot
   ! print*,this%parts(2)%pt%main%mcpot
    do n = 1, this%nparts
        call this%parts(n)%pt%part2d_t%collide()
    end do

end subroutine





end module
