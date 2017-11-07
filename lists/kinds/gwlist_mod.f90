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
module gwlist_mod
    use sublist_mod
    use stream_mod
    use collisionBGK_mod
    implicit none

    type, extends (ghostlist_t) :: gwlist_t
    contains
        procedure :: init01_sublist => init01_gwlist
    end type
contains

! Both stream and collision are initialized by default.
subroutine init01_gwlist(this,mlist,istinmlist,indinmlist,ist)
    class(gwlist_t),    target                  :: this
    class(mainlist_t),  pointer,    intent(in)  :: mlist
    integer,intent(in), optional                :: ist
    integer,intent(in)                          :: istinmlist,indinmlist
    type(streamwall_t)                          :: str
    type(wallcollision_t)                       :: col

    call this%ghostlist_t%init(mlist=mlist,istinmlist=istinmlist,indinmlist=indinmlist,ist=ist)
    call this%allocstream(source=str)
    call this%alloccollision(source=col)

end subroutine

end module
