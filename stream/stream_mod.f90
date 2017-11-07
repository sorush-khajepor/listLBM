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
module stream_mod
    use sublist_mod
    use ghostlist_mod
    implicit none

    type, extends(stream_t) :: streamfluid_t
    contains
        procedure :: stream1 => opposite_streamfluid
        procedure :: stream2 => reversepdfdirs_streamfluid
    end type

    type, extends(streamfluid_t) :: streamghostf_t
    contains
        procedure :: stream1 => opposite_streamghostf
        procedure :: stream2 => reversepdf_ghostcopy_streamghostf
    end type

    type, extends(streamfluid_t) :: streamwall_t
    contains
        procedure :: stream2 => blank_streamwall
    end type

    contains


    subroutine  opposite_streamfluid(this)
    class(streamfluid_t) :: this
    integer :: q,icell,dir
    real (kind=wp) :: tmp
    associate(ml=>this%sl%mlist,sl=>this%sl)

    do icell = sl%istinmlist,sl%indinmlist
        do q=1,nq_1/2
            dir = halfqs(q)
            tmp = ml%pdf(dir,icell)
            ml%pdf(dir,icell) = ml%pdf(opposite(dir),ml%neilist(dir,icell))
            ml%pdf(opposite(dir),ml%neilist(dir,icell)) = tmp
        end do
    enddo

    end associate
    end subroutine


    subroutine  reversepdfdirs_streamfluid(this)
    class(streamfluid_t) :: this
    integer :: q,icell,dir
    real (kind=wp) :: tmp
    associate(ml=>this%sl%mlist)

    do icell = this%sl%istinmlist,this%sl%indinmlist
        do q=1,nq_1/2
            dir = halfqs(q)
            tmp = ml%pdf(dir,icell)
            ml%pdf(dir,icell) = ml%pdf(opposite(dir),icell)
            ml%pdf(opposite(dir),icell) = tmp
        enddo
    enddo
    end associate
    end subroutine


    subroutine  blank_streamwall(this)
    class(streamwall_t) :: this

    end subroutine


    subroutine  opposite_streamghostf(this)
    class(streamghostf_t) :: this
    integer :: q,icell,dir
    associate(ml=>this%sl%mlist,sl=>this%sl)

    do icell = sl%istinmlist,sl%indinmlist
        do q=1,nq_1/2
            dir = halfqs(q)
            ml%pdf(dir,icell) = ml%pdf(opposite(dir),ml%neilist(dir,icell))
        end do
    enddo

    end associate
    end subroutine


    subroutine  reversepdf_ghostcopy_streamghostf(this)
    class(streamghostf_t) :: this
    real (kind=wp) :: tmp
    integer :: igcell,limit,body,icell,i,q
    associate(ml=>this%sl%mlist)

    call this%streamfluid_t%stream2

    ! Copys ghost fluid cells to their periodic correspondant
    select type (gfl=>this%sl)
    class is (ghostlist_t)

    do igcell = gfl%ist,gfl%ind
        limit = ubound(gfl%outlet(igcell)%dir,dim=1)
        body =  gfl%strbodycell(igcell)
        icell = gfl%getimcell(igcell)
        do i=1,limit
            q = gfl%outlet(igcell)%dir(i)
            ml%pdf(q,body) = ml%pdf(q,icell)
        end do
    end do

    end select

    end associate
    end subroutine


end module



