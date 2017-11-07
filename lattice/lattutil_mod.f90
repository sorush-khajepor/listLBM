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
module lattutil_mod
    use lattice_mod
    implicit none

contains

! Gets an integer vector and returns integer single-qvector which is
! with 0,1, and -1 components, like (4,-5,0) => (1,-1,0).
function vec2qvec(vec) result(singleqvec)
    integer,    intent(in)  :: vec(:)
    integer                 :: d,singleqvec(1:numd)

    singleqvec = 0
    do d = 1, numd
        if(vec(d).ne.0) singleqvec(d) = vec(d)/abs(vec(d))
    end do

end function

! This doesn't recognise q=0
! Converts an integer vector to its equivlant q vector
function vec2q(vec) result(qout)
    integer,    intent(in)  :: vec(:)
    integer                 :: ct,d,q,qout,singleqvec(1:numd)

    singleqvec = vec2qvec(vec)

    qout = -1
    do q = 1, nq_1
        ct = 0
        do d = 1, numd
            if (singleqvec(d).eq.qvec(q,d)) ct = ct + 1
        end do
        if (ct.eq.numd) qout = q
    end do

end function

! Gets the normal vector of a plane and finds q vectors crossing
! a plane toward the normal of it.
! They are used for Zou and He boundary condition 
subroutine findcrossplanedirs(normal,dirs)
    integer,intent(in)                      :: normal(:)
    integer,allocatable,    intent(inout)   :: dirs(:)
    integer                                 :: q,n,tmpd(1:nq_1),tmpn(1:numd)

    tmpn = vec2qvec(normal)
    n = 0
    do q = 1, nq_1
        if (sum(tmpn(1:numd)*qvec(q,1:numd)).gt.0) then
            n = n + 1
            tmpd(n) = q
        end if
    end do

    if (allocated(dirs)) deallocate(dirs)
    allocate(dirs(1:n))
    dirs(1:n) = tmpd(1:n)

end subroutine

! Gets the normal vector of a plane and finds q vectors tangent
! to the plane. Zero vector is one of them.
subroutine findtanplanedirs(normal,dirs)
    integer,intent(in)                      :: normal(:)
    integer,allocatable,    intent(inout)   :: dirs(:)
    integer                                 :: q,n,tmpd(1:nq_1),tmpn(1:numd)

    tmpn = vec2qvec(normal)
    n = 0
    do q = 0, nq_1
        if (sum(tmpn(1:numd)*qvec(q,1:numd)).eq.0) then
            n = n + 1
            tmpd(n) = q
        end if
    end do

    if (allocated(dirs)) deallocate(dirs)
    allocate(dirs(1:n))
    dirs(1:n) = tmpd(1:n)

end subroutine


end module
