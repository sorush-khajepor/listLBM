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
module outflowbc_mod
    use collisionbc_mod
    use lattutil_mod
    implicit none

    type,   extends(collisionbc_t)  ::  outflowbc_t
!		normal is the unit normal vector from boundary face toward outside
        integer,    allocatable :: normal(:)
!		Lattice directions (velocity vectors) toward inside domain
        integer,    allocatable :: inletdirs(:)
!		Lattice vector normal to boundary face toward outside
        integer                 :: qnormal
    contains
        procedure :: init01_outflowbc   => init01_outflowbc
        generic   :: init               => init01_outflowbc
        procedure :: updatepdf          => updatepdf_outflow
    end type

contains

subroutine init01_outflowbc(this,normal,fcol)
    class(outflowbc_t)              :: this
    class(collision_t)              :: fcol
    integer,            intent(in)  :: normal(:)

!	Assigning the bulk fluid collision which is applied to boundary
!	nodes too after having all their pdfs found (after updatepdf step). 
    call this%allocfluidcol(source=fcol)

    allocate(this%normal(1:numd))
    this%normal = normal

    this%qnormal = vec2q(vec=normal)

!	Finding inlet directions (velocity vectors), the vectors crossing
!	the boundary from outside to inside domain.
    call findcrossplanedirs(normal=-normal,dirs=this%inletdirs)

end subroutine

!	First order outflow boundary condition, for example,
! 	at right boundary pdf(i,x) = pdf(i,x-1).
subroutine updatepdf_outflow(this,icell)
    class(outflowbc_t)                  :: this
    integer,                intent(in)  :: icell
    integer                             :: n, nei
    associate(pdf=>this%sl%mlist%pdf,neilist=>this%sl%mlist%neilist)

    nei= neilist(opposite(this%qnormal),icell)
    do n = 1, ubound(this%inletdirs,dim=1)
        pdf(this%inletdirs(n),icell) = pdf(this%inletdirs(n),nei)
    end do

    end associate
end subroutine
end module
