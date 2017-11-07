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
module outflowcrossbc_mod
    use collisionbc_mod
    use outflowbc_mod
    implicit none

! 	A type for handling boundary condition at an edge where
!	there is outflow condition besides an arbitrary boundary
!	condition.
	type,   extends(collisionbc_t)  :: outflowcrossbc_t
!		Outflow boundary condition(s) at an edge/corner
        class(outflowbc_t),     allocatable :: ofbc(:)
!		The other boundary condition
        class(collisionbc_t),   allocatable :: otherbc
    contains
        procedure :: init2d_outflowcrossbc  => init2d_outflowcrossbc
        generic   :: init                   => init2d_outflowcrossbc
!		Updating distributioin functions at the boundary
        procedure :: updatepdf              => updatepdf_outflowcrossbc
!		Setting sublist of this collision		
        procedure :: setsublist             => setsublist_outflowcrossbc
    end type

contains

subroutine init2d_outflowcrossbc(this,outflow,otherbc)
    class(outflowcrossbc_t)             :: this
    class(outflowbc_t),     intent(in)  :: outflow
    class(collisionbc_t),   intent(in)  :: otherbc
    integer                             :: iof,nof

!	number of outflow boundaries
    nof = 1
    allocate(this%ofbc(1:nof))
    allocate(this%otherbc,source = otherbc)
    call this%allocfluidcol(otherbc%fcol)
    do iof = 1, nof
        call this%ofbc(iof)%init(outflow%normal,outflow%fcol)
    end do

end subroutine

subroutine setsublist_outflowcrossbc(this,sl)
    class(outflowcrossbc_t)         :: this
    class(sublist_t),       target  :: sl
    integer                         :: iof

    this%sl=>sl
    this%fcol%sl=>sl
    call this%otherbc%setsublist(sl)
    do iof =1, ubound(this%ofbc,dim=1)
        call this%ofbc(iof)%setsublist(sl)
    end do

end subroutine

!	Each boundary condition updates the distribution functions
!	which have responsibility of.
subroutine updatepdf_outflowcrossbc(this,icell)
    class(outflowcrossbc_t)             :: this
    integer,                intent(in)  :: icell
    integer                             :: iof

    do iof = 1, ubound(this%ofbc,dim=1)
        call this%ofbc(iof)%updatepdf(icell)
    end do
    call this%otherbc%updatepdf(icell)

end subroutine
end module
