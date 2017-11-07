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
! pbc (periodic boundary condition) is type which should be set
! before image processing. It helps setting cellkinds like
! ghost wall and ghost fluid.
! axis(1:3)  1=x, 2=y, 3=z
! dir(1:nbcsit)  : to check if a boundary (face,edge,or corner) is
! periodic. The cubevect points to all of them which are 8 boundaries
! in a 2D-domain and 26 in a 3D-domain.
module pbc_mod
    use lattice_mod
    implicit none

    type pbc_t
        logical, allocatable  :: axis(:)
        logical, allocatable  :: dir(:)
    contains
        procedure :: is_axis => is_axis_pbc
        procedure :: is_dir => is_dir_pbc
        procedure :: are_allaxes => are_allaxes_pbc
        procedure :: non_ofaxes  => non_ofaxes_pbc
        procedure :: init => init_pbc
    end type

contains 

    subroutine init_pbc(this,axis)
        class (pbc_t)     :: this
        logical,intent(in)  :: axis(:)
        integer :: q,d,a
        allocate(this%axis(1:numd),this%dir(1:nbcsit))
        this%axis(1:numd) = axis(1:numd)
        this%dir = .true.
        do d = 1,numd
            do q=1,nbcsit
                if (abs(cubevec(q,d)).eq.1.and..not.axis(d)) then
                    this%dir(q) = .false.
                endif
            enddo
        enddo
    end subroutine

    function is_axis_pbc(this,iaxis) result (val)
        class (pbc_t) :: this
        integer :: iaxis
        logical :: val
        val = this%axis(iaxis)
    end function

    function is_dir_pbc(this,idir) result (val)
        class (pbc_t) :: this
        integer :: idir
        logical :: val
        val = this%dir(idir)
    end function


    function are_allaxes_pbc(this) result (val)
        class (pbc_t) :: this
        integer :: idir,counter
        logical :: val
        counter = 0
        do idir = 1,numd
            if (this%is_axis(idir)) then
                counter = counter + 1
            end if
        end do

        if (counter.eq.numd) then
            val = .true.
        else
            val = .false.
        end if
    end function

    function non_ofaxes_pbc(this) result (val)
        class (pbc_t) :: this
        integer :: idir,counter
        logical :: val
        counter = 0
        do idir = 1,numd
            if (this%is_axis(idir)) then
                counter = counter + 1
            end if
        end do

        if (counter.eq.0) then
            val = .true.
        else
            val = .false.
        end if
    end function

end module pbc_mod
