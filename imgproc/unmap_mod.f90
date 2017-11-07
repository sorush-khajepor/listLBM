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
module unmap_mod
    use image2d_mod
    implicit none

    ! unmap_t converts the cell indexes (i,j,k) to its position in
    ! the main list.The cell(:) array stores all cells' position
    ! of an input geometry including the ones which will be removed
    ! during simulation like black-wall cells. Therefore, while it
    ! is very helpful for image processing, it should be deallocated
    ! afterwards.
    type unmap2d_t
        integer, allocatable  :: xyzlen(:)
        integer, allocatable  :: cell(:)
    contains
        procedure ::    init01_unmap2d => init01_unmap2d
        generic   ::    init => init01_unmap2d
        procedure ::    geticell => geticell_unmap2d
        procedure ::    setcellijk => setcellijk_unmap2d
        procedure ::    arr2vec => arr2vec_unmap2d
        procedure ::    dealloc => dealloc_unmap2d
    end type

    type, extends(unmap2d_t) :: unmap3d_t
    contains
        procedure ::    init01_unmap2d => init01_unmap3d
        procedure ::    arr2vec => arr2vec_unmap3d
    end type
contains

! The cell vector is initialized with 0, or bin node in the main
! list.
    subroutine init01_unmap2d(this,img)
        class(unmap2d_t) :: this
        class(imageobj2d_t),pointer,intent(in) :: img
        integer :: allnodes
        allocate(this%xyzlen(1:2))
        this%xyzlen(1) = img%dims(1)
        this%xyzlen(2) = img%dims(2)
        allnodes = (img%dims(1)+2)*(img%dims(2)+2)
        allocate(this%cell(0:allnodes-1))
        this%cell(0:allnodes-1) = 0
    end subroutine

    ! converts cell indexes (i,j,k) to its position
    ! in the main list
    function geticell_unmap2d(this,ijk) result(icell)
        class(unmap2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer :: icell
        icell = this%cell(this%arr2vec(ijk))
    end function

    subroutine setcellijk_unmap2d(this,ijk,icell)
        class(unmap2d_t) :: this
        integer, intent(in) :: icell,ijk(:)
        this%cell(this%arr2vec(ijk)) = icell
    end subroutine

    ! converts the 3D index (i,j,k) to a scalar value
    function arr2vec_unmap2d(this,ijk) result(out)
        class(unmap2d_t) :: this
        integer, Intent(in):: ijk(:)
        Integer :: out
        integer :: i,j,xlen
        out = ijk(2)*(this%xyzlen(1)+2) + ijk(1)
    end function

    subroutine dealloc_unmap2d(this)
        class(unmap2d_t) :: this
        deallocate(this%cell)
    end subroutine

!========Unmap3d_t routines==============


    subroutine init01_unmap3d(this,img)
        class(unmap3d_t) :: this
        class(imageobj2d_t),pointer,intent(in) :: img
        integer :: allnodes
        allocate(this%xyzlen(1:3))
        this%xyzlen(1) = img%dims(1)
        this%xyzlen(2) = img%dims(2)
        this%xyzlen(3) = img%dims(3)
        allnodes = (img%dims(1)+2)*(img%dims(2)+2)*(img%dims(3)+2)
        allocate(this%cell(0:allnodes-1))
    end subroutine

    function arr2vec_unmap3d(this,ijk) result(out)
        class(unmap3d_t) :: this
        integer, Intent(in):: ijk(:)
        Integer :: out
        associate(xyz=>this%xyzlen)
            out = ijk(3)*(xyz(1)+2)*(xyz(2)+2) + ijk(2)*(xyz(1)+2) + ijk(1)
        end associate
    end function

end module unmap_mod
