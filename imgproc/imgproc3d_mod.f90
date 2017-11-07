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
module imgproc3d_mod
    use imgproc2d_mod
    use presetghosts3d_mod
    use search3d_mod
    implicit none
    type, extends( imgproc2d_t) :: imgproc3d_t
    contains
        procedure :: init01             => init01_imgproc3d
        procedure :: searchneighbours   => searchneighbours_imgproc3d
        procedure :: mapcells           => mapcells_imageproc3d
    end type
contains

subroutine init01_imgproc3d(this,imageobj,part)
    class(imgproc3d_t)              :: this
    class(imageobj2d_t),    target  :: imageobj
    class(part_t),          target  :: part

    this%part => part
    this%imageobj => imageobj
    allocate(presetghosts3d_t::this%psg)
    allocate(search3d_t::this%search)
    allocate(unmap3d_t:: this%unmap)
    allocate(listlimits3d_t:: this%lislim)
    call this%unmap%init(imageobj)
    call this%lislim%init(imageobj%dims)

end subroutine



! Neighbours of a cell are set to be Bin=node(0) by default.
subroutine searchneighbours_imgproc3d(this)
    class (imgproc3d_t) :: this
    integer :: i,j,k,q,icell
    integer :: ijk(1:numd),anei(1:numd),nei(1:numd,1:nq_1)
    associate (img=> this%imageobj,neilist=>this%part%main%neilist, &
               xlen=>this%imageobj%dims(1),ylen=>this%imageobj%dims(2),   &
               zlen=>this%imageobj%dims(3))

    do k = 0, zlen+1
        do j = 0, ylen+1
            do i = 0, xlen+1
                ijk = (/i,j,k/)
                if (img%is_cellkind(ijk,"blackwall").or. &
                    img%is_cellkind(ijk,"blackghost")) cycle

                icell = this%mapijk2icell(ijk)


                nei = img%getnei(ijk)


                if (img%is_indomain(ijk)) then
                    do q=1,nq_1
                        neilist(q,icell)=this%mapijk2icell(nei(1:numd,q))
                    enddo
                ! Ghost  cells neighbour search:
                ! Ghosts don't know each other as neighbour.
                else
                    do q=1,nq_1
                        anei = nei(1:numd,q)
                        if (img%is_acell(anei).and.img%isnt_g_bg(anei)) &
                           neilist(q,icell) = this%mapijk2icell(anei)
                    enddo
                end if
            enddo
        enddo
    enddo

    end associate
end subroutine


subroutine mapcells_imageproc3d(this)
    class(imgproc3d_t) :: this
    integer :: n,i,j,k, icell,ijk(1:numd)
    associate(img=>this%imageobj, part => this%part,unmap=>this%unmap)

    do n = 1, part%nslists
        icell = part%slists(n)%pt%istinmlist
        do k = 0, img%dims(3)+1
            do j = 0, img%dims(2)+1
                do i = 0, img%dims(1)+1
                    ijk = (/i,j,k/)
                    if (img%is_cellkind(ijk,part%slists(n)%pt%name)) then
                        call part%main%setrcell(icell,ijk)
                        call unmap%setcellijk(ijk=ijk,icell=icell)
                        part%main%islist(icell) = n
                        icell = icell + 1
                    endif
                enddo
            enddo
        enddo
    enddo

    end associate
end subroutine


end module
