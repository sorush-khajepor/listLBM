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
module imgproc2d_mod
    use kinds_f
    use image2d_mod,   only : imageobj2d_t
    use part_mod,      only : part_t
    use lattice_mod
    use pbc_mod
    use search2d_mod
    use presetghosts2d_mod
    use unmap_mod
    use listlimits_mod
    implicit none

    type imgproc2d_t
        class (imageobj2d_t),       pointer     :: imageobj
        class (part_t),             pointer     :: part
        class (search2d_t),         allocatable :: search
        class (presetghosts2d_t),   allocatable :: psg
        class (unmap2d_t),          allocatable :: unmap
        class (listlimits2d_t),     allocatable :: lislim
        type (pbc_t)                            :: periodic
    contains
        procedure :: init01             => init01_imgproc2d
        generic   :: init               => init01
        procedure :: is_pbc             => is_pbc_imgproc2d
        procedure :: setpbc             => setpbc_imgproc2d
        procedure :: allocatepart       => allocatepart_imgproc2d
        procedure :: searchneighbours   => searchneighbours_imgproc2d
        procedure :: mapijk2icell       => mapijk2icell_imgproc2d
        procedure :: mapcells           => mapcells_imageproc2d
        procedure :: destruct           => destruct_imgproc2d
        procedure :: run                => run_imgproc2d
        procedure :: writesum           => writesummery_imgproc2d
    end type




!===================procedures================================
contains

subroutine destruct_imgproc2d(this)
    class(imgproc2d_t) :: this
    deallocate(this%unmap,this%search,this%psg)
end subroutine

subroutine init01_imgproc2d(this,imageobj,part)
    class(imgproc2d_t)              :: this
    class(imageobj2d_t),    target  :: imageobj
    class(part_t),          target  :: part

    this%part => part
    this%imageobj => imageobj
    allocate(this%search,this%psg,this%unmap,this%lislim)
    call this%unmap%init(imageobj)
    call this%lislim%init(imageobj%dims)

end subroutine

subroutine setpbc_imgproc2d (this,axis)
    class (imgproc2d_t) :: this
    logical             :: axis(:)

    call this%periodic%init(axis=axis)

end subroutine

! run all relevant procedures after initialization.
! except mapping procedures which can be replace each other
! the others should stay in the blow order. Inside mapping procedures
! can be manipulated.
subroutine run_imgproc2d(this)
    class(imgproc2d_t),target :: this
    write(*,'(A40)',advance='no')" Image processing is started (5 steps)  "
    call this%search%run(this%imageobj,this%periodic,this%lislim)
    write(*,'(A2)',advance='no')"1 "
    call this%allocatepart()
    write(*,'(A2)',advance='no')"2 "
    call this%mapcells
    write(*,'(A2)',advance='no')"3 "
    call this%searchneighbours()
    write(*,'(A2)',advance='no')"4 "
    call this%psg%run(this%part,this%imageobj,this%unmap,this%periodic)
    call this%writesum
    write(*,'(A12)')"5      done!"
end subroutine

! here assigning cellkinds to imageobj nodes is finished
! part are initialized.
! Note that if the list is empty ist=1, ind=0 which means
! the loops over that list won't be running.
! The order of assigning iend should be the same as
! initialization of slists(:) in the part class.
! for D2Q9 there are 8 Boundary lists.
subroutine allocatepart_imgproc2d (this)
    class (imgproc2d_t) :: this
    integer :: mxcelf,n
    integer,allocatable :: iend(:)
    associate(ll => this%lislim, img=>this%imageobj,part => this%part)

    call part%init(ll,img%dims)

    end associate
end subroutine


! A handy function which receives i,j,k and uses unmap to give
! icell or index of cell in main list.
integer function mapijk2icell_imgproc2d(this,ijk) result(icell)
    class(imgproc2d_t) :: this
    integer, intent(in) :: ijk(:)
    icell = this%unmap%geticell(ijk=ijk)
end function

! Neighbours of a cell are set to be Bin=node(0) by default.
subroutine searchneighbours_imgproc2d(this)
    class (imgproc2d_t) :: this
    integer :: i,j,q,icell
    integer :: ijk(1:numd),anei(1:numd),nei(1:numd,1:nq_1)
    associate (img=> this%imageobj,neilist=>this%part%main%neilist, &
               xlen=>this%imageobj%dims(1),ylen=>this%imageobj%dims(2))
    do j = 0, ylen+1
        do i = 0, xlen+1
            ijk = (/i,j/)
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
    end associate
end subroutine


subroutine mapcells_imageproc2d(this)
    class(imgproc2d_t) :: this
    integer :: n,i,j, icell,ijk(1:numd)
    associate(img=>this%imageobj, part => this%part,unmap=>this%unmap)
    do n = 1, part%nslists
        icell = part%slists(n)%pt%istinmlist
        do j = 0, img%dims(2)+1
            do i = 0, img%dims(1)+1
                ijk = (/i,j/)
                if (img%is_cellkind(ijk,part%slists(n)%pt%name)) then
                    call part%main%setrcell(icell,ijk)
                    call unmap%setcellijk(ijk=ijk,icell=icell)
                    part%main%islist(icell) = n
                    icell = icell + 1
                endif
            enddo
        enddo
    enddo
    end associate
end subroutine

function is_pbc_imgproc2d(this,axis) result (val)
    class(imgproc2d_t) :: this
    integer,intent(in) :: axis
    logical :: val
    val = this%periodic%is_axis(axis)
end function

subroutine writesummery_imgproc2d(this)
    class(imgproc2d_t)  :: this
    integer             :: funit,n
    associate (ll=>this%lislim)

    open(newunit = funit, file="imgproc.rep",status="replace")
    write(funit,*) "The number of cells are :"
    write(funit,*) "total w/o image-processing =", ll%total
    write(funit,*) "total w/o image-processing & w/o ghosts =", ll%indomain
    write(funit,*) "computational =", ll%total-ll%blackwall-ll%blackghost
    write(funit,*) "black =", ll%blackwall
    write(funit,*) "black ghost =", ll%blackghost
    do n=1,this%part%nslists
        write(funit,*) this%part%slists(n)%pt%name//' = ',this%part%slists(n)%pt%getnumocells()
    end do
    write(funit,*) "Porosity =", 1._wp- (ll%blackwall+ll%wall)/real(ll%indomain,kind=wp)
    close(funit)

    end associate
end subroutine
end module
