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
module image2d_mod
    use kinds_f
    use lattice_mod
    implicit none

    ! During image processing each cell is assigned a kind which is
    ! defined here.
    type cellkind2d_t
        integer(kind=i1p) :: fluid=1, wall=2, ghostw=3, ghostf=4
        integer(kind=i1p) :: blackwall=5, ghost=6, blackghost=7
        ! influid=internal fluid cells without boundary ones.
        ! notacell = the cell doesn't exist.
        integer(kind=i1p) :: influid=8,notacell=9
        ! boundary cell of fluid kind (their value starts with 100).
        integer(kind=i1p), allocatable :: bcf(:)
    contains
        procedure :: init => init_cellkind2d
    end type

    ! Imageobj is the 2d (or 3d) object which includes a 3D array
    ! (image). It saves the cell kind of each cell.
    ! Such an object helps mapping the 3D box of simulation
    ! onto various lists with different collision types.
    ! dim = dimensions of the simulation domain in the form of
    ! (x,y,z).
    type imageobj2d_t
        integer(kind=i1p),allocatable       :: image(:,:,:)
        integer,          allocatable       :: dims(:)
        class (cellkind2d_t),allocatable    :: cellkind
    contains
        procedure :: init01 => init01_imageobj2d
        generic   :: init => init01
        procedure :: getimage => getimage_imageobj2d
        procedure :: setimage => setimage_imageobj2d
        procedure :: isnt_w_bw => isnt_wall_or_blackwall_imageobj2d
        procedure :: isnt_g_bg=>isnt_ghost_or_blackghost_imageobj2d
        procedure :: is_cellkind => is_cellkind_imageobj2d
        procedure :: is_indomain => is_indomain_imageobj2d
        procedure :: is_bcell => is_boundarycell_imageobj2d
        procedure :: is_acell => is_acell_imageobj2d
        procedure :: getpnei => getperiodicneibour_imageobj2d
        procedure :: getnei => getneibour_imageobj2d
        procedure :: getpbodycell => getperiodicbodycell_imageobj2d
        procedure :: getghostsit => getghostsituation_imageobj2d
        procedure :: getbcsit => getbcellsituation_imageobj2d
        procedure :: setcellas => setcellas_imageobj2d
        procedure :: getbclimits => getbclimits_imageobj2d
        procedure :: destruct => destruct_imageobj2d

    end type
contains


    !   Different lattices have different boundary cell
    !   locations
    subroutine init_cellkind2d(this)
        class(cellkind2d_t) :: this
        integer :: q
        allocate(this%bcf(1:nbcsit))
        do q=1,nbcsit
            this%bcf(q) = int(99+q,i1p)
        enddo
    end subroutine

    ! Destructor of Imageobj
    subroutine destruct_imageobj2d(this)
        class(imageobj2d_t):: this
        deallocate(this%image)
    end subroutine

    ! Initialization of Imageobj
    subroutine init01_imageobj2d(this,dims)
        class(imageobj2d_t) :: this
        integer, intent(in) :: dims(:)
        allocate(this%dims(1:numd))
        this%dims(1:numd) = dims(1:numd)
        ! ghost cells are added to the domain therefore
        ! Indexes start from 0 to xlen+1.
        allocate(this%image(0:dims(1)+1,0:dims(2)+1,1:1))
        allocate(this%cellkind)
        call this%cellkind%init()
        ! Image is initialized with fluid
        this%image = this%cellkind%fluid
    end subroutine

    ! Returning the image or cellkind of cell (ijk)
    function getimage_imageobj2d(this,ijk) result(img)
        class(imageobj2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer :: img
        if(this%is_acell(ijk)) then
            img = int(this%image(ijk(1),ijk(2),1),i4p)
        else
            img = int(this%cellkind%notacell,i4p)
        endif
    end function

    ! Set image of cell ijk as cellkind {fluid,influid,wall,...}
    subroutine setimage_imageobj2d(this,ijk,cellkind)
        class(imageobj2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer(i1p) :: cellkind
        this%image(ijk(1),ijk(2),1) = cellkind
    end subroutine

    function isnt_wall_or_blackwall_imageobj2d (this,ijk) result (is)
        class (imageobj2d_t) :: this
        integer, intent(in)::ijk(:)
        logical :: is
        if (this%getimage(ijk) .eq. this%cellkind%wall.or.   &
            this%getimage(ijk) .eq. this%cellkind%blackwall) then
            is =.false.
        else
            is =.true.
        end if
    end function

    function isnt_ghost_or_blackghost_imageobj2d (this,ijk) result (is)
        class (imageobj2d_t) :: this
        integer, intent(in)  ::  ijk(:)
        logical :: is
        if (this%getimage(ijk) .eq. this%cellkind%ghost      .or. &
            this%getimage(ijk) .eq. this%cellkind%blackghost .or. &
            this%getimage(ijk) .eq. this%cellkind%ghostw     .or. &
            this%getimage(ijk) .eq. this%cellkind%ghostf)    then
            is =.false.
        else
            is =.true.
        end if
    end function

    ! In a periodic domain each ghost cell is representive of a cell
    ! (bodycell) in the domain. This function finds the bodycell of the ghost
    ! cell.
    function getperiodicbodycell_imageobj2d(this,ijk) result(bodycell)
        class(imageobj2d_t) :: this
        integer :: bodycell(1:numd),ijk(:)
        integer :: x,y
        associate(xlen=>this%dims(1),ylen=>this%dims(2))

            if(this%is_indomain(ijk)) then
                print*,"the input cell is=",ijk
                stop "Error in getperiodicbodycell_imageobj2d : the input cell is not ghost."
            endif

            x=ijk(1)
            y=ijk(2)
            if (ijk(1).eq.xlen+1) then
                x = 1
            elseif (ijk(1).eq.0) then
                x = xlen
            endif
            if (ijk(2).eq.ylen+1) then
                y = 1
            elseif (ijk(2).eq.0) then
                y = ylen
            endif

            bodycell(1:2) = (/x,y/)
        end associate
    end function


    subroutine getbclimits_imageobj2d(this,dir,l1xyz,l2xyz)
        class(imageobj2d_t)  :: this
        integer, intent(in)  :: dir
        integer, intent(out) :: l1xyz(1:numd),l2xyz(1:numd)
        logical :: rev(1:numd)
        integer :: irev(1:numd),cb(1:numd)
        associate(xlen=>this%dims(1),ylen=>this%dims(2))
        cb = cubevec(dir,1:numd)
        rev = cb
        rev = .not.rev
        irev = rev
        l1xyz(1) = irev(1)      + abs(cb(1)) * (cb(1)*(xlen-1)+xlen+1)/2
        l1xyz(2) = irev(2)      + abs(cb(2)) * (cb(2)*(ylen-1)+ylen+1)/2
        l2xyz(1) = irev(1)*xlen + abs(cb(1)) * (cb(1)*(xlen-1)+xlen+1)/2
        l2xyz(2) = irev(2)*ylen + abs(cb(2)) * (cb(2)*(ylen-1)+ylen+1)/2
        end associate
    end subroutine

    ! Based on cube vector directions, the function returens
    ! the boundary direction the ghost cell is placed. It should be
    ! noted qvec is used for finding LB interoutlet neighbours.
    ! But cubevec is used for geometry analysis of the cells.
    function getghostsituation_imageobj2d(this,ijk) result(sit)
        class(imageobj2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer :: i,j,sit,dir
        logical :: loc(-1:1,1:numd)
        associate(xlen=>this%dims(1),ylen=>this%dims(2))

            if(this%is_indomain(ijk)) then
                print*,"the input cell is=",ijk
                stop "Error in getghostsituation_imageobj2d : the input cell is not ghost."
            endif

            i=ijk(1)
            j=ijk(2)
            !   loc(cubevector,xyz)
            loc( 0,1) = i.ge.1.and.i.le.xlen
            loc(-1,1) = i.eq.0
            loc( 1,1) = i.eq.xlen+1

            loc( 0,2) = j.ge.1.and.j.le.ylen
            loc(-1,2) = j.eq.0
            loc( 1,2) = j.eq.ylen+1


            do dir = 1,nbcsit
                if (loc(cubevec(dir,1),1).and.loc(cubevec(dir,2),2)) then
                    sit = dir
                    exit
                endif
            enddo

        end associate
    end function


    ! Based on cube vector directions, the function returens
    ! the boundary situation the ghost cell is placed.
    function getbcellsituation_imageobj2d(this,ijk) result(sit)
        class(imageobj2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer :: i,j,sit,dir
        logical :: loc(-1:1,1:nbcsit)
        associate(xlen=>this%dims(1),ylen=>this%dims(2))

            if(.not.this%is_bcell(ijk)) then
                print*,"The input cell is=",ijk
                stop "Error in getbcellsituation_imageobj2d: The input cell is not a boundary cell"
            endif

            i=ijk(1)
            j=ijk(2)
            !   loc(cubevector,xyz)
            loc( 0,1) = i.ge.2.and.i.le.xlen-1
            loc(-1,1) = i.eq.1
            loc( 1,1) = i.eq.xlen

            loc( 0,2) = j.ge.2.and.j.le.ylen-1
            loc(-1,2) = j.eq.1
            loc( 1,2) = j.eq.ylen

            do dir = 1,nbcsit
                if (loc(cubevec(dir,1),1).and.loc(cubevec(dir,2),2)) then
                    sit = dir
                    exit
                endif
            enddo

        end associate
    end function


    ! Returns the array of neighbours around cell ijk. A neighbour
    ! can be outside of simulation domain and ghost halo.
    function getneibour_imageobj2d(this,ijk) result(nei)
        class(imageobj2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer :: nei(1:numd,1:nq_1)
        integer :: q
        do q=1,nq_1
            nei(1:numd,q)=ijk(1:numd)+ qvec(q,1:numd)
        end do
    end function


    ! Returns the array of neighbours around cell ijk. If a
    ! neighbour is outside of domain, a periodic cell is given.
    function getperiodicneibour_imageobj2d(this,ijk) result(pnei)
        class(imageobj2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer :: pnei(1:numd,1:nq_1),dir(-1:1,1:numd)
        integer :: q,d
        dir(0,1) = ijk(1)
        dir(0,2) = ijk(2)
        !   right direction +x
        dir(1,1) =  mod(ijk(1),this%dims(1))+1
        !   top direction +y
        dir(1,2) =  mod(ijk(2),this%dims(2))+1
        !   left direction -x
        dir(-1,1) =  mod(ijk(1)+this%dims(1)-2,this%dims(1))+1
        !   bottom direction -y
        dir(-1,2) =  mod(ijk(2)+this%dims(2)-2,this%dims(2))+1

        do q=1,nq_1
            do d=1,numd
                pnei(d,q) = dir(qvec(q,d),d)
            enddo
        enddo

    end function


    ! Returns if the cell ijk is in simulation domain or not. The
    ! ghost halo is not considered a part of simulation domain.
    function is_indomain_imageobj2d (this,ijk) result (is)
        class (imageobj2d_t) :: this
        integer, intent(in) :: ijk(:)
        integer :: i,j
        logical :: is
        i = ijk(1)
        j = ijk(2)
        if (i.ge.1 .and.i.le.this%dims(1).and.j.ge.1.and.j.le.this%dims(2)) then
            is = .true.
        else
            is = .false.
        endif
    end function

    ! Returns if the cell ijk is a boundary cell.
    function is_boundarycell_imageobj2d (this,ijk) result (is)
        class (imageobj2d_t) :: this
        integer,intent(in)   ::  ijk(:)
        logical :: is
        if (ijk(1).eq.1.or.ijk(1).eq.this%dims(1).or.  &
            ijk(2).eq.1.or.ijk(2).eq.this%dims(2))      then
            is = .true.
        else
            is = .false.
        endif
    end function

    ! Sets cell ijk as the input cellkind. This subroutine
    ! gets the cellkind as a character variable.
    subroutine setcellas_imageobj2d(this,ijk,cellkind)
        class (imageobj2d_t) :: this
        integer, intent(in)::ijk(:)
        character(len=*) :: cellkind
        logical :: found
        integer :: q
        found = .true.
        if (cellkind.eq."fluid") then
            call this%setimage(ijk,this%cellkind%fluid)
        elseif(cellkind.eq."wall") then
            call this%setimage(ijk,this%cellkind%wall)
        elseif(cellkind.eq."blackwall") then
            call this%setimage(ijk,this%cellkind%blackwall)
        elseif(cellkind.eq."blackghost") then
            call this%setimage(ijk,this%cellkind%blackghost)
        elseif(cellkind.eq."ghostf") then
            call this%setimage(ijk,this%cellkind%ghostf)
        elseif(cellkind.eq."ghostw") then
            call this%setimage(ijk,this%cellkind%ghostw)
        elseif(cellkind.eq."ghost") then
            call this%setimage(ijk,this%cellkind%ghost)
        elseif(cellkind.eq."influid") then
            call this%setimage(ijk,this%cellkind%influid)
        else
            found = .false.
        endif
        ! Cheking all boundary cell situations {rig,toprig,lef,...}
        do q = 1,nbcsit
            if(cellkind.eq.bcsname(q)%mem) then
                call this%setimage(ijk,this%cellkind%bcf(q))
                found = .true.
            endif
        enddo
        if (.not.found) then
            print *,"The requested cellkind is =",cellkind
            stop "Error: in  setcellas_imageobj2d cellkind is not recognized."
        endif
    end subroutine

    ! Returns if the cell exists in simulation domain or ghost halo.
    function is_acell_imageobj2d(this,ijk) result(is)
        class (imageobj2d_t) :: this
        integer, intent(in)::ijk(:)
        logical :: is
        if((ijk(1).ge.0.and.ijk(1).le.this%dims(1)+1).and. &
            (ijk(2).ge.0.and.ijk(2).le.this%dims(2)+1)) then
            is = .true.
        else
            is = .false.
        endif
    end function

    ! Returns if cell ijk has the input cellkind or not. The cellkind
    ! is a character variable.
    ! If ghost cellkind asked as input it considered as any of
    ! {ghost,ghostw and ghostf} cellkinds but not blackghost.
    ! if fluid cellkind asked as input it considered as any of
    ! {fluid,influid and all boundary fluid nodes}.
    function is_cellkind_imageobj2d(this,ijk,cellkind) result(is)
        class (imageobj2d_t) :: this
        integer, intent(in)::ijk(:)
        character(len=*) :: cellkind
        logical :: is,found
        integer :: q
        is = .false.
        found = .true.
        if (cellkind.eq."fluid") then
            if(this%getimage(ijk) .eq. this%cellkind%fluid)   is = .true.
            if(this%getimage(ijk) .eq. this%cellkind%influid)  is = .true.
            do q=1,nbcsit
                if(this%getimage(ijk) .eq. this%cellkind%bcf(q))    is = .true.
            enddo
        elseif(cellkind.eq."wall") then
            if(this%getimage(ijk) .eq. this%cellkind%wall) is = .true.
        elseif(cellkind.eq."blackwall") then
            if(this%getimage(ijk) .eq. this%cellkind%blackwall) is = .true.
        elseif(cellkind.eq."blackghost") then
            if(this%getimage(ijk) .eq. this%cellkind%blackghost) is = .true.
        elseif(cellkind.eq."ghostf") then
            if(this%getimage(ijk) .eq. this%cellkind%ghostf) is = .true.
        elseif(cellkind.eq."ghostw") then
            if(this%getimage(ijk) .eq. this%cellkind%ghostw) is = .true.
        elseif(cellkind.eq."ghost") then
            if(this%getimage(ijk) .eq. this%cellkind%ghost)  is = .true.
            if(this%getimage(ijk) .eq. this%cellkind%ghostf) is = .true.
            if(this%getimage(ijk) .eq. this%cellkind%ghostw) is = .true.
        elseif(cellkind.eq."influid") then
            if(this%getimage(ijk) .eq. this%cellkind%influid)   is = .true.
        elseif(cellkind.eq."notacell") then
            if(this%getimage(ijk) .eq. this%cellkind%notacell)   is = .true.
        else
            found = .false.
        endif

        do q=1,nbcsit
            if(cellkind.eq.bcsname(q)%mem) then
                found = .true.
                if (this%getimage(ijk) .eq. this%cellkind%bcf(q)) is = .true.
            endif
        enddo

        if(.not.found) then
            print *,"The requested cellkind is =",cellkind
            stop "Error: in  is_cellkind_imageobj2d cellkind is not recognized."
        endif

    end function
end module
