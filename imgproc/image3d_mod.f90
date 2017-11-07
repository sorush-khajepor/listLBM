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
module image3d_mod
    use image2d_mod
    implicit none

    type,extends(imageobj2d_t) :: imageobj3d_t
    contains
        procedure :: init01 => init01_imageobj3d
        procedure :: getimage => getimage_imageobj3d
        procedure :: getpbodycell=>getperiodicbodycell_imageobj3d
        procedure :: is_acell=>is_acell_imageobj3d
        procedure :: setimage=>setimage_imageobj3d
        procedure :: getpnei=>getperiodicneibour_imageobj3d
        procedure :: is_bcell=>is_boundarycell_imageobj3d
        procedure :: getbcsit=>getbcellsituation_imageobj3d
        procedure :: getghostsit=>getghostsituation_imageobj3d
        procedure :: getbclimits => getbclimits_imageobj3d
        procedure :: is_indomain => is_indomain_imageobj3d
    end type

contains

! Initialization of Imageobj
subroutine init01_imageobj3d(this,dims)
    class(imageobj3d_t) :: this
    integer, intent(in) :: dims(:)
    allocate(this%dims(1:numd))
    this%dims(1:numd) = dims(1:numd)
    ! ghost cells are added to the domain therefore
    ! Indexes start from 0 to xlen+1.
    allocate(this%image(0:dims(1)+1,0:dims(2)+1,0:dims(3)+1))
    allocate(this%cellkind)
    call this%cellkind%init()
    ! Image is initialized with fluid
    this%image = this%cellkind%fluid
end subroutine


function getimage_imageobj3d(this,ijk) result(img)
    class(imageobj3d_t) :: this
    integer, intent(in) :: ijk(:)
    integer :: img
    if(this%is_acell(ijk)) then
        img = int(this%image(ijk(1),ijk(2),ijk(3)),i4p)
    else
        img = int(this%cellkind%notacell,i4p)
    endif
end function

subroutine setimage_imageobj3d(this,ijk,cellkind)
    class(imageobj3d_t) :: this
    integer, intent(in) :: ijk(:)
    integer(i1p) :: cellkind
    this%image(ijk(1),ijk(2),ijk(3)) = cellkind
end subroutine

function getperiodicbodycell_imageobj3d(this,ijk) result(bodycell)
    class(imageobj3d_t) :: this
    integer :: bodycell(1:numd),ijk(:)
    integer :: x,y,z
    associate(xlen=>this%dims(1),ylen=>this%dims(2),zlen=>this%dims(3))

    if(this%is_indomain(ijk)) then
        print*,"the input cell is=",ijk
        stop "Error in getperiodicbodycell_imageobj3d : the input cell is not ghost."
    endif

    x=ijk(1)
    y=ijk(2)
    z=ijk(3)
    if (x.eq.xlen+1) then
        x = 1
    elseif (x.eq.0) then
        x = xlen
    endif

    if (y.eq.ylen+1) then
        y = 1
    elseif (y.eq.0) then
        y = ylen
    endif

    if (z.eq.zlen+1) then
        z = 1
    elseif (z.eq.0) then
        z = zlen
    endif

    bodycell(1:numd) = (/x,y,z/)
    end associate
end function

function is_acell_imageobj3d(this,ijk) result(is)
    class (imageobj3d_t) :: this
    integer, intent(in)::ijk(:)
    logical :: is

    if((ijk(1).ge.0.and.ijk(1).le.this%dims(1)+1).and. &
       (ijk(2).ge.0.and.ijk(2).le.this%dims(2)+1).and. &
       (ijk(3).ge.0.and.ijk(3).le.this%dims(3)+1)) then
        is = .true.
    else
        is = .false.
    endif
end function

function is_boundarycell_imageobj3d (this,ijk) result (is)
    class (imageobj3d_t) :: this
    integer,intent(in)   ::  ijk(:)
    logical :: is
    if (ijk(1).eq.1.or.ijk(1).eq.this%dims(1).or. &
        ijk(2).eq.1.or.ijk(2).eq.this%dims(2).or. &
        ijk(3).eq.1.or.ijk(3).eq.this%dims(3))    then
        is = .true.
    else
        is = .false.
    endif
end function

function getperiodicneibour_imageobj3d(this,ijk) result(pnei)
    class(imageobj3d_t) :: this
    integer,intent(in) :: ijk(:)
    integer :: dir(-1:1,1:numd),pnei(1:numd,1:nq_1)
    integer :: q,d
    dir(0,1) = ijk(1)
    dir(0,2) = ijk(2)
    dir(0,3) = ijk(3)
!   right direction +x
    dir(1,1) =  mod(ijk(1),this%dims(1))+1
!   top direction +y
    dir(1,2) =  mod(ijk(2),this%dims(2))+1
!   front direction +z
    dir(1,3) =  mod(ijk(3),this%dims(3))+1
!   left direction -x
    dir(-1,1) =  mod(ijk(1)+this%dims(1)-2,this%dims(1))+1
!   bottom direction -y
    dir(-1,2) =  mod(ijk(2)+this%dims(2)-2,this%dims(2))+1
!   back direction -z
    dir(-1,3) =  mod(ijk(3)+this%dims(3)-2,this%dims(3))+1

    do q=1,nq_1
        do d=1,numd
            pnei(d,q) = dir(qvec(q,d),d)
        enddo
    enddo

end function

! cubevector is used here. As qvec can be misleading.
function getbcellsituation_imageobj3d(this,ijk) result(sit)
    class(imageobj3d_t) :: this
    integer, intent(in) :: ijk(:)
    integer :: i,j,k,sit,dir
    logical :: loc(-1:1,1:nbcsit)
    associate(xlen=>this%dims(1),ylen=>this%dims(2),zlen=>this%dims(3))

    if(.not.this%is_bcell(ijk)) then
        print*,"The input cell is=",ijk
        stop "Error in getbcellsituation_imageobj3d: The input cell is not a boundary cell"
    endif

    i=ijk(1)
    j=ijk(2)
    k=ijk(3)
!   loc(cubevector,xyz)
    loc( 0,1) = i.ge.2.and.i.le.xlen-1
    loc(-1,1) = i.eq.1
    loc( 1,1) = i.eq.xlen

    loc( 0,2) = j.ge.2.and.j.le.ylen-1
    loc(-1,2) = j.eq.1
    loc( 1,2) = j.eq.ylen

    loc( 0,3) = k.ge.2.and.k.le.zlen-1
    loc(-1,3) = k.eq.1
    loc( 1,3) = k.eq.zlen


    do dir = 1,nbcsit
        if (loc(cubevec(dir,1),1).and.loc(cubevec(dir,2),2).and.loc(cubevec(dir,3),3)) then
            sit = dir
            exit
        endif
    enddo

    end associate
end function


subroutine getbclimits_imageobj3d(this,dir,l1xyz,l2xyz)
        class(imageobj3d_t)  :: this
        integer, intent(in)  :: dir
        integer, intent(out) :: l1xyz(1:numd),l2xyz(1:numd)
        logical :: rev(1:numd)
        integer :: irev(1:numd),cb(1:numd)
        associate(xlen=>this%dims(1),ylen=>this%dims(2),zlen=>this%dims(3))
        cb = cubevec(dir,1:numd)
        rev = cb
        rev = .not.rev
        irev = rev
        l1xyz(1) = irev(1)      + abs(cb(1)) * (cb(1)*(xlen-1)+xlen+1)/2
        l1xyz(2) = irev(2)      + abs(cb(2)) * (cb(2)*(ylen-1)+ylen+1)/2
        l1xyz(3) = irev(3)      + abs(cb(3)) * (cb(3)*(zlen-1)+zlen+1)/2
        l2xyz(1) = irev(1)*xlen + abs(cb(1)) * (cb(1)*(xlen-1)+xlen+1)/2
        l2xyz(2) = irev(2)*ylen + abs(cb(2)) * (cb(2)*(ylen-1)+ylen+1)/2
        l2xyz(3) = irev(3)*zlen + abs(cb(3)) * (cb(3)*(zlen-1)+zlen+1)/2
        end associate
    end subroutine

function getghostsituation_imageobj3d(this,ijk) result(sit)
    class(imageobj3d_t) :: this
    integer, intent(in) :: ijk(:)
    integer :: i,j,k,sit,dir
    logical :: loc(-1:1,1:numd)
    associate(xlen=>this%dims(1),ylen=>this%dims(2),zlen=>this%dims(3))

    if(this%is_indomain(ijk)) then
        print*,"the input cell is=",ijk
        stop "Error in getghostsituation_imageobj3d : the input cell is not ghost."
    endif
    i=ijk(1)
    j=ijk(2)
    k=ijk(3)
!   loc(cubevector,xyz)
    loc( 0,1) = i.ge.1.and.i.le.xlen
    loc(-1,1) = i.eq.0
    loc( 1,1) = i.eq.xlen+1

    loc( 0,2) = j.ge.1.and.j.le.ylen
    loc(-1,2) = j.eq.0
    loc( 1,2) = j.eq.ylen+1

    loc( 0,3) = k.ge.1.and.k.le.zlen
    loc(-1,3) = k.eq.0
    loc( 1,3) = k.eq.zlen+1


    do dir = 1,nbcsit
        if (loc(cubevec(dir,1),1).and.loc(cubevec(dir,2),2).and.loc(cubevec(dir,3),3)) then
            sit = dir
            exit
        endif
    enddo

    end associate
end function

! Returns if the cell ijk is in simulation domain or not. The
! ghost halo is not considered a part of simulation domain.
function is_indomain_imageobj3d (this,ijk) result (is)
    class (imageobj3d_t) :: this
    integer, intent(in) :: ijk(:)
    logical :: is
    if (ijk(1).ge.1.and.ijk(1).le.this%dims(1).and.&
        ijk(2).ge.1.and.ijk(2).le.this%dims(2).and.&
        ijk(3).ge.1.and.ijk(3).le.this%dims(3)) then
        is = .true.
    else
        is = .false.
    endif
end function

end module
