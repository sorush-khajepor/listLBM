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
module part_mod
    use kinds_f
    use mainlist_mod
    use sublist_mod
    use ghostlist_mod
    use gflist_mod
    use gwlist_mod
    use walllist_mod
    use fluidlist_mod
    use stream_mod
    use listlimits_mod
    implicit none

    type, abstract :: part_t
        integer :: nslists,partid
        class(listlimits2d_t), pointer :: lislim
        class(mainlist_t),   pointer :: main       => null()
        class(sublistpt_t),  pointer :: slists(:)  => null()
        class(sublistpt_t),  pointer :: fluids(:)  => null()
        class(sublistpt_t),  pointer :: walls(:)   => null()
        class(sublist_t),    pointer :: influid    => null()
        class(sublist_t),    pointer :: bcfluid(:) => null()
        class(ghostlist_t),  pointer :: ghostf     => null()
        class(ghostlist_t),  pointer :: ghostw     => null()
    contains
        procedure :: init01             => init01_part
        generic   :: init               => init01
        procedure :: stream             => stream_part
        procedure :: collide            => collide_part
        procedure :: setfluidcol        => setfluidcol_part
        procedure :: setinfluidcol      => setinfluidcol_part
        procedure :: setbcfluidcol      => setbcfluidcol_part
        procedure :: checkcollisions    => checkcollisions_part
        procedure :: allocsublist       => allocatesublist_part
        procedure :: pointsublist       => pointsublist_part
        procedure :: allocmainlist      => allocmainlist_part
        procedure :: setlists           => setlists_part
        procedure :: copy             => copy_part

        !procedure :: assign             => assign_part
        !generic   :: assignment(=)      => assign
    end type

    type :: partpt_t
        class(part_t),pointer :: pt
    end type

    type,extends(part_t) :: part2d_t
    contains
        procedure :: init01             => init01_part2d
        procedure :: setlists           => setlists_part2d
        procedure :: allocmainlist      => allocmainlist_part2d
    end type

    type,extends(part2d_t) :: part3d_t
    contains
        procedure :: allocmainlist      => allocmainlist_part3d
    end type

contains


!==========================part procedures ==================================
subroutine copy_part(this,rhs)
    class(part_t)   :: this
    class(part_t),intent(in)    :: rhs

    call this%init(rhs%lislim,(/rhs%main%xlen,rhs%main%ylen,rhs%main%zlen/))
    call this%main%copy (   rhs%main)
    call this%ghostw%copy(  rhs%ghostw)
    call this%ghostf%copy(  rhs%ghostf)


end subroutine



subroutine init01_part(this,lislim,dims)
    class(part_t) :: this
    class(listlimits2d_t)               :: lislim
    integer, intent(in)                 :: dims(1:numd)
    stop "Error : init01_part is an abstract function."
end subroutine

subroutine setlists_part(this,lislim,dims)
    class(part_t) :: this
    class(listlimits2d_t),target               :: lislim
    integer, intent(in)                 :: dims(1:numd)
    stop "Error : setlists_part is an abstract function."
end subroutine

subroutine allocmainlist_part(this)
    class(part_t)   :: this
    stop "Error : allocmainlist_part is an abstract function."
end subroutine

subroutine stream_part(this)
    class(part_t) :: this
    integer :: n

    do n = 1, this%nslists
        call this%slists(n)%pt%stream1
    end do

    do n = 1, this%nslists
        call this%slists(n)%pt%stream2
    end do

end subroutine

subroutine collide_part(this)
    class(part_t) :: this
    integer :: n
    do n = 1, this%nslists
        call this%slists(n)%pt%collide
    end do
end subroutine


subroutine setfluidcol_part(this,col)
    class(part_t)       :: this
    class(collision_t)  :: col
    integer :: n
    do n = 0, nbcsit
        call this%fluids(n)%pt%alloccollision(source=col)
    end do
    ! Finalize the input col
end subroutine

subroutine setbcfluidcol_part(this,sit,col)
    class(part_t)       :: this
    class(collision_t)  :: col
    integer, intent(in) :: sit

    call this%fluids(sit)%pt%alloccollision(source=col)
    ! Finalize the input col
end subroutine

subroutine setinfluidcol_part(this,col)
    class(part_t)             :: this
    class(collision_t),target :: col

    call this%fluids(0)%pt%setcollision(col=col)
end subroutine


subroutine checkcollisions_part(this)
    class(part_t)   :: this
    integer         :: n

    do n = 1, this%nslists
        if (.not.associated(this%slists(n)%pt%col)) then
            print*,"Error: The collision for "//this%slists(n)%pt%name//' boundary list is not set.'
            stop
        endif
    end do
    Print*,'All collsions of sublists are assigned.'

end subroutine

! sl=sublist
subroutine allocatesublist_part(this,slid,source,slname)
    class(part_t)                   :: this
    integer,intent(in)              :: slid
    class(sublist_t),   intent(in)  :: source
    character(len=*)                :: slname
    integer                         :: l

    allocate(this%slists(slid)%pt,source=source)
    this%slists(slid)%pt%id = slid
    l = len(trim(slname))
    allocate(character(len=l)::this%slists(slid)%pt%name)
    this%slists(slid)%pt%name=trim(slname)
end subroutine

subroutine pointsublist_part(this,slid,source,slname)
    class(part_t)                           :: this
    integer,intent(in)                      :: slid
    class(sublist_t),   intent(in), target  :: source
    character(len=*)                        :: slname
    integer                                 :: l

    this%slists(slid)%pt => source
    this%slists(slid)%pt%id = slid
    l = len(trim(slname))
    allocate(character(len=l)::this%slists(slid)%pt%name)
    this%slists(slid)%pt%name=trim(slname)
end subroutine

!==========part2D procedures ===================

subroutine allocmainlist_part2d(this)
    class(part2d_t)   :: this
    allocate(mainlist2d_t::this%main)
end subroutine


! The order of assigning slist(:) must be in accordance with
! allocatepart_imgproc2d


subroutine init01_part2d(this,lislim,dims)
    class(part2d_t)                     :: this
    class(listlimits2d_t)               :: lislim
    integer, intent(in)                 :: dims(1:numd)

    call this%setlists(lislim,dims)

end subroutine


subroutine setlists_part2d(this,lislim,dims)
    class(part2d_t)                     :: this
    class(listlimits2d_t),target        :: lislim
    integer, intent(in)                 :: dims(1:numd)
    type(fluidlist_t)                   :: fluidlist
    type(walllist_t)                    :: walllist
    integer i,q,n


    this%nslists = nbcsit+4

    this%lislim => lislim

    allocate(this%slists(1:this%nslists))
    allocate(this%fluids(0:nbcsit))
    allocate(this%walls(1:1))
    allocate(gflist_t::this%ghostf)
    allocate(gwlist_t::this%ghostw)
    call this%allocmainlist()

    n = 1
    call this%allocsublist(slid=n,source=fluidlist,slname="influid")
    this%fluids(n-1)%pt => this%slists(n)%pt

    do q=1,nbcsit
        n = q + 1
        call this%allocsublist(slid=n,source=fluidlist,slname=bcsname(q)%mem)
        this%fluids(q)%pt => this%slists(n)%pt
    enddo

    n = n+1
    call this%allocsublist(slid=n,source=walllist,slname="wall")

    this%walls(1)%pt => this%slists(n)%pt

    n = n+1
    call this%pointsublist(slid=n,source=this%ghostw,slname="ghostw")

    n = n+1

    call this%pointsublist(slid=n,source=this%ghostf,slname="ghostf")

    ! All sublists have the address of each other
    do i=1,this%nslists
        this%slists(i)%pt%partslists=>this%slists
    enddo


    block
    integer,allocatable ::  iend(:)
    integer             ::  mxcelf

    allocate(iend(0:this%nslists))
    lislim%fluid = lislim%indomain - lislim%blackwall - lislim%wall
    mxcelf     = lislim%fluid
    n=0
    iend(n) = 0
    n=n+1
    iend(n) = lislim%influid
    do n=2,nbcsit+1
        iend(n) = iend(n-1) + lislim%bcfluid(n-1)
    end do
    iend(n) = iend(n-1)  + lislim%wall
    n=n+1
    iend(n) = iend(n-1) + lislim%ghostwall
    n=n+1
    iend(n) = iend(n-1) + lislim%ghostfluid

    call this%main%init(ist=0,ind=iend(this%nslists),xyzlen=dims)
    do n = 1, this%nslists
        call this%slists(n)%pt%init(mlist=this%main,ist=1,  &
                             istinmlist=iend(n-1)+1,     &
                             indinmlist=iend(n))
    enddo
    end block
end subroutine


!=============Part3d=================
subroutine allocmainlist_part3d(this)
    class(part3d_t)   :: this
    allocate(mainlist3d_t::this%main)
end subroutine


end module
