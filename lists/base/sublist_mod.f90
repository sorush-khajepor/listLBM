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
! Due to mutual dependency of sublist and collision they should be
! placed in the same module.
module sublist_mod
    use list_mod
    use mainlist_mod
    implicit none

    type,  extends(list_t) :: sublist_t
        !pointer to array of sublists (fluid,wall,ghostwall,...) in the same part
        class(sublistpt_t),      pointer :: partslists(:)   => null()
        class(mainlist_t),       pointer :: mlist           => null()
        class(collision_t),      pointer :: col             => null()
        class(stream_t),         pointer :: str             => null()
        integer            :: id=0
        !character (len=20) :: name
        character (len=:),allocatable :: name
!       starting ith-cell  in the main list, ending ith-cell in main list
        integer :: istinmlist, indinmlist
    contains
        procedure :: init01_sublist => init01_sublist
        generic   :: init           => init01_sublist
        procedure :: getijk     => getijk_sublist
        procedure :: setbounds  => setbounds_sublist
        procedure :: setmlist   => setmlist_sublist
        procedure :: getrcell   => getrcell_sublist
        procedure :: getimcell  => getimcell_sublist !  i cell in mainlist
        procedure :: getiscell  => getiscell_sublist !  i cell in sublist
        procedure :: getpdf     => getpdf_sublist
        procedure :: initpdf    => initpdf_sublist
        procedure :: getneilist => getneilist_sublist

        procedure :: alloccollision => alloccollision_sublist
        procedure :: setcollision   => setcollision_sublist
        procedure :: collide    => collide_sublist
        procedure :: updatepot  => updatepot_sublist
        procedure :: getrho     => getrho_sublist
        procedure :: getvel     => getvel_sublist
        procedure :: getfeq     => getfeq_sublist
        procedure :: getaverho  => getaverho_sublist
        procedure :: allocstream    => allocstream_sublist
        procedure :: stream1        => stream1_sublist
        procedure :: stream2        => stream2_sublist
    end type

    type sublistpt_t
        class(sublist_t),pointer :: pt => null()
    end type

    type, extends(sublist_t) :: sublist2d_t
    end type

    type, extends(sublist2d_t) :: sublist3d_t
    end type

!   getvel = actual velocity of fluid
!   getvel1stmoment = velocity obtained from first moment
!   For example in pseudo-potential model they are different from each other.
    type, abstract  :: collision_t
        class(sublist_t), pointer  :: sl
    contains
        procedure(collide_collision),       deferred  :: collide
        procedure(updatepot_collision),     deferred  :: updatepot
        procedure(cellcollide_collision),   deferred  :: cellcollide
        procedure(getfeq_collision),        deferred  :: getfeq
        procedure(getvel_collision),        deferred  :: getvel
        procedure(getvel_collision),        deferred  :: getvel1stmom
        procedure(get1stmoment_collision),  deferred  :: get1stmoment
        procedure(getrho_collision),        deferred  :: getrho
        procedure(getaverho_collision),     deferred  :: getaverho
        procedure(setsublist),              deferred  :: setsublist
    end type

    type, abstract :: stream_t
        class(sublist_t), pointer  :: sl
    contains
        procedure(opposite_stream),   deferred :: stream1
        procedure(reversepdf_stream), deferred :: stream2
    end type

    interface stream1_sub
    subroutine opposite_stream(this)
        import :: stream_t
        class(stream_t) :: this
    end subroutine
    end interface

    interface stream2_sub
    subroutine reversepdf_stream(this)
        import :: stream_t
        class(stream_t) :: this
    end subroutine
    end interface


    interface collision_subs
    subroutine collide_collision(this)
        import              :: collision_t
        class(collision_t)  :: this
    end subroutine



    subroutine cellcollide_collision(this,icell)
        import                  :: collision_t
        class(collision_t)      :: this
        integer,    intent(in)  :: icell
    end subroutine

    subroutine setsublist(this,sl)
        import                      :: collision_t, sublist_t
        class(collision_t)          :: this
        class(sublist_t),   target  :: sl
    end subroutine
    end interface


    interface collision_subs2
    subroutine updatepot_collision(this)
        import              :: collision_t
        class(collision_t)  :: this
    end subroutine
    end interface
    interface collision_funcs01
    pure function getfeq_collision(this,rho,u) result(feq)
        import :: collision_t,wp,numd,nq_1
        class(collision_t), intent(in)  :: this
        real(kind=wp),      intent(in)  :: rho
        real(kind=wp),      intent(in)  :: u(1:numd)
        real(kind=wp)                   :: feq(0:nq_1)
    end function

    pure function getvel_collision(this,icell,rho) result(u)
        import :: collision_t,wp,numd
        class(collision_t),intent(in) :: this
        integer,intent(in)          :: icell
        real(kind=wp),intent(in)    :: rho
        real(kind=wp)               :: u(1:numd)
    end function

    pure function get1stmoment_collision(this,icell) result(fm)
        import :: collision_t,wp,numd
        class(collision_t),intent(in) :: this
        real(kind=wp)      :: fm(1:numd)
        integer,intent(in) :: icell
    end function

    elemental function getaverho_collision(this) result(averho)
        import :: collision_t,wp
        class(collision_t),intent(in) :: this
        real(kind=wp) :: averho
    end function
    end interface

    interface collision_funcs02
    elemental function getrho_collision(this,icell) result(rho)
        import :: collision_t,wp
        class(collision_t),intent(in) :: this
        real(kind=wp) :: rho
        integer, intent(in) :: icell
    end function
    end interface

contains

subroutine init01_sublist(this,mlist,istinmlist,indinmlist,ist)
    class(sublist_t),target :: this
    class(mainlist_t), pointer, intent(in) :: mlist
    integer,intent(in),optional            :: ist
    integer,intent(in)                     :: istinmlist,indinmlist
    if (present(ist)) then
        call this%setbounds(istinmlist=istinmlist,indinmlist=indinmlist,ist=ist)
    else
        call this%setbounds(istinmlist=istinmlist,indinmlist=indinmlist)
    endif
    call this%setmlist(mlist)
end subroutine


subroutine setcollision_sublist(this,col)
    class(sublist_t),       target  :: this
    class(collision_t),     target  :: col
    this%col => col
    call this%col%setsublist(this)
end subroutine

subroutine alloccollision_sublist(this,source)
    class(sublist_t),       target  :: this
    class(collision_t)              :: source

    allocate(this%col,source=source)
    call this%col%setsublist(this)

end subroutine

subroutine allocstream_sublist(this,source)
    class(sublist_t),       target  :: this
    class(stream_t)                 :: source
    allocate(this%str,source=source)
    this%str%sl => this
end subroutine

function getpdf_sublist(this,iscell) result (pdf)
    class(sublist_t)        :: this
    integer, intent(in)     :: iscell
    real(kind=wp),pointer   :: pdf(:)
    integer                 :: imcell
    imcell = this%getimcell(iscell)
    pdf(0:nq_1) => this%mlist%pdf(0:nq_1,imcell)
end function

function getneilist_sublist(this,iscell) result(nei)
    class(sublist_t)    :: this
    integer, intent(in) :: iscell
    integer, pointer    :: nei(:)
    integer             :: imcell
    imcell = this%getimcell(iscell)
    nei(1:nq_1) => this%mlist%neilist(1:nq_1,imcell)
end function

! istinmlist = the first index of sublist inside mainlist
! indinmlist = the last index of sublist inside mainlist
! ist = the first index of the sublist (usually ist = 1)
subroutine setbounds_sublist(this,istinmlist,indinmlist,ist)
    class(sublist_t) :: this
    integer,intent(in) :: istinmlist,indinmlist
    integer,intent(in), optional :: ist
    this%istinmlist = istinmlist
    this%indinmlist = indinmlist
    if (present(ist)) then
        this%ist = ist
    else
        this%ist = 1
    endif
    this%ind = indinmlist - istinmlist + this%ist
end subroutine

function getijk_sublist(this,icell) result(ijk)
    class (sublist_t) :: this
    integer :: ijk(1:numd),imcell
    integer,intent(in) :: icell
    imcell = this%getimcell(icell)
    ijk = this%mlist%getijk(imcell)
end function

function getimcell_sublist(this,iscell)
    class (sublist_t) :: this
    integer :: getimcell_sublist
    integer,intent(in) :: iscell
    getimcell_sublist = this%istinmlist + iscell - this%ist
end function

elemental function getiscell_sublist(this,imcell)
    class (sublist_t), intent(in) :: this
    integer :: getiscell_sublist
    integer,intent(in) :: imcell
    getiscell_sublist = -this%istinmlist + imcell + this%ist
end function


function getrcell_sublist(this,iscell)
    class(sublist_t) :: this
    integer,intent(in) :: iscell
    integer :: imcell, getrcell_sublist
    imcell = this%getimcell(iscell)
    getrcell_sublist = this%mlist%rcell(imcell)
end function

subroutine setmlist_sublist(this,mlist)
    class(sublist_t) :: this
    class (mainlist_t), target :: mlist
    this%mlist => mlist
end subroutine

function getfeq_sublist(this,rho,u) result(feq)
    class(sublist_t),intent(in) :: this
    real(kind=wp),   intent(in) :: rho
    real(kind=wp)               :: u(1:numd)
    real(kind=wp)               :: feq(0:nq_1)

    feq(0:nq_1) = this%col%getfeq(rho,u)

end function

function getvel_sublist(this,iscell) result(u)
    class(sublist_t),intent(in) :: this
    integer,intent(in)          :: iscell
!    real(kind=wp),intent(in)    :: rho
    real(kind=wp)               :: u(1:numd)
    integer                     :: imcell
    imcell = this%getimcell(iscell)
    u(1:numd) = this%col%getvel(icell=imcell,rho=this%col%getrho(imcell))

end function

function getaverho_sublist(this) result(averho)
    class(sublist_t),intent(in) :: this
    real(kind=wp) :: averho

    averho = this%col%getaverho()

end function

function getrho_sublist(this,iscell) result(rho)
    class(sublist_t),intent(in) :: this
    real(kind=wp)       :: rho
    integer, intent(in) :: iscell
    integer             :: imcell
    imcell = this%getimcell(iscell)
    rho = this%col%getrho(imcell)
end function

subroutine collide_sublist(this)
    class(sublist_t) :: this

    call this%col%collide

end subroutine

subroutine updatepot_sublist(this)
    class(sublist_t) :: this

    call this%col%updatepot

end subroutine

subroutine initpdf_sublist(this,iscell,rho,u)
    class(sublist_t) :: this
    real(kind=wp),  intent(in) :: rho,u(1:numd)
    integer,        intent(in) :: iscell
    integer                    :: imcell
    imcell = this%getimcell(iscell)
    this%mlist%pdf(0:nq_1,imcell) = this%col%getfeq(rho,u)
end subroutine


subroutine stream1_sublist(this)
    class(sublist_t) :: this

    call this%str%stream1

end subroutine

subroutine stream2_sublist(this)
    class(sublist_t) :: this

    call this%str%stream2

end subroutine



end module
