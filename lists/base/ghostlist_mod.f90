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
module ghostlist_mod
    use kinds_f
    use sublist_mod
    implicit none

    type :: outletdir_t
        integer,allocatable :: dir(:)
    end type
!====================== ghost node =====================
!
!           * (away)
! boundary
!    cells  * (closest2gcell)          * gcell          * strbodycell
!
!           * (away)
!====================== ghost node =====================


    type,extends(sublist_t) :: ghostlist_t
        integer,            pointer :: strbodycell(:)       => null()
        integer,            pointer :: potbodycell(:)       => null()
        integer,            pointer :: closest2gcell(:)     => null()
        integer,            pointer :: linkPE(:)            => null()
        type(outletdir_t),  pointer :: outlet(:)            => null()
    contains
        procedure :: init01_sublist     =>  init01_ghostlist
        procedure :: setstrbodycell     =>  setstrbodycell_ghostlist
        procedure :: copy             =>  copy_ghostlist
        !procedure :: assign             =>  assign_ghostlist
        !generic   :: assignment(=)      =>  assign
    end type

    type, extends (ghostlist_t) :: ghostlist2d_t
    end type

    type, extends (ghostlist2d_t) :: ghostlist3d_t
    end type

contains


!==========ghostlist procedures=======

subroutine copy_ghostlist(this,rhs)
    class(ghostlist_t) :: this
    class(ghostlist_t), intent(in ) :: rhs

   ! this%ist = rhs%ist
   ! this%ind = rhs%ind
    this%istinmlist = rhs%istinmlist
    this%indinmlist = rhs%indinmlist

    this%id = rhs%id
    this%name = rhs%name

    this%strbodycell => rhs%strbodycell
    this%potbodycell => rhs%potbodycell
    this%closest2gcell => rhs%closest2gcell
    this%linkPE => rhs%linkPE
    this%outlet => rhs%outlet

end subroutine


subroutine init01_ghostlist(this,mlist,istinmlist,indinmlist,ist)
    class(ghostlist_t), target                  :: this
    class(mainlist_t),  pointer,    intent(in)  :: mlist
    integer,intent(in), optional                :: ist
    integer,intent(in)                          :: istinmlist,indinmlist

    if (present(ist)) then
        call this%setbounds(istinmlist=istinmlist,indinmlist=indinmlist,ist=ist)
    else
        call this%setbounds(istinmlist=istinmlist,indinmlist=indinmlist)
    endif
    call this%setmlist(mlist)
    allocate(this%strbodycell(this%ist:this%ind),this%closest2gcell(this%ist:this%ind), &
        this%linkPE(this%ist:this%ind), &
        this%outlet(this%ist:this%ind))
    this%strbodycell = 0
    this%closest2gcell = 0
    this%linkPE = 0

!   potbody cell points to strbodycell or closest2gcell
end subroutine


subroutine setstrbodycell_ghostlist(this,strbodycell)
    class(ghostlist_t) :: this
    integer,pointer,intent(in) :: strbodycell(:)
    this%strbodycell => strbodycell
end subroutine


end module
