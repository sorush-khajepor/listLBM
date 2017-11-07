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
module mainlist_mod
    use kinds_f
    use list_mod
    implicit none

    type, abstract, extends(list_t) :: mainlist_t
!       Domain lengths
        integer :: xlen=-1, ylen=-1, zlen=-1
!		ID of the current part
        integer :: partid=0
!		Particle distribution function
        real(kind=wp), pointer :: pdf(:,:)      => null()
!		Neighbour list of a node
        integer,       pointer :: neilist(:,:)  => null()
! 		Real cell position in a 3D matrix
        integer,       pointer :: rcell(:)      => null() 
!		Index of the sublist
        integer,       pointer :: islist(:)
! 		(pseudo)potential for internal interactions of a component
		real(kind=wp), pointer :: pot(:)        => null() 
! 		Multi-component (pseudo)potential
! 		Multicomponent potential has values for interactions of the
! 		focused compoent with the others. Therefore, different potentials
! 		can be defined for different interactions.        
        real(kind=wp), pointer :: mcpot(:,:)    => null() 
    contains
        procedure :: init01_mainlist    =>  init01_mainlist
        generic   :: init               =>  init01_mainlist
        procedure :: setPDF             =>  setpdf_mainlist
        procedure :: setneilist         =>  setneilist_mainlist
        procedure :: setrcell           =>  setrcell_mainlist
        procedure :: setdims            =>  setdimensions_mainlist
        procedure :: allocpot           =>  allocpot_mainlist
        procedure :: allocmcpot         =>  allocmcpot_mainlist
        procedure :: copy               =>  copy_mainlist
     end type

    type,extends(mainlist_t) :: mainlist2d_t
    contains
!		Set dimensions
        procedure :: setdims => setdimensions_mainlist2d
!		Set real position of a cell
        procedure :: setrcell =>   setrcell_mainlist2d
!		Get x,y,z index (location) of a node
        procedure :: getijk     =>       getijk_mainlist2d
    end type

    type,extends(mainlist2d_t) :: mainlist3d_t
    contains
        procedure :: setdims  => setdimensions_mainlist3d
        procedure :: setrcell => setrcell_mainlist3d
        procedure :: getijk   => getijk_mainlist3d
    end type

contains

!==========mainlist procedures=======

! Copy current mainlist on another one 
subroutine copy_mainlist(this,rhs)
    class(mainlist_t) :: this
    class(mainlist_t), intent(in)  :: rhs

    this%ist=rhs%ist
    this%ind=rhs%ind
    if (associated(this%pdf))       deallocate(this%pdf)
    if (associated(this%neilist))   deallocate(this%neilist)
    if (associated(this%rcell))     deallocate(this%rcell)
    if (associated(this%islist))    deallocate(this%islist)
    allocate (this%pdf(0:nq_1,this%ist:this%ind))
    this%neilist =>  rhs%neilist
    this%rcell   =>  rhs%rcell
    this%islist  =>  rhs%islist
    this%pdf = 0._wp
    call this%setdims(xyzlen=(/rhs%xlen,rhs%ylen,rhs%zlen/))

end subroutine

subroutine setdimensions_mainlist(this,xyzlen)
    class(mainlist_t) :: this
    integer,intent(in):: xyzlen(:)
    stop "setdimensions_mainlist is an abstract subroutine."
end subroutine

subroutine init01_mainlist(this,ist,ind,xyzlen)
    class(mainlist_t)  :: this
    integer,intent(in) :: ist, ind
    integer,intent(in) :: xyzlen(:)
    this%ist=ist
    this%ind=ind
    allocate (this%pdf(0:nq_1,ist:ind))
    allocate (this%neilist(1:nq_1,ist:ind))
    allocate (this%rcell(ist:ind))
    allocate (this%islist(ist:ind))
    this%pdf = 0._wp
    this%neilist = 0
    this%rcell = 0
    this%islist = 0
    call this%setdims(xyzlen)
end subroutine


subroutine allocpot_mainlist(this)
    class(mainlist_t) :: this

    allocate(this%pot(this%ist:this%ind))
    this%pot = 0._wp

end subroutine


subroutine allocmcpot_mainlist(this,ncomps)
    class(mainlist_t)   :: this
    integer, intent(in) :: ncomps
    if (associated(this%mcpot)) &
        stop 'Error: in allocmcpot_mainlist mcpot is associated!'
    allocate(this%mcpot(1:(ncomps-1),this%ist:this%ind))
    this%mcpot = 0._wp

end subroutine


subroutine setpdf_mainlist(this,pdf)
    class(mainlist_t) :: this
    real(kind=wp),pointer,intent(in) :: pdf(:,:)
    this%pdf => pdf
end subroutine

subroutine setneilist_mainlist(this,neilist)
    class(mainlist_t) :: this
    integer, pointer, intent(in) :: neilist(:,:)
    this%neilist => neilist
end subroutine

subroutine setrcell_mainlist(this,icell,ijk)
    class(mainlist_t) :: this
    integer,intent(in) :: icell
    integer, intent (in) :: ijk(:)
    stop "Error : setrcell_mainlist is an abstract function."
end subroutine

!======mainlist2D procedures===========

subroutine setdimensions_mainlist2d(this,xyzlen)
    class(mainlist2d_t) :: this
    integer,intent(in):: xyzlen(:)
    this%xlen = xyzlen(1)
    this%ylen = xyzlen(2)
    this%zlen = 1
end subroutine

function getijk_mainlist2d(this,icell) result(ijk)
    class (mainlist2d_t) :: this
    integer,intent(in) :: icell
    integer :: ijk(1:numd)
    ijk(1) = mod(this%rcell(icell),this%xlen+2)
    ijk(2) = this%rcell(icell)/(this%xlen+2)
end function

subroutine setrcell_mainlist2d(this,icell,ijk)
    class(mainlist2d_t) :: this
    integer,intent(in) :: icell
    integer, intent (in) :: ijk(:)
    this%rcell(icell) =  ijk(2)*(this%xlen+2) + ijk(1)
end subroutine

!======mainlist3D procedures===========

subroutine setdimensions_mainlist3d(this,xyzlen)
    class(mainlist3d_t) :: this
    integer,intent(in):: xyzlen(:)
    this%xlen = xyzlen(1)
    this%ylen = xyzlen(2)
    this%zlen = xyzlen(3)
end subroutine

function getijk_mainlist3d(this,icell) result(ijk)
    class (mainlist3d_t) :: this
    integer,intent(in) :: icell
    integer :: ijk(1:numd),r
    r = mod(this%rcell(icell),((this%xlen+2)*(this%ylen+2)))
    ijk(1) = mod(r,this%xlen+2)
    ijk(2) = r/(this%xlen+2)
    ijk(3) = this%rcell(icell)/((this%xlen+2)*(this%ylen+2))
end function


subroutine setrcell_mainlist3d(this,icell,ijk)
    class(mainlist3d_t) :: this
    integer,intent(in) :: icell
    integer, intent (in) :: ijk(:)
    this%rcell(icell) = ijk(3) *(this%xlen+2)*(this%ylen+2) + ijk(2)*(this%xlen+2) + ijk(1)
end subroutine



end module
