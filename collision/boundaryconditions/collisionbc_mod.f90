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
module collisionbc_mod
    use kinds_f
    use sublist_mod
    implicit none

!	collisionbc_t is the base collision type for boundary conditions.
!	The main responsibility is to find the missing pdfs of the boundary
! 	nodes using the updatepdf procedure. collisionbc_t contains a pointer,
! 	fcol, to the bulk fluid collision for running collision on boundary 
!	nodes and calculating densities and velocities. Therefore, collisionbc_t
! 	works well without the need for knowing what collision is defined for 
!	the fluid.

    type,   extends(collision_t),abstract    :: collisionbc_t
        class(collision_t), pointer  :: fcol
    contains
        procedure(updatepdf_collisionbc), deferred :: updatepdf
        procedure :: allocfluidcol => allocfluidcol_collisionbc
        procedure :: collide        => collide_collisionbc
        procedure :: getrho         => getrho_collisionbc
        procedure :: getvel         => getvel_collisionbc
        procedure :: getfeq         => getfeq_collisionbc
        procedure :: getaverho      => getaverho_collisionbc
        procedure :: cellcollide    => cellcollide_collisionbc
        procedure :: get1stmoment   => get1stmoment_collisionbc
        procedure :: getvel1stmom   => getvel1stmom_collisionbc
        procedure :: setsublist     => setsublist_collisionbc
        procedure :: updatepot      => updatepot_collisionbc
    end type

interface updatepdf_sub
subroutine updatepdf_collisionbc(this,icell)
    import                              :: collisionbc_t
    class(collisionbc_t)                :: this
    integer,                intent(in)  :: icell
end subroutine
end interface

contains

subroutine setsublist_collisionbc(this,sl)
    class(collisionbc_t)            :: this
    class(sublist_t),       target  :: sl
    this%sl=>sl
    this%fcol%sl=>sl
end subroutine

subroutine allocfluidcol_collisionbc(this,source)
    class(collisionbc_t),       target  :: this
    class(collision_t)                  :: source

    allocate(this%fcol,source=source)
    if (associated(this%sl)) this%fcol%sl => this%sl

end subroutine

subroutine collide_collisionbc(this)
    class(collisionbc_t)  :: this
    integer :: icell
    associate(sl=>this%sl)

    do icell = sl%istinmlist,sl%indinmlist
        call this%updatepdf(icell)
        call this%cellcollide(icell)
    end do

    end associate
end subroutine


subroutine updatepot_collisionbc(this)
    class(collisionbc_t)  :: this
    integer :: icell
    call this%fcol%updatepot()
end subroutine


subroutine cellcollide_collisionbc(this,icell)
    class(collisionbc_t)      :: this
    integer,    intent(in)  :: icell
    call this%fcol%cellcollide(icell)
end subroutine

pure function getfeq_collisionbc(this,rho,u) result(feq)
    class(collisionbc_t),   intent(in)  :: this
    real(kind=wp),          intent(in)  :: rho
    real(kind=wp),          intent(in)  :: u(1:numd)
    real(kind=wp)                       :: feq(0:nq_1)
    feq(0:nq_1) = this%fcol%getfeq(rho,u)
end function

pure function getvel_collisionbc(this,icell,rho) result(u)
    class(collisionbc_t),   intent(in)  :: this
    integer,                intent(in)  :: icell
    real(kind=wp),          intent(in)  :: rho
    real(kind=wp)                       :: u(1:numd)
    u(1:numd) = this%fcol%getvel(icell,rho)
end function

pure function getvel1stmom_collisionbc(this,icell,rho) result(u)
    class(collisionbc_t),intent(in) :: this
    integer,intent(in)          :: icell
    real(kind=wp),intent(in)    :: rho
    real(kind=wp)               :: u(1:numd)
    u(1:numd) = this%fcol%getvel1stmom(icell,rho)
end function

pure function get1stmoment_collisionbc(this,icell) result(fm)
    class(collisionbc_t),intent(in) :: this
    real(kind=wp)      :: fm(1:numd)
    integer,intent(in) :: icell
    fm(1:numd) = this%fcol%get1stmoment(icell)
end function

elemental function getaverho_collisionbc(this) result(averho)
    class(collisionbc_t),intent(in) :: this
    real(kind=wp) :: averho
    averho = this%fcol%getaverho()
end function

elemental function getrho_collisionbc(this,icell) result(rho)
    class(collisionbc_t),intent(in) :: this
    real(kind=wp) :: rho
    integer, intent(in) :: icell
    rho = this%fcol%getrho(icell)
end function


end module
