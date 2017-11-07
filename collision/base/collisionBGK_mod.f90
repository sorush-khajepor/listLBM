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

! This is the base type (class) for different BGK collisions.
module collisionBGK_mod
    use kinds_f
    use lattice_mod
    use sublist_mod, only : sublist_t,collision_t
    implicit none

!	Collision BGK, defined from collision_t type(calss) which is defined in
!	sublist_mod.f90.
    type, extends(collision_t) :: collisionBGK_t
        real(kind=wp) :: omega,tau
    contains
        procedure :: init01_BGK     => init01_BGK
        generic   :: init           => init01_BGK
!		Calculating and getting equilibrium distribution function
        procedure :: getfeq         => getfeq_BGK
!		Getting equilibrium distribution in a specific direction 
        procedure :: getfeqdir      => getfeqdir_BGK
!		Calculating and getting velocity of a node
        procedure :: getvel         => getvel_BGK
!		Calculating and getting velocity from the first moment
        procedure :: getvel1stmom   => getvel_BGK
!		Calculating and getting first moment (sum(fi.ei))
        procedure :: get1stmoment   => get1stmoment_BGK
!		Getting the density of a node
        procedure :: getrho         => getrho_BGK
!		Calculating and getting average density of this sublist
        procedure :: getaverho      => getaverho_BGK
!		Run collision on its sublist (each collision has a sublist.)
        procedure :: collide        => collide_bgk
!		Run collision on a node
        procedure :: cellcollide    => cellcollide_bgk
!		Updating (Shan-Chen) potential value at all nodes
        procedure :: updatepot      => updatepot_blank
!		Set a sublist for this collision
        procedure :: setsublist     => setsublist_BGK
    end type

!	A collision without any calculation
    type, extends(collisionBGK_t) :: blankcollision_t
    contains
        procedure :: collide => collide_blankcollision
    end type

!	A wall collision for pseudopotential LBM. A density is
!	defined for the wall. 
    type, extends(blankcollision_t) :: wallcollision_t
        real(kind=wp)   :: rho_w=0._wp
    contains
        procedure :: getrho         => getrho_wallcollision
        procedure :: getvel         => getvel_wallcollision
    end type

contains


subroutine setsublist_BGK(this,sl)
    class(collisionBGK_t)           :: this
    class(sublist_t),       target  :: sl

    this%sl => sl

end subroutine


subroutine collide_BGK(this)
    class(collisionBGK_t) :: this
    integer :: icell
    associate(sl=>this%sl)

    do icell = sl%istinmlist,sl%indinmlist
        call this%cellcollide(icell)
    end do

    end associate
end subroutine


subroutine updatepot_blank(this)
    class(collisionBGK_t) :: this

end subroutine

subroutine cellcollide_BGK(this,icell)
    class(collisionBGK_t)               :: this
    integer,                intent(in)  :: icell
    real(kind = wp)                     :: rho,u(1:numd),feq(0:nq_1)
    integer :: q
    associate(omega => this%omega,pdf=>this%sl%mlist%pdf)

    rho = this%getrho(icell)
    u   = this%getvel(icell,rho)
    feq = this%getfeq(rho,u)
    do q = 0, nq_1
         pdf(q,icell) = (1._wp-omega)* pdf(q,icell) + omega* feq(q)
    end do

    end associate
end subroutine


subroutine collide_blankcollision(this)
    class(blankcollision_t) :: this

end subroutine

subroutine init01_BGK(this,tau)
    class(collisionBGK_t) :: this
    real(kind=wp),intent(in) :: tau
    this%tau = tau
    this%omega = 1.0_wp/tau
end subroutine

! This functions of collision are using mainlist index
! otherwise an internal function should be called to
! convert sublist index to mainlist index which seems
! to be unecessarily time consuming.
elemental function getrho_BGK(this,icell) result(rho)
    class(collisionBGK_t),intent(in) :: this
    real(kind=wp) :: rho
    integer, intent(in) :: icell
    rho = sum(this%sl%mlist%pdf(0:nq_1,icell))
end function



! Calculate the average density of the sublist
elemental function getaverho_BGK(this) result(averho)
    class(collisionBGK_t),intent(in) :: this
    real(kind=wp) :: averho
    associate (ml=>this%sl%mlist,sl=>this%sl)
    averho = sum(ml%pdf(0:nq_1,sl%istinmlist:sl%indinmlist))
    averho = averho / sl%getnumocells()
    end associate
end function

! rho should not be zero, I don't use "if condition block"
! as subroutine is calling in them main loop.
pure function getvel_BGK(this,icell,rho) result(u)
    class(collisionBGK_t),intent(in) :: this
    real(kind=wp),intent(in) :: rho
    real(kind=wp)      :: u(1:numd)
    integer,intent(in) :: icell
    integer :: d,q
    associate (ml=>this%sl%mlist)
    do d = 1, numd
        u(d) =  sum(cqvec(1:nq_1,d)*ml%pdf(1:nq_1,icell)) / rho
    end do
    end associate
end function


pure function get1stmoment_BGK(this,icell) result(fm)
    class(collisionBGK_t),intent(in) :: this
    real(kind=wp)      :: fm(1:numd)
    integer,intent(in) :: icell
    integer :: d,q
    associate (ml=>this%sl%mlist)
    do d = 1, numd
        fm(d) = sum(cqvec(1:nq_1,d)*ml%pdf(1:nq_1,icell))
    end do
    end associate
end function

pure function getfeq_BGK(this,rho,u) result(feq)
    class(collisionBGK_t),intent(in) :: this
    real(kind=wp), intent(in) :: rho
    real(kind=wp), intent(in) :: u(1:numd)
    real(kind=wp) :: feq(0:nq_1)
    real(kind=wp) :: uq,uSqr
    integer :: q
    uSqr = sum(u*u)
    do q = 0, nq_1
        uq = sum(u*cqvec(q,1:numd) )
        feq(q) = weight(q) * rho * (1.0_wp + 3.0_wp * uq/cSqr + 4.5_wp * uq * uq/cqad    &
                                    - 1.5_wp * uSqr / cSqr)
    end do
end function

pure function getfeqdir_BGK(this,q,rho,u) result(feq)
    class(collisionBGK_t),intent(in) :: this
    real(kind=wp), intent(in) :: rho
    real(kind=wp), intent(in) :: u(1:numd)
    real(kind=wp) :: feq
    real(kind=wp) :: uq,uSqr
    integer,intent(in) :: q
    uSqr = sum(u*u)
    uq = sum(u*cqvec(q,1:numd) )
    feq = weight(q) * rho * (1.0_wp + 3.0_wp * uq/cSqr + 4.5_wp * uq * uq/cqad    &
                                    - 1.5_wp * uSqr / cSqr)
end function

elemental function getrho_wallcollision(this,icell) result(rho)
    class(wallcollision_t), intent(in)  :: this
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: rho

    rho = this%rho_w

end function

pure function getvel_wallcollision(this,icell,rho) result(u)
    class(wallcollision_t), intent(in)  :: this
    real(kind=wp),          intent(in)  :: rho
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: u(1:numd)

    u = 0._wp

end function



end module
