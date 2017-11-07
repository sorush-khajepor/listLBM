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

module colBGKextForce_mod
    use lattice_mod
    use collisionBGK_mod
    implicit none

! 	A collision type for imposing external field force such as gravity.
    type, extends(collisionBGK_t)   :: colBGKextForce_t
!	External force allocatable in 2D and 3D
    real(kind=wp), allocatable  :: extforce(:)
    contains
        procedure :: init01_colBGKextForce      => init01_colBGKextForce
        generic   :: init                       => init01_colBGKextForce
!		Getting velocity of a node
        procedure :: getvel         => getvel_colBGKextForce
!		Getting velocity shifted by the force
        procedure :: getvelshifted  => getvelshifted_colBGKextForce
!		Getting force imposed on a node
        procedure :: getforce       => getforce_colBGKextForce
!		Doing collision step on a node
        procedure :: cellcollide    => cellcollide_colBGKextForce
!		Calculating and getting the source term due to external force
!	 	should be added to LB equation with Guo et al. method (PRE,2002.65(4))
        procedure :: getsource      => getsource_colBGKextForce
    end type

contains

subroutine init01_colBGKextForce(this,tau,extforce)
    class(colBGKextForce_t)             :: this
    real(kind=wp),          intent(in)  :: extforce(1:numd),tau

    call this%init(tau)
    allocate(this%extforce(1:numd))
    this%extforce = extforce

end subroutine

subroutine cellcollide_colBGKextForce(this,icell)
    class(colBGKextForce_t) :: this
    integer,    intent(in)  :: icell
    real(kind = wp)         :: rho,u(1:numd),force(1:numd),feq(0:nq_1),source(0:nq_1)
    integer                 :: q
    associate(omega => this%omega,pdf =>this%sl%mlist%pdf)

    rho = this%getrho(icell)
    force = this%getforce(icell)
    u   = this%getvelshifted(icell,rho,force)
    feq = this%getfeq(rho,u)
    source = this%getsource(u,force)
    do q = 0, nq_1
         pdf(q,icell) = (1._wp-omega)* pdf(q,icell) + omega * feq(q) + source(q)
    end do

    end associate
end subroutine


! Calculating actual velocity of a node under an external force
pure function getvel_colBGKextForce(this,icell,rho) result(u)
    class(colBGKextForce_t) ,intent(in) :: this
    real(kind=wp)           ,intent(in) :: rho
    integer                 ,intent(in) :: icell
    real(kind=wp)                       :: u(1:numd)

    u(1:numd) = this%getvel1stmom(icell,rho) + this%getforce(icell) / (2._wp*rho)

end function

! Shifted velocity due to external force
pure function getvelshifted_colBGKextForce(this,icell,rho,force) result(u)
    class(colBGKextForce_t) ,intent(in) :: this
    real(kind=wp)           ,intent(in) :: rho,force(1:numd)
    integer                 ,intent(in) :: icell
    real(kind=wp)                       :: u(1:numd)

    u(1:numd) = this%getvel1stmom(icell,rho) + force(1:numd) / (2._wp*rho)

end function


pure function getforce_colBGKextForce(this,icell) result(f)
    class(colBGKextForce_t) ,intent(in) :: this
    integer                 ,intent(in) :: icell
    real(kind=wp)                       :: f(1:numd)

    f = this%extforce(1:numd)

end function

! Adding external force with the method of Guo et al. (Physical Review E, 2002. 65(4): p. 046308)
pure function getsource_colBGKextForce (this,u,force) result(source)
    class(colBGKextForce_t),intent(in)  :: this
    real(kind=wp),          intent(in)  :: force(1:numd),u(1:numd)
    real(kind=wp)                       :: source(0:nq_1)
    integer                             :: q

    do q = 0, nq_1
         source(q) =    (1._wp-0.5*this%omega)*weight(q)* (   &
                        sum(force(1:numd)*(cqvec(q,1:numd)-u(1:numd)))/(csqr/3._wp) + &
                        sum(force(1:numd)*cqvec(q,1:numd))* sum(u(1:numd)*cqvec(q,1:numd))/(cqad/9._wp) )
     enddo
end function

end module
