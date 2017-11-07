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
! Note: ppot = pseudo-potential

module ppot_mod
    use potentialfunc_mod
    use colBGKextForce_mod
    use ghostlist_mod
    implicit none

!   BGKPP = BGK Pseudo-Potential
    type, extends(colBGKextForce_t):: colBGKPP_t
!		Parameters of the potential		
        class (potparam_t), pointer     :: potparam
!		Weights for the interaction force
        real(kind=wp),      allocatable :: weight(:)
!		fcqvec = fcc*lattice velocity vectors
        real(kind=wp),      allocatable :: fcqvec(:,:)
!   	fcc = lattice spacing for force interaction 
!	   	(fcc should be replaced by c)
        real(kind=wp)                   :: fcc
    contains
       	procedure :: init_colBGKPP   => init_colBGKPP
       	generic   :: init            => init_colBGKPP
!		Calculating and getting the potential of a node
       	procedure :: getpot          => getpot_colBGKPP
!		Getting the force imposed on a node
       	procedure :: getforce        => getforce_colBGKPP
!		Calculating single-component-multiphase force
       	procedure :: getSCMPforce    => getSCMPforce_colBGKPP
!		Calculating and updating a node potential
       	procedure :: updatenodepot   => updatenodepot_colBGKPP
!		Updating the potential of nodes in a list
       	procedure :: updateslpot     => updateslpot_colBGKPP
!		Updating the potentials in this collision/sublist
       	procedure :: updatepot       => updateslpot_colBGKPP
    end type

! 	Wall collision for pseudo-potential model
    type, extends(colBGKPP_t):: colBGKPPwall_t
!		Density is defined for the wall.		
        real(kind=wp) :: rho_wall = 0._wp
    contains
       	procedure :: init_colbgkppwall   => init_colbgkppwall
       	generic   :: init                => init_colbgkppwall
!		Get density of the wall
       	procedure :: getrho              => getrho_colbgkppwall
!		Get velocity of the wall
       	procedure :: getvel              => getvel_colbgkppwall
!		Run collision on a cell/node
       	procedure :: cellcollide         => cellcollide_colbgkppwall
    end type

! 	Ghost fluid collision for pseudo-potential model
!	Ghost fluid = ghost nodes which are placed out of the
!	domain boundary and resemble a fluid node
! 	It doesn't need initialization as it uses other lists'
! 	collisions.
    type, extends(colBGKPPwall_t):: colBGKPPgf_t
    contains
       procedure :: updatenodepot       => updatenodepot_colbgkppgf
       procedure :: getrho              => getrho_colBGKPPgf
       procedure :: getvel              => getvel_colBGKPPgf
    end type

!========================== Procedures ==========================

contains


subroutine init_colBGKPP(this,tau,fcc,potparam)
    class(colBGKPP_t)                           :: this
    class(potparam_t),  target                  :: potparam
    real(kind=wp),      intent(in)              :: tau
    real(kind=wp),      intent(in), optional    :: fcc


!	Initializing relaxation time
    call this%init(tau=tau)

    allocate(this%weight(1:nq_1))
    allocate(this%fcqvec(0:nq_1,1:numd))

    this%potparam => potparam

!	fcc=c=1.0 mostly in LBM
    if (present(fcc)) then
        this%fcc = fcc
    else
        this%fcc = 1._wp
    end if

    this%fcqvec = this%fcc*qvec

!	Weight of the interaction force is 3 times weight of 
!	LB equation. LB weight is defined in lattice_mod.f90.
    this%weight(1:nq_1) = weight(1:nq_1)*3._wp

end subroutine

pure function getforce_colBGKPP(this,icell) result(force)
    class(colBGKPP_t),  intent(in)  :: this
    integer,            intent(in)  :: icell
    real(kind=wp)                   :: force(1:numd)

    force = this%getSCMPforce(icell)

end function

! Calculating Shan-Chen force
pure function getSCMPforce_colBGKPP(this,icell) result(force)
    class(colBGKPP_t),  intent(in)  :: this
    integer,            intent(in)  :: icell
    real(kind=wp)                   :: force(1:numd)
    real(kind=wp)                   :: sumpot(1:numd)
    integer                         :: d,q
    associate (pot=>this%sl%mlist%pot,ml=>this%sl%mlist, &
               w=>this%weight, G=>this%potparam%G,neilist=>this%sl%mlist%neilist, &
               fcqvec=>this%fcqvec)

    sumpot = 0._wp
    do d = 1, numd
        do q = 1, nq_1
            sumpot(d) = sumpot(d) + w(q) * fcqvec(q,d) * pot(neilist(q,icell))
        end do
    end do

    force(1:numd) = (-G) * pot(icell) * sumpot(1:numd)

    end associate
end function

function getpot_colBGKPP(this,rho)  result(pot)
    class(colBGKPP_t), intent(in)               :: this
    real(kind=wp)                   :: pot
    real(kind=wp),      intent(in)  :: rho

    pot = this%potparam%getpot(rho)

end function

subroutine updateslpot_colBGKPP(this)
    class(colBGKPP_t)   :: this
    integer             :: icell
    associate(sl=>this%sl)

    do icell = sl%istinmlist,sl%indinmlist
        call this%updatenodepot(icell)
    end do

    end associate
end subroutine

subroutine updatenodepot_colBGKPP(this,icell)
    class(colBGKPP_t)               :: this
    integer,            intent(in)  :: icell
    associate(ml=>this%sl%mlist)

    ml%pot(icell) = this%getpot(this%getrho(icell))

    end associate
end subroutine
!============== colBGKPPwall ==================
subroutine init_colBGKPPwall(this,potparam,rho_wall)
    class(colBGKPPwall_t)                       :: this
    class(potparam_t),      target              :: potparam
    real(kind=wp),          intent(in)          :: rho_wall

    this%potparam => potparam
    this%rho_wall = rho_wall

end subroutine

elemental function getrho_colBGKppwall(this,icell) result(rho)
    class(colBGKppwall_t),  intent(in)  :: this
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: rho

    rho = this%rho_wall

end function

!	Velocity of the wall is set to zero.
pure function getvel_colBGKppwall(this,icell,rho) result(u)
    class(colBGKppwall_t),  intent(in)  :: this
    real(kind=wp),          intent(in)  :: rho
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: u(1:numd)

    u = 0._wp

end function

subroutine cellcollide_colBGKppwall(this,icell)
    class(colBGKppwall_t)               :: this
    integer,                intent(in)  :: icell
    ! Blank Collision
end subroutine


!==================== ghost fluid collision procedures ====================


subroutine init_colBGKPPgf(this,potparam )
    class(colBGKPPgf_t)                     :: this
    class(potparam_t),      target          :: potparam

    this%potparam => potparam

end subroutine


! Get density of the ghost node:
! Firstly, finding the body cell (the actual node the ghost is 
! representing of), then getting the density of that.
! ghost cell => body cell
! in periodic boundary condition ghost cell is representing 
! the node from the other side of the domain.
elemental function getrho_colBGKPPgf(this,icell) result(rho)
    class(colBGKPPgf_t),    intent(in)  :: this
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: rho
    integer                             :: isl,igcell,ibcell
    select type (sl=>this%sl)
    class is (ghostlist_t)
!		index of the ghost cell
        igcell = sl%getiscell(icell)
!		index of the body cell (actual cell)
        ibcell = sl%potbodycell(igcell)
!		index of the sublist of the body cell
        isl = sl%mlist%islist(ibcell)
!		Calling the collision of the body cell to get density
        rho = sl%partslists(isl)%pt%col%getrho(ibcell)
    end select
end function



pure function getvel_colBGKPPgf(this,icell,rho) result(u)
    class(colBGKPPgf_t),  intent(in)  :: this
    real(kind=wp),          intent(in)  :: rho
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: u(1:numd)
    integer                             :: isl,igcell,ibcell

    select type (sl=>this%sl)
    class is (ghostlist_t)
!		index of the ghost cell
        igcell = sl%getiscell(icell)
!		index of the body cell (actual cell)
        ibcell = sl%potbodycell(igcell)
!		index of the sublist of the body cell
        isl = sl%mlist%islist(ibcell)
!		Calling the collision of the body cell to get velocity
        u = sl%partslists(isl)%pt%col%getvel(ibcell,rho)
    end select

end function

subroutine updatenodepot_colbgkppgf(this,icell)
    class(colbgkppgf_t)               :: this
    integer,                intent(in)  :: icell
    integer :: ibcell,igcell,isl
    associate(ml=>this%sl%mlist)
    select type (sl=>this%sl)
    class is (ghostlist_t)
        igcell = sl%getiscell(icell)
        ibcell = sl%potbodycell(igcell)
        isl = sl%mlist%islist(ibcell)
        select type (col=> sl%partslists(isl)%pt%col)
        class is (colBGKpp_t)
            ml%pot(icell) = col%getpot(col%getrho(ibcell))
        class default
            stop'Error: ghost fluid node reperesents wrong node.'
        end select
    end select
    end associate

end subroutine


end module
