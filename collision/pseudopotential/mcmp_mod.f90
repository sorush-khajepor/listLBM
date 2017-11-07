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
module mcmp_mod
    use ppot_mod
    use partmcmp_mod
    implicit none

    type, extends(colBGKPP_t):: colBGKmcmp_t
        class(partmcmppt_t),    pointer :: parts(:)
!		Multi-component potential parameters
        class (potparam_t), pointer     :: mcpotparam(:) 
        integer,            allocatable :: othercomppotid(:,:),othercomp(:,:)
!		Number of parts or components
        integer                         :: nparts
    contains
       	procedure :: init_colBGKmcmp => init_colBGKmcmp
       	generic   :: init            => init_colBGKmcmp
       	procedure :: getforce        => getforce_colBGKmcmp
!		Get interaction force between components
       	procedure :: getmcmpforce    => getmcmpforce_colBGKmcmp
!		Update node potential
       	procedure :: updatenodepot   => updatenodepot_colBGKmcmp
!		update multi-component node potential
       	procedure :: updatemcmpnodepot => updatemcmpnodepot_colBGKmcmp
!		Get multi-component potential of a node
       	procedure :: getmcpot        => getmcpot_colBGKmcmp
!		Run collision on a node (or cell)
       	procedure :: cellcollide     => cellcollide_colBGkscmcmp
    end type

!   Wall collision for Multi-component model
    type, extends(colBGKmcmp_t):: colBGKmcmpwall_t
!		Density of the wall
        real(kind=wp) :: rho_wall = 0._wp
    contains
       procedure :: init_colbgkmcmpwall => init_colbgkmcmpwall
       generic   :: init                => init_colbgkmcmpwall
       procedure :: getrho              => getrho_colbgkmcmpwall
       procedure :: getvel              => getvel_colbgkmcmpwall
       procedure :: cellcollide         => cellcollide_colbgkmcmpwall
    end type

    type, extends(colBGKmcmp_t):: colBGKmcmpgf_t
    contains
       procedure :: updatemcmpnodepot   => updatemcmpnodepot_colBGKmcmpgf
       procedure :: cellcollide         => cellcollide_colBGKmcmpgf
       procedure :: getrho              => getrho_colBGKmcmpgf
       procedure :: getvel              => getvel_colBGKmcmpgf
       procedure :: getbcad             => getbodycelladdress_colBGKmcmpgf
       procedure :: get1stmoment        => get1stmoment_colBGKmcmpgf
    end type


contains

!===================== colBGKmcmpwall =============================
subroutine cellcollide_colbgkmcmpwall(this,icell)
    class(colBGKmcmpwall_t)               :: this
    integer,                intent(in)  :: icell
    ! Blank Collision
end subroutine

elemental function getrho_colbgkmcmpwall(this,icell) result(rho)
    class(colbgkmcmpwall_t),  intent(in)    :: this
    integer,                intent(in)      :: icell
    real(kind=wp)                           :: rho

    rho = this%rho_wall

end function

pure function getvel_colbgkmcmpwall(this,icell,rho) result(u)
    class(colbgkmcmpwall_t),  intent(in)    :: this
    real(kind=wp),          intent(in)      :: rho
    integer,                intent(in)      :: icell
    real(kind=wp)                           :: u(1:numd)

    u = 0._wp

end function


subroutine init_colBGKmcmpwall(this,mcpotparam,rho_wall)
    class(colbgkmcmpwall_t)                     :: this
    class(potparam_t),      target              :: mcpotparam(:)
    real(kind=wp),          intent(in)          :: rho_wall

    this%mcpotparam => mcpotparam
    this%rho_wall = rho_wall

end subroutine


!================== cell collide shan chen =====================================

! Shan-Chen multi-component force interaction
subroutine cellcollide_colBGkscmcmp(this,icell)
    class(colBGKmcmp_t) 	:: this
!	Index of the cell
    integer,    intent(in)  :: icell
    real(kind = wp)         :: rho,force(1:numd),feq(0:nq_1)
!	Velocity of the mixture
    real(kind = wp)         :: u(1:numd)
    integer                 :: q,d
!	Index of the other component
    integer                 :: other
!	Index of sublist of the other component
    integer                 :: islother
!	First moment of the other component
    real(kind=wp)           :: fmother(1:numd)
!	Density of the other component
    real(kind=wp)           :: rhother
    associate(omega => this%omega,pdf =>this%sl%mlist%pdf,slid=>this%sl%id)

    rho = this%getrho(icell)
    force = this%getforce(icell)
    if (this%sl%mlist%partid.eq.1) then
        other = 2
    else
        other = 1
    endif

    fmother = this%parts(other)%pt%slists(slid)%pt%col%get1stmoment(icell)
    rhother = this%parts(other)%pt%slists(slid)%pt%col%getrho(icell)
    u  = (this%get1stmoment(icell)+ fmother)/(rho+rhother) + force/rho

    if (rho+rhother.le.eps_wp) then
        u = 0._wp
    endif

    feq = this%getfeq(rho,u)

    do q = 0, nq_1
         pdf(q,icell) = (1._wp-omega) * pdf(q,icell) + omega * feq(q)
    end do

    end associate
end subroutine


!============================== colBGkmcmp ======================================

! Multi-component collision
subroutine cellcollide_colBGkmcmp(this,icell)
    class(colBGKmcmp_t) :: this
    integer,    intent(in)  :: icell
    real(kind = wp)         :: rho,u(1:numd),force(1:numd),feq(0:nq_1),source(0:nq_1)
    integer                 :: q
    integer                 :: other,d
    real(kind=wp)           :: fmother(1:numd),rhother
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

subroutine init_colBGKmcmp(this,tau,fcc,mcpotparam,parts,extforce)
    class(colBGKmcmp_t)                         :: this
    class(potparam_t),  target                  :: mcpotparam(:)
    real(kind=wp),      intent(in)              :: tau
    class(partmcmppt_t),    target              :: parts(:)
    real(kind=wp),      intent(in), optional    :: fcc,extforce(1:numd)
    integer                                     :: k,comp,otherid,other

    call this%init(tau=tau)
    allocate(this%extforce(1:numd))
    if (present(extforce)) then
        this%extforce = extforce
    else
        this%extforce = 0._wp
    end if
    allocate(this%weight(1:nq_1))
    allocate(this%fcqvec(0:nq_1,1:numd))
    this%mcpotparam => mcpotparam
    if (present(fcc)) then
        this%fcc = fcc
    else
        this%fcc = 1._wp
    end if
    this%fcqvec = this%fcc*qvec

    this%weight(1:nq_1) = weight(1:nq_1)*3._wp

    this%parts => parts
!	Number of parts
    this%nparts = ubound(parts,dim=1)-lbound(parts,dim=1)+1

    allocate(this%othercomppotid(1:this%nparts,1:this%nparts-1))
    allocate(this%othercomp(1:this%nparts,1:this%nparts-1))

    k = 0
    do comp =1, this%nparts
        do otherid = 1, this%nparts-1
            k = k + 1
            this%othercomppotid(comp,otherid) = (k-1)/this%nparts + 1
        end do
    enddo


!	Each component is interacting with n-1 other components.
!	We store these interactions in "othercomp" array.
    do comp =1, this%nparts
        otherid = 0
        do other = 1, this%nparts
            if (comp.eq.other) cycle
            otherid = otherid + 1
            this%othercomp(comp,otherid) = other
        end do
    enddo

end subroutine

! MCMP force + external force
pure function getforce_colBGKmcmp(this,icell) result(force)
    class(colBGKmcmp_t),  intent(in)  :: this
    integer,            intent(in)  :: icell
    real(kind=wp)                   :: force(1:numd)

    force = this%getmcmpforce(icell) + this%extforce*this%getrho(icell)

end function

pure function getmcmpforce_colBGKmcmp(this,icell) result(force)
    class(colBGKmcmp_t),  intent(in)  :: this
    integer,            intent(in)  :: icell
    real(kind=wp)                   :: force(1:numd)
    real(kind=wp)                   :: sumpot(1:numd),G
    integer                         :: d,q,comp,oc,ocpid,otherid
    associate (mcpot=>this%sl%mlist%mcpot,ml=>this%sl%mlist, &
               w=>this%weight,neilist=>this%sl%mlist%neilist, &
               fcqvec=>this%fcqvec,comp => this%sl%mlist%partid,ncomps=>this%nparts)

    force(1:numd) = 0._wp
    do otherid = 1, ncomps-1

        oc = this%othercomp(comp,otherid)
        ocpid = this%othercomppotid(comp,otherid)
        G=this%mcpotparam(otherid)%G

        sumpot = 0._wp
        do d = 1, numd
            do q = 1, nq_1
                sumpot(d) = sumpot(d) + w(q) * fcqvec(q,d) * this%parts(oc)%pt%main%mcpot(ocpid,neilist(q,icell))
            end do
        end do

        force(1:numd) = force(1:numd) + (-G) * mcpot(otherid,icell) * sumpot(1:numd)

    end do

    end associate

end function


subroutine updatenodepot_colBGKmcmp(this,icell)
    class(colBGKmcmp_t)             :: this
    integer,            intent(in)  :: icell

    call this%updatemcmpnodepot(icell)

end subroutine

subroutine updatemcmpnodepot_colBGKmcmp(this,icell)
    class(colBGKmcmp_t)             :: this
    integer,            intent(in)  :: icell
    integer                         :: otherid
    associate(ml=>this%sl%mlist,ncomps=>this%nparts)

    do otherid = 1, ncomps-1
        ml%mcpot(otherid,icell) = this%getmcpot(otherid,this%getrho(icell))
    end do

    end associate
end subroutine


function getmcpot_colBGKmcmp(this,otherid,rho)  result(mcpot)
    class(colBGKmcmp_t), intent(in)   :: this
    real(kind=wp)                   :: mcpot
    real(kind=wp),      intent(in)  :: rho
    integer,            intent(in)  :: otherid

    mcpot = this%mcpotparam(otherid)%getpot(rho)

end function


!====================== Ghost fluid collision =======================


subroutine updatemcmpnodepot_colBGKmcmpgf(this,icell)
    class(colbgkmcmpgf_t)               :: this
    integer,                intent(in)  :: icell
    integer :: ibcell,isl,otherid
    associate(ml=>this%sl%mlist,ncomps=>this%nparts)

        call this%getbcad(icell=icell,ibcell=ibcell,isl=isl)
        select type (col=> this%sl%partslists(isl)%pt%col)
        class is (colbgkmcmp_t)
            do otherid = 1, ncomps-1
                ml%mcpot(otherid,icell) = col%getmcpot(otherid,col%getrho(ibcell))
            end do
        class default
            stop'Error: ghost fluid node reperesents wrong node.'
        end select

    end associate
end subroutine


subroutine cellcollide_colBGKmcmpgf(this,icell)
    class(colBGKmcmpgf_t)               :: this
    integer,                intent(in)  :: icell
    ! Blank Collision
end subroutine


! icell (for ghost) and ibcell (for bodycell) are
! index of cell in the mainlist. isl= index of sublist.
elemental subroutine getbodycelladdress_colbgkmcmpgf(this,icell,ibcell,isl)
    class(colbgkmcmpgf_t),  intent(in)  :: this
    integer,                intent(out) :: ibcell,isl
    integer,                intent(in)  :: icell
    integer                             :: igcell

    select type (sl=>this%sl)
    class is (ghostlist_t)
        igcell = sl%getiscell(icell)
        ibcell = sl%potbodycell(igcell)
        isl = sl%mlist%islist(ibcell)
    end select

end subroutine

elemental function getrho_colBGKmcmpgf(this,icell) result(rho)
    class(colbgkmcmpgf_t),    intent(in)  :: this
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: rho
    integer                             :: isl,ibcell

    call this%getbcad(icell=icell,ibcell=ibcell,isl=isl)
    rho = this%sl%partslists(isl)%pt%col%getrho(ibcell)

end function

pure function getvel_colBGKmcmpgf(this,icell,rho) result(u)
    class(colbgkmcmpgf_t),  intent(in)  :: this
    real(kind=wp),          intent(in)  :: rho
    integer,                intent(in)  :: icell
    real(kind=wp)                       :: u(1:numd)
    integer                             :: isl,ibcell

    call this%getbcad(icell=icell,ibcell=ibcell,isl=isl)
    u = this%sl%partslists(isl)%pt%col%getvel(ibcell,rho)

end function

pure function get1stmoment_colBGKmcmpgf(this,icell) result(fm)
    class(colbgkmcmpgf_t),intent(in)    :: this
    integer,intent(in)                  :: icell
    real(kind=wp)                       :: fm(1:numd)
    integer                             :: isl,ibcell

    call this%getbcad(icell=icell,ibcell=ibcell,isl=isl)
    fm(1:numd) = this%sl%partslists(isl)%pt%col%get1stmoment(ibcell)

end function



end module
