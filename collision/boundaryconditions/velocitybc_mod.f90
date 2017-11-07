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
module velocitybc_mod
    use outflowbc_mod
    implicit none

!	A type for putting velocity condition on a boundary for D2Q9, D3Q15 and D3Q19. But D3Q27 is not supported!
!   For more detials please see papers below for more details:
!	(Zou, Q. and X. He, Physics of Fluids, 1997. 9(6): p. 1591-1598)
!	(Hecht, M., & Harting, J. Journal of Statistical Mechanics: Theory and Experiment, 2010(01), P01018)
    type,   extends(outflowbc_t)    ::  velocitybc_t
!		Velocity vectors tangent to the boundary
        integer,        allocatable ::  tandirs(:)
!		Velocity vectors crossing the boundary from inside toward outside
        integer,        allocatable ::  outletdirs(:)
!		An array needed for the boundary condition calculations (refer to papers above)
        integer,        allocatable ::  nxyz(:,:)
!		Tangent to boundary component of all velocity vectors (Q vectors)
        integer,        allocatable ::  qtan(:,:)
        real(kind=wp),  allocatable ::  u(:)
!		Constants depend on the lattice
        real(kind=wp)               ::  m1,m2,m3,m4
    contains
        procedure :: init01_velocitybc  => init01_velocitybc
        generic   :: init               => init01_velocitybc
!		Updating missing distribution functions of boundary nodes
        procedure :: updatepdf          => updatepdf_velocitybc
!		Calculating and getting boundary node density
        procedure :: getbcrho           => getbcrho_velocitybc
!		Setting Nxyz (see papers above)
        procedure :: setNxyz            => setNxyz_velocitybc
    end type

    contains

subroutine init01_velocitybc(this,u,normal,fcol)
    class(velocitybc_t)             :: this
    class(collision_t)              :: fcol
    integer,            intent(in)  :: normal(:)
    real(kind=wp),      intent(in)  :: u(:)

    allocate(this%u(1:numd))
    allocate(this%qtan(0:nq_1,1:numd))
    this%u = u
    this%qtan = 0
    call this%outflowbc_t%init(normal,fcol)
    call findtanplanedirs(normal=normal,dirs=this%tandirs)
    call findcrossplanedirs(normal=normal,dirs=this%outletdirs)

    call this%setNxyz()

    if (latname.eq.'D2Q9') then
        this%m1 = 2._wp
        this%m2 = 1._wp/6._wp
        this%m3 = 1._wp/3._wp
        this%m4 = 0.5_wp
    elseif (latname.eq.'D3Q19') then
        this%m1 = 1._wp
        this%m2 = 1._wp/6._wp
        this%m3 = 1._wp/3._wp
        this%m4 = 0.5_wp
    elseif (latname.eq.'D3Q15') then
        this%m1 = 2._wp
        this%m2 = 1._wp/12._wp
        this%m3 = 1._wp/6._wp
        this%m4 = 0.25_wp
    else
        stop 'Error in init01_velocitybc: lattice is not recognised.'
    end if

end subroutine

subroutine updatepdf_velocitybc(this,icell)
    class(velocitybc_t)                  :: this
    integer,                intent(in)  :: icell
    integer                             :: iq,j,q
    real(kind=wp)::rho,tmp
    associate(pdf=>this%sl%mlist%pdf)

    rho = this%getbcrho(icell)
    pdf(opposite(this%qnormal),icell) = pdf(this%qnormal,icell) - this%m1 * rho*sum(this%u*this%normal)/3._wp
    do q = 2, ubound(this%inletdirs,dim=1)
        iq = this%inletdirs(q)
        tmp = pdf(opposite(iq),icell) + this%m2 * rho * sum(qvec(iq,1:numd)*this%u(1:numd)) &
                                      + this%m3 * rho * sum(this%qtan(iq,1:numd)*this%u(1:numd))
        do j=1,ubound(this%nxyz,dim=1)
            tmp = tmp - this%m4 * pdf(abs(this%nxyz(j,q)),icell) * sign(1,this%nxyz(j,q))
        end do
        pdf(iq,icell) = tmp
    end do

    end associate
end subroutine

!	Initializing Nxyz vector
subroutine setNxyz_velocitybc(this)
    class(velocitybc_t)                 :: this
    integer                             :: q,val,n,niq,i,iq
    integer, allocatable                :: tmpnxyz(:,:)

    niq = ubound(this%inletdirs,dim=1)
    allocate(tmpnxyz(1:30,niq))

    do q = 0,nq_1
        this%qtan(q,1:numd) = qvec(q,1:numd) - this%normal(1:numd) *sum(qvec(q,1:numd)*this%normal(1:numd))
    end do

    do i = 1, niq
        iq = this%inletdirs(i)
        if (iq.eq.opposite(this%qnormal) ) then
            tmpnxyz(:,i) = 0
        else
        n = 0
        do q = 0,nq_1
            val = sum(this%qtan(iq,1:numd)*qvec(q,1:numd)) * (1-abs(sum(this%normal(1:numd)*qvec(q,1:numd))))
            if (val.ne.0) then
                n = n+1
                tmpnxyz(n,i) = q*val
            end if
        end do
        endif
    end do

    allocate(this%nxyz(1:n,1:niq))
    this%nxyz(1:n,1:niq) = tmpnxyz(1:n,1:niq)
    deallocate(tmpnxyz)

end subroutine

!	Calculating the density of the boundary node
function getbcrho_velocitybc(this,icell) result(rho)
    class(velocitybc_t)     :: this
    integer,    intent(in)  :: icell
    real(kind=wp)           :: rho
    integer                 :: q
    associate(pdf=>this%sl%mlist%pdf)

    rho = 0._wp
    do q=1,ubound(this%tandirs,dim=1)
        rho = rho + pdf(this%tandirs(q),icell)
    end do

    do q=1,ubound(this%outletdirs,dim=1)
        rho = rho + 2._wp*pdf(this%outletdirs(q),icell)
    end do

    rho = rho /(1._wp + sum(this%u(1:numd)*this%normal(1:numd)))

    end associate
end function

end module
