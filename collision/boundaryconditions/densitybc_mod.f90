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

module densitybc_mod
    use velocitybc_mod
    implicit none

!	This type puts a density constraint on the selected boundary.
    type,   extends(velocitybc_t)    ::  densitybc_t
        real(kind=wp)   ::  rho,signu
        integer         ::  iu
    contains
        procedure :: init01_densitybc  => init01_densitybc
        generic   :: init               => init01_densitybc
        procedure :: updatepdf          => updatepdf_densitybc
        procedure :: getbcvel           => getbcvel_densitybc
    end type

    contains

subroutine init01_densitybc(this,rho,normal,fcol)
    class(densitybc_t)              :: this
    class(collision_t)              :: fcol
    integer,            intent(in)  :: normal(:)
    real(kind=wp),      intent(in)  :: rho
    real(kind=wp)                   :: u(1:numd)
    integer                         :: tmp

    this%rho = rho
    tmp = sum((/1,2,3/)*normal)
    this%signu = sign(1._wp,real(tmp,kind=wp))
    this%iu = abs(tmp)

    u = 0._wp
    call this%velocitybc_t%init(u=u,normal=normal,fcol=fcol)

end subroutine

subroutine updatepdf_densitybc(this,icell)
    class(densitybc_t)                  :: this
    integer,                intent(in)  :: icell
    integer                             :: iq,j,q
    real(kind=wp)::u(1:numd),tmp
    associate(pdf=>this%sl%mlist%pdf)

    u = this%getbcvel(icell)
    pdf(opposite(this%qnormal),icell) = pdf(this%qnormal,icell) - this%m1 * this%rho*sum(u*this%normal)/3._wp
    do q = 2, ubound(this%inletdirs,dim=1)
        iq = this%inletdirs(q)
        tmp = pdf(opposite(iq),icell) + this%m2 * this%rho * sum(qvec(iq,1:numd)*u(1:numd)) &
                                      + this%m3 * this%rho * sum(this%qtan(iq,1:numd)*u(1:numd))
        do j=1,ubound(this%nxyz,dim=1)
            tmp = tmp - this%m4 * pdf(abs(this%nxyz(j,q)),icell) * sign(1,this%nxyz(j,q))
        end do
        pdf(iq,icell) = tmp
    end do

    end associate
end subroutine

function getbcvel_densitybc(this,icell) result(u)
    class(densitybc_t)     :: this
    integer,    intent(in)  :: icell
    real(kind=wp)           :: u(1:numd),tmp
    integer                 :: q
    associate(pdf=>this%sl%mlist%pdf)

    tmp = 0._wp
    do q=1,ubound(this%tandirs,dim=1)
        tmp = tmp + pdf(this%tandirs(q),icell)
    end do

    do q=1,ubound(this%outletdirs,dim=1)
        tmp = tmp + 2._wp*pdf(this%outletdirs(q),icell)
    end do

    u = 0._wp
    u(this%iu) = (tmp/ this%rho - 1._wp)*this%signu

    end associate
end function

end module
