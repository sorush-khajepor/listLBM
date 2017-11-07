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
module initfield_mod
    use kinds_f
    use lattice_mod
    use sublist_mod
    implicit none

    type,   abstract    :: initfield_t
        integer,    allocatable :: lowbound(:),upbound(:)
    contains    
        procedure :: run => run_initfield
        procedure :: is_infield => is_infield_initfield
        procedure(getvel_initfield),deferred    :: getvel
        procedure(getrho_initfield),deferred    :: getrho
    end type
    
    type, extends(initfield_t) :: initfielduniform_t
        real(kind=wp),  allocatable :: u(:)
        real(kind=wp)               :: rho
    contains
        procedure  :: init01_initfielduniform   => init01_initfielduniform
        generic    :: init                      => init01_initfielduniform
        procedure  :: getvel                    => getvel_initfielduniform
        procedure  :: getrho                    => getrho_initfielduniform
    end type

    type, extends(initfield_t) :: initbubble_t
        real(kind=wp),  allocatable :: u(:)
        real(kind=wp)               :: rho_in,rho_out
        integer                     :: interface_width,rad
        integer,        allocatable :: ijkref(:)
    contains
        procedure :: init01_initbubble      => init01_initbubble
        generic   :: init                   => init01_initbubble
        procedure :: getrho                 => getrho_initbubble
        procedure :: getvel                 => getvel_initbubble
    end type

    type, extends(initbubble_t) :: initflatiface_t
        real(kind=wp), allocatable :: normal(:)
    contains
    procedure :: getrho                 => getrho_initflatiface
    procedure :: init01_initflatiface   => init01_initflatiface
    generic   :: init                   => init01_initflatiface

    end type


interface getvelfield
function getvel_initfield(this,ijk) result(u)
    import initfield_t,wp,numd
    class(initfield_t)              :: this
    real(kind=wp)                   :: u(1:numd)
    integer,            intent(in)  :: ijk(1:numd)
end function
end interface

interface getrhofield
function getrho_initfield(this,ijk) result(rho)
    import initfield_t,wp,numd
    class(initfield_t)              :: this
    real(kind=wp)                   :: rho
    integer,            intent(in)  :: ijk(1:numd)
end function
end interface


contains


subroutine run_initfield(this,slist)
    class(initfield_t)              :: this
    class(sublist_t),   intent(in)  :: slist
    integer :: ijk(1:numd),iscell

    do iscell = slist%ist, slist%ind
        ijk = slist%getijk(iscell)
        if (this%is_infield(ijk)) &
            call slist%initpdf(iscell=iscell,rho=this%getrho(ijk),u=this%getvel(ijk))
    end do

end subroutine

function is_infield_initfield(this,ijk) result(is)
    class(initfield_t)              :: this
    integer,            intent(in)  :: ijk(1:numd)
    logical                         :: is
    integer                         :: val,d,intlog

    val = 0
    do d = 1, numd
        intlog = ijk(d).le.max(this%lowbound(d),this%upbound(d))
        val = val + intlog
        intlog = ijk(d).ge.min(this%lowbound(d),this%upbound(d))
        val = val + intlog
    end do

    if (val.eq.2*numd) then
        is = .true.
    else
        is = .false.
    end if

end function

!======================Init Uniform Field=================
subroutine init01_initfielduniform(this,rho,u,lowbound,upbound)
    class(initfielduniform_t)                                       :: this
    real(kind=wp),                          intent(in)              :: rho
    real(kind=wp),                          intent(in), optional    :: u(1:numd)
    integer,                                            optional    :: lowbound(1:numd),upbound(1:numd)

    allocate(this%u(1:numd),this%lowbound(1:numd),this%upbound(1:numd))
    if(present(u)) then
        this%u = u
    else
        this%u = 0._wp
    end if
    this%rho = rho

    if (present(lowbound)) then
        this%lowbound = lowbound
    else
        this%lowbound = -huge(0)
    endif

    if (present(upbound)) then
        this%upbound = upbound
    else
        this%upbound = huge(0)
    endif

end subroutine


function getvel_initfielduniform(this,ijk) result(u)
    class(initfielduniform_t)               ::this
    real(kind=wp)                           :: u(1:numd)
    integer,                    intent(in)  :: ijk(1:numd)

    u = this%u

end function

function getrho_initfielduniform(this,ijk) result(rho)
    class(initfielduniform_t)               :: this
    real(kind=wp)                           :: rho
    integer,                    intent(in)  :: ijk(1:numd)

    rho = this%rho

end function

!============================= Sphere Gradient===============

subroutine init01_initbubble(this,ijkref,rho_in,rho_out,rad,interface_width,u,lowbound,upbound)
    class(initbubble_t) :: this
    integer,                            intent(in)              :: ijkref(1:numd),interface_width,rad
    real(kind=wp),                      intent(in)              :: rho_in,rho_out
    real(kind=wp),                      intent(in), optional    :: u(1:numd)
    integer,                                        optional    :: lowbound(1:numd),upbound(1:numd)


    this%interface_width = interface_width
    this%rad=rad
    this%rho_in  = rho_in
    this%rho_out = rho_out
    allocate(this%u(1:numd),this%ijkref(1:numd),this%lowbound(1:numd),this%upbound(1:numd))
    this%ijkref = ijkref
    if(present(u)) then
        this%u = u
    else
        this%u = 0._wp
    end if
    if (present(lowbound)) then
        this%lowbound = lowbound
    else
        this%lowbound = -huge(0)
    endif

    if (present(upbound)) then
        this%upbound = upbound
    else
        this%upbound = huge(0)
    endif
end subroutine



function getrho_initflatiface(this,ijk) result(rho)
    class(initflatiface_t)             :: this
    integer,                intent(in)  :: ijk(1:numd)
    real(kind=wp)                       :: rho,dr
    associate(rho_in=>this%rho_in,rho_out=>this%rho_out)

    dr = abs(sum((ijk(1:numd)-this%ijkref(1:numd))*this%normal(1:numd)))

    rho = (rho_in + rho_out)/2.0_wp-(rho_in-rho_out)/2.0_wp*tanh(2.0_wp*(dr-this%rad)/ &
        this%interface_width)

    end associate
end function


function getrho_initbubble(this,ijk) result(rho)
    class(initbubble_t)             :: this
    integer,                intent(in)  :: ijk(1:numd)
    real(kind=wp)                       :: rho,dr,tmp
    integer                             :: d
    associate(rho_in=>this%rho_in,rho_out=>this%rho_out)

    tmp = 0._wp
    do d = 1, numd
       tmp = tmp + real((ijk(d)-this%ijkref(d))**2,wp)
    end do
    dr = sqrt(tmp)

    rho = (rho_in + rho_out)/2.0_wp-(rho_in-rho_out)/2.0_wp*tanh(2.0_wp*(dr-this%rad)/ &
        this%interface_width)

    end associate
end function

function getvel_initbubble(this,ijk) result(u)
    class(initbubble_t)                     ::this
    real(kind=wp)                           :: u(1:numd)
    integer,                    intent(in)  :: ijk(1:numd)

    u = this%u

end function

subroutine init01_initflatiface(this,ijkref,rho_in,rho_out,rad,interface_width,normal,u,lowbound,upbound)
    class(initflatiface_t) :: this
    integer,                            intent(in)              :: ijkref(1:numd),interface_width,rad,normal(1:numd)
    real(kind=wp),                      intent(in)              :: rho_in,rho_out
    real(kind=wp),                      intent(in), optional    :: u(1:numd)
    integer,                                        optional    :: lowbound(1:numd),upbound(1:numd)


    this%interface_width = interface_width
    this%rad=rad
    this%rho_in  = rho_in
    this%rho_out = rho_out
    allocate(this%u(1:numd),this%ijkref(1:numd),this%lowbound(1:numd),this%upbound(1:numd))
    allocate(this%normal(1:numd))
    this%normal = normal/sqrt(real(sum(normal*normal),wp))
    this%ijkref = ijkref
    if(present(u)) then
        this%u = u
    else
        this%u = 0._wp
    end if
    if (present(lowbound)) then
        this%lowbound = lowbound
    else
        this%lowbound = -huge(0)
    endif

    if (present(upbound)) then
        this%upbound = upbound
    else
        this%upbound = huge(0)
    endif
end subroutine


end module
