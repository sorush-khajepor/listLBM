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
program main
    use imgproc3d_mod
    use image3d_mod
    use readsolid3d_mod
    use writeoutput3d_mod
    use stream_mod
    use lattutil_mod
    use partppot_mod
    use ppot_mod
    use initfield_mod
    implicit none
    type(partppot3d_t),     target  :: part
    type(imageobj3d_t),     target  :: imgobj
    type(imgproc3d_t)               :: imgpro
    type(solidconstruct3d_t)        :: solid
    type(solidfile_t)               :: solidparam
    type(colBGKpp_t)                :: col
    type(blankcollision_t)          :: blankcol
    type(colBGKppwall_t)            :: colw
    type(colBGKppgf_t)              :: colgf
    type(output3d_t)                :: out


    integer :: xlen=50,ylen=50,zlen=50

    call initD3Q19lattice()

    call imgobj%init(dims=(/xlen,ylen,zlen/))

    call solid%init(imgobj)

!    block ! set all BCs as wall
!    integer :: i
!    do i = 5, 8
!        call solid%setbcaswall(dir = i)
!    end do
!    call solid%setbcaswall(dir = 1)
!    call solid%setbcaswall(dir = 3)
!    endblock

    call solid%box((/xlen/2-5,ylen/2-5,zlen/2-5/),(/xlen/2+5,ylen/2+5,zlen/2+5/))


    call imgpro%init(imageobj=imgobj,part=part)
    call imgpro%setpbc (axis=(/.true.,.true.,.true./))
    call imgpro%run()
    call imgpro%destruct()
    call imgobj%destruct()

    call out%init(part)


    block ! Collision Initialization
        type(potparam_t) :: pm
        pm%getpot => pot_SC
        pm%G = -1.5_wp
        pm%rho_0 = 1._wp
        call col%init(tau=1._wp,potparam=pm)
        call colw%init(rho_wall=2._wp,potparam=pm)

        if (.not.associated(part%main%pot)) &
            call part%main%allocpot()

        call part%setfluidcol(col)
        call part%slists(28)%pt%alloccollision(colw)
        call part%slists(29)%pt%alloccollision(colw)
        call part%slists(30)%pt%alloccollision(colgf)
    end block


    block ! Domain Density and Velocity Initialization

    integer :: icell,n
    type(initbubble_t) :: initbub

    call initbub%init(ijkref=(/25,25,25/),rho_in=0.25_wp,rho_out=1.2_wp,&
                      rad=5,interface_width=3, &
                      lowbound = (/1,1,1/),upbound=(/50,50,50/))
    do n = 1, part%nslists
        call initbub%run(slist = part%slists(n)%pt)
    end do

    end block

    block ! LBM Iterations
    integer :: tstep
    real :: t1,t2

    call out%writevtk(0)

    call cpu_time(t1)

    do tstep = 1, 10000
        call part%stream
        call part%collide
        if(mod(tstep,200).eq.0) then
            write(*,*) "tstep = ", tstep
            call out%writevtk(tstep)
        endif
    end do

    call cpu_time(t2)

    write(*,*) "run time =",t2-t1,"sec"
    end block



end program main
