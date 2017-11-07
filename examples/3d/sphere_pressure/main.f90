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

! This example simulates flow over a sphere placed in a channel
! Pressure boundaries are set at left and right of the domain.
! All the other boundaries are wall.
program main
    use imgproc3d_mod
    use image3d_mod
    use readsolid3d_mod
    use plotcellkinds_mod
    use writeoutput3d_mod
    use stream_mod
    use colBGKextForce_mod
    use lattutil_mod
    use densitybc_mod
    implicit none
    integer :: xlen=50,ylen=50,zlen=20
    type(part3d_t),target           :: part
    type(imgproc3d_t)               :: imgpro
    type(imageobj3d_t),target       :: imgobj
    type(solidconstruct3d_t)        :: solid
    type(output3d_t)                :: out
    type(collisionBGK_t)            :: col
    type(densitybc_t)               :: lefden,rigden
    type(blankcollision_t)          :: blankcol
    integer q

    call initD3Q19lattice()
    call imgobj%init(dims=(/xlen,ylen,zlen/))
    call solid%init(imgobj)
    call solid%sphere(centerxyz=(/20,20,10/),rad=8.0_wp)

!	Simulation domain is a box like a cell.
!	Set all surfaces, edges, corners which have y value as wall
!	q is counter for lattice velocity vectors
    do q=1,nbcsit
        if (cubevec(q,2).ne.0) &
        call solid%setbcaswall(dir=q)
    end do

    call imgpro%init(imageobj=imgobj,part=part)
    call imgpro%setpbc(axis=(/.true.,.true.,.true./))
    call imgpro%run()
    call imgpro%destruct()
    call imgobj%destruct()

    block ! Collision Initialization
    integer :: sit
    call col%init(tau=1._wp)
    call lefden%init(fcol=col,normal = (/-1,0,0/),rho = 1.05_wp)
    call rigden%init(fcol=col,normal = (/ 1,0,0/),rho = 1.0_wp)
    call part%setinfluidcol(col)

    do sit=1,nbcsit
        if (cubevec(sit,1).eq.-1) then
            call part%setbcfluidcol(sit,lefden)
        elseif (cubevec(sit,1).eq.1) then
            call part%setbcfluidcol(sit,rigden)
        elseif (sit.eq.5.or.sit.eq.6) then
            call part%setbcfluidcol(sit,col)
        else
            ! whatever we put here it is not running because they don't have cell
            call part%setbcfluidcol(sit,blankcol)
        endif
    end do
    end block

    call out%init(part)
    block ! Domain Density and Velocity initialization
    integer n,icell
    do n = 1, part%nslists
        do icell = part%slists(n)%pt%ist, part%slists(n)%pt%ind
            if (.not.associated(part%slists(n)%pt%col)) then
                print*,"Error: The collision for "//trim(part%slists(n)%pt%name)//' boundary list is not set.'
            endif
            call part%slists(n)%pt%initpdf(iscell=icell,rho=1._wp,u=(/0.0_wp,0.0_wp,0.0_wp/))
        end do
    end do

    end block

    block ! Stream and Collision
    integer     :: tstep
    real        :: t1,t2
    call out%writevtk(0)
    call cpu_time(t1)
    do tstep = 1, 5000
        call part%stream
        call part%collide
        if(mod(tstep,1000).eq.0) then
            write(*,*) tstep
            call out%writevtk(tstep)
        endif
    end do
    call cpu_time(t2)

    print *, "time=",t2-t1,"sec"
    end block




end program
