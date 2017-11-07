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
    use collisionBGK_mod
    use densitybc_mod
    implicit none
    type(part3d_t),         target  :: part
    type(imgproc3d_t)               :: imgpro
    type(imageobj3d_t),     target  :: imgobj
    type(solidconstruct3d_t)        :: solid
    type(solidfile_t),      target  :: sfile
    type(output3d_t)                :: out
    type(collisionBGK_t)            :: col
    type(densitybc_t)               :: lefden,rigden

    integer :: xlen=250,ylen=250,zlen=4

    call initD3Q19lattice()

    call sfile%init01(file="3D_porous_250x250x4.dat",origin=(/1,1,1/),dims=(/xlen,ylen,zlen/),resolution=1)
    call imgobj%init(dims=(/xlen,ylen,zlen/))
    call solid%init01(imgobj,sfile)
    call solid%readfile

    call imgpro%init(imageobj=imgobj,part=part)
    call imgpro%setpbc(axis=(/.false.,.true.,.true./))
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
        else
            ! set collision of other boundary sitations
            call part%setbcfluidcol(sit,col)
        endif
    end do

    call part%checkcollisions
    end block

    call out%init(part)

    block ! Domain Density and Velocity initialization
    integer n,icell
    do n = 1, part%nslists
        do icell = part%slists(n)%pt%ist, part%slists(n)%pt%ind
            call part%slists(n)%pt%initpdf(iscell=icell,rho=1._wp,u=(/0.0_wp,0.0_wp,0.0_wp/))
        end do
    end do
    end block

    block ! Stream and Collision
    integer     :: tstep
    real        :: t1,t2
    call out%writevtk(0)
    call cpu_time(t1)
    do tstep = 1, 10000
        call part%stream
        call part%collide
        if(mod(tstep,100).eq.0) then
            write(*,*) tstep
            call out%writevtk(tstep)
        endif
    end do
    call cpu_time(t2)

    print *, "time=",t2-t1,"sec"
    end block

end program
