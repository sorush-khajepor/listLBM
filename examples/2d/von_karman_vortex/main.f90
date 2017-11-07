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
    use plotcellkinds_mod
    use writeoutput3d_mod
    use stream_mod
    use colBGKextForce_mod
    use lattutil_mod
    use velocitybc_mod
    implicit none
    integer :: xlen=400,ylen=100
    type(part2d_t),target ::    part
    type(imgproc2d_t) :: imgpro
    type(imageobj2d_t),target  :: imgobj
    type(solidconstruct2d_t) ::  solid
    type(solidfile_t) :: solidparam
    type(plotcellkinds2d_t) :: pcellkinds
    type(output2d_t) :: out
    type(collisionBGK_t)            :: col
    type(outflowbc_t)               :: outbc
    type(velocitybc_t)              :: velbc
    type(blankcollision_t)          :: blankcol
    integer q

    call initD2Q9lattice()

    call imgobj%init(dims=(/xlen,ylen/))

    call solid%init(imgobj)

!	The cylinder is a tiny bit off-center.
    call solid%sphere(centerxyz=(/ xlen/4, (ylen/2)-1 /),rad=8.0_wp)
     do q=1,nbcsit
        if (cubevec(q,2).ne.0) &
        call solid%setbcaswall(dir=q)
    end do

    call imgpro%init(imageobj=imgobj,part=part)
    call imgpro%setpbc (axis=(/.true.,.true./))
    call imgpro%run()
    call imgpro%destruct()
    call imgobj%destruct()

    call out%init(part)


   block ! Initialization
    integer :: sit,n,icell
    call col%init(tau=.55_wp)
    call outbc%init(fcol=col,normal = (/1,0/))
    call velbc%init(fcol=col,normal = (/-1,0/),u=(/0.1_wp,0.0_wp/))
    call part%setinfluidcol(col)

    do sit=1,nbcsit
        if (cubevec(sit,1).eq.-1) then
            call part%setbcfluidcol(sit,velbc)
        elseif (cubevec(sit,1).eq.1) then
            call part%setbcfluidcol(sit,outbc)
        elseif(sit.eq.2.or.sit.eq.4) then
            call part%setbcfluidcol(sit,col)
        else
            ! whatever we put here it is not running because don't have cell
            call part%setbcfluidcol(sit,blankcol)
        endif
    end do


    do n = 1, part%nslists
        do icell = part%slists(n)%pt%ist, part%slists(n)%pt%ind
            if (.not.associated(part%slists(n)%pt%col)) then
                print*,"Error: The collision for "//part%slists(n)%pt%name//' boundary list is not set.'
            endif
            call part%slists(n)%pt%initpdf(iscell=icell,rho=1._wp,u=(/0.0_wp,0.0_wp,0.0_wp/))
        end do
    end do

    end block

    block ! collision
    integer :: tstep
    real :: t1,t2

    call out%writevtk(0)
    call cpu_time(t1)
    do tstep = 1, 6000
        call part%stream
        call part%collide
        if(mod(tstep,200).eq.0) then
            write(*,*) tstep
            call out%writevtk(tstep)
        endif
    end do
    call cpu_time(t2)

    print *, "time=",t2-t1,"sec"
    end block


end program main
