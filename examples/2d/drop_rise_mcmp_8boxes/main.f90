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

! Multi-component simulation of drop rising in the secondary liquid.
! Both liquids have almost the same density. Shan-Chen model is used.
! To make simulation more complex, 8 non-wettable boxes are placed
! in the container.
! Find the video of this simulation in my youtube channel mentioned above.

program main
    use imgproc2d_mod
    use image2d_mod
    use readsolid2d_mod
    use writeoutput2d_mod
    use stream_mod
    use lattutil_mod
    use multipart_mod
    use mcmp_mod
    use initfield_mod

    implicit none

    type(partmcmppt_t),     target  :: parts(1:2)
    type(multipart2d_t),    target  :: mtpart
    type(imageobj2d_t),     target  :: imgobj
    type(imgproc2d_t)               :: imgpro
    type(solidconstruct2d_t)        :: solid
    type(solidfile_t)               :: solidparam
    type(colBGKmcmp_t)              :: col1,col2
    type(blankcollision_t)          :: blankcol
    type(colBGKmcmpwall_t)          :: colw1,colw2
    type(colBGKmcmpgf_t)            :: colgf
    type(output2d_t)                :: out

    integer :: xlen=500,ylen=600

    call initD2Q9lattice()
!	Allocate two components via parts
    allocate(parts(1)%pt,parts(2)%pt)

    call imgobj%init(dims=(/xlen,ylen/))

    call solid%init(imgobj)

    block ! set wall cells
    integer  i,j
    ! set boundaries as wall
    call solid%setbcaswall(dir = 2)
    call solid%setbcaswall(dir = 4)
    do i=5,6
        call solid%setbcaswall(dir = i)
    end do
    ! add 8 boxes to domain
    j = 350
    do i=120,420,100
        call solid%box(p1xyz=(/i,j/),p2xyz=(/i+40,j+40/))
    enddo

    i = 160
    j = 250
    do i=90,390,100
        call solid%box(p1xyz=(/i,j/),p2xyz=(/i+40,j+40/))
    enddo
    endblock


    call imgpro%init(imageobj=imgobj,part=parts(1)%pt)
    call imgpro%setpbc (axis=(/.false.,.true./))
    call imgpro%run()
    call imgpro%destruct()
    call imgobj%destruct()

    call parts(2)%pt%copy( parts(1)%pt)

    call out%init(parts(2)%pt)

    call mtpart%init(parts)

    block ! Collision Initialization
        type(potparam_t),target :: pm(1:1)
        pm(1)%getpot => pot_givesrho
        pm(1)%G = 0.75_wp

        call col1%init(tau=1._wp,mcpotparam=pm,parts=parts,extforce=(/0.0_wp,-0.000_wp/))
        call colw1%init(rho_wall=01.0_wp,mcpotparam=pm)
        colgf%nparts = 2

        call parts(1)%pt%setfluidcol(col1)
        call parts(1)%pt%slists(10)%pt%alloccollision(colw1)
        call parts(1)%pt%slists(11)%pt%alloccollision(colw1)
        call parts(1)%pt%slists(12)%pt%alloccollision(colgf)

    end block

    block ! Collision Initialization
        type(potparam_t),target :: pm(1:1)
        pm(1)%getpot => pot_givesrho
        pm(1)%G = 0.75_wp

        call col2%init(tau=1._wp,mcpotparam=pm,parts=parts,extforce=(/0.0_wp,+0.001_wp/))
        call colw2%init(rho_wall=03.01_wp,mcpotparam=pm)
        ! colgf%nparts = 2
        call parts(2)%pt%setfluidcol(col2)
        call parts(2)%pt%slists(10)%pt%alloccollision(colw2)
        call parts(2)%pt%slists(11)%pt%alloccollision(colw2)
        call parts(2)%pt%slists(12)%pt%alloccollision(colgf)

    end block



    block ! Domain Density and Velocity Initialization First component

    integer :: icell,n
    type(initbubble_t) :: initbub
    type(initflatiface_t) :: inittop

    call initbub%init(ijkref=(/xlen/2,50/),rho_in=0.05_wp,rho_out=1.19_wp,&
                      rad=30,interface_width=3, &
                      lowbound = (/1,1/),upbound=(/500,500/))

    call inittop%init(ijkref=(/xlen/2,ylen/),rho_in=0.05_wp,rho_out=1.19_wp,&
                      rad=90,interface_width=3, &
                      lowbound = (/1,501/),upbound=(/500,600/),normal = (/0,-1/))
    do n = 1, parts(1)%pt%nslists
        call initbub%run(slist = parts(1)%pt%slists(n)%pt)
        call inittop%run(slist = parts(1)%pt%slists(n)%pt)
    end do

    end block

    block ! Domain Density and Velocity Initialization Second Component

    integer :: icell,n
    type(initbubble_t) :: initbub
    type(initflatiface_t) :: inittop

    call initbub%init(ijkref=(/xlen/2,50/),rho_in=1.19_wp,rho_out=.05_wp,&
                      rad=30,interface_width=3, &
                      lowbound = (/1,1/),upbound=(/500,500/))

    call inittop%init(ijkref=(/xlen/2,ylen/),rho_in=1.19_wp,rho_out=0.05_wp,&
                      rad=90,interface_width=3, &
                      lowbound = (/1,501/),upbound=(/500,600/),normal = (/0,-1/))
    do n = 1, parts(2)%pt%nslists
        call initbub%run(slist = parts(2)%pt%slists(n)%pt)
        call inittop%run(slist = parts(2)%pt%slists(n)%pt)
    end do

    end block

    block ! LBM Iterations
    integer :: tstep
    real :: t1,t2

    call out%writevtk(0)

    call cpu_time(t1)

    do tstep = 1, 15000
        call mtpart%stream
        call mtpart%collide
        if(mod(tstep,5).eq.0) then
            write(*,*) "tstep = ", tstep
            call out%writevtk(tstep)
        endif
    end do

    call cpu_time(t2)

    write(*,*) "run time =",t2-t1,"sec"
    end block



end program main
