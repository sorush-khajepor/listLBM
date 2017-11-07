!============================== License GPLv3 ===================================
!    listLBM is a lattice Boltzmann solver for multiphase flow in porous media
!    Copyright (C) 2017  Sorush Khajepor sk451@hw.ac.uk
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


! This example simulates liquid drop in equilibrium with its vapor 
! with Shan-Chen model[1]. The domain is 3D with the size of 50x50x50.
! Shan-Chen potential is 
! pot_SC = rho_0 *[1-exp(-rho/rho_0)] 
! where G=-1.5 and rho_0= 1.0. You can choose different potentials from
! potentialfunc_mod.f90 module. At the equilibrium there is a vapor
! bubble in the center of a liquid box with densities of 0.22 and 1.4 in
! in the lattice system.
!
! [1] Shan, X. and H. Chen, Physical Review E, 1993. 47(3): p. 1815-1819.


program main
    use imgproc3d_mod
    use image3d_mod
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
    type(colBGKpp_t)                :: col
    type(blankcollision_t)          :: blankcol
    type(colBGKppwall_t)            :: colw
    type(colBGKppgf_t)              :: colgf
    type(output3d_t)                :: out

!	Setting domain size
    integer :: xlen=50,ylen=50,zlen=50

!	Initializing D3Q9 lattice
    call initD3Q19lattice()

!	Initializing imgobj ( 3D array of the domain box )
    call imgobj%init(dims=(/xlen,ylen,zlen/))

!	Initializing imgpro (image processing)
    call imgpro%init(imageobj=imgobj,part=part)

!	Setting periodic boundaries in x,y,z direction
!	which is necessary to remove unwanted solid nodes.
    call imgpro%setpbc (axis=(/.true.,.true.,.true./))
    call imgpro%run()
    call imgpro%destruct()
    call imgobj%destruct()

!	Initialize part (the container of sublists)
    call out%init(part)

! 	Collision Initialization
    block 
!		Pseudopotential model initialization
        type(potparam_t) :: pm

!		Choosing potential function
        pm%getpot => pot_SC
        pm%G = -1.5_wp
        pm%rho_0 = 1._wp

!		Setting relaxation time
        call col%init(tau=1._wp,potparam=pm)

!		Setting density for the the wall which is used
!		for fluid-solid force interactions. 
        call colw%init(rho_wall=2._wp,potparam=pm)

!		Checking if potential is allocated.
        if (.not.associated(part%main%pot)) &
            call part%main%allocpot()

!		Allocating the collision of different lists.

        call part%setfluidcol(col)
        call part%slists(28)%pt%alloccollision(colw)
        call part%slists(29)%pt%alloccollision(colw)
        call part%slists(30)%pt%alloccollision(colgf)
    end block


! 	Domain Density and Velocity Initialization
    block 

    integer :: icell,n
    type(initbubble_t) :: initbub

    call initbub%init(ijkref=(/25,25,25/),rho_in=0.20_wp,rho_out=1.25_wp,&
                      rad=5,interface_width=3, &
                      lowbound = (/1,1,1/),upbound=(/50,50,50/))
    do n = 1, part%nslists
        call initbub%run(slist = part%slists(n)%pt)
    end do

    end block

! 	LBM Iterations, stream & collision
    block 
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
