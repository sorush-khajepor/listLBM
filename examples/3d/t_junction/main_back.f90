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
    implicit none
    integer :: xlen=60,ylen=50,zlen=10
    type(part3d_t),target      ::    part
    type(imgproc3d_t)          :: imgpro
    type(imageobj3d_t),target  :: imgobj
    type(solidconstruct3d_t)   ::  solid
    type(output3d_t) :: out

    call initD3Q19lattice()

    call imgobj%init(xlen=xlen,ylen=ylen,zlen=zlen)
    call solid%init01(imgobj)
    call solid%box((/1,1,1/),(/30,20,10/))
    call solid%box((/1,31,1/),(/30,50,10/))
    call solid%box((/41,11,1/),(/60,40,10/))

    call imgpro%init(imageobj=imgobj,part=part)
    call imgpro%setpbc (axis=(/.true.,.true.,.true./))
    call imgpro%run()
    call imgpro%destruct()
    call imgobj%destruct()
    call out%init(part)
    call out%writevtk(0)

end program
