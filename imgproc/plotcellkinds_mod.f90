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
module plotcellkinds_mod
    use part_mod
    use utility_mod, only: arrayofnames_t
    implicit none
    type plotcellkinds2d_t
    contains
        procedure :: plot => plot_plotcellkinds2d
    end type

    type plotcellkinds3d_t
    contains
        procedure :: plot => plot_plotcellkinds3d
    end type
contains

    subroutine plot_plotcellkinds2d(this,part,resolution,z)
    class(plotcellkinds2d_t) :: this
    class(part_t) :: part
    integer :: funit,icell,n,xlen,ylen,xres,yres
    integer :: resolution(:)
    integer, optional :: z

    xres = resolution(1)
    yres = resolution(2)

    xlen = part%main%xlen
    ylen = part%main%ylen
    call system("mkdir -p plotcellkinds")

    do n=1,part%nslists
    open(newunit=funit,file="./plotcellkinds/"//part%slists(n)%pt%name//'.dat',status="replace")
    write(funit,*) "# x,y,sublist_ID"
    do icell =part%slists(n)%pt%ist,part%slists(n)%pt%ind
        write(funit,*) part%slists(n)%pt%getijk(icell), part%slists(n)%pt%id
    end do
    close(funit)
    enddo

    open(newunit=funit,file="./plotcellkinds/grid_gp.cfg")
    write(funit,*) "set terminal png size", xres,",",yres
    write(funit,*)'set output "./plotcellkinds/grid.png"'
    write(funit,*)'set key outside'
    write(funit,*)'set xr[-1:',xlen+3,']'
    write(funit,*)'set yr[-1:',ylen+3,']'
    write(funit,*)'p  "./plotcellkinds/'//trim(part%slists(1)%pt%name)//'.dat"',' u 1:2:3 title', &
                  '"'//trim(part%slists(1)%pt%name)//'"',' w p, \'
    do n=2,part%nslists-1
        if(part%slists(n)%pt%ist.le.part%slists(n)%pt%ind) &
        write(funit,*)'"./plotcellkinds/'//part%slists(n)%pt%name//'.dat"',' u 1:2:3 title',   &
                      '"'//part%slists(n)%pt%name//'"',' w p, \'
    enddo
    if(part%slists(n)%pt%ist.le.part%slists(n)%pt%ind) &
    write(funit,*) '"./plotcellkinds/'//part%slists(n)%pt%name//'.dat"',' u 1:2:3 title', &
            '"'//part%slists(n)%pt%name//'"',' w p'
    close (funit)
    call system("gnuplot ./plotcellkinds/grid_gp.cfg")
    do n=1,part%nslists
        call system('rm ./plotcellkinds/'//part%slists(n)%pt%name//'.dat')
    enddo

    end subroutine



    subroutine plot_plotcellkinds3d(this,part,resolution,z)
    class(plotcellkinds3d_t) :: this
    class(part_t) :: part
    integer :: funit,icell,n,xlen,ylen,xres,yres,zres,zlen
    integer :: resolution(:)
    integer :: ijk(1:numd)
    integer, optional :: z
    character(len=10) :: zc
    logical :: avail(1:part%nslists)
    type(arrayofnames_t) :: fnames(1:part%nslists)

    avail = .false.
    write(zc,'(i3)') z


    xres = resolution(1)
    yres = resolution(2)
    zres = resolution(3)

    xlen = part%main%xlen
    ylen = part%main%ylen
    zlen = part%main%zlen
    call system("mkdir -p plotcellkinds")


    do n=1,part%nslists
        fnames(n)%mem = "./plotcellkinds/"//part%slists(n)%pt%name//'.dat'
    end do

    do n=1,part%nslists
    open(newunit=funit,file=fnames(n)%mem,status="replace")
    write(funit,*) "# x,y,sublist_ID"
    do icell =part%slists(n)%pt%ist,part%slists(n)%pt%ind
        ijk = part%slists(n)%pt%getijk(icell)
        if (ijk(3).eq.z) then
        avail(n) = .true.
        write(funit,*) ijk(1),ijk(2), part%slists(n)%pt%id
        endif
    end do
    close(funit)
    enddo

    open(newunit=funit,file="./plotcellkinds/grid_gp.cfg",status="replace")
    write(funit,*) "set terminal png size", xres,",",yres
    write(funit,*)'set output "./plotcellkinds/grid'//trim(zc)//'.png"'
    write(funit,*)'set key outside'
    write(funit,*)'set key spacing 6'
    write(funit,*)'set xr[-1:',xlen+3,']'
    write(funit,*)'set yr[-1:',ylen+3,']'
    write(funit,'(A2)',advance='no')'p '
   ! if(avail(1)) &
    write(funit,*)'"'//fnames(1)%mem//'"',' u 1:2:3 title', &
                  '"'//trim(part%slists(1)%pt%name)//'"',' w p pointsize 7, \'
    do n=2,part%nslists-1
    !    if(avail(n)) &
        write(funit,*)'"'//fnames(n)%mem//'"',' u 1:2:3 title',   &
                      '"'//part%slists(n)%pt%name//'"',' w p pointsize 7 , \'
    enddo
    !if(avail(n)) &
    write(funit,*) '"'//fnames(n)%mem//'"',' u 1:2:3 title', &
            '"'//part%slists(n)%pt%name//'"',' w p pointsize 7'
    close (funit)
    call system("gnuplot ./plotcellkinds/grid_gp.cfg")
    do n=1,part%nslists
        call system('rm '//fnames(n)%mem)
    enddo

    end subroutine

end module plotcellkinds_mod
