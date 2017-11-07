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
module readsolid3d_mod
    use readsolid2d_mod
    implicit none

    type,extends(solidconstruct2d_t) :: solidconstruct3d_t
    contains
        procedure :: box            => box_solidconstruct3d
        procedure :: sphere         => sphere_solidconstruct3d
        procedure :: setbcaswall    => setbcaswall_solidconstruct3d
        procedure :: readfile       => readfile_solidconstruct3d
    end type

contains

! p1xyz = point 1 with argument x,y,z
! p2xyz = point 2 with argument x,y,z
subroutine box_solidconstruct3d (this,p1xyz,p2xyz)
    class(solidconstruct3d_t)               :: this
    integer,                    intent(in)  :: p1xyz(:),p2xyz(:)
    integer                                 :: i,j,k,rig,top,bot,lef,fro,bac
    associate(img=>this%imageobj)

    rig = max(p1xyz(1),p2xyz(1))
    lef = min(p1xyz(1),p2xyz(1))
    top = max(p1xyz(2),p2xyz(2))
    bot = min(p1xyz(2),p2xyz(2))
    fro = max(p1xyz(3),p2xyz(3))
    bac = min(p1xyz(3),p2xyz(3))

    do i = lef, rig
        do j = bot, top
            do k = bac,fro
                if(img%is_indomain((/i,j,k/))) &
                    call img%setcellas((/i,j,k/),"wall")
            enddo
        end do
    end do
    end associate
end subroutine

! centerxyz = center of sphere with argument x,y,z
! rad = radius of the sphere
subroutine sphere_solidconstruct3d (this,centerxyz,rad)
    class(solidconstruct3d_t)               :: this
    integer,                    intent(in)  :: centerxyz(:)
    real(kind=wp),              intent(in)  :: rad
    integer                                 :: i,j,k
    real(kind=wp)                           :: r
    associate(img=>this%imageobj,c=>centerxyz)

    do i = c(1)-int(rad)-1, c(1)+int(rad)+1
    do j = c(2)-int(rad)-1, c(2)+int(rad)+1
    do k = c(3)-int(rad)-1, c(3)+int(rad)+1
        if(img%is_indomain((/i,j,k/))) then
            r = sqrt(real((c(1)-i)**2+(c(2)-j)**2+(c(3)-k)**2,kind=wp))
            if (r.le.rad) then
                call img%setcellas((/i,j,k/),"wall")
            endif
        endif
    enddo
    enddo
    enddo

    end associate
end subroutine

!    ========================================================
!     Reading solid objects from file
!
!       A) Increase resolution of solid
!
!           q subcell in y direction
!           p subcell in x direction
!           length in x-dir = p*n
!           length in y-dir = q*n
!    ========================================================
subroutine readfile_solidconstruct3d (this)
    class(solidconstruct3d_t)   :: this
    Integer                     :: i,j,k,funit,p,q,r,z,y,x,n,cc,iswall


    if (.not.this%is_readfile) &
        stop "Error: in readfile_solidconstruct2d sfile which contains geometry of porous is not entered! "
    associate(img=>this%imageobj, ss=>this%sfile)
    open(newunit=funit, file= trim(ss%file),status='old',action='read')
    write(*,'(A40)',advance='no'),"Reading solid file is started ...      "
    cc = 0
    n = ss%resolution


    do r = 1, ss%dims(3)
    do q = 1, ss%dims(2)
    do p = 1, ss%dims(1)
        read(funit,*) iswall
        do k = 1, n
        do j = 1, n
        do i = 1, n
            x=(p-1)*n+i+ss%origin(1)-1
            y=(q-1)*n+j+ss%origin(2)-1
            z=(r-1)*n+k+ss%origin(3)-1
            if (iswall .eq. 1) then
                cc = cc + 1
                call img%setcellas((/x,y,z/),"wall")
            end if
        enddo;        enddo;        enddo
    enddo;    enddo;    enddo
    close (funit)
    end associate
    write(*,*) "done."
end subroutine

! dir (direction) is based on cubevec.
subroutine setbcaswall_solidconstruct3d (this,dir)
    class(solidconstruct3d_t)               :: this
    integer,                    intent(in)  :: dir
    Integer                                 :: i,j,k,l1(1:numd),l2(1:numd)
    associate(img=>this%imageobj)

    ! dir = bcsitname%getno('top')
    call img%getbclimits(dir=dir,l1xyz=l1,l2xyz=l2)
    do i=l1(1),l2(1)
        do j=l1(2),l2(2)
            do k=l1(3),l2(3)
                call img%setcellas((/i,j,k/),"wall")
            enddo
        enddo
    end do

    end associate
end subroutine


end module
