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
module readsolid2d_mod
    use kinds_f
    use image2d_mod
    implicit none
    ! origin is the origin of the solid in the simulation domain.
    ! dims includes dimensions of the solid in x, y, and, z direction.
    ! resolution=n magnify each point of solid into n point in x and n point in y and n point in z direction
    ! therefore resolution = 2, each wall point of solid would be scaled to 4 wall point in 2D and 8 wall point in 3D LBM
    type solidfile_t
        character (len=500)  :: file = ''
        integer,allocatable  :: origin(:),dims(:)
        integer              :: resolution = 0
    contains
        procedure :: init01 => init01_solidfile
        generic   :: init   => init01
    end type

!   Duplicates in mirror way. nmx is the number of duplicationin x-dir
!   which for solid itself is zero. So, nmx =1, nmy=2 creates 2*3 domain
!   of copies of that solid
    type duplicateparam2d_t
        integer :: nmx = 0, nmy =0
    end type

    type solidconstruct2d_t
        class (imageobj2d_t),         pointer :: imageobj
        class (solidfile_t),          pointer :: sfile
        class (duplicateparam2d_t),   pointer :: dupparam
        logical  :: is_duplicate = .false.
        logical  :: is_readfile  = .false.
    contains
        procedure :: init01         => init01_solidconstruct2d
        generic   :: init           => init01
        procedure :: readfile       => readfile_solidconstruct2d
        procedure :: duplicate      => mirror_duplicate_solidconstruct2d
        procedure :: carve          => carve_solidconstruct2d
        procedure :: sphere         => sphere_solidconstruct2d
        procedure :: box            => box_solidconstruct2d
        procedure :: setbcaswall    => setbcaswall_solidconstruct2d
    end type

contains

!=============solidfile_t====================
subroutine init01_solidfile(this,file,origin,dims,resolution)
    class(solidfile_t) :: this
    character (len=*), intent(in) :: file
    integer, intent(in) ::  origin(:),dims(:),resolution
    allocate(this%origin(1:numd),this%dims(1:numd))
    this%file = trim(file)
    this%origin = origin
    this%dims = dims
    this%resolution = resolution
end subroutine

!============ solidconstruct2d_t============================

subroutine init01_solidconstruct2d(this,imageobj,sfile,dupparam)
    class(solidconstruct2d_t) :: this
    class(imageobj2d_t),               pointer,intent(in) :: imageobj
    type(solidfile_t),       optional, pointer,intent(in) :: sfile
    type(duplicateparam2d_t),optional, pointer,intent(in) :: dupparam

    this%imageobj  => imageobj
    if (present(dupparam)) then
        this%dupparam=>dupparam
        this%is_duplicate = .true.
    end if
    if (present(sfile)) then
        this%sfile => sfile
        this%is_readfile = .true.
    end if

end subroutine


!    ========================================================
!     Reading solid objects from file
!    ========================================================
subroutine readfile_solidconstruct2d (this)
    class(solidconstruct2d_t)   :: this
    Integer                     :: i,j,funit,p,q,y,x,n,cc
    Real(kind = wp)             :: tmp_real

    if (.not.this%is_readfile) &
        stop "Error: in readfile_solidconstruct2d sfile which contains geometry of porous is not entered! "
    associate(img=>this%imageobj, ss=>this%sfile)
    open(newunit=funit, file= trim(ss%file),status='old',action='read')
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!  Increase resolution of solid   !!
!   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! q subcell in y direction
    ! p subcell in x direction
    ! length in x-dir = p*n
    ! length in y-dir = q*n
    write(*,'(A40)',advance='no'),"Reading solid file is started ...      "
    cc = 0
    n = ss%resolution
    do q = 1, ss%dims(2)
        do p = 1, ss%dims(1)
            read(funit,103) tmp_real
            do i = 1, n
                do j = 1, n
                    x=(p-1)*n+i+ss%origin(1)-1
                    y=(q-1)*n+j+ss%origin(2)-1

                    if (tmp_real .GE. 0.90_wp) then
                        cc = cc + 1
                        call img%setcellas((/x,y/),"wall")
                    end if
                end do
            end do
        end do
    end do
103 format (G17.9)
    close (funit)
    end associate
    write(*,*) "done."
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Dublicating solid in mirror way !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! mx counter of mirrors in x-dir
    ! my counter of mirrors in y-dir
    ! first of all I go to fill first row in x direction then
    ! I replicate all that row in y direction
    ! nmx = n mirror in x-dir
    ! nmy = n mirror in y-dir
    ! bpx = boundary point in x-dir
    ! bpy = boundary point in y_dir    each block starts from bpx+1 --> nbpx
    ! nbpx = next bpx
!**********************************************************
!           mx=0                        1
!    !-----------------!       !-----------------!
!    1                bpx     bpx+1             nbpx
!**********************************************************
! x-dir replications

subroutine mirror_duplicate_solidconstruct2d (this)
    implicit none
    class(solidconstruct2d_t) :: this
    Integer              :: i,j,n,mx,my,bpx,bpy,nbpx,nbpy
    associate(img=>this%imageobj, sf=>this%sfile,image=>this%imageobj%image, dp=>this%dupparam)

    if (.not.this%is_duplicate) &
        stop "Error: in mirror_duplicate2d_t dupparam which contains duplication parameters of porous is not entered! "

    n = sf%resolution
    bpx = sf%dims(1)*n + sf%origin(1) - 1
    do mx = 1, dp%nmx
        nbpx = bpx + sf%dims(1)*n
        do i = bpx + 1, nbpx
            do j = sf%origin(2), n*sf%dims(2) + sf%origin(2) -1
                image(i,j,1) = image( 2*bpx-i+1,j,1)
            end do
        end do
        bpx = nbpx
    end do
! y-dir replication
    bpy = sf%dims(2)*n + sf%origin(2) - 1
    do my = 1, dp%nmy
        nbpy = bpy + sf%dims(2)*n
        do j = bpy + 1, nbpy
            do i = sf%origin(1), (dp%nmx+1)*(n*sf%dims(1)) + sf%origin(1) -1
                image(i,j,1) = image(i,2*bpy-j+1,1)
            end do
        end do
        bpy = nbpy
    end do
    end associate
end subroutine

!----------------------------------------------------------------------
!   Defining Solid Objects inside Domain
!----------------------------------------------------------------------
! centerxyz = center of sphere with argument x,y,z
! rad = radius of the sphere
subroutine sphere_solidconstruct2d (this,centerxyz,rad)
    class(solidconstruct2d_t) :: this
    integer, intent(in)       :: centerxyz(:)
    real(kind=wp),intent(in)  :: rad
    integer :: i,j
    real(kind=wp) :: r2
    associate(img=>this%imageobj,c=>centerxyz)
    do i = c(1)-int(rad)-1, c(1)+int(rad)+1
    do j = c(2)-int(rad)-1, c(2)+int(rad)+1
        if(img%is_indomain((/i,j/))) then
            r2 = real((c(1)-i)**2+(c(2)-j)**2,kind=wp)
            if (r2.le.rad**2) then
                call img%setcellas((/i,j/),"wall")
            endif
        endif
    enddo
    enddo
    end associate
end subroutine


! p1xyz = point 1 with argument x,y,z
! p2xyz = point 2 with argument x,y,z
subroutine box_solidconstruct2d (this,p1xyz,p2xyz)
    class(solidconstruct2d_t) :: this
    integer, intent(in) :: p1xyz(:),p2xyz(:)
    integer :: i,j,rig,top,bot,lef
    associate(img=>this%imageobj)

    rig = max(p1xyz(1),p2xyz(1))
    lef = min(p1xyz(1),p2xyz(1))
    top = max(p1xyz(2),p2xyz(2))
    bot = min(p1xyz(2),p2xyz(2))

    do i = lef, rig
        do j = bot, top
             if(img%is_indomain((/i,j/))) &
                call img%setcellas((/i,j/),"wall")
        end do
    end do
    end associate
end subroutine



!========================================================
! Change the type of the exposed layer of wall to fluid
!========================================================
subroutine carve_solidconstruct2d (this,stage)
    class(solidconstruct2d_t) :: this
    integer,intent(in)        :: stage
    Integer                   :: i,j,rig,lef,top,bot,summ,ct
    Integer(kind = i1p)       :: tmp_image(1:this%imageobj%dims(1),1:this%imageobj%dims(2))
    associate(img=>this%imageobj,image=>this%imageobj%image)
    do ct = 1, stage
    tmp_image(1:img%dims(1),1:img%dims(2)) = image(1:img%dims(1),1:img%dims(2),1)
    do i=1,img%dims(1)
        rig = mod(i,img%dims(1))+1
        lef = mod(i+img%dims(1)-2,img%dims(1))+1
        do j=1,img%dims(2)
            if (.not.img%is_cellkind((/i,j/),"wall")) cycle
            if (j.eq.1.or.j.eq.img%dims(2).or.i.eq.1.or.i.eq.img%dims(1)) cycle ! keep top lef rig bottom boundaries
            top = mod(j,img%dims(2))+1
            bot = mod(j+img%dims(2)-2,img%dims(2))+1
            summ = int(image(i,top,1),i4p)  + int(image(i,bot,1),i4p)  +int(image(rig,j,1),i4p)  +int(image(lef,j,1),i4p) + &
                   int(image(rig,top,1),i4p)+int(image(lef,top,1),i4p)+int(image(rig,bot,1),i4p)+int(image(lef,bot,1),i4p)
            if (summ.Ne.8*int(img%cellkind%wall,i4p)) then
                tmp_image(i,j)= img%cellkind%fluid ! make wall cell fluid if it is in contact with fluids.
            end if
        enddo
    enddo
    image(1:img%dims(1),1:img%dims(2),1) = tmp_image(1:img%dims(1),1:img%dims(2))
    enddo
end associate
end subroutine


! Direction is based on cubevec. Therefore, corners, edges and faces
! can be refered (without corners or edges)
subroutine setbcaswall_solidconstruct2d (this,dir)
    class(solidconstruct2d_t) :: this
    integer,intent(in) :: dir
    Integer                   :: i,j,k,l1(1:numd),l2(1:numd)
    logical :: rev(1:numd)
    integer :: irev(1:numd),cb(1:numd)
    associate(img=>this%imageobj)
    call img%getbclimits(dir=dir,l1xyz=l1,l2xyz=l2)
    do i=l1(1),l2(1)
        do j=l1(2),l2(2)
            call img%setcellas((/i,j/),"wall")
        enddo
    end do
    end associate
end subroutine

end module
