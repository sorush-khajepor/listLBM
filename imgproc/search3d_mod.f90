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
module search3d_mod
    use search2d_mod
    implicit none

    type, extends(search2d_t) :: search3d_t
    contains
        procedure :: blackwall => blackwall_search3d
        procedure :: ghost     => ghost_search3d
        procedure :: ghostw    => ghostw_search3d
        procedure :: ghostf    => ghostf_search3d
        procedure :: fluids   => fluids_search3d
        procedure :: delvainwall => delvainwall_search3d
    end type

contains

subroutine blackwall_search3d(this,imageobj,lislim)
    class(search3d_t)      :: this
    class(imageobj2d_t)    :: imageobj
    class(listlimits2d_t) :: lislim
    Integer        :: i,j,k,q,counter
    Integer        :: blackwno, wallno,pnei(1:numd,1:nq_1),ijk(1:numd)
    associate (wall =>imageobj%cellkind%wall,blackwall =>imageobj%cellkind%blackwall,  &
                image=>imageobj%image, xlen=>imageobj%dims(1),ylen=>imageobj%dims(2),zlen=>imageobj%dims(3))
    blackwno = 0
    wallno = 0
    do k = 1, zlen
        do j = 1, ylen
            do i = 1, xlen
                ijk = (/i,j,k/)
                pnei = imageobj%getpnei(ijk)
                counter = 0
                do q = 1, nq_1
                    if (imageobj%getimage(pnei(1:numd,q)).EQ.wall .or.  &
                        imageobj%getimage(pnei(1:numd,q)).EQ.blackwall) &
                        counter = counter + 1
                enddo
                if (counter.EQ.nq_1 ) then
                    call imageobj%setcellas(ijk,"blackwall")
                    blackwno = blackwno + 1
                end if
                if (imageobj%is_cellkind(ijk,"wall")) wallno = wallno + 1
            enddo
        enddo
    enddo
    lislim%blackwall = blackwno
    lislim%wall = wallno
    end associate
end subroutine


! finding number of ghost cells.
! As each ghost cell has some neighbours in domain, one of them
! -at least- should be fluid (not wall and not black_wall).
subroutine ghost_search3d (this,imageobj,lislim)
    class (search3d_t)     :: this
    class(imageobj2d_t)  :: imageobj
    class(listlimits2d_t)  :: lislim
    integer :: nei(1:numd,1:nq_1),ijk(1:numd)
    integer :: i,j,k,q,gno
    associate (img => imageobj, xlen=>imageobj%dims(1), ylen=>imageobj%dims(2),zlen=>imageobj%dims(3))

!   Initialization of all ghost cells with blackghost cellkind.
    do k = 0, zlen+1
        do j=0, ylen+1
            do i=0, xlen+1
                ijk = (/i,j,k/)
                if (.not.img%is_indomain(ijk)) call img%setcellas(ijk,"blackghost")
            enddo
        enddo
    enddo
    gno = 0
    do k = 0, zlen+1
        do j = 0, ylen+1
            do i=0, xlen+1
                ijk = (/i,j,k/)
                if (img%is_indomain(ijk)) cycle
                nei = img%getnei(ijk)
                do q = 1,nq_1
                    if (img%is_cellkind(nei(1:numd,q),"fluid")) then
                        call img%setcellas(ijk,"ghost")
                        gno = gno+1
                        exit
                    endif
                enddo
            end do
        end do
    enddo
    lislim%ghost = gno
    lislim%blackghost = lislim%total - lislim%indomain - gno

    end associate
end subroutine

! assigin ghostw based on the strain of the  periodic
! corresponding node
subroutine ghostw_search3d(this,imageobj,periodic,lislim)
    class(search3d_t)        :: this
    class(imageobj2d_t)    :: imageobj
    class(listlimits2d_t)    :: lislim
    class(pbc_t)         ::  periodic
    integer :: i,j,k,gwno,sit,pb(1:numd),ijk(1:numd)
    associate(xlen=>imageobj%dims(1),ylen=>imageobj%dims(2),zlen=>imageobj%dims(3),img=>imageobj)
    gwno = 0
    do k = 0, zlen+1
        do j = 0, ylen+1
            do i = 0, xlen+1
                ijk = (/i,j,k/)
                if (img%is_cellkind(ijk,"ghost")) then
                    pb = img%getpbodycell(ijk)
                    sit = img%getghostsit(ijk)
                    if (periodic%is_dir(sit) .and. &
                        img%is_cellkind(pb,"wall"))  then
                        call img%setcellas(ijk,"ghostw")
                        gwno = gwno + 1
                    endif
                endif
            enddo
        enddo
    enddo
    end associate
    lislim%ghostwall = gwno
end subroutine

! If the ghost cell is not gostwall, it is ghostfluid.
subroutine ghostf_search3d(this,imageobj,lislim)
    class(search3d_t)     :: this
    class(imageobj2d_t) :: imageobj
    class(listlimits2d_t) :: lislim
    integer :: i,j,k,gfno,ijk(1:numd)
    associate(xlen=>imageobj%dims(1),ylen=>imageobj%dims(2),zlen=>imageobj%dims(3),img=>imageobj)
    gfno = 0
    do k = 0, zlen+1
        do j = 0, ylen+1
            do i = 0, xlen+1
                ijk = (/i,j,k/)
                if (img%is_cellkind(ijk,"ghost").and.&
               .not.img%is_cellkind(ijk,"ghostw")) then
                        call img%setcellas(ijk,"ghostf")
                        gfno = gfno + 1
                endif
            enddo
        enddo
    enddo
    lislim%ghostfluid = gfno
    end associate
end subroutine


! The subroutine turns the useless wall cells  on the
! boundary to black wall after finding ghost-wall and
! ghost-fluid cells.
subroutine delvainwall_search3d(this,imageobj,lislim)
    class(search3d_t)     :: this
    class(imageobj2d_t)   :: imageobj
    class(listlimits2d_t) :: lislim
    integer :: i,j,k,bno,fno,q,nei(1:numd,1:nq_1),ijk(1:numd)
    associate(xlen=>imageobj%dims(1),ylen=>imageobj%dims(2),zlen=>imageobj%dims(3),img=>imageobj)
    bno = 0
    do k = 1, zlen
        do i= 1, xlen
            do j = 1, ylen
                ijk = (/i,j,k/)
                if (img%is_cellkind(ijk,"wall")) then
                    nei = img%getnei(ijk)
                    fno = 0
                    do q = 1,nq_1
                        if (img%is_cellkind(nei(1:numd,q),"fluid")) fno = fno+1
                    enddo
                    if (fno.eq.0) then
                        call img%setcellas(ijk,"blackwall")
                        bno=bno+1
                    endif
                endif
            end do
        end do
    enddo
    lislim%wall = lislim%wall - bno
    lislim%blackwall = lislim%blackwall + bno
    end associate
end subroutine

! Fluid cells including influid (internal fluid cells), and
! boundary fluid cells are found and are put in place of
! "fluid" kind but note that is_fluid function is still working.
subroutine fluids_search3d(this,imageobj,lislim)
    class(search3d_t) :: this
    class(imageobj2d_t)   :: imageobj
    class(listlimits2d_t) :: lislim
    integer :: i,j,k,sit,ijk(1:numd)
    associate(xlen=>imageobj%dims(1),ylen=>imageobj%dims(2),zlen=>imageobj%dims(3),img=>imageobj,ll=>lislim)
    ll%bcfluid = 0
    ll%influid = 0
!   Set internal cells as "influid" kind.
    do k = 2, zlen-1
        do i=2,xlen-1
            do j=2,ylen-1
                ijk = (/i,j,k/)
                if (img%is_cellkind(ijk,"fluid")) then
                    call img%setcellas(ijk,"influid")
                    ll%influid = ll%influid + 1
                endif
            end do
        end do
    end do
!   search faces (sits)
    do k=1,zlen
        do j=1,ylen
            do i =1,xlen
                ijk = (/i,j,k/)
                if ((.not.img%is_bcell(ijk)).or.(.not.img%is_cellkind(ijk,"fluid"))) cycle
                sit = img%getbcsit(ijk)
                call img%setcellas(ijk,bcsname(sit)%mem)
                ll%bcfluid(sit) = ll%bcfluid(sit) + 1
            enddo
        end do
    enddo
    end associate
end subroutine

end module
