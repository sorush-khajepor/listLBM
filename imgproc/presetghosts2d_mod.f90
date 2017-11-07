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
module presetghosts2d_mod
    use image2d_mod
    use part_mod
    use unmap_mod
    use pbc_mod
    implicit none

    type presetghosts2d_t
        class(part_t),          pointer :: part
        class(imageobj2d_t),    pointer :: imageobj
        class(unmap2d_t),       pointer :: unmap
        class(pbc_t),           pointer :: periodic
    contains
        procedure :: init               =>  init_presetghost2d
        procedure :: run                =>  run_presetghost2d
        procedure :: mapijk2icell       =>  mapijk2icell_presetghost2d
        procedure :: initgfoutletdirs   =>  initghostfoutletdirs_presetghost2d
        procedure :: setclosest2gcell   =>  setclosest2gcell_presetghost2d
        procedure :: setstrbodycell     =>  setstrbodycell_presetghost2d
        procedure :: setpotbodycell     =>  setpotbodycell_presetghost2d
    end type

contains

    subroutine run_presetghost2d(this,part,imageobj,unmap,periodic)
        class(presetghosts2d_t) :: this
        class(part_t),          pointer,    intent(in) :: part
        class(imageobj2d_t),    pointer,    intent(in) :: imageobj
        class(unmap2d_t),       pointer,    intent(in) :: unmap
        class(pbc_t),           pointer,    intent(in) :: periodic

        call this%init(part,imageobj,unmap,periodic)
        call this%setstrbodycell()
        call this%initgfoutletdirs()
        call this%setclosest2gcell()
        call this%setpotbodycell()

    end subroutine


    subroutine setpotbodycell_presetghost2d(this)
        class(presetghosts2d_t) :: this
        integer                 :: d,icell,ijk(1:numd)
        associate(lgf=>this%part%ghostf)

        if(this%periodic%are_allaxes()) then
            lgf%potbodycell => lgf%strbodycell
        elseif (this%periodic%non_ofaxes()) then
            lgf%potbodycell => lgf%closest2gcell
        else

        allocate(lgf%potbodycell(lgf%ist:lgf%ind))
        ! all potbodycells are streambodycells which are periodic
        lgf%potbodycell = lgf%strbodycell
        ! The none periodic nodes are pointing to closest-to-gcells
        ! corners and edges points to closest-to-gcell in the case of
        ! periodic and non-periodic intersection.

        do icell = lgf%ist, lgf%ind
            ijk = lgf%getijk(icell)
            do d = 1, numd
                if (.not.this%periodic%is_axis(d)) then
                    if (ijk(d).eq.0.or.ijk(d).eq.this%imageobj%dims(d)+1) then
                        lgf%potbodycell(icell) = lgf%closest2gcell(icell)
                    endif
                endif
            end do
        end do

        end if

        end associate
    end subroutine

    subroutine init_presetghost2d(this,part,imageobj,unmap,periodic)
        class(presetghosts2d_t) :: this
        class(part_t),          pointer,    intent(in) :: part
        class(imageobj2d_t),    pointer,    intent(in) :: imageobj
        class(unmap2d_t),       pointer,    intent(in) :: unmap
        class(pbc_t),           pointer,    intent(in) :: periodic


        this%imageobj => imageobj
        this%part => part
        this%unmap => unmap
        this%periodic => periodic

    end subroutine

    ! A handy function which receives i,j,k and uses unmap to give
    ! icell or index of cell in main list.
    integer function mapijk2icell_presetghost2d(this,ijk) result(icell)
        class(presetghosts2d_t)             :: this
        integer,                intent(in)  :: ijk(:)

        icell = this%unmap%geticell(ijk=ijk)

    end function

    ! strbodycell is an interior node (not ghost) which ghost is reperesentive of.
    ! It is used in the stream step.
    ! for periodic BC, it is in current domain. But for parallel extension strbodycell
    ! will be in another domain. strbodycell is identified by the place in mainlist.
    subroutine setstrbodycell_presetghost2d (this)
        class(presetghosts2d_t) ::  this
        integer                 ::  i,j,igcell,linkimcell,ghostimcell,ijk(1:numd)
        associate(xlen=>this%imageobj%dims(1),ylen=>this%imageobj%dims(2),img=>this%imageobj, &
            lgf=>this%part%ghostf, lgw=>this%part%ghostw)

            do j = 0, ylen+1
                do i = 0, xlen+1
                    ijk = (/i,j/)
                    if (img%is_indomain(ijk)) cycle
                    if (img%is_cellkind( ijk, "ghostw")) then
                        linkimcell = this%mapijk2icell(img%getpbodycell(ijk))
                        ghostimcell = this%mapijk2icell(ijk)
                        igcell = lgw%getiscell(ghostimcell)
                        lgw%strbodycell(igcell) = linkimcell
                    elseif(img%is_cellkind( ijk, "ghostf")) then
                        linkimcell = this%mapijk2icell(img%getpbodycell(ijk))
                        ghostimcell = this%mapijk2icell(ijk)
                        igcell = lgf%getiscell(ghostimcell)
                        lgf%strbodycell(igcell) = linkimcell
                    endif
                enddo
            enddo

        end associate
    end subroutine

    subroutine setclosest2gcell_presetghost2d(this)
        class(presetghosts2d_t)         :: this
        integer                         :: igcell,dir,opp
        integer,                pointer :: neilist(:)
        associate(gfl=>this%part%ghostf,gwl=>this%part%ghostw,img=>this%imageobj)

            do igcell = gfl%ist,gfl%ind
                dir = img%getghostsit(gfl%getijk(igcell))
                opp = opposite(dir)
                neilist => gfl%getneilist(igcell)
                gfl%closest2gcell(igcell) = neilist(opp)
            end do

            do igcell = gwl%ist,gwl%ind
                dir = img%getghostsit(gwl%getijk(igcell))
                opp = opposite(dir)
                neilist => gwl%getneilist(igcell)
                gwl%closest2gcell(igcell) = neilist(opp)
            end do

        end associate
    end subroutine

    ! outlet directions (outletdirs): each fluid ghost node has several actdirs. These are the
    ! directions that streams come into ghost from fluid cells inside of domain. And ghost most
    ! deliver them to its link (strbodycell).
    subroutine initghostfoutletdirs_presetghost2d(this)
        class(presetghosts2d_t) :: this
        integer                 :: n,q,igcell,ghostimcell,i,j
        integer                 :: dir(1:nq_1),nei(1:numd,1:nq_1),ijk(1:numd)
        associate(img=>this%imageobj,lgf=>this%part%ghostf, &
            xlen=>this%imageobj%dims(1),ylen=>this%imageobj%dims(2))

            do j = 0, ylen+1
                do i = 0, xlen+1
                    ijk = (/i,j/)
                    if (img%is_indomain(ijk)) cycle
                    if(img%is_cellkind( ijk, "ghostf")) then
                        nei = img%getnei(ijk)
                        n = 0
                        do q=1,nq_1
                            if (img%is_cellkind(nei(1:numd,q),"fluid")) then
                                n = n+1
                                dir(n) = opposite(q)
                            endif
                        enddo
                        ghostimcell = this%mapijk2icell(ijk)
                        igcell = lgf%getiscell(ghostimcell)
                        if (n.eq.0) then
                            print*, "ghost fluid cell has no fluid neighbour i,j,k,=",i,j
                            stop "Error in initghostfoutletdirs_presetghosts."
                        endif
                        allocate(lgf%outlet(igcell)%dir(1:n))
                        lgf%outlet(igcell)%dir(1:n) = dir(1:n)
                    endif
                enddo
            enddo

        end associate
    end subroutine

end module
