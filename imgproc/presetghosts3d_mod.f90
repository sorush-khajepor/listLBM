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
module presetghosts3d_mod
    use presetghosts2d_mod
    implicit none


    type, extends(presetghosts2d_t) :: presetghosts3d_t
    contains
        procedure :: initgfoutletdirs=> initghostfoutletdirs_presetghosts3d
        procedure :: setstrbodycell=> setstrbodycell_presetghosts3d
    end type

contains

    ! strbodycell is an interior node (not ghost) which ghost is reperesentive of.
    ! for periodic BC, it is in current domain. But for parallel extension strbodycell
    ! will be in another domain. strbodycell is identified by the place in mainlist.
    subroutine setstrbodycell_presetghosts3d (this)
        class(presetghosts3d_t)::this
        integer :: i,j,k,igcell,linkimcell,ghostimcell,ijk(1:numd)
        associate(xlen=>this%imageobj%dims(1),ylen=>this%imageobj%dims(2),zlen=>this%imageobj%dims(3), &
                  img=>this%imageobj, lgf=>this%part%ghostf, lgw=>this%part%ghostw)
        do k = 0, zlen+1
            do j = 0, ylen+1
                do i = 0, xlen+1
                    ijk = (/i,j,k/)
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
        enddo
        end associate
    end subroutine

    ! outlet directions (outletdirs): each fluid ghost node has several actdirs. These are the
    ! directions that streams come into ghost from inside of domain. And ghost most
    ! deliver them to its link (strbodycell).
    subroutine initghostfoutletdirs_presetghosts3d(this)
        class(presetghosts3d_t) :: this
        integer :: n,q,igcell,ghostimcell,i,j,k
        integer :: dir(1:nq_1),nei(1:numd,1:nq_1),ijk(1:numd)
        associate(img=>this%imageobj,lgf=>this%part%ghostf, &
            xlen=>this%imageobj%dims(1),ylen=>this%imageobj%dims(2), &
            zlen=>this%imageobj%dims(3))
        do k = 0, zlen+1
            do j = 0, ylen+1
                do i = 0, xlen+1
                    ijk = (/i,j,k/)
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
                            print*, "ghost fluid cell has no fluid neighbour i,j,k,=",i,j,k
                            stop "Error in initghostfoutletdirs_presetghosts."
                        endif
                        allocate(lgf%outlet(igcell)%dir(1:n))
                        lgf%outlet(igcell)%dir(1:n) = dir(1:n)
                    endif
                enddo
            enddo
        enddo
        end associate
    end subroutine

end module
