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

! This module uses FortranVTK library to write 3D unstructured vtk files which can
! be visualized in ParaView. 
module writeoutput3d_mod
    use writeoutput2d_mod
    implicit none

    type, extends(output2d_t) :: output3d_t
    contains
        procedure :: writeVTK       => writeVTK_output3d
        procedure :: does_form_a_cell => does_form_a_cell_output3d
    end type

contains


subroutine writeVTK_output3d(this,tstep)
    use Lib_VTK_IO
    class(output3d_t) :: this
    integer,Intent(IN):: tstep
    real   (kind = R4P), allocatable, dimension(:) :: X,Y,Z
    real   (kind = R8P), allocatable, dimension(:) :: rho,vxx,vyy,vzz
    integer(kind = I1P), allocatable, dimension(:) :: cell_type
    integer,             allocatable, dimension(:) :: offset,connect,icell2inode
    real(kind = R8P) :: rho_,u_(1:numd)
    integer(kind=I4P):: Nc, ii,icell,nn, e_io,inode,nmembers,member,q,n,isnode
    integer :: ivtknode(this%part%main%ist:this%part%main%ind),ijk(1:numd)
    associate(ml=>this%part%main,part=>this%part,&
              xlen=>this%part%main%xlen,ylen=>this%part%main%ylen,&
              zlen=>this%part%main%zlen,neilist=>this%part%main%neilist)

    write(*,'(A25)',advance='no') " Writing VTK file ...     "

    ! 3-dimensions of the domain, All nodes
    NN= 0
    do n=lbound(part%fluids,dim=1),ubound(part%fluids,dim=1)
        NN = NN + part%fluids(n)%pt%getnumocells()
    end do

    do n=lbound(part%walls,dim=1),ubound(part%walls,dim=1)
        NN = NN + part%walls(n)%pt%getnumocells()
    end do


    !number of cells (vtk cells)
    Nc = 1
    allocate (X(1:nn),Y(1:nn),Z(1:nn))
    allocate (vxx(1:nn),vyy(1:nn),vzz(1:nn),rho(1:nn))
    allocate (icell2inode(1:nn))

    ivtknode = 0
    ii = 0
    do n=lbound(part%fluids,dim=1),ubound(part%fluids,dim=1)
        do inode = part%fluids(n)%pt%istinmlist,part%fluids(n)%pt%indinmlist
            ii = ii + 1
            ivtknode(inode) = ii
            isnode = part%fluids(n)%pt%getiscell(inode)
            rho_ = part%fluids(n)%pt%getrho(isnode)
            if (rho_.gt.eps_wp) then
                u_ =  part%fluids(n)%pt%getvel(isnode)
            else
                u_ = 0._wp
            endif
            rho(ii) = rho_
            vxx(ii) = u_(1)
            vyy(ii) = u_(2)
            vzz(ii) = u_(3)
            ijk = ml%getijk(inode)
            X(ii)=int(ijk(1),i2p)
            Y(ii)=int(ijk(2),i2p)
            Z(ii)=int(ijk(3),i2p)
        end do
    enddo

    do n=lbound(part%walls,dim=1),ubound(part%walls,dim=1)
        do inode = part%walls(n)%pt%istinmlist,part%walls(n)%pt%indinmlist
        ii = ii + 1
        ivtknode(inode) = ii
        isnode = part%walls(n)%pt%getiscell(inode)
        rho_ =   part%walls(n)%pt%getrho(isnode)
        u_ =     part%walls(n)%pt%getvel(isnode)
        rho(ii) = rho_
        vxx(ii) = u_(1)
        vyy(ii) = u_(2)
        vzz(ii) = u_(3)
        ijk = ml%getijk(inode)
        X(ii)=int(ijk(1),i2p)
        Y(ii)=int(ijk(2),i2p)
        Z(ii)=int(ijk(3),i2p)
    end do
    end do

    ! number of members (nodes) of a cell
    nmembers = 8

    icell2inode = 0
    icell = 0

 ! 1st setp: Finding nodes which can form a cell:
 !      by forming a cell I mean that node is placed in left bottom
 !      back corner of a cell and other neighbours (rig, rig-up, up,
 !      ...) can be find with the aid of that node.
 !      Therefore, i=xlen and j=ylen are not able to form a cell.
 !      Moreover, wall nodes may have black wall in neighbour (rig,
 !      rig-up,up), so these nodes cann't form a cell.
 !      every node form a cell right and up of itself. So, the last
 !      xlen, ylen or zlen cells are not eligible.

    do n=lbound(part%fluids,dim=1),ubound(part%fluids,dim=1)
        do inode = part%fluids(n)%pt%istinmlist,part%fluids(n)%pt%indinmlist
            ijk = ml%getijk(inode)
            if (ijk(1).eq.xlen .or. ijk(2).eq.ylen.or. ijk(3).eq.zlen) cycle
            if (this%does_form_a_cell(inode)) cycle
            ! Counting cells which are eligible
            icell = icell + 1
            ! A funtion to relate number of cell to number of forming node.
            icell2inode(icell) = inode
        enddo
    enddo
    ! Neighbours of a node in d3q19 or d3q15 some of diagonal neighbours.
    ! For the sake of generalization, neighbours are found with the aid
    ! of only unique vectors of neilist.
    do n=lbound(part%walls,dim=1),ubound(part%walls,dim=1)
        do inode = part%walls(n)%pt%istinmlist,part%walls(n)%pt%indinmlist
            ijk = ml%getijk(inode)
            if (ijk(1).eq.xlen .or. ijk(2).eq.ylen.or. ijk(3).eq.zlen) cycle
            if (this%does_form_a_cell(inode)) cycle
            ! Counting cells which are eligible
            icell = icell + 1
            ! A funtion to relate number of cell to number of forming node.
            icell2inode(icell) = inode
        enddo
    end do
    ! Number of cells
    Nc = icell
    Allocate (cell_type(1:Nc),connect(1:Nc*nmembers),offset(1:Nc))

    cell_type =  11

    do icell = 1, Nc
        offset(icell) = nmembers*icell
        do Q = 1, nmembers
            inode = icell2inode(icell)
            if (Q.eq.1) member = ivtknode(inode)
            if (Q.eq.2) member = ivtknode(neilist(1,inode))
            if (Q.eq.3) member = ivtknode(neilist(3,inode))
            if (Q.eq.4) member = ivtknode(neilist(1,neilist(3,inode)))
            if (Q.eq.5) member = ivtknode(neilist(5,inode))
            if (Q.eq.6) member = ivtknode(neilist(1,neilist(5,inode)))
            if (Q.eq.7) member = ivtknode(neilist(3,neilist(5,inode)))
            if (Q.eq.8) member = ivtknode(neilist(1,neilist(3,neilist(5,inode))))

            ! number of points starts from zero in connects, but
            ! coordination indexes (for x,y,z) in our code from 1.
            member = member-1
            connect( (icell-1)*nmembers + Q) = member
        enddo
    enddo
    E_IO = VTK_INI_XML(output_format = 'binary',filename ='vtk/XML_UNST-binary'//trim(this%getfilenum(tstep))//'.vtu' &
    ,mesh_topology='UnstructuredGrid')
    E_IO = VTK_FLD_XML(fld_action='open')
    ! E_IO = VTK_FLD_XML(fld=0._R8P,fname='TIME')
    E_IO = VTK_FLD_XML(fld=tstep,fname='Time Step')
    E_IO = VTK_FLD_XML(fld_action='close')
    E_IO = VTK_GEO_XML(NN = Nn, NC = Nc, X = x, Y = y, Z = z)
    E_IO = VTK_CON_XML(NC = Nc, connect = connect, offset = offset, cell_type = cell_type )
    E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'open')
    E_IO = VTK_VAR_XML(NC_NN = Nn, varname = 'Density', var = real(rho(1:nn),R4P))
    E_IO = VTK_VAR_XML(NC_NN = Nn, varname = 'Velocity', varX=real(vxx(1:nn),R4P),varY=real(vyy(1:nn),R4P),varZ=real(vzz,R4P))
    E_IO = VTK_DAT_XML(var_location = 'node', var_block_action = 'close')
    E_IO = VTK_GEO_XML()
    E_IO = VTK_END_XML()

    Deallocate (X,Y,Z)
    Deallocate(cell_type,connect,offset,icell2inode)
    Deallocate (vxx,vyy,vzz,rho)
    write(*,*) "done."
    end associate
end subroutine

! Finds the nodes which have black-node nearby and cannot form a vtk
! cell.
function does_form_a_cell_output3d(this,inode) result(is)
    class(output3d_t) :: this
    integer,intent(in) :: inode
    logical :: is
    associate(neilist=>this%part%main%neilist)
    if ( neilist(1,inode).eq.0  .or.  &                     !(+1, 0, 0)
         neilist(3,inode).eq.0  .or.  &                     !( 0,+1, 0)
         neilist(5,inode).eq.0  .or.  &                     !( 0, 0,+1)
         neilist(1,neilist(3,inode)).eq.0 .or.  &           !(+1,+1, 0)
         neilist(1,neilist(5,inode)).eq.0 .or.  &           !(+1, 0,+1)
         neilist(3,neilist(5,inode)).eq.0 .or.  &           !( 0,+1,+1)
         neilist(1,neilist(3,neilist(5,inode))).eq.0) then  !(+1,+1,+1)
         is = .true.
     else
         is = .false.
     endif

    end associate
end function
end module
