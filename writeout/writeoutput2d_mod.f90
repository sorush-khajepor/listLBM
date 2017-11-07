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

! This module uses FortranVTK library to write 2D unstructured vtk files which can
! be visualized in ParaView. 
module writeoutput2d_mod
    use kinds_f
    use part_mod
    implicit none

    type output2d_t
        class(part_t),   pointer  :: part
        integer :: filecounter = 0
        integer  :: tstep_old = -1
    contains
        procedure :: init01 => init01_output2d
        generic   :: init => init01
        procedure :: getfilenum => getfilenumber_output2d
        procedure :: writeVTK => writeVTK_output2d
    end type


contains

subroutine init01_output2d(this,part)
    class(output2d_t) :: this
    class(part_t),   pointer,intent(in) :: part
    this%part=>part
end subroutine


function getfilenumber_output2d(this,tstep) result(filenumber)
    class(output2d_t) :: this
    integer,intent(in) :: tstep
    character (LEN=6) :: filenumber
    Integer:: i,tmp_fnum

    if (tstep.eq.this%tstep_old) this%filecounter = this%filecounter - 1
    filenumber = "000000"
    tmp_fnum = this%filecounter
    do i =5, 0, -1
        write(filenumber((6-i):(6-i)),'(i1)') tmp_fnum/(10**i)
        tmp_fnum=mod(tmp_fnum,10**i)
    end do
    this%filecounter = this%filecounter + 1
    this%tstep_old = tstep
end function

! A subroutine to save data in a vtk file in unstructured cube cells
! Note that it needs to save x,y,z (i,j,k), which are usually less
! than 1000, in 16 bits integer.
! Note that here cell is a vtk standard cell and node is the a data
! cell in the mainlist.
subroutine writeVTK_output2d(this,tstep)
    use Lib_VTK_IO
    use, intrinsic:: ISO_FORTRAN_ENV, only: stdout=>OUTPUT_UNIT, stderr=>ERROR_UNIT
    class(output2d_t) :: this
    integer,Intent(IN):: tstep
    real (kind = R4P), allocatable, dimension(:) :: X,Y,Z
    real   (kind = R8P), allocatable, dimension(:) :: rho,vxx,vyy,vzz
    integer(kind = I1P), allocatable, dimension(:) :: cell_type
    integer,             allocatable, dimension(:) :: offset,connect,icell2inode
    real(kind = R8P) :: rho_,u_(1:numd)
    integer(kind=I4P):: Nc, ii,icell,nn, e_io,inode,nmembers,member,q,n,isnode
    integer :: ivtknode(this%part%main%ist:this%part%main%ind),ijk(1:numd)

    associate(ml=>this%part%main,part=>this%part,&
              xlen=>this%part%main%xlen,ylen=>this%part%main%ylen,&
              neilist=>this%part%main%neilist)

    write(*,'(A25)',advance='no') " Writing VTK file ...     "

    ! 3-dimensions of the domain, All points
    NN= 0
    do n=lbound(part%fluids,dim=1),ubound(part%fluids,dim=1)
        NN = NN + part%fluids(n)%pt%getnumocells()
    enddo

    do n=lbound(part%walls,dim=1),ubound(part%walls,dim=1)
        NN = NN + part%walls(n)%pt%getnumocells()
    end do

    !number of cells
    Nc = 1
    Allocate (X(1:nn),Y(1:nn),Z(1:nn))
    Allocate (vxx(1:nn),vyy(1:nn),vzz(1:nn),rho(1:nn))
    allocate (icell2inode(1:nn))
    vzz=0.0_wp
    ivtknode=0
    ii = 0

    do n=lbound(part%fluids,dim=1),ubound(part%fluids,dim=1)
        do inode = part%fluids(n)%pt%istinmlist,part%fluids(n)%pt%indinmlist
            ii = ii + 1
            ! index of points starts from zero in the vtk file.
            ivtknode(inode)=ii-1
            isnode = part%fluids(n)%pt%getiscell(inode)
            rho_ = part%fluids(n)%pt%getrho(isnode)
            if (rho_.gt.eps_wp) then
                u_ =part%fluids(n)%pt%getvel(isnode)
            else
                u_ = 0._wp
            endif
            rho(ii) = rho_
            vxx(ii) = u_(1)
            vyy(ii) = u_(2)
            ijk = ml%getijk(inode)
            X(ii)=int(ijk(1),i2p)
            Y(ii)=int(ijk(2),i2p)
            Z(ii)=1_i2p
        end do
    enddo

    do n=lbound(part%walls,dim=1),ubound(part%walls,dim=1)
        do inode = part%walls(n)%pt%istinmlist,part%walls(n)%pt%indinmlist
        ii = ii + 1
        ivtknode(inode)=ii-1
        isnode = part%walls(n)%pt%getiscell(inode)
        rho_ =   part%walls(n)%pt%getrho(isnode)
        u_ =     part%walls(n)%pt%getvel(isnode)
        rho(ii) = rho_
        vxx(ii) = u_(1)
        vyy(ii) = u_(2)
        ijk = ml%getijk(inode)
        X(ii)=int(ijk(1),i2p)
        Y(ii)=int(ijk(2),i2p)
        Z(ii)=1_i2p
        end do
    end do

  !  NN = mxcelfw  ! 3 dimensions of the domain, All points
    nmembers = 4  ! number of members (points) of a cell

    icell2inode = 0
    icell = 0

 ! 1st setp: Finding nodes which can form a cell:
 !      by forming a cell I mean that node is placed in left bottom a cell
 !      and other neighbours (rig, rig-up, up) can be find with the aid of that node.
 !      Therefore, i=xlen and j=ylen are not able to form a cell.
 !      Moreover, wall nodes may have black wall in neighbour (rig,rig-up,up), so these
 !      nodes cann't form a cell.
 !      every node form a cell right and up of itself. So, the last xlen or ylen cells are not able.
 !      rig = 1, up = 2, rig-up=5

    do n=lbound(part%fluids,dim=1),ubound(part%fluids,dim=1)
        do inode = part%fluids(n)%pt%istinmlist,part%fluids(n)%pt%indinmlist
            ijk = ml%getijk(inode)
            if (ijk(1).eq.xlen .or. ijk(2).eq.ylen) cycle
            ! Counting cells which are eligible.
            icell = icell + 1
            ! A funtion to relate number of cell to number of forming node.
            icell2inode(icell) = inode
        enddo
    enddo

    do n=lbound(part%walls,dim=1),ubound(part%walls,dim=1)
        do inode = part%walls(n)%pt%istinmlist,part%walls(n)%pt%indinmlist
            ijk = ml%getijk(inode)
            if (ijk(1).eq.xlen .or. ijk(2).eq.ylen) cycle
            if ( neilist(1,inode).eq.0  .or.  &
                 neilist(5,inode).eq.0  .or.  &
                 neilist(2,inode).eq.0) cycle
            ! Counting cells which are eligible.
            icell = icell + 1
            ! A funtion to relate number of cell to number of forming node.
            icell2inode(icell) = inode
        enddo
    end do


    ! Number of cells
    Nc = icell
    Allocate (cell_type(1:Nc),connect(1:Nc*nmembers),offset(1:Nc))
    cell_type = 8
    do icell = 1, Nc
        offset(icell) = nmembers*icell
        do Q = 1, nmembers
            inode = icell2inode(icell)
            if (Q.eq.1) member = ivtknode(inode)
            if (Q.eq.2) member = ivtknode(neilist(1,inode))
            if (Q.eq.3) member = ivtknode(neilist(2,inode))
            if (Q.eq.4) member = ivtknode(neilist(5,inode))

            connect( (icell-1)*nmembers + Q) = member

           ! Q = 0 node itself
           ! Q = 1 Right neighbor of node neilist 1
           ! Q = 2 Up neilist 2
           ! Q = 3 Right-Up neilist 5
        enddo
    enddo

   E_IO = VTK_INI_XML(output_format = 'binary',filename ='vtk/XML_UNST-binary'//trim(this%getfilenum(tstep))//'.vtu' &
   ,mesh_topology='UnstructuredGrid')
   E_IO = VTK_FLD_XML(fld_action='open')
!  E_IO = VTK_FLD_XML(fld=0._R8P,fname='TIME')
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


end module
