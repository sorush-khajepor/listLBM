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

!========================================================
!   Lattice type
!========================================================
module lattice_mod
    use kinds_f
    use utility_mod, only: arrayofnames_t
    implicit none

    character (len=:), protected, allocatable :: latname
    ! ID number of lattice
    integer, protected :: latid
    ! numd = number of dimensions,  numq = number of Qs (like D2Q9)
	! Note: in the whole code Q or q means lattice velocity vector(s)
    integer, protected  :: nd_1, nq_1
    integer, protected  :: numd, numq
    ! number of boundary cell situations (corners and edges)
    integer, protected :: nbcsit
    ! micro-velocity (c in LBM)
    real(kind = wp) :: cvel,csqr,cqad
    ! opposits of Qs
    integer,            protected, allocatable, dimension(:)    :: opposite, halfqs
    ! Q vectors    = all lattice velocity vectors (including zero)
	! cube vectors = all lattice velocity vectors without zero
    integer,            protected, allocatable, dimension(:,:)  :: qvec,cubevec
    ! weight functions
    real(kind = wp),    protected, allocatable, dimension(:)    :: weight
    ! cvel(c)*Qvectors
    real(kind = wp),    protected, allocatable, dimension(:,:)  :: cqvec
    ! boundary cell situation name
    type(arrayofnames_t),allocatable :: bcsname(:)

    private :: numq,memalloc,setmicrovel,nd_1
!====================procedures==========================
contains

subroutine initD2Q9lattice(cvel_)
    real(kind = wp), optional :: cvel_
    latname = "D2Q9"
    latid = 1
    numd = 2
    numq = 9
    nq_1 = numq-1
    nd_1 = numd-1
    nbcsit = 8
    call memalloc()

    if (present(cvel_)) then
        call setmicrovel( cvel_)
    else
        call setmicrovel( 1.0_wp)
    endif

    qvec(0:8,1:2)  =    reshape((/0,+1, 0,-1, 0,   +1,-1,-1,+1, &          ! x = 0
                                  0, 0,+1, 0,-1,   +1,+1,-1,-1/), (/9,2/)) ! y = 1
    cubevec(1:8,1:2)  = reshape((/  +1, 0,-1, 0,   +1,-1,-1,+1, &            ! x = 0
                                     0,+1, 0,-1,   +1,+1,-1,-1/), (/8,2/))   ! y = 1
    cqvec = cvel*real(qvec,wp)

    weight(0)   = 4.0_wp/9.0_wp
    weight(1:4) = 1.0_wp/9.0_wp
    weight(5:8) = 1.0_wp/36.0_wp
!   opposite directions of {0 1 2 3 4 5 6 7 8}
    opposite(0:8) =       (/0,3,4,1,2,7,8,5,6/)

!	Half of lattice vectors which are in the same side if a lattice
!   halfed with a staight line.		
    halfqs(1:nq_1/2) = (/1,2,5,6/)

!   Setting name for boundary nodes (mem=member)
    bcsname(1)%mem = "rig"
    bcsname(2)%mem = "top"
    bcsname(3)%mem = "lef"
    bcsname(4)%mem = "bot"
    bcsname(5)%mem = "rigtop"
    bcsname(6)%mem = "leftop"
    bcsname(7)%mem = "lefbot"
    bcsname(8)%mem = "rigbot"
end

subroutine initD3Q27lattice(cvel_)
    real(kind = wp), optional :: cvel_
    latname = "D3Q27"
    latid = 4
    numd = 3
    numq = 27
    nbcsit = 26
    nq_1 = numq-1
    nd_1 = numd-1
    call memalloc()
    if (present(cvel_)) then
        call setmicrovel( cvel_)
    else
        call setmicrovel( 1.0_wp)
    endif
                              !0  1  2  3  4  5  6    7  8  9  10   11 12 13 14    15 16 17 18   19 20 21 22     23 24 25 26
    qvec(0:26,1:3) = reshape((/0,+1,-1, 0, 0, 0, 0,   +1,-1,+1,-1,  +1,-1,+1,-1,    0, 0, 0, 0,   +1,-1,+1,-1,    +1,-1,+1,-1,   &      !x 1
                               0, 0, 0,+1,-1, 0, 0,   +1,-1,-1,+1,   0, 0, 0, 0,   +1,-1,+1,-1,   +1,+1,-1,-1,    +1,+1,-1,-1,   &      !y 2
                               0, 0, 0, 0, 0,+1,-1,    0, 0, 0, 0,  +1,-1,-1,+1,   +1,-1,-1,+1,   +1,+1,+1,+1,    -1,-1,-1,-1 /),&    !z 3
                               (/27,3/))

  cubevec(1:26,1:3) = reshape((/+1,-1, 0, 0, 0, 0,   +1,-1,+1,-1,  +1,-1,+1,-1,     0, 0, 0, 0,   +1,-1,+1,-1,    +1,-1,+1,-1,   &      !x 1
                                 0, 0,+1,-1, 0, 0,   +1,-1,-1,+1,   0, 0, 0, 0,    +1,-1,+1,-1,   +1,+1,-1,-1,    +1,+1,-1,-1,   &      !y 2
                                 0, 0, 0, 0,+1,-1,    0, 0, 0, 0,  +1,-1,-1,+1,    +1,-1,-1,+1,   +1,+1,+1,+1,    -1,-1,-1,-1 /),&     !z 3
                               (/26,3/))

    cqvec = cvel*real(qvec,wp)

    weight(0)    =  8.0_wp/27.0_wp
    weight(1:6)  =  2.0_wp/27.0_wp
    weight(7:18) =  1.0_wp/54.0_wp
    weight(19:26) = 1.0_wp/216.0_wp
                     ! 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    opposite(0:26) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,26,25,24,23,22,21,20,19/)
    halfqs(1:nq_1/2) = (/1,3,5,7,9,11,13,15,17,19,20,21,22/)

    block;
    character(len=3) :: xn,yn,zn
    integer :: q
    do q=1,nbcsit
        xn="";yn="";zn=""
        if (cubevec(q,1).eq.1) then
            xn = "rig"
        elseif(cubevec(q,1).eq.-1) then
            xn = "lef"
        endif
        if (cubevec(q,2).eq.1) then
            yn = "top"
        elseif(cubevec(q,2).eq.-1) then
            yn = "bot"
        endif
        if (cubevec(q,3).eq.1) then
            zn = "fro"
        elseif(cubevec(q,3).eq.-1) then
            zn = "bac"
        endif
         bcsname(q)%mem = trim(xn)//trim(yn)//trim(zn)
    enddo
    endblock

end subroutine
subroutine initD3Q19lattice(cvel_)
    real(kind = wp), optional :: cvel_
    latname = "D3Q19"
    latid = 3
    numd = 3
    numq = 19
    nq_1 = numq-1
    nd_1 = numd-1
    nbcsit = 26
    call memalloc()
    if (present(cvel_)) then
        call setmicrovel( cvel_)
    else
        call setmicrovel( 1.0_wp)
    endif
    !                          0  1  2  3  4  5  6     7  8  9  10   11 12 13 14     15 16 17 18
    qvec(0:18,1:3) = reshape((/0,+1,-1, 0, 0, 0, 0,   +1,-1,+1,-1,  +1,-1,+1,-1,     0, 0, 0, 0,   &      !x 1
                               0, 0, 0,+1,-1, 0, 0,   +1,-1,-1,+1,   0, 0, 0, 0,    +1,-1,+1,-1,   &      !y 2
                               0, 0, 0, 0, 0,+1,-1,    0, 0, 0, 0,  +1,-1,-1,+1,    +1,-1,-1,+1/), &      !z 3
                               (/19,3/))
!                                                                                                 19 20 21 22     23 24 25 26
  cubevec(1:26,1:3) = reshape((/+1,-1, 0, 0, 0, 0,   +1,-1,+1,-1,  +1,-1,+1,-1,     0, 0, 0, 0,   +1,-1,+1,-1,    +1,-1,+1,-1,   &      !x 1
                                 0, 0,+1,-1, 0, 0,   +1,-1,-1,+1,   0, 0, 0, 0,    +1,-1,+1,-1,   +1,+1,-1,-1,    +1,+1,-1,-1,   &      !y 2
                                 0, 0, 0, 0,+1,-1,    0, 0, 0, 0,  +1,-1,-1,+1,    +1,-1,-1,+1,   +1,+1,+1,+1,    -1,-1,-1,-1 /), &     !z 3
                               (/26,3/))

    cqvec = cvel*real(qvec,wp)

    weight(0)    = 1.0_wp/3.0_wp
    weight(1:6)  = 1.0_wp/18.0_wp
    weight(7:18) = 1.0_wp/36.0_wp
                     ! 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26
    opposite(0:26) = (/0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17,26,25,24,23,22,21,20,19/)
    halfqs(1:nq_1/2) = (/1,3,5,7,9,11,13,15,17/)

    block;
    character(len=3) :: xn,yn,zn
    integer :: q
    do q=1,nbcsit
        xn="";yn="";zn=""
        if (cubevec(q,1).eq.1) then
            xn = "rig"
        elseif(cubevec(q,1).eq.-1) then
            xn = "lef"
        endif
        if (cubevec(q,2).eq.1) then
            yn = "top"
        elseif(cubevec(q,2).eq.-1) then
            yn = "bot"
        endif
        if (cubevec(q,3).eq.1) then
            zn = "fro"
        elseif(cubevec(q,3).eq.-1) then
            zn = "bac"
        endif
         bcsname(q)%mem = trim(xn)//trim(yn)//trim(zn)
    enddo
    endblock

end subroutine



subroutine initD3Q15lattice(cvel_)
    real(kind = wp), optional :: cvel_
    latname = "D3Q15"
    latid = 2
    numd = 3
    numq = 15
    nq_1 = numq-1
    nd_1 = numd-1
    nbcsit = 26
    call memalloc()
    if (present(cvel_)) then
        call setmicrovel( cvel_)
    else
        call setmicrovel( 1.0_wp)
    endif
    !                          0  1  2  3  4  5  6     7  8  9  10   11 12 13 14
    qvec(0:14,1:3) = reshape((/0,+1,-1, 0, 0, 0, 0,    +1,-1,+1,-1,  +1,-1,+1,-1,   &      !x 1
                               0, 0, 0,+1,-1, 0, 0,    +1,+1,-1,-1,  +1,+1,-1,-1,   &      !y 2
                               0, 0, 0, 0, 0,+1,-1,    +1,+1,+1,+1,  -1,-1,-1,-1 /), &     !z 3
                               (/15,3/))
!                                                                                                 19 20 21 22     23 24 25 26
  cubevec(1:26,1:3) = reshape((/+1,-1, 0, 0, 0, 0,     +1,-1,+1,-1,  +1,-1,+1,-1,   +1,-1,+1,-1,  +1,-1,+1,-1,     0, 0, 0, 0,    &      !x 1
                                 0, 0,+1,-1, 0, 0,     +1,+1,-1,-1,  +1,+1,-1,-1,   +1,-1,-1,+1,   0, 0, 0, 0,    +1,-1,+1,-1,    &      !y 2
                                 0, 0, 0, 0,+1,-1,     +1,+1,+1,+1,  -1,-1,-1,-1,    0, 0, 0, 0,  +1,-1,-1,+1,    +1,-1,-1,+1 /), &      !z 3
                               (/26,3/))

    cqvec = cvel*real(qvec,wp)

    weight(0)    = 16.0_wp/72.0_wp
    weight(1:6)  = 8.0_wp/72.0_wp
    weight(7:14) = 1.0_wp/72.0_wp
                     ! 0 1 2 3 4 5 6   7  8  9  10 11 12 13 14   15 16 17 18 19 20 21 22 23 24 25 26
    opposite(0:26) = (/0,2,1,4,3,6,5,  14,13,12,11,10, 9, 8, 7,  16,15,18,17,20,19,22,21,24,23,26,25/)
    halfqs(1:nq_1/2) = (/1,3,5,7,8,9,10/)
    block;
    character(len=3) :: xn,yn,zn
    integer :: q
    do q=1,nbcsit
        xn="";yn="";zn=""
        if (cubevec(q,1).eq.1) then
            xn = "rig"
        elseif(cubevec(q,1).eq.-1) then
            xn = "lef"
        endif
        if (cubevec(q,2).eq.1) then
            yn = "top"
        elseif(cubevec(q,2).eq.-1) then
            yn = "bot"
        endif
        if (cubevec(q,3).eq.1) then
            zn = "fro"
        elseif(cubevec(q,3).eq.-1) then
            zn = "bac"
        endif
         bcsname(q)%mem = trim(xn)//trim(yn)//trim(zn)
    enddo
    endblock

end subroutine

! Memory allocation for lattice components
subroutine memalloc()
    allocate(opposite(0:nbcsit),weight(0:nq_1))
    allocate(qvec(0:nq_1,1:numd),cqvec(0:nq_1,1:numd))
    allocate (bcsname(1:nbcsit))
    allocate(cubevec(1:nbcsit,1:numd))
    allocate(halfqs(1:nq_1/2))
end

subroutine setmicrovel(cvel_)
    real (kind = wp), intent(in) :: cvel_
    cvel = cvel_
    csqr = cvel **2
    cqad = cvel **4
end    

end module lattice_mod

