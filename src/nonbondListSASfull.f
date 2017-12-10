c * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * 
c nonbonded list calculation
c
c Yuri Vorobjev 2002, 2005 
c
c READ in:  molecPDBxyz, etc.,
c
c REULT: 
c         nbpairList() - sotred neighbours List             
c         startnbPairL(ia) - start in nbpairList() array
c         nnbPairL(ia)   - number of neighbours for aton ia
c         nbpairListD2() - sorted D2 by increasing
c All physical neighbours 12,13,14,are included in the list
c
c * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * 
c
        subroutine nonbondListFull(ncallPL,
     &          natom,atomXYZ,rcut,
     &          nbpairList,startnbPairL,nnbPairL,
     &          nbpairListD2,nnbpLT,nnbpLMX)
c
c        implicit none
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        real atomXYZ(*)
        real rcut
        integer ncallPL
        integer natom
        integer nnbpLMX,nnbpLT
cRESULT:
        integer nbpairList(*)
        integer startnbPairL(*),nnbPairL(*)
        real nbpairListD2(*)
c local variables
c linked list ! can be made global
        integer ncellMAX
        parameter (ncellMAX = 60)    !~150 A size molecule
        integer ncell3MAX
        parameter (ncell3MAX=ncellMAX*ncellMAX*ncellMAX)
c
        integer headat(ncell3MAX)
        integer natNonBondListLocMAX
        parameter (natNonBondListLocMAX=2*natomMAX)
        integer linkList(natNonBondListLocMAX)
c
        real cellh
        real rcut2
        real big,low,dij2
        real xyzi(3),xyzj(3)
        real xMAX(3),xMIN(3)
        integer nsiz(3),ncell,ncellMX
        integer i,j,k 
        integer i3,j3
        integer ia,ja,iam
        integer ix,iy,iz
        integer ixm,iym,izm
        integer ix0,iy0,iz0
        integer ni3,ni2
        integer deltaN
        integer i1s,i1f,i2s,i2f,i3s,i3f
        integer nnbPs 
        logical commute
        real bb
        integer nbb
        integer kanalp
        logical CONTROL0,CONTROL,CONTROL1
c
        kanalp=kanalRunOut
c
        CONTROL0 = .true. ! is on usually
        CONTROL  = .false. ! .false.
        CONTROL1 = .false.
c
        cellh = 2.5     !cell size
c
        rcut2 = rcut**2     ! range for pairList       
c
        if(CONTROL0)then
        write(kanalp,*)'nonbondListFL:  rcut2:',rcut2
        write(kanalp,*)'                natom:',natom
        end if
c
c initialize
c
        if(natom .gt. natomMAX)then
        write(kanalp,*)'ERROR!:nonbondListFL: natomLoc >natomMAXglob:',
     &  natom,natomMAX
        stop 
        end if
c init linkedList 
        do i = 1,natNonBondListLocMAX
        linkList(i) = 0
        end do 
c
        do i =1,ncell3MAX
        headat(i) = 0        
        end do! i
c
c init  pair List
c
        nnbPs = 0
c
        do i=1,natom
        startnbPairL(i)=0
        nnbPairL(i)=0
        end do!i
c
        do i =1,nnbpLMX
        nbpairList(i) = 0
        nbpairListD2(i) = 0.0
        end do!i
c
c define cell coords
c
        big = 10.0e10
        low = -big
        do k=1,3
        xyzi(k)=0.0
        xMAX(k)=low
        xMIN(k)=big
        end do !k
c        
        do ia=1,natom
        i3=3*ia-3
        do k=1,3
        if(xMAX(k) .lt. atomXYZ(i3+k)) xMAX(k)=atomXYZ(i3+k)
        if(xMIN(k) .gt. atomXYZ(i3+k)) xMIN(k)=atomXYZ(i3+k)
        end do!k
c
        if(CONTROL)then
        write(kanalp,*)'at i, xyz:',ia,(atomXYZ(i3+k),k=1,3)
        end if!C
c
        end do!ia
c molecule size
        if(CONTROL0)then
        write(kanalp,*)'nonbondListFL:  MolSize;' 
        write(kanalp,*)'nonbondListFL: W2: xMAX:',xMAX
        write(kanalp,*)'nonbondListFL: W2: xMIN:',xMIN
        end if
c
        ncellMX = 1
        do k=1,3
        nsiz(k) = (xMAX(k) - xMIN(k))/cellh+1 
        ncellMX = ncellMX*nsiz(k)
        end do!k
c
        if(CONTROL0)then
        if(nsiz(1) .lt. 1 .or. nsiz(2) .lt. 1 .or. nsiz(3) .lt. 1 )then
        write(kanalp,*)'nonbondListSASfull:ERROR in atomXYZ !! '
c
        write(kanalp,*)'xMAX(3):',xMAX
        write(kanalp,*)'xMIN(3):',xMIN
c
        do k=1,natom
        write(kanalp,'(a6,i5,6x,3f8.3)')
     &  'ATOM  ',k, atomXYZ(3*k-2),atomXYZ(3*k-1),atomXYZ(3*k) 
        end do !k
        stop
        end if !nsiz(1) .lt. 1
c
        do k=1,3
        if(nsiz(k) .gt. ncellMAX)then
        write(kanalp, *)'ERROR! nonbondListSASfull:ncellMAX is low !!'
        write(kanalp,*)'nonbondListFL: W3: xMAX:',xMAX
        write(kanalp,*)'nonbondListFL: W3: xMIN:',xMIN
        write(kanalp,*)'moleculeize is too large to be in the BOX!'
        stop
        end if!
        end do !k
        end if !CONTROL0
c
c distribute atoms over cells
c make linked list of atoms in cells
c headat(n) - head(incellN)
c linkList(ia) - linkedList
        ixm=1
        iym=1
        izm=1
        do ia = 1,natom
c calculate cell numb
        i3=3*ia-3
        xyzi(1)=atomXYZ(i3+1)-xMIN(1)
        xyzi(2)=atomXYZ(i3+2)-xMIN(2)
        xyzi(3)=atomXYZ(i3+3)-xMIN(3)
        ix = xyzi(1)/cellh+1
        iy = xyzi(2)/cellh+1
        iz = xyzi(3)/cellh+1
        if(ixm .lt. ix)ixm = ix
        if(iym .lt. iy)iym = iy
        if(izm .lt. iz)izm = iz
c cell number
        ncell = ix + (iy-1)*nsiz(1) + (iz-1)*nsiz(1)*nsiz(2)
        if(ncell .gt. ncell3MAX)then
        write(kanalp,*)'ERROR!:nonbondList: ncell3MAX is low !!'
        stop
        end if!
c make linked list        
        linkList(ia) =  headat(ncell) 
        headat(ncell) = ia
        end do !ia
c CONTROL
        if(CONTROL1)then
        write(kanalp,*)'nonbondListFL: W4:headat(): ncellMX:',ncellMX
        call print10(ncellMX,headat)
        write(kanalp,*)'nonbondListFL: W5:linkList():'
        call print10(natom,linkList)
        end if
c make neighbour list
        deltaN = int(rcut/cellh)+1
c
c loop over Atoms
c
        do ia = 1,natom         
c
        startnbPairL(ia) = nnbPs + 1
c        
        i3=3*ia-3
        xyzi(1)=atomXYZ(i3+1)-xMIN(1)
        xyzi(2)=atomXYZ(i3+2)-xMIN(2)
        xyzi(3)=atomXYZ(i3+3)-xMIN(3)
c cell coord of central atom ia
        ix0 = xyzi(1)/cellh+1
        iy0 = xyzi(2)/cellh+1
        iz0 = xyzi(3)/cellh+1
c
c take all atom in the cells in the vicinity +-deltaN
         i3s = iz0 - deltaN
         if(i3s .le. 1 )i3s = 1
         i3f = iz0 + deltaN
         if(i3f .gt. izm ) i3f = izm
c
         i2s = iy0 - deltaN
         if(i2s .le. 1 )i2s = 1
         i2f = iy0 + deltaN
         if(i2f .gt. iym ) i2f = iym
c
         i1s = ix0 - deltaN
         if(i1s .le. 1 )i1s = 1
         i1f = ix0 + deltaN
         if(i1f .gt. ixm ) i1f = ixm
c
c loop over lattice neighbours
        do iz = i3s,i3f
         ni3 = (iz-1)*nsiz(1)*nsiz(2)
        do iy = i2s,i2f      
         ni2 = (iy-1)*nsiz(1)
        do ix = i1s,i1f      
c cell number
         ncell = ix + ni2 + ni3
c extract ja atoms from cell ncell
         ja = headat(ncell)
c
100      if(ja .gt. 0 )then         
        j3 = 3*(ja - 1)                                 
        xyzj(1)=atomXYZ(j3+1)-xMIN(1)
        xyzj(2)=atomXYZ(j3+2)-xMIN(2)
        xyzj(3)=atomXYZ(j3+3)-xMIN(3)
        dij2 = (xyzi(1)-xyzj(1))**2+
     &         (xyzi(2)-xyzj(2))**2 + (xyzi(3)-xyzj(3))**2  
c
        if( ia .ne. ja )then     
c
         if(CONTROL1)then
         write(kanalp,*)'nonbondListFL: W7:ia,ja,ncell,dij2:',
     &        ia,ja,ncell,dij2 
         end if      
c  
c nbpairList: includes all 12,13,14 bonded pairs
       if( dij2 .lt. rcut2)then
        nnbPairL(ia) = nnbPairL(ia) + 1
        nnbPs = nnbPs + 1
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondList:W9: ia,ja,nnbPs, dij2:',
     &  ia,ja,nnbPs, dij2
        end if !C1
c
        if(CONTROL0)then
        if(nnbPs .gt. nnbpLMX) then
        write(kanalp,*)
     &	'ERROR:nonbondList: nnbpLMX is low!!=',nnbpLMX
        write(kanalp,*)
     &  'increase param: nnbpLMAX:in nbPairSAS03.h:',nnbPs
c
        write(kanalPStat,*)mError,
     &  'nonbondList: nnbpLMX is low!=',nnbpLMX,
     &  '  increase param: nnbpLMAX in nbPairSAS03.h'
c
        stop
        end if
        end if !CONTROL0
c
c nbpairList (): is symmerical pairList
c
        nbpairList(nnbPs) = ja
        nbpairListD2(nnbPs) = dij2    ! list of ij dist2
c
        end if ! pair
c
	end if !ia#ja   
c
        else
        goto 101  ! goto new cell
        end if ! ja >0
c 
c next ja
102    ja = linkList(ja)
        goto 100
c
101     continue
        end do !iz
        end do !iy
        end do !ix
c        
c list for atom ia is complete
c
c sorting of the list for ia by i-j distances
c
300     commute = .false.
c
        do j = startnbPairL(ia),nnbPs-1
        if(nbpairListD2(j) .gt. nbpairListD2(j+1) )then
        commute = .true.
        bb = nbpairListD2(j)
        nbpairListD2(j) = nbpairListD2(j+1)
        nbpairListD2(j+1) = bb
        nbb = nbpairList(j)
        nbpairList(j) = nbpairList(j+1)
        nbpairList(j+1) = nbb
        end if
        end do !j
c
        if(commute)goto 300
c
c after this line the portion of nbpairListDj2() and nbpairList() 
c related to the ia atom are sorted
c
        end do !ia
c
        nnbpLT = nnbPs
        ncallPL = ncallPL + 1
c
        if(CONTROL0)then
        write(kanalp,'(a22,i8,a42,3i8)')
     &  'nonbondListFL:  ncall:',ncallPL,
     &  'numbPairTot: nnbPs:',nnbpLT
        end if !C
c
        if(CONTROL)then
        write(kanalp,*)'nonbondListFL: RESULT: '
        write(kanalp,*)'startnbPairL():'
        call print10(natom,startnbPairL)
        write(kanalp,*)'nnbPairL():'
        call print10(natom,nnbPairL)
        write(kanalp,*)'nbpairList():nnbPs:',nnbPs 
        call print10(nnbPs,nbPairList)
        write(kanalp,*)'nbpairListD2():nnbPs:',nnbPs
        call print10r(kanalp,nnbPs,nbPairListD2)
        end if !CONTROL
c
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB

            return
            end
