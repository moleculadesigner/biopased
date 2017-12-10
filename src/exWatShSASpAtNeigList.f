c * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * 
c ProtSASatom-exWatShell neighbor List for restr_exWatShSASpotential
c
c Yuri Vorobjev 2009
c
c READ in: natom,atomXYZ(*) - proteinSASatXYZ, rcut - rCutoff
c          natomLg,atomXYZLg(*) - ligand/or waterSH atoms 
c RESULT: 
c         nbpairList() - sotred neighbours List of ProtSASat
c                        to LigAtoms
c         startnbPairL(ia) - start in nbpairList() array
c                        for LigAtom= ia
c         nnbPairL(ia)   - number of neighbours for atom ia
c         nbpairListD2() - sorted D2 by increasing
c         nnbpLT - total length of array nbpairList()
c         nnbpLMX - MAX allowed length of array nbpairList()
c         indxSASatNeigh(iLg)=ia : indx ia of neighbor SAS atom to LigAtom=iLg
c         d2SASatNeigh(iLg) = D2 to neighbor atom ia
c
c * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * 
c
        subroutine exWatShSASpAtNeigList(ncallPL,
     &      natom,atomXYZ,rcut,
     &      natomLg,atomXYZLg,
     &      nbpairList,startnbPairL,nnbPairL,
     &      nbpairListD2,indxSASatNeigh,d2SASatNeigh,nnbpLT,nnbpLMX)
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
	real atomXYZLg(*)
	integer natomLg
        integer indxSASatNeigh(*)   !*=ia: index of nearest SASatom to ia
        real d2SASatNeigh(*)        !*=ia :d2 from ia to nearest SASatom
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
        real bb,scale1,rcut1
        integer nbb
        integer kanalp
        logical CONTROL0,CONTROL,CONTROL1
c
        kanalp=kanalRunOut
c
        CONTROL0 = .false. ! is on usually
        CONTROL  = .false. ! .false.
        CONTROL1 = .false.
c
        cellh = 2.5     !cell size
c
        rcut1=rcut
        scale1 = 1.5    ! increment for rcut1
c
1000    continue
        rcut2 = rcut1**2     ! range for pairList       
c
        if(CONTROL1)then
	write(*,*)'expLPotLgSASneigList:Start!:'
        write(kanalp,*)'expLPotLgSASneigList:Start! rcut2:',rcut2
        write(kanalp,*)'                natom:',natom
        end if !C1
c
        if(CONTROL1)then
        write(kanalp,*)'ProtSAS atoms: XYZ:'
	do ia = 1,natom
	write(kanalp,*)'ATOMp ',ia, ia, (atomXYZ(3*ia-3+k),k=1,3)
	end do!ia
c
        write(kanalp,*)'ProtSAS atoms: XYZ:'
        do ia = 1,natomLg
        write(kanalp,*)'ATOMlg',ia, ia, (atomXYZLg(3*ia-3+k),k=1,3)
        end do!ia
c
	end if!C1
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
c init  pair List for natomLg
c
        nnbPs = 0
c
        do i=1,natomLg
        startnbPairL(i)=0
        nnbPairL(i)=0
        end do!i
c
        do i =1,nnbpLMX
        nbpairList(i) = 0
        nbpairListD2(i) = 0.0
        end do!i
c
c define cell coords for proteinSAS atoms
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
        if(CONTROL1)then
        write(kanalp,*)'at i, xyz:',ia,(atomXYZ(i3+k),k=1,3)
        end if!C
c
        end do!ia
c molecule size
        if(CONTROL0)then
        write(kanalp,*)'expLPotLgSASneigList:Molsize:'
        write(kanalp,*)'expLPotLgSASL: W2: xMAX:',xMAX
        write(kanalp,*)'expLPotLgSASL: W2: xMIN:',xMIN
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
       write(kanalp,*)'expLPotLgSAS:ERROR in atomXYZ !! '
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
        write(kanalp, *)'ERROR! expLPotLgSASL:ncellMAX is low !!'
        write(kanalp,*)'expLPotLgSASL: W3: xMAX:',xMAX
        write(kanalp,*)'expLPotLgSASL: W3: xMIN:',xMIN
        write(kanalp,*)
     &	'expLPotLgSAS: molecule size is too large to be in the BOX!'
        stop
        end if!
        end do !k
        end if !CONTROL0
c
c distribute natoms atomXYZ(*) over cells
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
        write(kanalp,*)'ERROR!:expLPotLgSASL:ncell3MAX is low!!'
        stop
        end if!
c make linked list        
        linkList(ia) =  headat(ncell) 
        headat(ncell) = ia
        end do !ia
c CONTROL
        if(CONTROL1)then
        write(kanalp,*)'expLPotLgSASL: W4:headat(): ncellMX:',ncellMX
        call print10(ncellMX,headat)
        write(kanalp,*)'expLPotLgSASL: W5:linkList():'
        call print10(natom,linkList)
        end if !C1
c make neighbour list
        deltaN = int(rcut/cellh)+1
c
c loop over Lig Atoms natomLg, atomXYZLg(*)
c
        do ia = 1,natomLg         
c
        startnbPairL(ia) = nnbPs + 1
c        
        i3=3*ia-3
        xyzi(1)=atomXYZLg(i3+1)-xMIN(1)
        xyzi(2)=atomXYZLg(i3+2)-xMIN(2)
        xyzi(3)=atomXYZLg(i3+3)-xMIN(3)
c cell coord of central atom ia
        ix0 = xyzi(1)/cellh+1
        iy0 = xyzi(2)/cellh+1
        iz0 = xyzi(3)/cellh+1
c
c take all SAS atoms in the cells in the vicinity +-deltaN
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
         ja = headat(ncell)      ! ja - atom of ProtSAS
c
100      if(ja .gt. 0 )then         
        j3 = 3*(ja - 1)                                 
        xyzj(1)=atomXYZ(j3+1)-xMIN(1)
        xyzj(2)=atomXYZ(j3+2)-xMIN(2)
        xyzj(3)=atomXYZ(j3+3)-xMIN(3)
        dij2 = (xyzi(1)-xyzj(1))**2+
     &         (xyzi(2)-xyzj(2))**2 + (xyzi(3)-xyzj(3))**2  
c
cx        if( ia .ne. ja )then     
c
         if(CONTROL1)then
         write(kanalp,*)'expLPotLgSASneigList: W7:ia,ja,ncell,dij2:',
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
        write(kanalp,*)'ERROR:nonbondList: nnbpLMX is low!!'
        write(kanalp,*)
     &  'increase param: nnbpLMAX:in nbPairSAS03.h:',nnbPs
c
        write(kanalPStat,*)mError,
     &  'nonbondList: nnbpLMX is low!',
     &  '  increase param: nnbpLMAX in nbPairSAS03.h'
c
        stop
        end if !nnbPs .gt. nnbpLMX
        end if !CONTROL0
c
c nbpairList (): is symmerical pairList
c
        nbpairList(nnbPs) = ja
        nbpairListD2(nnbPs) = dij2    ! list of ij dist2
c
        end if ! pair
c
cx	end if !ia#ja   
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
        if(nnbPairL(ia).eq.0) then
        rcut1 = scale1*rcut1
        goto 1000               ! no neigbor SAS atoms are found for the ia
c                               ! increase rcut1 and repeate
        end if !nnbPairL(ia).eq.0
c
        indxSASatNeigh(ia) = nbpairList(startnbPairL(ia))
        d2SASatNeigh(ia) = nbpairListD2(startnbPairL(ia))
c
        if(CONTROL)then
        write(kanalp,*)'ia: indxSASatNeigh(ia) d2SASatNeigh(ia)=',
     &    ia,indxSASatNeigh(ia),d2SASatNeigh(ia)
        end if !C
c
        end do !ia
c
        nnbpLT = nnbPs
        ncallPL = ncallPL + 1
c
        if(CONTROL0)then
        write(kanalp,'(a22,i8,a42,3i8)')
     &  'expLPotLgSASneigList:  ncall:',ncallPL,
     &  'numbPairTot: nnbPs:',nnbpLT
        end if !C
c
        if(CONTROL)then
        write(kanalp,*)'expLPotLgSASneigList: RESULT: '
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
