c * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * 
c nonbonded list calculation
c
c Yuri Vorobjev 2002 
c
c READ in:  moveAtomList, molecPDBxyz, etc.,
c           vbondpair 12, 13, 14 info
c makeVdW: 0/1 VDW(coul) R1 range pairList, 
c makeCL : 0/1 COUL R2 range List   
c 12, 13,14 neigbours are excluded from VDW R1 and COUL R2 pairLists
c
c makeS  : 0/1 Solvation GS model pairList, includes all 12,13,14 
c
c  nbPairListV(ia) is calculated for the central ia atoms from the moveAtomList()
c over all atomXYZ
c
c RESULT: nbpairListV()-for VanderWaals force calculation and COUL in R1 range
c         startnbPairLV(ia) - ks, point to the startPosition of 
c                             the sequense for atom ia 
c         in the nbpairListV(), the nbpairListV(ks),(ks+1),.., 
c                                          - atom N of neighbours 
c         nnbPairLV(ia) - number of neighbours
c 12,13,14 - pairs are excluded, 
c ia < ja - nonsymmetrical pairList V,C
c
c ---------------------------------------------------------------------
c         nbpairListS() - neighbours for GS solvation model             
c         startnbPairLS()
c         nnbPairLS()   - 
c All physical neighbours 12,13,14,are included in the list
c
c --------------------------------------------------------------------
c         nbpairListC()-for Coulombic force calculation in R1-R2 range
c         startnbPairLC()
c         nnbPairLC() - 
c 12,13,14 - pairs are excluded
c 
c * * * * * * * * * * * * * * * * * * * *  * * * * * * * * * * * * * * * 
c
        subroutine nonbondListVCS(ncallPL,rcutV,rcutC,rcutS,
     &          atomXYZ,atomQ,
     &          rbuffV,rbuffC,rbuffS,
     &          makeVdW,makeCL,makeS,
     &          natom,atomNumb,atomName,resName,chName,resNumb,
     &          nres,resNameRes,atomBlockName,
     &          atomNameEx,startAtInRes,
     &          nmoveatom,moveAtomList,moveFlag,
     &          pair12List,startPairL12,nPairL12,
     &          pair13List,startPairL13,nPairL13,
     &          pair14List,startPairL14,nPairL14,
     &          nbpairListV,startnbPairLV,nnbPairLV,nnbpLVT,nnbpLVMAX,
     &          nbpairListC,startnbPairLC,nnbPairLC,nnbpLCT,nnbpLCMAX,
     &          nbpairListS,startnbPairLS,nnbPairLS,nnbpLST,nnbpLSMAX)
c
c        implicit none
        include 'xyzPDBsize.h'
c           
        real atomXYZ(*)
        real atomQ(*)
        real rcutV,rcutC,rcutS
        real rbuffV,rbuffC,rbuffS
        integer ncallPL
        integer makeVdW,makeCL,makeS
        integer natom
        character*4 atomName(*)
        character*4 resName(*)
        character*1 chName(*)
        integer atomNumb(*),resnumb(*)
        integer nres
        character*4 resNameRes(*)
c        character*1 chNameRes(*) 
        character*4 atomBlockName(*)
        character*8 atomNameEx(*)
        integer startAtInRes(*)
        integer nmoveatom,moveAtomList(*)
        integer moveFlag(*)
        integer nnbpLVMAX,nnbpLVT
        integer nnbpLCMAX,nnbpLCT
        integer nnbpLSMAX,nnbpLST
        integer pair12List(*),startPairL12(*),nPairL12(*)
        integer pair13List(*),startPairL13(*),nPairL13(*)
        integer pair14List(*),startPairL14(*),nPairL14(*)
cRESULT:
        integer nbpairListV(*)
        integer startnbPairLV(*),nnbPairLV(*)
        integer nbpairListC(*)
        integer startnbPairLC(*),nnbPairLC(*)
        integer nbpairListS(*)
        integer startnbPairLS(*),nnbPairLS(*)
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
c local variables
c linked list ! can be made global
        integer ncellMAX
        parameter (ncellMAX = 60)    !~150A size molecule
        integer ncell3MAX
        parameter (ncell3MAX=ncellMAX*ncellMAX*ncellMAX)
        integer natomLMAX
        parameter (natomLMAX = natomMAX)
        integer headat(ncell3MAX)
        integer linkList(natomLMAX)
c
        real cellh
        real rcutV2,rcutC2p,rcutS2
        real rcutV2m,rcutVMIN
        real big,low,dij2
        real xyzi(3),xyzj(3)
        real xMAX(3),xMIN(3)
        real r122max,r132max,r142max
        integer nsiz(3),ncell,ncellMX
        integer i,j,k 
        integer i3,j3
        integer ia,ja,iam
        integer ix,iy,iz
        integer ixm,iym,izm
        integer ix0,iy0,iz0
        integer ni3,ni2
        integer deltaN,deltaNV,deltaNC,deltaNS
        integer i1s,i1f,i2s,i2f,i3s,i3f
        integer nnbPsV,nnbPsC,nnbPsS 
        logical OPT_removePairInHCycBase
        integer kanalp
        logical CONTROL0,CONTROL,CONTROL1
        logical CONTROLS
        logical find12,find13,find14
c
        if(makeVDW .eq. 0 .and. makeS .eq. 0 
     &                  .and. makeCL .eq. 0 )return
c
c check compatibility
        if(makeCL .eq. 1) then
         if(makeVDW .eq. 0) makeVDW = 1
        end if !makeCL
c
        kanalp=kanalRunOut 
c
        CONTROL0 = .true.  ! is on usually
c C,C1,CS : regular = .false.
        CONTROL  = .false.
        CONTROL1 = .false.
        CONTROLS = .false.
        OPT_removePairInHCycBase = .true. ! atoms donot interact VDW if belongs to
c                                           the same block
c
        cellh = 2.5     !cell size
        r122max = 2.55**2  !max 12bond
        r132max = 3.55**2  !max 13dist
        r142max = 4.50**2  !max 14 dist
c
        
        rcutV2 = (rcutV + rbuffV)**2     ! range for List1
        rcutV2m = (rcutV - rbuffC)**2    ! range for List2
        rcutC2p = (rcutC + rbuffC)**2    ! range for List2
        rcutS2 = (rcutS + rbuffS)**2     ! range for SolvationGSList
c
        if(CONTROL)then
        write(kanalp,*)'nonbondListVCS: makeVdW,makeCL,makeS:',
     &  makeVdW,makeCL,makeS,':rcutV2,rcutV2m,rcutC2p,rcutS2:',
     &  rcutV2,rcutV2m,rcutC2p,rcutS2
        end if
c
c initialize
        if(natom .gt. natomLMAX)then
        write(kanalp,*)'ERROR!:nonbondListVCS:natomLMAX is low!!'
        stop 
        end if
c init linkedList 
        do i = 1,natomLMAX
        linkList(i) = 0
        end do 
        do i =1,ncell3MAX
        headat(i) = 0        
        end do! i
c
c init VDW/COUL  R1 pair List
        if(makeVDW .eq. 1)then
        nnbPsV = 0
c
        do i=1,natom
        startnbPairLV(i)=0
        nnbPairLV(i)=0
        end do!i
c
        do i =1,nnbpLVMAX
        nbpairListV(i) = 0
        end do!i
c
        end if !makeVDW
c
c init COULR2 pairList 
        if (makeCL .eq. 1) then
        nnbPsC = 0
c
        do i=1,natom
        startnbPairLC(i)=0
        nnbPairLC(i)=0
        end do!i
c
        do i =1,nnbpLCMAX
        nbpairListC(i) = 0
        end do!i
        end if  !makeCL
c
c init GS solvation model pair List
        if (makeS .eq. 1 )then
        nnbPsS = 0
c
        do i=1,natom
        startnbPairLS(i)=0
        nnbPairLS(i)=0
        end do!i
c
        do i =1,nnbpLSMAX
        nbpairListS(i) = 0
        end do!i
        end if ! SL
c
        if(CONTROL)then
        write(kanalp,*)'nonbondListVCS: initializationFinish:'
        end if !C

c define cell coords
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
c        if(CONTROL)then
c        write(kanalp,*)'at i, xyz:',ia,(atomXYZ(i3+k),k=1,3)
c        end if!C
        end do!ia
c molecule size
        if(CONTROL)then
        write(kanalp,*)'nonbondListVCS: MolSize;' 
        write(kanalp,*)'nonbondListVCS:W2: xMAX:',xMAX
        write(kanalp,*)'nonbondListVCS:W2: xMIN:',xMIN
        end if
c
        ncellMX = 1
        do k=1,3
        nsiz(k) = (xMAX(k) - xMIN(k))/cellh+1 
        ncellMX = ncellMX*nsiz(k)
        end do!k
c
        if(CONTROL0)then
c
        if(nsiz(1) .lt. 1 .or. nsiz(2) .lt. 1 .or. nsiz(3) .lt. 1 )then
        write(kanalp,*)'ERROR1:nonbondListVCS: in atomXYZ !! '
        do k=1,natom
        write(kanalp,'(a6,i5,6x,3f8.3)')
     &  'ATOM  ',k, atomXYZ(3*k-2),atomXYZ(3*k-1),atomXYZ(3*k) 
        end do !k
c
        write(kanalPStat,*)mError, 
     &  ' molecule XYZ box is too large ... ',
     &  ' molecule is unstable ...'
c
        stop
        end if
c
        do k=1,3
        if(nsiz(k) .gt. ncellMAX)then
        write(kanalp, *)'ERROR! nonbondList:ncellMAX is low !!'
        write(kanalp,*)'nonbondListVCS:W3: xMAX:',xMAX
        write(kanalp,*)'nonbondListVCS:W3: xMIN:',xMIN
        write(kanalp,*)'moleculeSize is too large to be in the BOX!',
     &  ' increase ncellMAX in the nonbondListVCS_BS.f !'
c
        write(kanalPStat,*)mError, 
     &  ' moleculeSize is too large to be in the BOX!',
     &  'increase ncellMAX in the nonbondListVCS_BS.f !'
c
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
        write(kanalp,*)'ERROR2!:nonbondList: ncell3MAX is low !!'
        write(kanalp,*)'... increase ncell3MAX in nonbondListVCS_BS.f'
        stop
        end if!
c make linked list        
        linkList(ia) =  headat(ncell) 
        headat(ncell) = ia
        end do !ia
c CONTROL
        if(CONTROL1)then
        write(kanalp,*)'nonbondListVCS:W4:headat(): ncellMX:',ncellMX
        call print10(ncellMX,headat)
        write(kanalp,*)'nonbondListVCS:W5:linkList():'
        call print10(natom,linkList)
        end if
c
c make neighbour list
        rcutVMIN = 6.0           !minimal rcutVdW
c
        if( rcutV .le. rcutVMIN) rcutV = rcutVMIN
        if( rcutC .le. rcutV ) rcutC = rcutV
        if( rcutS .le. rcutV ) rcutS = rcutV
        deltaNS = int((rcutS+rbuffS)/cellh)+1
        deltaNV = int((rcutV+rbuffV)/cellh)+1
        deltaNC = int((rcutC+rbuffC)/cellh)+1
c
        deltaN = 1
        if(makeVdW .eq. 1 ) deltaN = deltaNV 
        if(makeS .eq. 1 )then
        if(deltaN .lt. deltaNS) deltaN = deltaNS
        end if
        if(makeCL .eq. 1 )then
        if(deltaN .lt. deltaNC) deltaN = deltaNC
        end if 
c
c loop over movingAtomList
c
c nnbPairLV(ia) is calculated for the central ia from the moveAtomList()
c over all atomXYZ 
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondListVCS: moveAtomList: iam ia:'
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        write(kanalp,'(25x,2i5)') iam, ia
        end do !
        end if !CONTROL
c
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
c
        if(makeVdW .eq. 1 )then 
        startnbPairLV(ia) = nnbPsV + 1
        end if
c
        if(makeCL .eq. 1 )then  
        startnbPairLC(ia) = nnbPsC + 1
        end if
c
        if( makeS .eq. 1)then 
        startnbPairLS(ia) = nnbPsS + 1
        end if ! makeS
c        
        i3=3*ia-3
        xyzi(1)=atomXYZ(i3+1)-xMIN(1)
        xyzi(2)=atomXYZ(i3+2)-xMIN(2)
        xyzi(3)=atomXYZ(i3+3)-xMIN(3)
c cell coord if central atom
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
         write(kanalp,*)'nonbondListVCS:W7:ia,ja,ncell,dij2:',
     &        ia,ja,ncell,dij2 
         end if      
c  
c nbpairListS: includes all 12,13,14 pairs
       if( makeS .eq. 1. and. dij2 .lt. rcutS2)then
        nnbPairLS(ia) = nnbPairLS(ia) + 1
        nnbPsS = nnbPsS + 1
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondListS:W9: ia,ja,nnbPsS, dij2:',
     &  ia,ja,nnbPsS, dij2
        end if
c
        if(CONTROL0)then
        if(nnbPsS .gt. nnbpLSMAX) then
        write(kanalp,*)'ERROR:nonbondListCVS: nnbpLSMAX is low!!'
        write(kanalp,*)'increase param: nnbpLSMAX:',nnbPsS
        stop
        end if
        end if !CONTROL0
c
        nbpairListS(nnbPsS) = ja
        end if ! makeS .eq. 1 ! pairS
c
c exclude 12, 13,14 neigb for VDW and COUL pairList
c
        if(dij2 .le. r122max) then
c
        call compAtomPairList(ja,ia,
     &             pair12List,startPairL12,npairL12,find12)
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondListVCS:W8:find12:',find12
        end if
c
        if(find12) goto 102
        end if !12
c
       if ( dij2 .le. r132max) then
c 13 pair
        call compAtomPairList(ja,ia,
     &             pair13List,startPairL13,npairL13,find13)
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondListVCS:W8:find13:',find13
        end if
c
        if(find13) goto 102
        end if !13
c
       if( dij2 .le. r142max) then 
c 14 pair
        call compAtomPairList(ja,ia,
     &             pair14List,startPairL14,npairL14,find14)
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondListVCS:W8:find14:',find14
        end if
c
        if(find14) goto 102
        end if !14
c
c exclude ja if atom ja belongs to the same resNumb and 
c                           to the same HetCyclicbaseBlock in the RESidue
c
        if(OPT_removePairInHCycBase)then
        if((resNumb(ia) .eq. resNumb(ja)) .and.
     &     (atomBlockName(ia)(1:1) .eq. 'B' .and. 
     &      atomBlockName(ja)(1:1) .eq. 'B') .and.
     &      (atomBlockName(ia) .eq. atomBlockName(ja)))then
          goto 102             ! skip atom ja
         end if ! atomBlockName
        end if !OPT_removePairInHCycBase
c
c from this point atom ja should be included in the pairListV and C
c include in VdW pairList
        if( ( makeVDW .eq. 1 .and. dij2 .lt. rcutV2) .and.
cx     &      ((moveFlag(ja) .eq. 1 .and. ia .lt. ja) .or.                !lastmod
     &      ((moveFlag(ja) .eq. moveFlag(ia) .and. ia .lt. ja) .or.
     &       (moveFlag(ja) .eq. 0 )) ) then

        nnbPairLV(ia) = nnbPairLV(ia) + 1
        nnbPsV = nnbPsV + 1
c
        if(nnbPsV .gt. nnbpLVMAX) then
        write(kanalp,*)'ERROR:nonbondListVCS.f: nnbpLVMAX is low!!'
        write(kanalp,*)'increase param: nnbpLVMAX:',nnbPsV
        stop
        end if
c
        nbpairListV(nnbPsV) = ja
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondListV:W10: ia,ja,nnbPsV, dij2:',
     &  ia,ja,nnbPsV, dij2
        end if !C
c
        end if ! rcutV2
c
c make twin-Range R2 
c include in CL(R1-R2) list
       if( (makeCL .eq. 1 .and. 
     &     (dij2 .gt. rcutV2m .and. dij2 .lt. rcutC2p)) .and. 
cx     &      ((moveFlag(ja) .eq. 1 .and. ia .lt. ja) .or.             !lastmod
     &     ((moveFlag(ja) .eq. moveFlag(ia) .and. ia .lt. ja) .or.
     &      (moveFlag(ja) .eq. 0 )) ) then
c
        nnbPairLC(ia) = nnbPairLC(ia) + 1
        nnbPsC = nnbPsC + 1
c
        if(nnbPsC .gt. nnbpLCMAX) then
        write(kanalp,*)'ERROR:nonbondListVCS.f: nnbpLCMAX is low!!'
        write(kanalp,*)'increase param: nnbpLCMAX:',nnbPsC 
        stop
        end if
c
        nbpairListC(nnbPsC) = ja
c
        if(CONTROL1)then
        write(kanalp,*)'nonbondListC:W11: ia,ja,nnbPsC, dij2:',
     &  ia,ja,nnbPsC, dij2
        end if !Control
c
        end if !CL
c
	end if !ia#ja   
c
        else
        goto 101
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
        end do !iam
c
        if(makeVdW .eq. 1) nnbpLVT = nnbPsV
        if(makeS .eq. 1) nnbpLST = nnbPsS
        if(makeCL .eq. 1) nnbpLCT = nnbPsC
        ncallPL = ncallPL + 1
c
        if(CONTROL0)then
        write(kanalp,'(a22,i8,a42,3i8)')
     &  'nonbondListVCS: ncall:',ncallPL,
     &  'numbPairTot:nnbPsV,nnbPsC,nnbPsS:'
     &  ,nnbpLVT,nnbpLST,nnbpLCT
        end if !C
c
	if(CONTROL1)then
        write(kanalp,*)'nonbondListVCS:RESULT: '
        write(kanalp,*)'startnbPairLV():'
        call print10(natom,startnbPairLV)
        write(kanalp,*)'nnbPairLV():'
        call print10(natom,nnbPairLV)
        write(kanalp,*)'nbpairListV():nnbPsV:',nnbPsV 
        call print10(nnbPsV,nbPairListV)
c
        write(kanalp,*)'startnbPairLC():'
        call print10(natom,startnbPairLC)
        write(kanalp,*)'nnbPairLC():'
        call print10(natom,nnbPairLC)
        write(kanalp,*)'nbpairListC():nnbPsC:',nnbPsC 
        call print10(nnbPsC,nbPairListC)
        end if !CONTROL1
c
        if(CONTROLS)then
        write(kanalp,*)'startnbPairLS():'
        call print10(natom,startnbPairLS)
        write(kanalp,*)'nnbPairLS():'
        call print10(natom,nnbPairLS)
        write(kanalp,*)'nbpairListS():nnbPsS:',nnbPsS 
        call print10(nnbPsS,nbPairListS)
        end if !CONTROLS
c
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB

            return
            end
c
c compare atom ja with neihgbours of atom ia
	subroutine compAtomPairList(ja,ia,
     &             pairList,startPairL,npairL,find)
c
c pairList() - consequtive List of neighbours for ia=1,...,natom
c startPairL(ia) = startPos in pairList() for neighbours of atom ia
c nPairL(ia) = numbers of neighbours in the pairList() file
c find = .true. if Ja in the neighbour List of atom ia
c
        implicit none
        integer ia,ja
        integer pairList(*)
        integer startPairL(*)
        integer npairL(*)
        logical find
c
        integer i,natom,np12s

        find = .false.
c
        if(startPairL(ia) .gt. 0 .and. 
     &                 npairL(ia) .gt. 0) then  
    
        do i = startPairL(ia),startPairL(ia)+npairL(ia)-1
c
        if( pairList(i) .eq. ja ) then
        find = .true.
        return
        end if
c
        end do !
        end if ! nonzero Plist
c
        return
        end
c
