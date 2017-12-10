c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*                                                                           *
c*   getHbXHYtriplet                                                         *
c* Defines  the the X-H ...Y  H-bohd triplets                                *
c*                                                                           *
c*     Yury Vorobjev 2004                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine getHbXHYtriplet(atomXYZ,ffAtomName,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListV,startnbPairLV,nnbPairLV,startPairL12,
     &           rcutHb128,nYHXtriplet,nYHXtripletList,
     &           Hb128ParamList)
c
c look for X-H ... Y triplet in nbpairListV and pair14List
c
c  INput: natom,nmoveatom,moveAtomList
c         atomXYZ,ffAtomName,
c         nbpairListV,startnbPairLV,nnbPairLV,
c  OUt :
c        rcutHb128 - cutofDist
c        nYHXtriplet, number of triplets
c       nYHXtripletList(1,2,3) = (Y,H,X)
c
	implicit none
c
        real atomXYZ(*) 
        character*2 ffAtomName(*)
        integer natom
        integer nmoveatom,moveAtomList(*)
        integer nbpairListV(*)
        integer startnbPairLV(*),nnbPairLV(*)
        integer startPairL12(*)
        real rcutHb128
        integer nYHXtriplet,nYHXtripletList(*)
        real Hb128ParamList(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer i,i3,ia,iam
        integer j,ja,j3
        integer hbt(3)
        real rij(3),dij2,dij1
        real rcutHb2
        real hbrm12,hbem12,A12,B12
        real xi(3)
        integer nt3,nt4
	integer kanalp
        logical OPT_do14    
        logical hb128status
        logical CONTROL,CONTROL0
c
        OPT_do14 = .true. 
c
        kanalp = kanalRunOut
        CONTROL0 = .false.  !print
        CONTROL = .false.   !print
c
        rcutHb2 = rcutHb128**2
c
        if(CONTROL0)then
        write(kanalp,*)'* * * * * * * * * * * * * * * *' 
        write(kanalp,*)'Update  HbXHYtriplet List:'
        write(kanalp,*)'* * * * * * * * * * * * * * * *' 
        end if
c initialize
c
	if(CONTROL0)then
        write(kanalp,*)'getHbXHYtriplet: nmoveatom:',nmoveatom
        write(kanalp,*)'moveAtomList:'
        write(kanalp,*)(moveAtomList(i),i=1,nmoveatom)
        end if !C0
c
        nYHXtriplet = 0
c
c loop over movingAtomList
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3 = 3*ia-3
        xi(1) = atomXYZ(i3+1)
        xi(2) = atomXYZ(i3+2)
        xi(3) = atomXYZ(i3+3)
c
c collect full list of neighbours = nbPairLV + nbPairL14
c
c loop over shortR1 neighbours        
c
	do j = startnbPairLV(ia),startnbPairLV(ia)+nnbPairLV(ia)-1
        ja = nbpairListV(j)
c atoms ia,ja can interact via Hb128
        j3 = 3*ja-3
        rij(1) = atomXYZ(j3+1) - xi(1)
        rij(2) = atomXYZ(j3+2) - xi(2)
        rij(3) = atomXYZ(j3+3) - xi(3)
        dij2 = rij(1)**2 + rij(2)**2 + rij(3)**2 
c
        if( dij2 .lt. rcutHb2)then  
         dij1 = sqrt(dij2)
c
        call getHb128ijPar(ia,ja,hbrm12,hbem12,A12,B12,hb128status)
c
        if(hb128status)then
        nYHXtriplet = nYHXtriplet + 1
        nt3 = nYHXtriplet*3 -3         
c ia,ja atoms makes Hb128
c define which atom is H
        hbt(1) = ia     ! Y
        hbt(2) = ja     ! Hx
        if (ffAtomName(ia)(1:1) .eq. 'H') then
        hbt(2) = ia     ! Hx
        hbt(1) = ja     ! Y 
        end if
c define atom X
        hbt(3) = startPairL12(hbt(2))
c        
        nYHXtripletList(nt3+1)=hbt(1)
        nYHXtripletList(nt3+2)=hbt(2)
        nYHXtripletList(nt3+3)=hbt(3)
c        
        nt4 =  nYHXtriplet*4 -4
        Hb128ParamList(nt4+1) = hbrm12
        Hb128ParamList(nt4+2) = hbem12
        Hb128ParamList(nt4+3) = A12
        Hb128ParamList(nt4+4) = B12
c
        if(CONTROL)then
        write(kanalp,'(a23,i4,3i5,a18,2f7.4)')
     &      'getHbXHYtrip: nt,1,2,3:',nYHXtriplet,hbt,
     &                 ' Hb128Par: rm,em:',hbrm12,hbem12
        end if !C
c
        end if ! hb128status         
        end if ! dij2 .lt. rcutHb128
c
        end do ! j
        end do !iam
c
	return
	end
