c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*                                                                           *
c*  getHbTripletAll                                                          *
c* Defines  the hBond Triplets X-H...Y                                       *
c*                                                                           *
c*     Yury Vorobjev 2004                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine getHbTripletAll(atomXYZ,
     &           natom,nmoveatom,moveAtomList,
     &           nbpairListV,startnbPairLV,nnbPairLV,
     &           rcutHb128,nHBbondHxYmx,
     &           nHBbondHxY,hB128List,hB128TypeList,
     &           hB128PotParList)
c
	implicit none
c
        real atomXYZ(*) 
        integer natom
        integer nmoveatom,moveAtomList(*)
        integer nbpairListV(*)
        integer startnbPairLV(*),nnbPairLV(*)
        real rcutHb128
        integer nHBbondHxYmx,nHBbondHxY
        integer hB128List(*)
        real hB128PotParList(*)
        integer hB128TypeList(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer i,i2,i3,i4,ia,iam
        integer j,ja,j3,k
        integer tijk(3)
        integer kanalp
        real rij(3),xi(3)
        real dij2,dij1
        real A12,B12
        integer nHb3,nHb4
        integer hB128Type
        real hbrm,hbem,rcutHb128p2
        logical hb128status
        logical CONTROL,CONTROL0
c
        kanalp = kanalRunOut
        CONTROL0 = .false.  !print
        CONTROL = .false.   !print
c
        if(CONTROL0)then
        write(kanalp,*)'* * * * * * * * * * * * * * * *' 
        write(kanalp,*)'Update getHbTripletAll:'
        write(kanalp,*)'* * * * * * * * * * * * * * * *' 
        end if
c initialize
c
c        write(kanalp,*)'nHBbondHxYmx :',nHBbondHxYmx
c
        rcutHb128p2 = rcutHb128**2
        nHBbondHxY = 0
	do i = 1, nHBbondHxYmx
        i3=3*i-2
        hB128List(i3) = 0
        hB128List(i3+1) = 0
        hB128List(i3+2) = 0
c
        i4 = 4*i-3
        hB128PotParList(i4) = 0.0
        hB128PotParList(i4+1) = 0.0 
        hB128PotParList(i4+2) = 0.0
        hB128PotParList(i4+3) = 0.0
        end do
c
	if(CONTROL0)then
        write(kanalp,*)'getHbondTriplets: nmoveatom:',nmoveatom
        write(kanalp,*)'moveAtomList:'
        write(kanalp,*)(moveAtomList(i),i=1,nmoveatom)
        end if !C0
c
c loop over movingAtomList
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3 = 3*ia-3
        xi(1) = atomXYZ(i3+1)
        xi(2) = atomXYZ(i3+2)
        xi(3) = atomXYZ(i3+3)
c
c loop over shortR1 neighbours        
	do j = startnbPairLV(ia),startnbPairLV(ia)+nnbPairLV(ia)-1
        ja = nbpairListV(j)
c atoms ia,ja interacts
        j3 = 3*ja-3
        rij(1) = atomXYZ(j3+1) - xi(1)
        rij(2) = atomXYZ(j3+2) - xi(2)
        rij(3) = atomXYZ(j3+3) - xi(3)
        dij2 = rij(1)**2 + rij(2)**2 + rij(3)**2 
c
c Hbond128
        if( dij2 .lt. rcutHb128p2)then
        dij1 = sqrt(dij2)
        call getHb128ijkTripPar
     &          (ia,ja,tijk,hbrm,hbem,A12,B12,hb128status,hB128Type)
c
        if(hb128status)then
        nHBbondHxY = nHBbondHxY + 1
c
        if(nHBbondHxY .gt. nHBbondHxYmx)then
        write(kanalp,*)
     & 'getHb128TrpletAll: ERROR! nHBbondHxYMAX is small !'
        stop
        end if
c
        hB128TypeList(nHBbondHxY) = hB128Type 
        nHb3=3*nHBbondHxY-2
        hB128List(nHb3)=tijk(1)
        hB128List(nHb3+1)=tijk(2)
        hB128List(nHb3+2)=tijk(3)
c
        nHb4=4*nHBbondHxY-3
        hB128PotParList(nHb4) = A12
        hB128PotParList(nHb4+1) = B12 
        hB128PotParList(nHb4+2) = hbrm
        hB128PotParList(nHb4+3) = hbem
c
        end if
        end if ! dij2 .lt. rcutHb128p2
c
        end do !iam
        end do !j
c
        if(CONTROL0)then
        write(kanalp,*)'getHb128TrpletAll: nHBbondHxY=',nHBbondHxY
        if(nHBbondHxY .gt. 0)then
        do  i = 1,nHBbondHxY
        nHb3=3*i-3
        nHb4=4*i-4
        write(kanalp,*)'Y,H,X =',(hB128List(nHb3+k),k=1,3)
        write(kanalp,*)
     &  'A12,B12,rm,em: ',(hB128PotParList(nHb4+k),k=1,4) 
        end do !i
        end if 
        end if !C
c
c	stop
c
	return
	end
