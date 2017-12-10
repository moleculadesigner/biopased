c calculate forces for BornPolarization E term
c
	subroutine getBornPolzEforce(atomXYZ,
     &               bornPolzEng,bornPolzEngAt,bornPolzForceAt)
c
c calculates bornPolzEng,bornPolzForceAt(*) : enForce.h
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'dataSASdr.h'
        include 'nbPairSAS03.h'
        include 'solvGBorn.h'
        include 'coulEnPar.h'
        include 'movingAtom.h'
c
        real atomXYZ(*)
        real bornPolzEng,bornPolzEngAt(*)
        real bornPolzForceAt(*)
c local
        real epsSOL
        real forceSum(3)
        real forceSumM(3),forceSumP(3)
        real scalP(3),scalM(3)
        integer iTypeCorr
        integer ia,ia3,k
c
        logical CONTROL
        integer kanalp
c
        CONTROL=.false. 
c
cx        kanalp = 6
        kanalp = kanalRunOut
c
cx        write(kanalp,*)'Start :getBornPolzEforce:'
        
c
        epsSOL = epsSol_OPT
c
        do k=1,3*natomMAX
        bornPolzForceAt(k)=0.0
        end do
c
        do k=1,natomMAX
        bornPolzEngAt(k) = 0.0
        end do!k
        bornPolzEng = 0.0 
c
        iTypeCorr = 0
cx        iTypeCorr = 1    ! 0=equal correction, 1=proportional correction
c
        do k=1,3
        forceSum(k)=0.0
        forceSumP(k)=0.0
        forceSumM(k)=0.0
        scalP(k)=0.0
        scalM(k)=0.0
        end do!k
c
cx        sckal = 332.0 ! scaled in call selfEPolEngGBorn
c
        do ia = 1,natom        
        ia3=3*ia-3
c
	call selfEPolEngGBorn(atomQ(ia),atBornRad(ia),
     &    atBornRadDr(ia3+1),bornPolzEngAt(ia),
     &    bornPolzForceAt(ia3+1))        
c
        bornPolzEng = bornPolzEng + bornPolzEngAt(ia)
        end do!iam 
c
        if(CONTROL)then
        do ia = 1,natom
        ia3 = 3*(ia-1)
        write(kanalp,7003)
     &  'ATOM  ',ia,atomName(ia),resName(ia),resNumb(ia),
     &  (atomXYZ(ia3+k),k=1,3),
     &  (bornPolzForceAt(ia3+k),k=1,3),bornPolzEngAt(ia),
     &   atomQ(ia)
c
        do k=1,3
        forceSum(k) = forceSum(k)+bornPolzForceAt(ia3+k)/natom
        end do!k
c
        end do!ia
c
        write(kanalp,*)'bornPolzForceAtSUM: ',(forceSum(k),k=1,3)
        end if !C
c
c force correction by zeroSum condition
c
        if(iTypeCorr .eq. 1)then 
        do ia=1,natom
        ia3=3*ia-3
        do k=1,3
        if(bornPolzForceAt(ia3+k) .gt. 0.0)
     &  forceSumP(k)=forceSumP(k) + bornPolzForceAt(ia3+k)
        if(bornPolzForceAt(ia3+k) .lt. 0.0)
     &  forceSumM(k)=forceSumM(k) + bornPolzForceAt(ia3+k)
        end do!k
        end do !ia
c
        do k=1,3
        if(ForceSumP(k) .gt. 0.0)
     &        scalP(k)=ForceSUM(k)*0.5/ForceSumP(k)
        if(ForceSumM(k) .lt. 0.0)
     &        scalM(k)=ForceSUM(k)*0.5/ForceSumM(k)
        end do !k
        end if !iTypeCorr .eq. 1
c
        do ia=1,natom         
        ia3=3*ia-3
        do k=1,3
        if(iTypeCorr .eq. 1)then
        if(bornPolzForceAt(ia3+k) .gt. 0.0)
     &  bornPolzForceAt(ia3+k) = bornPolzForceAt(ia3+k)*(1.0 - scalP(k))
        if(bornPolzForceAt(ia3+k) .lt. 0.0)
     &  bornPolzForceAt(ia3+k) = bornPolzForceAt(ia3+k)*(1.0 - scalM(k))
c
        else
        bornPolzForceAt(ia3+k)=bornPolzForceAt(ia3+k) - ForceSUM(k)
        end if !iTypeCorr
        end do!k
        end do !ia
c
        if(CONTROL)then
        ForceSUM(1)=0.0
        ForceSUM(2)=0.0 
        ForceSUM(3)=0.0 
        do ia=1,natom
        ia3=3*ia-3
        do k =1,3
        ForceSUM(k)=ForceSUM(k)+bornPolzForceAt(ia3+k)
        end do!k
c
        write(kanalp,7003)
     &  'ATOM  ',ia,atomName(ia),resName(ia),resNumb(ia),
     &  (atomXYZ(ia3+k),k=1,3),
     &  (bornPolzForceAt(ia3+k),k=1,3),bornPolzEngAt(ia),
     &   atomQ(ia)
c
        end do !ia
        write(kanalp,*)'scalP(k):',scalP
        write(kanalp,*)'scalM(k):',scalM 
        write(kanalp,*)
     &  'getBornPolzEforce corrected: ',(ForceSUM(k),k=1,3) 
        end if !C
c
	return
c
7003   format(a4,2x,i5,1x,a4,a4,1x,i4,4x,3f8.3,1x,3f8.3,f7.2,1x,f7.4) ! PDB 
        end
