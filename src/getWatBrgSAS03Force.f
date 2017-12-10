c calculate forces from WaterBridges  
c solvateMoL03  solvation model ElectrFieldWaterShell
c
	subroutine getWatBrgForceSAS03(atomXYZs,
     &                   watBrgEnergySoL,watBrgForceAt)
c
c calculates watBrgEnergySoL,watBrgForceAt(*) : enForce.h
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
        include 'solvWBrg01.h'
        include 'coulEnPar.h'
c
        real watBrgEnergySoL
        real watBrgForceAt(*)
        real atomXYZs(*)
c local
        real forceSum(3)
        real forceSumM(3),forceSumP(3)
        real scalP(3),scalM(3)
        integer iTypeCorr
        integer natInListMAX
        integer natInList
        parameter (natInListMAX = 3000)
        integer atWbrgNeibList(natInListMAX)
c
        real watBrgForceSUM(3)
        real sckal,cReps
        real rcut(2),rsc,rcut2
        real dij1,dij2
        real rij(3)
        real ecoul,fi(3)
        real ecoulSC,fiSC(3)
        real ixyz(3)
        real qi,qj
        integer k,iw,ia,ia3,ip
        integer iaw,iawb,iawb3
        integer jn,ja,ja3
        integer natInOneWbrg
c
        real scaleEFS_OPT
c
        logical CONTROL
        integer kanalp
c
        CONTROL=.false.  
        kanalp =  kanalRunOut 
c
        scaleEFS_OPT = 0.80  !empirical scaling coeff
c
        do k=1,3*natomMAX
        watBrgForceAt(k) = 0.0
        end do
        watBrgEnergySoL = 0.0 
c
        iTypeCorr = 0
cx        iTypeCorr = 1    ! 0=equal correction, 1=proportiona correction
c
        do k=1,3
        watBrgForceSUM(k) = 0.0
        forceSum(k)=0.0
        forceSumP(k)=0.0
        forceSumM(k)=0.0
        scalP(k)=0.0
        scalM(k)=0.0
        end do!k
c
        natInOneWbrg = 3
        ip =1
        sckal = 332.0
        cReps = cReps_OPT
        rcut2 = RcutEFSMAX**2
        rcut(1) = 1.0/RcutEFSMAX
        rcut(2) = rcut(1)**2
        rsc = 1.2 !linear interpolation rij<rsc
c
cx        write(kanalp,*)'Start getWatBrgSAS03Force:'
cx        write(kanalp,*)
cx     & 'getWatBrgForce:cReps,rcut:',cReps,rcut
c
        if(nWBrgNow .gt. 0)then 
        do iw = 1,nWBrgNow
c
        scaleWBrgEng(iw) = 0.0
        if(wBrgEpotQ(iw) .lt. wbrgEpotQMAX_OPT)then ! over strong WatBrg
        scaleWBrgEng(iw) = 1.0 - wbrgEpotQMAX_OPT/wBrgEpotQ(iw)
c
        ia = dot_IATNUM(ideFMX(iw),ip)    ! SAS atom for iProbe=1
c loop over neighbour atoms
cx        write(kanalp,*)
cx     & 'getWatBrgSAS03Force:iw,ideFMX(iw),ia:',iw,ideFMX(iw),ia
cx        write(kanalp,*)'nnbPairLSAS03(ia) :',nnbPairLSAS03(ia)
c
        if(nnbPairLSAS03(ia) .ge. 1) then
c
        do iaw =1,natInOneWbrg
        iawb=(iw-1)*natInOneWbrg+iaw
        iawb3=3*(iawb-1)
        qi=wBrgAtomQ(iawb)*sckal
c
        do k=1,3
        ixyz(k) = wBrgWmoLXYZ(iawb3+k)
        end do!k
c
cx        write(kanalp,*)'iawb,ixyz:',iawb,ixyz
c
c make atWbrgNeibList(*)
        natInList=1
        atWbrgNeibList(natInList) = ia
        do jn = startnbPairLSAS03(ia),
     &          startnbPairLSAS03(ia) + nnbPairLSAS03(ia)-1
           ja = nbpairListSAS03(jn)   ! atom j
         natInList=natInList + 1
         if(natInList .gt. natInListMAX)then
         write(kanalp,*)'ERROR!:natInListMAX low! getWatBrgSAS03Force.f'
         write(kanalPStat,*)mError,
     &    ' natInListMAX is low! (MAX number of WatBridges)',
     &    ' in getWatBrgSAS03Force.f '
         stop
         end if 
         atWbrgNeibList(natInList) = ja 
         end do !jn
c
         do jn = 1,natInList
           ja = atWbrgNeibList(jn)
           ja3 = ja*3-3
c
        qj = atomQ(ja)
c
        dij2=0.0
        do k=1,3
        rij(k) = atomXYZs(ja3+k)- ixyz(k) 
        dij2=dij2 + rij(k)**2
        end do!k
c
cx        write(kanalp,*)' jn,ja:',jn,ja
cx        write(kanalp,*)' dij2, rij:',dij2, rij
c
        if(dij2 .lt. rcut2)then
        dij1 = sqrt(dij2)
c
c force on atom j = -fi
	call coulenforceij02SC
     &  (cReps,rcut(1),rcut(2),dij2,dij1,rij,qi,qj,rsc,ecoul,fi)
c
cx        write(kanalp,*)' qi,qj;',qi/sckal,qj
cx        write(kanalp,*)' ecoul,fi:',ecoul,fi
c
        ecoulSC = ecoul*scaleWBrgEng(iw)
        watBrgEnergySoL = watBrgEnergySoL + ecoulSC
        do k=1,3
        fiSC(k)= fi(k)*scaleWBrgEng(iw)
        watBrgForceAt(ja3+k) = watBrgForceAt(ja3+k) - fiSC(k)
        watBrgForceSUM(k) = watBrgForceSUM(k) - fiSC(k)
        end do !k
        end if ! dij2 .lt. rcut(2)
c 
         end do! jn
         end do! iaw
         end if !nnbPairLSAS03(ia) .ge. 1
c
         end if !wBrgEpotQ(iw) .lt. wbrgEpotQMAX_OPT)then ! over strong WatBrg
         end do !iw
         end if !nWBrgNow .gt. 0
c
        if(CONTROL)then
        write(kanalp,*)
     &  'getWatBrgSAS03Force: watBrgEnergySoL:',watBrgEnergySoL
        do ia = 1,natom
        ia3 = 3*(ia-1)
        write(kanalp,7003)
     &  'ATOM  ',ia,atomName(ia),resName(ia),resNumb(ia),
     &  (atomXYZs(ia3+k),k=1,3),(watBrgForceAt(ia3+k),k=1,3),
     &   atomQ(ia)
        end do!ia
        write(kanalp,*)'watBrgForceSUM: ',(watBrgForceSUM(k),k=1,3)
        end if !C
c
c force correction by zeroSum condition
c
        do ia=1,natom
        ia3=3*ia-3
        do k=1,3
        forceSum(k)=forceSum(k) + watBrgForceAt(ia3+k)
        if(watBrgForceAt(ia3+k) .gt. 0.0)
     &  forceSumP(k)=forceSumP(k) + watBrgForceAt(ia3+k)
        if(watBrgForceAt(ia3+k) .lt. 0.0)
     &  forceSumM(k)=forceSumM(k) + watBrgForceAt(ia3+k)
        end do!k
        end do !ia
c
        do k=1,3
        if(ForceSumP(k) .gt. 0.0)
     &        scalP(k)=ForceSUM(k)*0.5/ForceSumP(k)
        if(ForceSumM(k) .lt. 0.0)
     &        scalM(k)=ForceSUM(k)*0.5/ForceSumM(k)
        ForceSUM(k) = ForceSUM(k)/natom          
        end do !k
c
        do ia=1,natom         
        ia3=3*ia-3
        do k=1,3
        if(iTypeCorr .eq. 1)then
        if(watBrgForceAt(ia3+k) .gt. 0.0)
     &  watBrgForceAt(ia3+k) = watBrgForceAt(ia3+k)*(1.0 - scalP(k))
        if(watBrgForceAt(ia3+k) .lt. 0.0)
     &  watBrgForceAt(ia3+k) = watBrgForceAt(ia3+k)*(1.0 - scalM(k))
c
        else
        watBrgForceAt(ia3+k)=watBrgForceAt(ia3+k) - ForceSUM(k)
        end if !iTypeCorr
c scale
        watBrgForceAt(ia3+k)=watBrgForceAt(ia3+k)*scaleEFS_OPT
        end do!k
        end do !ia
c
       watBrgEnergySoL = scaleEFS_OPT*watBrgEnergySoL
c
        if(CONTROL)then
        do k =1,3      
        ForceSUM(k)=0.0
        do ia=1,natom
        ia3=3*ia-3
        ForceSUM(k)=ForceSUM(k)+watBrgForceAt(ia3+k)
        end do !ia
        end do !k        
        write(kanalp,*)'scalP(k):',scalP
        write(kanalp,*)'scalM(k):',scalM 
        write(kanalp,*)'watBrgForceSUM:corrected: ',(ForceSUM(k),k=1,3) 
        end if !C
c
	return
c
7003   format(a4,2x,i5,1x,a4,a4,1x,i4,4x,3f8.3,1x,3f8.3,f8.5) ! PDB 
        end
