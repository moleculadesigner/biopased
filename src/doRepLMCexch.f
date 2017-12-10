c do replica i,i+1 exchange xyz, vel       
c
	subroutine doRepLMCexch(nMCcycle)  
c
c nMCcycle - the current number of the MC replica exchange cycle     
c          = nMDcycle the number of the MD replica mixingCycle, 
c	implicit none
        include 'xyzPDBsize.h' 
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'replicaEx.h'
c        include 'mdRunPar.h'
        include 'randomGen.h'
c
        integer nMCcycle
c local
        real tempT1,tempT2
        real tempi,tempj
        real scKb,scKb1
        real scaleVij
        real eDelta,expDelta
        real eDeltaij,eDeltaji,dBetaij
        real expDeltaij,expDeltaji
        real eDeltaTcorr,deltaTemp,ePotTempCoeff
        integer nRepLicMin
        integer i,k,na3,nRpL
        integer pShift
        integer eAvShifti,eAvShiftj
        integer xShifti,xShiftj
        logical OPT_minEngSnap,OPT_eTempCorr
        integer nEngTermRef
        real eDelta0,eDelta1
        logical doRepLexchij,doRepLexchji
        logical CONTROL
        integer kanalp
c
        scKb = 0.002    !scaling kBotzmann kcal/mol
        scKb1 = 1.0/scKb
        ePotTempCoeff = 0.75E-03    !empirical coeff ePot(T) = epot0 + epC*Temp*natom  
c
cx       kanalp = 6
       kanalp = kanalRunOut
       CONTROL = .true.
c
       OPT_minEngSnap = .true.
cx      OPT_minEngSnap = .false.
c
       nEngTermRef=2
       OPT_eTempCorr = .true.
c init replicas XYZ
       na3 = 3*natom
       nRpL = nRepLic
       nRepLicMin=2
       pShift = 2*nRepLic*(nMCcycle - 1)
c copy to newReplicas XYZV
        do k=1,na3*nRepLic
        atomXYZRepLnew(k) = atomXYZRepL(k)
        atomVelRepLnew(k) = atomVelRepL(k)
        end do !i
c
        if(CONTROL)then
        write(kanalp,*)'doRepLMCexch: START * * * * * * * * *'
        write(kanalp,*)'doRepLMCexch: nMCexCycle=',nMCcycle
        write(kanalp,*)'doRepLMCexch: nMDsimCycle=',nMDsimCycle
        write(kanalp,*)'OPT_minEngSnap,OPT_eTempCorr :',
     &  OPT_minEngSnap,OPT_eTempCorr
        end if !C
c Starts exchange from highT to lowT i-->j
c i = nRpL
c j = nRpL-1
c
100    if(nRpL .lt. nRepLicMin) goto 101
c replics indexes
c i = nRpL
c j = nRpL-1
       tempi = tempRepLic(nRpL)
       tempj = tempRepLic(nRpL-1)
       deltaTemp = tempi - tempj 
c energyDiff
       eAvShifti = engTermMAX*nRepLic*(nMDsimCycle-1)
     &           + engTermMAX*(nRpL-1) 
       eAvShiftj = engTermMAX*nRepLic*(nMDsimCycle-1)
     &           + engTermMAX*(nRpL-2)
c calculate energyDecrement
       eDelta0 = mdEngAvRepLtra(eAvShifti+nEngTermRef) 
     &          - mdEngAvRepLtra(eAvShiftj+nEngTermRef)
c
       if(OPT_minEngSnap)then
       eDelta0 = mdEngRepL0MixMin(nEngTermRef,nRpL) - 
     &          mdEngRepL0MixMin(nEngTermRef,nRpL-1) 
       end if !OPT_minEngSnap
c
        eDeltaTcorr = 0.0
       if(OPT_eTempCorr)then
        eDeltaTcorr = deltaTemp*ePotTempCoeff*natom 
       end if!OPT_eTempCorr
        eDelta1 = eDelta0 - eDeltaTcorr
c
       dBetaij = 1.0/tempRepLic(nRpL-1)-1.0/tempRepLic(nRpL)
       eDelta = abs(dBetaij)*eDelta1
       eDeltaij = eDelta*scKb1
       eDeltaji = -eDeltaij
c
c calculate Probability of exchange i->j (replace j by i)
       if( eDeltaij .le. 0.0)then
       expDeltaij = 1.0
       else
       expDeltaij = exp(-eDeltaij)
       end if
c calculate Probability of exchange j->i (replace i by j)
       if( eDeltaji .le. 0.0)then
       expDeltaji = 1.0
       else
       expDeltaji = exp(-eDeltaji)
       end if
c
       repLijExProbabilityTra(pShift+2*nRpL-1) = expDeltaij
       repLijExProbabilityTra(pShift+2*nRpL) = expDeltaji 
c
       if(CONTROL)then
       write(kanalp,*)'doRepLMCexch: nR->nR-1, nRpL=',nRpL
       write(kanalp,*)'doRepLMCexch: tempi->tempj:',tempi,tempj
       write(kanalp,*)'doRepLMCexch: eDelt0,eDelt1,eCorrT (kcal/mol):',
     & eDelta0,eDelta1,eDeltaTcorr
       write(kanalp,*)'doRepLMCexch: ei-ej=beta*Deltaij:',eDeltaij 
       write(kanalp,*)'doRepLMCexch: eAvShifti,eAvShiftj:',
     & eAvShifti,eAvShiftj
       end if!C 
c Metropolis Criterium
       doRepLexchij = .false. 
       if(expDeltaij .ge. 1.0) then 
        doRepLexchij = .true.
       else
       call RANDOM(randV01,randSeedIG)
       if(expDeltaij .ge. randV01)doRepLexchij = .true.
       end if ! expDelta .ge. 1.0
c
       doRepLexchji = .false.
       if(expDeltaji .ge. 1.0) then
        doRepLexchji = .true.
       else
       call RANDOM(randV01,randSeedIG)
       if(expDeltaji .ge. randV01)doRepLexchji = .true.
       end if ! expDelta .ge. 1.0
c
       if(CONTROL)then
       write(kanalp,*)
     & 'doRepLMCexch: expDeltaij,expDeltaji:',expDeltaij,expDeltaji
       write(kanalp,*)'doRepLMCexch: doRepLexchij =',doRepLexchij
       write(kanalp,*)'doRepLMCexch: doRepLexchji =',doRepLexchji
       end if!C
c
c do replicaExchange
cx       xShiftj = na3*(nRpL - 1)
cx       xShifti = na3*nRpL
       xShifti = na3*(nRpL - 1)
       xShiftj = na3*(nRpL - 2)
       if(CONTROL)then
       if(xShifti + na3. gt. 3*natomMAX*nRepLicMAX) then
       write(kanalp,*)
     & 'doRepLMCexch:ERROR! xShifti + na3 LARGE!:',(xShifti+na3)
       stop
       end if
       end if !C
       scaleVij = sqrt(tempi/tempj)
c
       if(CONTROL)then
       write(kanalp,*)'doRepLMCexch: xShifti,xShiftj,scaleVij:',
     & xShifti,xShiftj,scaleVij
       write(kanalp,*)'doRepLMCexch: nTimeStepRepL(nRpL):',
     & nTimeStepRepL(nRpL)
       end if !C
c
       if(doRepLexchij)then
c store replicaExchangeMatrix
c init to Zero in the initReplica
       repLijTranzitionTra(pShift+2*nRpL-1) = 1        
c
        do k=1,na3 
        atomXYZRepLnew(xShiftj+k) = atomXYZRepL(xShifti+k)
        atomVelRepLnew(xShiftj+k) = atomVelRepL(xShifti+k)/scaleVij
        end do!k 
        end if !doRepLexchij
c
        if(doRepLexchji)then
c store replicaExchangeMatrix
c init to Zero in the initReplica
       repLijTranzitionTra(pShift+2*nRpL) = 1
c
        do k=1,na3
        atomXYZRepLnew(xShifti+k) = atomXYZRepL(xShiftj+k)
        atomVelRepLnew(xShifti+k) = atomVelRepL(xShiftj+k)*scaleVij
        end do!k
        end if !doRepLexchji 
c
        nRpL = nRpL - 1
        goto 100
101     continue
c
c copy newReplicas XYZV to the newAccepted XYZVstate
        do k=1,na3*nRepLic
        atomXYZRepL(k) = atomXYZRepLnew(k)
        atomVelRepL(k) = atomVelRepLnew(k)
        end do !i
c
       if(CONTROL)then
       write(kanalp,*)'doRepLMCexch: nMCcycleL:',nMCcycle
       write(kanalp,*)'doRepLMCexch: nMDstep DONE - - - - - - -'
       end if !C
c
       return
       end
