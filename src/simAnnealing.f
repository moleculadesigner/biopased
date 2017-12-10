c simAnnealing  protocol
c
c Yuri Vorobjev   2002
c
        subroutine simAnnealing(nSAstep,nSAparMX,SAProtcol)
c
c        implicit none
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'mdRunPar.h'
        include 'vdw12Par.h'
        include 'engXYZtra.h'
        include 'optionPar.h'
        include 'enForce.h'
        include 'hbond128.h'
c
        integer nSAstep,nSAparMX
        real SAProtcol(*)
c
        integer nsas,shiftSAP
        real tempTg1,tempTg2
        real aSC,sc1,sc2
        real tempInSA
c
        integer i,k
        integer kanalp
        logical CONTROL,CONTROL1
        kanalp = kanalRunOut
        CONTROL = .true.
        CONTROL1 = .false.
c
        if(CONTROL)then
        write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * '
        write(kanalp,*) 'simAnnealing: Start'
        write(kanalp,*)'nStep:',nSAstep,' nSAparMX:',nSAparMX,
     &   ' SAprotcol[nStep,Tmd,aSCore,hB128,Coul]:'
        do k=1,nSAstep
        write(kanalp,*)'nS:',k,
     &                 (SAprotcol((k-1)*nSAparMX+i),i=1,nSAparMX)
        end do!k
        write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * '
        end if
c
c run SimAnnealing Protocol
c
       nsas = 0
c
100    nsas = nsas + 1
c run SA step = nsas
       if(nsas .gt. nSAstep ) goto 101
c
       ntime0  = ntime
c
c aug.19, 2004: assign VariableForceField SAprotol parameters to FField:
c
       shiftSAP = nSAparMX*(nsas-1)
       ntimeMX = SAProtcol(shiftSAP + 1) 
       tempTg2 =  SAProtcol(shiftSAP + 2) 
c vdw SoftCore
       aSC =  SAProtcol(shiftSAP + 3)   ! aSC : 1 -->0 (rigid-->soft)
       if(aSC .lt. 0.0)aSC = 0.0        ! maxSoft
       if(aSC .gt. 1.0)aSC = 1.0        ! maxRigid
       sc1 = 0.65        ! MINscPar= rSC/rMIN
       sc2 = 0.975       ! MAXscPar
       aSoftCore = sc1*aSC + sc2*(1.0 - aSC)  ! vdw12Par.h 
c  CoulForces
c       fEngWF(6) = SAProtcol(shiftSAP + 6)
c       fEngWF(7) = fEngWF(6)             
c hB128 forces    
       fEngWF(10) =  SAProtcol(shiftSAP + 4)      
       fEngWF(11) =  SAProtcol(shiftSAP + 5) 
c Update hB128WfList(i)
       if(CONTROL1)then
       write(kanalp,*)'hb128: hb128TypeList: 1=(BackB),2=(BS,SS)',
     & ' fEngWF(10),(11):',fEngWF(10),fEngWF(11), 
     & ' hB128WfScaleALL=',hB128WfScaleALL
       end if !C
c
       do i = 1,nHBbondHxY
       if(hB128TypeList(i) .eq. 1)
     & hB128WfList(i)=fEngWF(10)*hB128WfScaleALL
       if(hB128TypeList(i) .eq. 2)
     & hB128WfList(i)=fEngWF(11)*hB128WfScaleALL
c
       if(CONTROL1)then
       write(kanalp,*)'hb128:ih:',i,' hb128TypeList: hB128WfList:',
     &  hB128TypeList(i),hB128WfList(i)
       end if
       end do !i
c 
      if(CONTROL1)stop
c
       optra = 0         ! open *.tra files 
       if(nsas .eq. 1 .and. (.not. OPT_doMD))then
       tempInSA = 50  ! default in inputMDSApar                              
       tempTg = tempInSA
c
       call initMDStart(tempTg)
c
       optra = 1        
       end if !nsas = 1       
c
       atype = 1         ! 0- free MD, 1- NTV
       wtra  = 1         ! write *.tra
       cltra = 0         ! close *.tra files
c
       tempTg1 = tempTg 
c
       if(CONTROL)then
       write(kanalp,*)'simAnnealing: step:',nsas
       write(kanalp,'(a22,i6,a13,2f8.2)')
     & 'MD step to be done:',ntimeMX,' T1->T2(K):',tempTg1,tempTg2
       end if
c
	tempTg1 = tempTg
c
	call mdRun(ntimeMX,ntime0,ntime,ntimeR1,ntimeR2,
     &             ntimeF1,ntimeF2,ntimeF3,deltat, 
     &       tempTg1,tempTg2,tauTRF,atype,optra,wtra,nwtra,cltra)
c
        tempTg = tempTg2
c
       if(CONTROL)then
       write(kanalp,*)'simAnnealing: step:',nsas, ' FINISH!'
c average engPot
       engPOTENTTraAv = 0.0
       engPOTENTTraSd = 0.0 
       engVDWR1TraAv = 0.0
       engVDWR1TraSd = 0.0 
       engCOULTraAv = 0.0
       engCOULTraSd = 0.0 
       hBHxYeng128TraAv = 0.0
       hBHxYeng128TraSd = 0.0
       engMolSolAv = 0.0
       engMolSolSd = 0.0 
       tempT0traAv = 0.0
       tempT0traSd = 0.0
c      
       do i = 1,nStepTra 
       engPOTENTTraAv = engPOTENTTraAv + engPOTENTTra(i)
       engVDWR1TraAv = engVDWR1TraAv + engVDWR1Tra(i)
       engCOULTraAv = engCOULTraAv + engCOULTra(i)
       tempT0traAv = tempT0traAv + tempT0tra(i)
       hBHxYeng128TraAv = hBHxYeng128TraAv +  hBHxYeng128Tra(i)
       engMolSolAv = engMolSolAv + eMolSolvTra(i)
       end do    
c
       engPOTENTTraAv = engPOTENTTraAv/nStepTra
       engVDWR1TraAv = engVDWR1TraAv/nStepTra
       engCOULTraAv = engCOULTraAv/nStepTra 
       hBHxYeng128TraAv = hBHxYeng128TraAv/nStepTra
       engMolSolAv = engMolSolAv/nStepTra
       tempT0traAv = tempT0traAv/nStepTra
c
       write(kanalp,*)'Average Eng over nSnapTra:',nStepTra
       write(kanalp,*)'engPOTENTTraAv: ',engPOTENTTraAv
       write(kanalp,*)'engVDWR1TraAv :',engVDWR1TraAv
       write(kanalp,*)'engCOULTraAv  :',engCOULTraAv
       write(kanalp,*)'hBHxYeng128TraAv:',hBHxYeng128TraAv
       write(kanalp,*)'engMolSolAv:',engMolSolAv 
       write(kanalp,*)'tempT0traAv:',tempT0traAv 
c
       end if !CONTROL
c
        goto 100
101     continue
c
        if(CONTROL)then
        write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * '
        write(kanalp,*) 'ff-variable simAnnealing: Finish'
        write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * '
        end if
c
         return
	 end
