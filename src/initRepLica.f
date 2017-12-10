c init Replica's T(i), etc.
c
	subroutine initReplicaSet
c
c	implicit none
        include 'xyzPDBsize.h' 
        include 'xyzPDBinfo.h'
        include 'replicaEx.h'
        include 'mdRunPar.h'
        include 'mdAtomXYZvel.h'
        include 'randomGen.h'
        include 'enForce.h'
        include 'engXYZtra.h'
        include 'ligInfo.h'
        include 'compGeoEF.h'
        include 'hbond128.h'
        include 'pair1234array.h'
        include 'optionPar.h' 
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
c local
        character*7 eMdRepLf,mMdRepLf
        character*4 nRePC,nRePC4
        character*(charLenMAX) filePdbN
cx        integer kanalPDBres 
        real tempTg1,tempTg2
        integer i,k,k3,j,ich
        integer na3,nReP,xShift
        logical doWRpdb
        integer kanalp
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .true.
c
        if(CONTROL)then
        write(kanalp,*)'initRepLica: START * * * * * * * * * *'
        end if
c
        eMdRepLf = 'eMdRepL'
        mMdRepLf = 'mMdRepL'
c
c replica temperatures
        tempRepLic(1) = temp0RepL       ! lowest T
        nStepTraRepL(1) = 0
        nRecPdbRepL(1) = 0
        nTimeStepRepL(1) =0
        nSnapMixMdRun = 0
c
        do i = 1,2*nRepLicMAX*nMCexCycleMAX
        repLijExProbabilityTra(i) = 0.0
        repLijTranzitionTra(i) =0
        end do
c
        do i = 2,nRepLicMAX
        nStepTraRepL(i) = 0
        nRecPdbRepL(i) = 0
        nTimeStepRepL(i) = 0
        tempRepLic(i) = tempRepLic(i-1)
     &                  *(1.0 + deltaT0RepL/temp0RepL)
        end do !i
c
c init replicas XYZ
       ntimeMX = ntimeREMDinit
       na3 = 3*natom
       nReP = 0
cinit randomGeneratorSeed
       randSeedIG = ntimeMX*na3*nRepLic  
c
100    nReP = nReP + 1
       if(nReP .gt. nRepLic) goto 101
c
       call RANDOM(randV01,randSeedIG) 
c
          ntime0 = nTimeStepRepL(nReP)
          tempTg1 = tempRepLic(nReP)
          tempTg2 = tempRepLic(nReP)
c
       atype = 1         ! 0- free MD, 1- NTV
       wtra  = 0         ! write *.tra
       cltra = 0         ! close *.tra files
c
       if(CONTROL)then
       write(kanalp,*)'initRepLica: mdRun START nRepL:',nReP,
     & ' - - - - - - - - - - '
       end if
c
        call mdRun(ntimeMX,ntime0,ntime,ntimeR1,ntimeR2,
     &             ntimeF1,ntimeF2,ntimeF3,deltat,
     &       tempTg1,tempTg2,tauTRF,atype,optra,wtra,nwtra,cltra)
c
       xShift = na3*(nReP - 1)
       do k = 1,na3
       atomXYZRepL(xShift+k) = atomXYZ0(k) 
       atomVelRepL(xShift+k) = atomVel0(k)
       end do !k  
c
       tempTg = tempTg2 
       nTimeStepRepL(nReP) = ntime
c define fileName:
         call conv_numb_to_char(nReP,nRePC)
c
         if(nReP .le. 9) then
         nRePC4 = '000'//nRePC(1:1)
          else
           if(nReP .le. 99 )then
            nRePC4 = '00'//nRePC(1:2)
             else
              if(nReP .le. 999 )then
               nRePC4 = '0'//nRePC(1:3)
               else
                if(nReP .le. 9999 )then
                 nRePC4 = nRePC(1:4)
                else
                 write(kanalp,*)
     &           'initReplica: ERROR!:number of nReP > 9999 '
                 stop
                end if
               end if
              end if
             end if
c
        kanalmdEngRepLtraf(nRep) = 50 + nRep
        mdEngRepLtraf(nRep) = eMdRepLf//nRePC4//'.tra'
        mMdRepLPDBf(nRep)   = mMdRepLf//nRePC4 
c write 
        kanalPDBres = 34
        filePdbN = mMdRepLf//nRePC4//'.0000.pdb' 
        open(unit=kanalPDBres, file=filePdbN, form='formatted',
     &       status='unknown')
c
           write(kanalPDBres,'(a31,f10.4 )')
     &     "REMARK: Md result : MdTime(ps):", ntime*deltat
          write(kanalPDBres,'(a9,i7)')'$nstep:  ', ntime
          write(kanalPDBres,'(a9,f8.1)')'$ tgT(K):',tempTg
          write(kanalPDBres,'(a8,2x,a6)')'$ENERGY: ',' :Kcal'
          write(kanalPDBres,'(a10,f14.5)')'eVbondDef:',eVbondDef
          write(kanalPDBres,'(a10,f14.5)')'eVangDef :',eVangDef
          write(kanalPDBres,'(a10,f14.5)')'eImpDef  :',eImpDef
          write(kanalPDBres,'(a10,f14.5)')'eTorsDef :',eTorsDef
          write(kanalPDBres,'(a10,f14.5)')'engVDWR1 :',engVDWR1
          write(kanalPDBres,'(a12,f14.5)')'hBHxYeng128:',hBHxYeng128
          write(kanalPDBres,'(a10,f14.5)')'engCOULR1:',engCOULR1
           write(kanalPDBres,'(a10,f14.5)')'engCOULR2:',engCOULR2
          write(kanalPDBres,'(a10,f14.5)')'restr1Eng:',restr1Eng
          write(kanalPDBres,'(a10,f14.5)')'eRstHW1ML:',restr1MHWEng
          write(kanalPDBres,'(a10,f14.5)')'eRstA2Eng:',restrDistA2Eng
          write(kanalPDBres,'(a10,f14.5)')'eCompcEng:',compactGeoEn
          write(kanalPDBres,'(a10,f14.5)')'eGeoDef  :',eGeoDef
          write(kanalPDBres,'(a10,f14.5)')'engCOUL  :',engCOUL
          write(kanalPDBres,'(a10,f14.5)')'engSolv  :',molSolEn
          write(kanalPDBres,'(a10,f14.5)')'engPOTENT:',engPOTENT
          write(kanalPDBres,'(a10,f14.5)')'kinEng   :',kinEng
          write(kanalPDBres,'(a10,f14.5)')'Temp(K)  :',tempT0
          write(kanalPDBres,'(a10,f14.5)')'engTotal :',engTotal
c
          if(vLigFlag .eq. 1) then
          write(kanalPDBres,'(a8)')'########'
          write(kanalPDBres,'(a8)')'$ENELIG: '
          write(kanalPDBres,'(a12,f14.5)')'eVbondDefLG:',eVbondDefLg
          write(kanalPDBres,'(a12,f14.5)')'eVangDefLG :',eVangDefLg
          write(kanalPDBres,'(a12,f14.5)')'eImpDefLG  :',eImpDefLg
          write(kanalPDBres,'(a12,f14.5)')'eTorsDefLG :',eTorsDefLg
          write(kanalPDBres,'(a12,f14.5)')'engVDWR1LG :',engVDWR1Lg
          write(kanalPDBres,'(a14,f14.5)')'hBHxYeng128Lg:',hBHxYeng128Lg
          write(kanalPDBres,'(a12,f14.5)')'engCOULR1LG:',engCOULR1Lg
          write(kanalPDBres,'(a12,f14.5)')'engCOULR2LG:',engCOULR2Lg
          write(kanalPDBres,'(a12,f14.5)')'restr1EngLG:',restr1EngLg
          write(kanalPDBres,'(a12,f14.5)')'eRstHW1MLLG:',restr1MHWEngLg
          write(kanalPDBres,'(a12,f14.5)')'eGeoDefLG  :',eGeoDefLg
          write(kanalPDBres,'(a12,f14.5)')'engCOULLG  :',engCOULLg
          write(kanalPDBres,'(a12,f14.5)')'engSolvLG  :',molSolEnLg
          write(kanalPDBres,'(a12,f14.5)')'engPOTENTLG:',engPOTENTLg
          end if !Lig
c
          write(kanalPDBres,'(a12)')'REMARK: PDB: '
c
           do ich=1,nChain
           do k = startAtInCha(ich),stopAtInCha(ich)
           k3 = 3*k-3
c
           doWRpdb=.true.
           if(OPT_doLigDock .eq. 1)then
           doWRpdb=.false.
           if(resName(k) .eq. resName(atomInLigList(1))) doWRpdb=.true.
           end if !OPT_doLigDock .eq. 1
c
            if(doWRpdb)then
           write(kanalPDBres,7071)"ATOM  ",
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atomXYZ0(k3+j),j=1,3), atomQ0(k)
           end if !doWRpdb
           end do !k
c TERminal chainLine
           if( (resName(stopAtInCha(ich))(1:3) .ne. 'HOH' .and.
     &        resName(stopAtInCha(ich))(1:3) .ne. 'ION') .or.
     &        ich .eq. nChain)
     &      write(kanalPDBres,7071)"TER   ", stopAtInCha(ich)
           end do !ich
c END      line
           write(kanalPDBres,7071)"END   "
c
           close(kanalPDBres)
c
c
       if(CONTROL)then
       write(kanalp,*)'initialized XYZ, Vel, at T(K) = ',
     &                 tempRepLic(nReP) 
       write(kanalp,*)'file: mdEngRepLtraf(nReP) =',mdEngRepLtraf(nReP)
       write(kanalp,*)
     & 'kanalmdEngRepLtraf(nReP)=',kanalmdEngRepLtraf(nReP)
       write(kanalp,*)'fileHead: mMdRepLPDBf(nReP) = ',mMdRepLPDBf(nReP)
       write(kanalp,*)'initRepLica: initEND nRepL:',nReP 
       end if !C
       goto 100
101    continue
c
       if(CONTROL)then
       write(kanalp,*)'initRepLica: FINISH!! * * * * * * * * * * '
       write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * '
       end if
c
       return
c
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5) ! PDB rq
7072    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,1x,3f8.3) ! PDB vel
c
       end
