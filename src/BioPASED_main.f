c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c   BioPolymer Analysis Structure Energy Dynamics                   * 
c                                                                   *
c  Yury Vorobjev 2002                                               *
c                2003                                               *
c                2004                                               *
c                2005                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	program bioPASED       
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'filedat.h'
        include 'vdw12Par.h'
        include 'loopInfo.h'
        include 'movingAtom.h'
        include 'enForce.h'
        include 'restrainInfo.h'
        include 'restrain1MHW.h'
        include 'restrain1RHW.h'
        include 'mdAtomXYZvel.h'
        include 'mdRunPar.h'
        include 'simAnneal.h'
        include 'optionPar.h'
        include 'solvGSarray.h'
        include 'optimiz.h'
        include 'xyzPDBopt.h'
        include 'ligInfo.h'
        include 'restrDistA2.h'
        include 'replicaEx.h'
        include 'rigBody01.h'
        include 'dataSASdr.h'
        include 'solvate01.h'
c
        integer iLoop
c
        integer i
        integer ia,ia3,k
        integer iwx,iwe
        integer kanalp
        logical CONTROL
c
c define kanalUse.h
        kanalSTD6 = 6
        kanalPStat = kanalSTD6
        kanalRunOut = 10
c        
        kanalp = kanalRunOut
c
        kanalInPdb = 14
        kanalInMdynPar = 15
        kanalInMoveRes = 16
        kanalInSAprot = 17
        kanalInRestrALLType = 18
c
        kanal_htg = 21
        kanal_AmbZmConv = 22
        kanal_ResTop01 = 24
        kanal_ResTop02 = 25
        kanal_ResTop03 = 26
        kanal_solvAtType = 27
        kanal_solvGSdat = 28
        kanal_watBox0 = 29
c
c Results:
        kanalXYZtra = 31
        kanalEngTra = 32
        kanalPDBres = 33
        kanalPDBfin = 34
        kanalpHB128 = 35
        kanalwWatBrg = 36
        kanalModXYZ = 37
        engTablkanaL1= 38
        engTablkanaL2 = 39
        kanalModFreq = 41
c        
        call getArgLine
c
c init Messages
        iMolTopoMX = 1
        iFFbuildMX = 1
        iengOptimizMX = 14
        imdStartRunMX = 3
        iSolvModMX = 1
        imdSArunMX = 85
        iCPUMX = 100
c
        mError = 'Error:'
        messageStatus = ' run BioPASED ...'
        mStatus = 'Status:'
c
        iStatus = 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
        CONTROL = .true.
c
        write(kanalp,*)'BioPASED (by Yuri Vorobjev): START: '
        write(kanalp,*)
c
        call inputMDSApar
c
        call initMolecTopSeq01
c
        iStatus = iStatus + iMolTopoMX
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' made molec topology ...'
c
c read molec.pdb for complete sequence of Res., defines XYZ atomList(*)
c
c store atomXYZinitial
        call initSwap01(natom,atomXYZ,atomXYZin)
c
        call initLoopMoveAtList(nMoveLoop,moveLoopList,
     &                OPT_LoopMD)
c
c init readRestr1ResHW:
        i=0
        nRestr1RHWSegMX = nRestr1RHWSegMAX
        call readRestr1ResHW
     &  (i,restr1RHWFile,nRestr1RHWSegMX,nRestr1RHWSeg,resEndRestr1RHW)
c
c counterIon shell
        if(OPT_coIonShell)then
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' made Ion solvation shell ..'
c
        call ionShellSAS03(atomXYZ)
c
        end if !OPT_coIonShell
c
c solvationExplicitWaterShell:
        if(OPT_SolvateExWat)then 
        write(kanalp,*)'DO solvateMoL01 :'
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' start Water solvation shell ..'
c
        call  solvateMoL01
c
        iStatus = iStatus + iSolvModMX-1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' Water solvation shell is done ...'
        end if !OPT_SolvateExWat
c PWeX: define atomPWatTag() = 0/1 (solute/solvent)
        call initPWatAtomTag
c
        call initMolecTopSeq02 
c
        call initSwap01(natom,atomXYZ,atomXYZin)
c
        iStatus = iStatus + iMolTopoMX
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' made solvated molec topology ..'
c
c  solvate03 : Water Bridges
        if(OPT_SolvEFS)then
        write(kanalp,*)'DO solvateMoL03 :' 
        call solvateMoL03(atomXYZ)
        write(kanalp,*)'STOP solvateMoL03 :' 

        call getWatBrgForceSAS03(atomXYZ,
     &                   watBrgEnergySoL,watBrgForceAt)
c
        call getHPolvForceSAS03(atomXYZ,
     &                   hpSoLEnergy,hpSoLForceAt)
c
        call writewBrgWmoLXYZ 
        if(OPT_writeSASdotXYZ) call writeDotSASXYZ(0)
c  
        iStatus = iStatus + iSolvModMX
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' WatBridges is done ...'
c        
        end if !OPT_SolvEFS
c
c BornRad
        if(OPT_SolvGBorn)then
        call initBornRadSAS(atomXYZ)
c
         call getBornPolzEforce(atomXYZ,
     &               bornPolzEng,bornPolzEngAt,bornPolzForceAt)
c
        iStatus = iStatus + iSolvModMX
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' init GenBorn solvation model ..'
c
        end if !OPT_SolvGBorn
c
c initLigand: ligInfo.h data
        nAtomInLig = 0
        ncallNBPLLig = 0
        if(OPT_LigRes)then
        nLigdMX = nLigdMAX
        call readLigInf(ligResFile,nLigdMX,nLigd,resStEndLig)
c
        call initLigAtList
c
        call initLigFlag1234
c
c write molAtXYZ In the vicinity of Lig 
        if(OPT_WRmolAtIntLig)
     &       call wrAtXYZinLigLoc
c
        end if !OPT_LigRes
c
c init FField and SolvationGSmod for the current PDBXYZ
c
        call initFFieldParam
c
        iStatus = iStatus + iFFbuildMX
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' init ForceField parameters ..'
c
c init restraints
c define data restrainInfo.h { resr1AtomList() etc.} for harmonically restrained atoms
        nRestr1Seg = 0
        nRestr1atom = 0
c
        if(OPT_harmAt1PosRst)then
        nRestr1SegMX = nRestr1SegMAX
c yv180209
        call readRestr1AtInf09
     &             (restr1File,nRestr1SegMX,nRestr1Seg,resEndRestr1,
     &              resResSegKhConst,resRest1AtLine)
c
        call initRestr1AtList09(nRestr1Seg,resEndRestr1,
     &                 resResSegKhConst,resRest1AtLine,
     &                 nRestr1atom,restr1atomList,refAtomPos,
     &                 restr1AtConst,restr1AtConstList)
        end if !OPT_harm1APosRst for mainMoleculeAtoms   
c
c init targetAtomPositional restraints:
        if(OPT_TargAtPos)then
c
        call readTargAtPos(targAtPosFile,natomMX,nresMX,
     &             natomTargAtPos,atNumbTargAtPos,
     &             atNameTargAtPos,resNameTargAtPos,
     &             resNumbTargAtPos,atXYZTargAtPos)
c
        call initTargAtPosList(natomTargAtPos,atNumbTargAtPos,
     &             atNameTargAtPos,resNameTargAtPos,
     &             resNumbTargAtPos,atXYZTargAtPos,
     &             nRestr1atom,restr1atomList,refAtomPos,
     &             targAtPosHConst,restr1AtConstList)
        end if!OPT_targAtPos
c end init targetAtomPositional restraints
c
c define data restrain1MHW.h [restr1MHWatList() etc. for MoleculeCMRestraints
        nRestr1MHWSeg = 0
        nRestr1MHWatom = 0
c
        if(OPT_HW1MolPosRst)then
        nRestr1MHWSegMX = nRestr1MHWSegMAX
c
        call readRestr1MolHW
     &  (restr1MHWFile,nRestr1MHWSegMX,nRestr1MHWSeg,resEndRestr1MHW)
c
        call initRestr1MHWAtList(nRestr1MHWSeg,resEndRestr1MHW,
     &                 nRestr1MHWatom,restr1MHWatList,ref1MHWPos,
     &                 atMs1MHWatList)
c
        end if ! OPT_HW1MolPosRs
c
c restraintsRES1pos: - HarmWallRestraints for RES c.m.
c
        if(OPT_HW1ResPosRst)then
c read restr1RES file: restrRes1.inp for RESidues of MAIN molecule!!
c solvationShell water are included in the resEndRestr1RHW(*) data
        if(.not.OPT_SolvateExWat)then
        nRestr1RHWSeg = 0
        nRestr1RHWatom = 0
        end if !.not.OPT_SolvateExWat
c
        nRestr1RHWSegMX = nRestr1RHWSegMAX
        i=1
        call readRestr1ResHW
     &  (i,restr1RHWFile,nRestr1RHWSegMX,nRestr1RHWSeg,resEndRestr1RHW)
c
        end if !OPT_HW1ResPosRst
c
cx        if(OPT_SolvateExWat .or. OPT_HW1ResPosRst)then
        if(OPT_HW1ResPosRst)then
c make rest1HWRES atoms list for  water in SolvShell
        call initRestr1RHWAtList(atomXYZ,nRestr1RHWSeg,resEndRestr1RHW,
     &        nRestr1ResidHW,restr1ResidHWList,
     &        nRestr1RHWAtom,restr1RHWatomList,
     &        restr1RefResidPosHW,restr1ResidMass)

        end if ! OPT_HW1ResPosRs
c
        nRestrDistA2 = 0
        if(OPT_restrDistA2)then
        nRestrDistA2MX = nRestrDistA2MAX
c
cx        write(kanalp,*)'BioPASED_main: restrA2File:',restrA2File
c
         call readRestrDistA2
     &           (restrA2File,nRestrDistA2MX,nRestrDistA2,
     &                       restrDistA2List,restrDistA2DistHK)  
c
        write(kanalp,*)'BioPASED_main: nRestrDistA2:',nRestrDistA2
        write(kanalp,*)'BioPASED_main: restrDistA2List:',
     &  (restrDistA2List(2*i-1), restrDistA2List(2*i),
     &   restrDistA2DistHK(2*i-1),restrDistA2DistHK(2*i),
     &   i=1,nRestrDistA2)
        write(kanalp,*)'- - - - - - - - - - - - - - - - - - - - - -'
c
         end if !OPT_restrDistA2
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c initRigidBody distribution:
         if(OPT_RigidBodyMd .or. OPT_RigidBodyCADist) then
c
         nRigBodyMX = nRigBodyMAX
         call  readRigBodyRes
     &             (rigBodyFile,nRigBodyMX,nRigBody,rigBodyStEndRes)
c
         call initRigBodyAtList(nRigBody,rigBodyStEndRes,
     &        nAtRigBody,nAtRigBodySeg,startAtInRGB,stopAtInRGB,
     &        atRigBodyFlag)
c
         end if !OPT_RigidBodyMd
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
          if(OPT_RigidBodyCADist) then
       call initRigBodyCA2List(cAcARigBodyHK,nRigBody,
     &                          startAtInRGB,stopAtInRGB,
     &                          atomName,atomXYZ,
     &                          nRestrDistA2MX,nRestrDistA2,
     &                          restrDistA2List,restrDistA2DistHK)
          end if ! OPT_RigidBodyCADist
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c init SolvationGSmod for the current PDBXYZ  
c        
         if(OPT_SolvGS)then
          call initSolvatGSmod
         iStatus = iStatus + iSolvModMX
         write(kanalPStat,*)mStatus,iStatus,messageStatus
     &   ,' init Gauss Shell solvation model ..'
         end if !OPT_SolvGS 

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	if(OPT_engCalc .and. .not. OPT_mdRestart)then
c calculate energy and forces for given conformation XYZ
        fcall = 0
        ncallNBPL = 0
c make atomXYZopt() = atomXYZ()
        call initSwap01(natom,atomXYZ,atomXYZopt)
	write(kanalp,*)'engCalculation: OPT_engCalc: Starts:'
c
        iwe = 1
        iwx = 1
        call engOptFunct(moveAtomXYZ,engPOTENT,atForceTot,iwe,iwx)
c
        iStatus = iStatus + 1             
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' initialXYZ energy calculation done ...'
c
c control hb128
        call initSwap01(natom,atomXYZ,atomXYZ0)
        call printHB128EngForces
c
        if(OPT_SolvEFS) call writewBrgWmoLXYZ  
        if(OPT_writeSASdotXYZ) call writeDotSASXYZ(1)
c
        end if !OPT_engCalc
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        if(OPT_engOptim .and. .not. OPT_mdRestart)then
c energyOptimization 
        fcall = 0
        ncallNBPL = 0
c init atomXYZopt:  initial point
        call initSwap01(natom,atomXYZ,atomXYZopt)
c
        write(kanalp,*)'flatcher_Optim: Starts: nOptIterMax:',nOptIter
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' start energyOptimization ...'
c
        nOptXYZ = 3*nmoveatom
        call flatcher_Optim(nOptXYZ,moveAtomXYZ,
     &                       engEps,nOptIter,engPotMin)
c
        call initSwap01(nmoveAtom,moveAtomXYZ,atomXYZmin)
c
        call initSwap01(natom,atomXYZopt,atomXYZ)
c control hb128
        call initSwap01(natom,atomXYZopt,atomXYZ0)
        call printHB128EngForces
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' energyOpimization is done ...'
c
        if(OPT_SolvEFS) call writewBrgWmoLXYZ 
        if(OPT_writeSASdotXYZ) call writeDotSASXYZ(1)
c
        write(kanalp,*)'engOptimization Final: engPotMin:',engPotMin
c
        end if!OPT_engOptim  
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c init MD run
        nRecPdb = 0      !number of MdResNNNN.pdb
        ntime0 = 0
        ntime = 0
        optra = 1   !open file : './job/xyzMd.tra'
        wtra = 1
        cltra = 0
c do MolDyn
         if(OPT_doMD)then
c
        if(OPT_mdRestart)then
        write(kanalp,*)'********************************************'
        write(kanalp,*)'MD RESTART: from mdXYZVin.pdb !'
        write(kanalp,*)'********************************************'
        end if ! mdRestart
c
        call initMDStart(tempT0)
c
        iStatus = iStatus + imdStartRunMX
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' start molDyn run ...'
c
        if(OPT_RigidBodyMd)then
        call initRigBodyMdStart
        end if !OPT_RigBodyMd
c
        call mdRun(ntimeMX,ntime0,ntime,ntimeR1,ntimeR2,
     &             ntimeF1,ntimeF2,ntimeF3,deltat,
     &             tempT0,tempTg,tauTRF,atype,optra,wtra,nwtra,cltra)
c
c get new atomXYZ(*) after MD
        call initSwap01(natom,atomXYZ0,atomXYZ)
c
        call printHB128EngForces 
c
        if(OPT_SolvEFS) call writewBrgWmoLXYZ 
c
        if(OPT_writeSASdotXYZ) then
        if(.not.OPT_SolvateExWat)call writeDotSASXYZ(1)
        if(OPT_SolvateExWat)call writeDotSASXYZ(0)
        end if!OPT_writeSASdotXYZ
c
        iStatus = iStatus + 1                
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' eqvilibration mDyn is done ...'
c        
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c simulated Annealing or arbitrary heating/cooling schedule
c
        if(OPT_MDSA) then
c
        if(nSAstep .gt. nSAstepMAX)then
        write(kanalp,*)'ERROR!!: nSAstep .gt. nSAstepMAX!!'
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' simulated annealing starts ...'
c
        write(kanalPStat,*)mError,
     &    ' (-sa saProtocol file) has too many stages ! '
        stop
        end if!nSAstep
c
	call simAnnealing(nSAstep,nSAparMX,SAProtcol)
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' simulated annealing is done ...'
c
        end if !MDSA
c
c get new atomXYZ(*) after MD
        call initSwap01(natom,atomXYZ0,atomXYZ)
        call printHB128EngForces 
c
        if(OPT_SolvEFS) call writewBrgWmoLXYZ 
c
        if(OPT_writeSASdotXYZ) then
	if(.not.OPT_SolvateExWat)call writeDotSASXYZ(1)
	if(OPT_SolvateExWat)call writeDotSASXYZ(0)
	end if!OPT_writeSASdotXYZ
c
        end if !do_MD
c
        if(OPT_EssModeAnalys)then
c
         write(kanalp,*)'start Essential mdModeAnalys01'

         call  mdModeAnalys01

         write(kanalp,*)'END mdModeAnalys01'
c
        if(OPT_freeEnergy) then
        call mdFreeEnergy01
        end if ! mdfreeEng
c
        end if! EssMode
c 
        write(kanalp,*)' bioPASED is successivelly finished !'
c
        iStatus = iStatus + 1
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' bioPASED is successivelly finished !' 
c
7071   format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,2x,f6.3,f8.3) ! PDB
c
         stop
	 end
