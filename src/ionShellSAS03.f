c put mobile IONs around SoluteMolecule 
c to the IonAccessibleSurface points of MAX ElectPotetial
c Yury Vorobjev 2005
c
	subroutine ionShellSAS03(atomXYZ)
c
c READ in:  molecPDBxyz, etc.,          
c
        include 'xyzPDBsize.h'
        real atomXYZ(*)
        include 'xyzPDBinfo.h'
        include 'coulEnPar.h'
        include 'dataIonAS.h'
        include 'ionShell.h'  
        include 'optionPar.h'
        include 'movingAtom.h'
        include "filedat.h"
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer natomMOL,nresMOL
        real big,engIon
        logical positiveIonShell
        real  ePotMIN
        real totalMolecQ
        integer iw,iw3,i,i3
        integer n3,ia3,im3
        integer id1,id3,k,id3mx
        integer ndotIonMX
        integer natomM0,nresM0
        logical ePotMXfound
        integer kanalp,kanalXYZion
        logical CONTROL,CONTROL0
c
        kanalp=kanalRunOut 
        CONTROL  = .true. 
        CONTROL0 = .true. 
c OPT_ variables:
        eIonBindMIN_OPT = -1.2 !min IonBindEpot
c
c get total ProteinQ
        big = 1.0e10
        totalMolecQ = 0.0
        do i=1,natom
        totalMolecQ=totalMolecQ+atomQ(i)
        end do
c
        positiveIonShell=.false.
        if(totalMolecQ .lt. 0.0)positiveIonShell=.true.
        nCounterIonMX=abs(totalMolecQ)
c
        write(kanalp,*)'ionShellSAS03: totalMolecQ:',totalMolecQ
        write(kanalp,*)'ionShellSAS03: nCounterIonMX:',nCounterIonMX
        if(positiveIonShell)then
        write(kanalp,*)'positive CounterIons added ...'
        else
        write(kanalp,*)'NO positive CounterIons will be added ...' 
        end if
c
        if(positiveIonShell)then
        write(kanalp,*)'nCounterIons (positive) will be added ...' 
        iStatus = iStatus + 1               
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' add counterIon shell ',nCounterIonMX, ' ions...'
        else
        write(kanalp,*)'NO positive CounterIons will be added ...'
        end if
c
c
        if(nCounterIonMX .eq. 0)return
        if(.not.positiveIonShell)return
        
c calculate IonAS
        ndotIonMX = ndotIasMAX
        do id1=1,ndotIonMX
        dotIon_ePot_occFlag(id1)=0
        end do !id
c
	call getAtomRadSAS03(natom,atomName,atomRad) 
c
        dotIonden = 4.0  !inputMDpar.f
        ionMolRad_OPT = 3.4
        dprobeIon = ionMolRad_OPT  ! = wDiam + ionRad ~ 2.8 + 0.6 = 3.4
c
        nIonAScall = 0
c dotIon_IonAS and dotIon_ePOT/eField
cx        if(nIonAScall .eq. 0 )then
c 
               call  surf_SAS04(atomXYZ,atomRad,
     &                 natom,dprobeIon,dotIonden,
     &                 dotIonXYZ(1),dotIonnrm(1),dotIonarea(1),
     &                 dotIon_IATNUM(1),ndotIon,ndotIonMX,
     &                 nsurfIonAt,nsurfIonAtList(1),
     &                 atSurfIonAr(1),atSurfIonNrm(1),
     &                 atSurfIonXYZ(1),atSurfIonProj(1),
     &                 bindSiteIonAt01(1),
     &                 dotIon_jATNUM(1),
     &                 head_dotIonNum,linkListDotIonAt,
     &                 atDistToIonAS,atIonASexp,
     &                 dotIon_eField,dotIon_ePot)
c
          nIonAScall = nIonAScall + 1
c
          write(kanalp,*)
     &    'ionShellSAS03: surf_SAS04 call: nIonAScall=',nIonAScall
c
cx          end if !
c
c 2) calculate MAX of ePot on the IonAS from MOLECatoms
c 3) put ion to the MAX of ePot if ePotMAX > ePotPOROG
c 4) correct potential/efield on IonAS taking the new IOn
c 5) goto 2
c 
c put ONLY Positive  univalence solvatedIons [Na,K,..]
c
               nIonTypeMX=1
               nIonType = 1
               zarIon(nIonType)=1.0
c
               nIonNow(nIonType) = 0
               iw = 1
c
1001           ePotdotMX(iw) = big
               iDotePotdotMX(iw)=0
c
               if(CONTROL)then
               write(kanalp,*)
     &         'ionShellSAS03: ndotIon:',ndotIon,' tryNew ion:',iw
               end if!C
c
               ePotMXfound = .false.
c find MIN ePot for ion on IonAS
c
            do id1=1,ndotIon
            if(dotIon_ePot_occFlag(id1) .eq. 0 )then
c         
            engIon = dotIon_ePot(id1)*zarIon(nIonType)
            if(engIon .lt. ePotdotMX(iw))then
            ePotdotMX(iw)=engIon                               
            iDotePotdotMX(iw)=id1
            ePotMXfound=.true.
            end if
c
            end if !occFlag
            end do!id1
c
	       if( .not. ePotMXfound)goto 1002
c               
               if(CONTROL)then
               write(kanalp,*)
     &         'ionShellSAS03 nIonNow(nIonType)=',nIonNow(nIonType),
     &         ' iDotePotdotMX(iw)=',iDotePotdotMX(iw),
     &         ' ePotdotMX(iw)=',ePotdotMX(iw)
               end if !C
c
c put the implicit water molecule
               if(ePotdotMX(iw) .lt. eIonBindMIN_OPT    
     &                  .and. (iw .le. nCounterIonMX)  
     &                   .and. (iw .le. nIonIASMAX) )then
c newIon accepted:
               ionIasEpot(iw) = ePotdotMX(iw)
               iw3=iw*3-3
               id3mx=3*iDotePotdotMX(iw)-3
c
               do k=1,3
               ionIASXYZ(iw3+k)=dotIonXYZ(id3mx+k)
               end do !k
               ionAtomQ(iw) = zarIon(nIonType)
c
               if(CONTROL)then
               write(kanalp,7003)
     &   'IONAS ',iDotePotdotMX(iw),
     &   'ATN','iAS',dotIon_IATNUM(iDotePotdotMX(iw)),
     &   (ionIASXYZ(iw3+k),k=1,3),ePotdotMX(iw)
               end if !C
c
               nIonNow(nIonType) = iw
c
c correct ePot and calculate dotIon_ePot_occFlag(*) due to newIon=iw
c
           call getIonEpotSAS03(ionIASXYZ(iw3+1),zarIon(nIonType),
     &                   ndotIon,dotIonXYZ(1),dotIon_ePot,
     &                                        dotIon_ePot_occFlag)
c
           if(CONTROL)then
           write(kanalp,*)'ionShellSAS03: eLpot on IonAS is corrected:'
           end if
c
           iw = iw+1 
           goto 1001     !gotoNextCandidate 
c
           else
           goto 1002  ! stop WB search
           end if  ! ePotdotMX(iw) .gt. eIonBindMIN_OPT
c
1002           continue
c
c
c add ION to the main atomXYZ(*) and correct MolTopo files:
c
         do i=1,nIonTypeMAX
         ionResName(i)='ION '
         ionAtomName(i)='NA+ '
         ionffatomName(i)='IB'
         end do !
c 
          nresM0=nres
          natomM0=natom
          do iw=1,nIonNow(nIonType)
          iw3=3*iw-3
          natom = natom + 1
          n3 = 3*natom-3
          nres=nres+1
          do k=1,3
          atomXYZ(n3+k)=ionIASXYZ(iw3+k)
          end do!k   
          atomQ(natom) = ionAtomQ(iw)
          atomName(natom)=ionAtomName(nIonType) 
          resName(natom)= ionResName(nIonType)
          atomNumb(natom) = natom
          ffAtomName(natom)=ionffatomName(nIonType)
          resNumb(natom)=nres
          resNameRes(nres)=ionResName(nIonType)
          chName(natom)=' '
          chNameRes(nres)=' '
          resListSeq(nres)=resNameRes(nres)
          resNameEx(natom)=resName(natom)
          atomNameEx(natom)=atomName(natom)
          atHvyNbInRes(nres)=1
c
          startAtInRes(nres) = natom
          stopAtInRes(nres) = natom
          startAtInRes(nres+1)=natom+1
c
          nmoveatom = nmoveatom + 1
          moveAtomList(nmoveatom)=natom 
          moveFlag(moveAtomList(nmoveatom)) = 1
          realAtomFlag(moveAtomList(nmoveatom)) = 1
c
          im3 = nmoveatom*3-3 
          ia3 =  3*moveAtomList(nmoveatom)-3
          moveAtomXYZ(im3+1) = atomXYZ(ia3+1)
          moveAtomXYZ(im3+2) = atomXYZ(ia3+2)
          moveAtomXYZ(im3+3) = atomXYZ(ia3+3)
c
c correct ionEnergies: ionIasEpot(iw)
          ionIasEpot(iw)=dotIon_ePot(iDotePotdotMX(iw))*zarIon(nIonType) 
c
          end do!iw
c
c update nChain and startAtInCha(ich),stopAtInCha(ich)
c coIonShell = new chain
          nChain = nChain + 1
          startAtInCha(nChain) = stopAtInCha(nChain-1) + 1
          stopAtInCha(nChain) = natom
          stopAtInMol(nMolec) = natom
c END correct TOPology
c
c
c write file ionSheLLXYZ.pdb
        write(kanalp,*)'write coIonSHeLLXYZ: '
c
        fileIonASXYZ = 'coIonSheLXYZin.pdb'
        if(nMolNameLett .ge. 1)
     &  fileIonASXYZ = 
     &   molNameArgLine(1:nMolNameLett)//'.coIonSheLXYZin.pdb'
c
        if(OPT_resFileDir)then
        fileIonASXYZ=resultFileDir(1:nResFileDirLett)
     &                                 //'coIonSheLXYZin.pdb'
        if(nMolNameLett .ge. 1)then
        fileIonASXYZ=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.coIonSheLXYZin.pdb'
        end if
        end if
c
        kanalXYZion = kanalwWatBrg
c
        open(unit=kanalXYZion,file=fileIonASXYZ, form='formatted',
     &       status='unknown')
c 
        if(nIonNow(nIonType) .ge. 1)then
        write(kanalXYZion,*)
     &  '#coIonShell.pdb:                XYZ         ePotIon(kcal/mol)',
     &  '   ionZ'
          do iw=1,nIonNow(nIonType)
          iw3=3*iw-3
          write(kanalXYZion,7003)
     &   'ATOM  ',iw,ionAtomName(nIonType),ionResName(nIonType),
     &   (nresM0+iw),(ionIASXYZ(iw3+k),k=1,3),
     &   ionIasEpot(iw),ionAtomQ(iw)
               end do !iw
c
           write(kanalXYZion,'(a4)')'TER '
           write(kanalXYZion,'(a4)')'END '
c
          end if !nIonNow(nIonType .ge. 1
c
          close(kanalXYZion)
c
          if(CONTROL0)then
          write(kanalp,*)'ionShellSAS03 final ionShell XYZ  ePotIon  Q:'
          do iw=1,nIonNow(nIonType)
          iw3=3*iw-3
          write(kanalp,7003)
     &   'ATOMi ',iw,'ion','ShL',nIonType,
     &   (ionIASXYZ(iw3+k),k=1,3),ionIasEpot(iw),ionAtomQ(iw)
           end do !iw
           end if!C
c
          write(kanalp,*)'ionShellSAS03: Finish!'
c
        return
7003   format(a4,1x,i6,2x,a4,a4,1x,i4,4x,3f8.3,1x,2f8.3) ! PDB
        end
