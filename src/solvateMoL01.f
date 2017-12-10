c generates water solvent around SoluteMolecule 
c Yury Vorobjev 2005
c
	subroutine solvateMoL01
c
c READ in:  molecPDBxyz, etc.,          
c
c OUT: PDB file of molec.pdb + watShell
c     name: fileWatSheLL = mol.watSheLLXYZ0.pdb : default
c           fileWatSheLL = molName.watSheLLXYZ0.pdb : if there ismolecule name 
c
c dSOLVSheLL: solvationSheLL thickness
c
c        implicit none
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBcrd.h'
        include 'xyzPDBinfo.h'
        include 'filedat.h'
        include 'solvate01.h'
        include 'movingAtom.h'
        include 'pair1234array.h'
        include 'restrain1RHW.h'
	include 'enForce_PW.h'
c local
        integer nWatMolAtAdd01
        integer nWatMolAdd01
        integer natMOLandWat
        integer natomMOL,nresMOL
c
c soluteMOL+Water Heavy atoms: noH
        integer natomMOLnoH
        integer natMOLandWatnoH
        integer natMOLandWatnoHMAX
        parameter (natMOLandWatnoHMAX=2*natomMAX/3+nWaterMolMAX)
        real atomXYZMolWatnoH(3*natMOLandWatnoHMAX)
        character*4 atomNamenoH(natMOLandWatnoHMAX)
        character*4 resNamenoH(natMOLandWatnoHMAX)
        integer resNumbnoH(natMOLandWatnoHMAX)
c WatBOX0 translation
        integer nWat126Load(3)
        real    molecMinMaxBox(6)
        real    molWatMinMaxBox(6) 
        real molWatBoxSize(3)
c
        real dSHELLMIN,dSHELL
c
c neihbourwPairList
        integer ncallwPL        
        real rCUTwPL
c local WpairList dataStructure
        integer natomLocMAX
        parameter (natomLocMAX=natMOLandWatnoHMAX)
        integer nnbwPLMAX
        parameter (nnbwPLMAX = natomLocMAX*100) ! rcutMAX < 5.5 A
        integer nnbWPairL(natomLocMAX)
        integer startnbWPairL(natomLocMAX)
        integer nbWPairList(nnbwPLMAX)
        real nbWPairListD2(nnbwPLMAX)
        integer nnbWpLT,nnbWpLMX
c
        integer ihsRadMAX
        parameter (ihsRadMAX=5)
        integer ihsRad
        real hsRad(ihsRadMAX)
c local
        character*4 watResName
        integer kanalWATbox0
        integer nWatAtinFile
        character*6 head
        character*80 line
c
        real watBox216size,watBox216size2
        real dAtMolWatMin,dAtMolWatMax
        real dAtMolWatMin2,dAtMolWatMax2
        real big,low,d
        real xyzi(3),xyzj(3)
        real xMAX(3),xMIN(3)
        real wXYZ(3)
        real watBoxShift(3)
c
        logical loadWat
        integer nsiz(3),ncell,ncellMX
        integer i,j,k,k3,n3,ir
        integer i1,i2,i3
        integer ia,ia3,iam,im3
        integer ja,jn
        integer wa1,wa3,ich
        integer ix,iy,iz
        integer ixm,iym,izm
        integer ix0,iy0,iz0
        integer iwM,iwat,iw3,iwa,iw
        logical readAtom 
        logical CONTROL0,CONTROL
        integer kanalp
        integer  kanalXYZwSheLL
c
        kanalp=kanalRunOut 
c
        CONTROL0 = .true.  ! is on usually
        CONTROL  = .true. 
c
        hsRad(1) = 1.4   ! WATer
        hsRad(2) = 1.75  ! C
        hsRad(3) = 1.40  ! O
        hsRad(4) = 1.40  ! N
        hsRad(5) = 1.80  ! P,S,
c init/correct
c
        write(kanalp,*)'solvate01: natom:',natom
        write(kanalp,*)'solvate: nres:',nres
        write(kanalp,*)'solvate: nChain:',nChain
        write(kanalp,*)'solvate: nMolec:',nMolec
c
        watBox216size = 18.6206 !watBox for spc216.amb.dat file
        watBox216size2 = watBox216size/2.0
c
        dAtMolWatMin = 2.8     ! minDistWat-Atom
        dSHELLMIN =  3.3     !minWaterSheLL
        dSHELL = dSOLVSheLL
        if(dSHELL .lt. dSHELLMIN)dSHELL=dSHELLMIN
c
        dAtMolWatMax = dSHELL       
        dAtMolWatMin2=dAtMolWatMin**2
        dAtMolWatMax2=dAtMolWatMax**2
c
c init:
c define Protein or Water atom Tag
         do ia=1,natom
	 atomPWatTag(ia)=1       !1 = Solute molecule
	 end do !i
	 do ia=natom+1,natomMAX
	 atomPWatTag(ia)=0
	 end do!ia
c
        if(CONTROL)then
        write(kanalp,*)'solvateMoL01: make hydrSHELL(A):',dSHELL,
     &  ' fileWATBOX0:',fileWATBOX0
        write(kanalp,*)'nChain:',nChain 
        end if
c read in water mol in water box216 :
c solvate.h 
        kanalWATbox0=kanal_watBox0
c
        open(unit=kanalWATbox0,file=fileWATBOX0,form='formatted',
     &      status='old' )
c read *.pdb file ~/dat/spc216.amb.dat
c       
        nWatAtinFile=0
        k=0
100     read(kanalWATbox0,'(a80)',end=101)line
        readAtom=.false.
        if(line(1:4).eq.'ATOM' )then   !ATOM
           readAtom=.true.
        end if !
c
        line(76:80) = ' Orig'
        write(kanalp,*)line
c
        if(readAtom)then
        k=k+1
        k3=k*3-3
        read(line,7071)head,
     &     watAtNumbInBox0(k), watAtNameInBox0(k),
     &     watResNameInBox0(k),watchNameInBox0(k),
     &     watResNumbInBox0(k),(watXYZinBox0(k3+i),i=1,3)
           if(chName(k) .eq. ' ')chName(k)='W' 
        end if !readAtom
        goto 100
101     continue
c
        close(kanalWATbox0)
c
        nWatAtinFile=k
        if(nWatAtInBox0Max .ne. nWatAtinFile)then
        write(kanalp,*)' solvateMoL01: ERROR! read  spc216.amb.dat:',
     &  ' nWatAtinFile.ne.nWatAtInBox0Max:',
     &  nWatAtinFile,nWatAtInBox0Max
c
        write(kanalPStat,*)mError, 
     &  ' solvateMoL01: ERROR! read  spc216.amb.dat:',
     &  ' nWatAtinFile.ne.nWatAtInBox0Max:'
c
        stop
        end if
c end of file
cinit
c define moleculeSize
        big = 10.0e10
        low = -big
        do k=1,3
        xyzi(k)=0.0
        xMAX(k)=low
        xMIN(k)=big
        end do !k
c
        natomMOL=natom
	natomSMOL=natom
        nresMOL = nres
	nresSMOL = nres
c
        do ia=1,natomMOL
        i3=3*ia-3
        do k=1,3
        if(xMAX(k) .lt. atomXYZ(i3+k)) xMAX(k)=atomXYZ(i3+k)
        if(xMIN(k) .gt. atomXYZ(i3+k)) xMIN(k)=atomXYZ(i3+k)
        end do!k
c
cx        if(CONTROL)then
cx        write(kanalp,*)'at i, xyz:',ia,(atomXYZ(i3+k),k=1,3)
cx        end if!C
c
        end do!ia
c
        do k=1,3 
        molecMinMaxBox(k)=xMIN(k)
        molecMinMaxBox(3+k)=xMAX(k)
        molWatMinMaxBox(k)=xMIN(k)-dSHELL
        molWatMinMaxBox(3+k)=xMAX(k)+dSHELL
        molWatBoxSize(k)=molWatMinMaxBox(3+k)-molWatMinMaxBox(k)  
        nWat126Load(k)=molWatBoxSize(k)/watBox216size +1
        end do !k
c molecule size
        if(CONTROL)then
        write(kanalp,*)'solvateMoL01 :  MolSize:'
        write(kanalp,*)'solvateMoL01 :  xMAX:',xMAX
        write(kanalp,*)'solvateMoL01 :  xMIN:',xMIN
        write(kanalp,*)'solvateMoL01 : molWatBoxSize():',molWatBoxSize
        write(kanalp,*)'solvateMoL01 : nWat126Load()  :',nWat126Load
        end if
c load atomXYZ               
        natMOLandWat=natomMOL
        nWatMolAtAdd01 = 0
        nWatMolAdd01 = 0   ! total number of water MOLec in BOX
        
c load watBox216
        do i1=1,nWat126Load(1)
        do i2=1,nWat126Load(2)
        do i3=1,nWat126Load(3)
        watBoxShift(1)=i1*watBox216size - watBox216size2
     &                 + molWatMinMaxBox(1)
        watBoxShift(2)=i2*watBox216size - watBox216size2
     &                 + molWatMinMaxBox(2)
        watBoxShift(3)=i3*watBox216size - watBox216size2
     &                 + molWatMinMaxBox(3)
c
c load water mol from the BOX i1,i2,i3
        do iwM=1,nWatInBox0MAX
        loadWat = .true.
cendUpdate
c we USe the special order O,H1,H2 atoms for water !
cx check the O atom = 1 first atom of WAterMol    
c iwat - watAtNumb in WatBox0
cx        iwat=natIn1watMol*(iwM-1)+1      
cx        iw3=3*iwat-3    
cx iwat - number of water atom in BOX0
c
        if(i1.eq.nWat126Load(1) .or. i2 .eq. nWat126Load(2) 
     &     .or. i3.eq.nWat126Load(3) )then
c lastWaterBox check inBox condition for O of waterMol
        do k = 1,3
        wXYZ(k) = watXYZinBox0(iw3+k) +  watBoxShift(k)
        if(watAtNameInBox0(iwat)(1:1) .eq. 'O' .and.  
     &    (wXYZ(k) .gt. molWatMinMaxBox(3+k)) ) loadWat=.false.
         end do!k
         end if !i1.eq.nWat126Load(1)
c
c load water MOL to SOLVation BOX  
        if(loadWat)then
        nWatMolAdd01=nWatMolAdd01+1
        do iwa = 1,natIn1watMol
        iwat=natIn1watMol*iwM-natIn1watMol+iwa
        iw3=3*iwat-3
        nWatMolAtAdd01 = nWatMolAtAdd01 + 1
        natMOLandWat = natMOLandWat + 1
c name
        if(nWatMolAtAdd01 .ge. nAtomWatMAX)then
        write(kanalp,*)
     &  'ERROR:solvateMoL01: nWatMolAtAdd01 is tooLarge>nAtomWatMAX:',
     &   nWatMolAtAdd01,' increase nWaterMolMAX in solvate01.h !'
c
        write(kanalPStat,*)mError,
     &  'ERROR:solvateMoL01: nWatMolAtAdd01 is tooLarge>nAtomWatMAX:'
c
        stop
c
        end if !nWatMolAtAdd01 .ge. nAtomWatMAX
        atomNameWat(nWatMolAtAdd01)=watAtNameInBox0(iwat)
        resNameWat(nWatMolAtAdd01)=watResNameInBox0(iwat)
        resNumbWat(nWatMolAtAdd01)=nWatMolAdd01
c XYZ
        n3=3*nWatMolAtAdd01-3
        do k =1,3
         wXYZ(k) = watXYZinBox0(iw3+k) +  watBoxShift(k)
         atomXYZWat(n3+k) = wXYZ(k)   
        end do !k
c
        if(CONTROL)then
        ia=nWatMolAtAdd01
        ia3=ia*3-3
        write(kanalp,7071)'atomAD',ia,atomNameWat(ia),
     &        resNameWat(ia),'W',resNumbWat(ia),
     &        (atomXYZWat(k+ia3),k=1,3)
        end if !C
c
        end do !iwa
        end if ! loadWaw(iw)
c
        end do !iwM 
c
        end do!i1
        end do!i2
        end do!i3 
c ENDsolvationBOX
c
        if(CONTROL)then
        write(kanalp,*)'solvateMoL01:ALL waterInsolBOX:'
        do ia=1,nWatMolAtAdd01
        ia3=ia*3-3
        write(kanalp,7071)'ATOMW ',ia,atomNameWat(ia),
     &        resNameWat(ia),'W',resNumbWat(ia),
     &        (atomXYZWat(k+ia3),k=1,3)
        end do!ia
        end if !C
c
        watResName = 'HOH '
        watAtName(1)='O   '
        watAtName(2)='H1  '
        watAtName(3)='H2  '
c init  waterMolTag(*)
         do i=1,nWaterMolMAX 
         waterMolTag(i)=0      !all exCluded
         end do
c
c define noH MOL+WATer atoms subset
c
         natomMOLnoH=0
         natMOLandWatnoH=0
c
         write(kanalp,*)
     &  'solvateMoL01:noH:MolecAtoms:natMOLandWatnoH:'
c
         do ia=1,natomMOL
         if(atomName(ia)(1:1) .ne. 'H')then
         natomMOLnoH=natomMOLnoH+1
         natMOLandWatnoH=natMOLandWatnoH+1
         n3=natMOLandWatnoH*3-3 
         ia3=3*ia-3
c
         if(natMOLandWatnoH .ge. natMOLandWatnoHMAX)then
         write(kanalp,*)
     &   'solvateMoL01: ERROR! natMOLandWatnoHMAX is Small:',
     &   natMOLandWatnoHMAX, ' increase param:natMOLandWatnoHMAX !'
c
         write(kanalPStat,*)mError,
     &  'ERROR:solvateMoL01:natMOLandWatnoHMAX is SmaLL:',
     &   natMOLandWatnoHMAX, ' increase param:natMOLandWatnoHMAX !'
c
         stop
         end if
c
         atomNamenoH(natMOLandWatnoH)=atomName(ia)
         resNamenoH(natMOLandWatnoH)=resName(ia)
         resNumbnoH(natMOLandWatnoH)=resNumb(ia)
         atomXYZMolWatnoH(n3+1)=atomXYZ(ia3+1)
         atomXYZMolWatnoH(n3+2)=atomXYZ(ia3+2)
         atomXYZMolWatnoH(n3+3)=atomXYZ(ia3+3)
c
        if(CONTROL)then
        write(kanalp,*)'solvateMoL01: natomMOLnoH,natMOLandWatnoH:',
     &  natomMOLnoH,natMOLandWatnoH
        write(kanalp,7071)'ATOMM ',natomMOLnoH,atomName(ia),
     &        resName(ia),'M',resNumb(ia),
     &        (atomXYZ(k+ia3),k=1,3)
        end if !C
c
        end if !H
        end do !ia
c
c collect O atoms of WATers 
          write(kanalp,*)
     &  'solvateMoL01:noH:    WaterInsolBOX:natMOLandWatnoH:'
     &  ,natomMOLnoH,natMOLandWatnoH
c
         do ia=1,nWatMolAtAdd01
         if(atomNameWat(ia)(1:1) .ne. 'H')then
         natMOLandWatnoH=natMOLandWatnoH+1
         n3=natMOLandWatnoH*3-3
         ia3=3*ia-3
         atomNamenoH(natMOLandWatnoH)=atomNameWat(ia)
         resNamenoH(natMOLandWatnoH)=resNameWat(ia)
         resNumbnoH(natMOLandWatnoH)=resNumbWat(ia)
         atomXYZMolWatnoH(n3+1)=atomXYZWat(ia3+1)
         atomXYZMolWatnoH(n3+2)=atomXYZWat(ia3+2)
         atomXYZMolWatnoH(n3+3)=atomXYZWat(ia3+3)
c
        if(natMOLandWatnoH .ge. natMOLandWatnoHMAX)then
         write(kanalp,*)
     &   'solvateMoL01: ERROR! natMOLandWatnoHMAX is Small:',
     &   natMOLandWatnoHMAX, ' increase param:natMOLandWatnoHMAX !'
c
         write(kanalPStat,*)mError,
     &   'solvateMoL01: ERROR! natMOLandWatnoHMAX is Small:',
     &   natMOLandWatnoHMAX, ' increase param:natMOLandWatnoHMAX !'
c
         stop
         end if
c
        if(CONTROL)then
        write(kanalp,7071)'ATOMWw',natMOLandWatnoH,atomNameWat(ia),
     &        resNameWat(ia),'W',resNumbWat(ia),
     &        (atomXYZWat(k+ia3),k=1,3)
        end if !C
c
        end if !.ne. 'H'
        end do !ia
c
        ncallwPL=0
        nnbWpLMX = nnbWpLMAX
cx        rCUTwPL = 6.0
        rCUTwPL = dSOLVSheLL + 0.5
c
         if(CONTROL)then
        write(kanalp,*)'solvateMoL01: start nonbondListFull:'
        end if !C
c
                 call nonbondListFull(ncallWPL,
     &          natMOLandWatnoH,atomXYZMolWatnoH,rCUTwPL,
     &          nbWPairList,startnbWPairL,nnbWPairL,
     &          nbWPairListD2,nnbWpLT,nnbWpLMX)
c
c do loop over all Water Oxygens
        do ia = natomMOLnoH+1,natMOLandWatnoH
        if(nnbWPairL(ia) .ge. 1) then
c take the first closest MOLatol
        do jn = startnbWPairL(ia),startnbWPairL(ia) + nnbWPairL(ia)-1
        ja = nbWpairList(jn)
c
        if(ja .le. natomMOLnoH)then
c ja = soluteMOL atom
c define dAtMolWatMin2 for atom ja-WAter
         ihsRad=ihsRadMAX
         if(atomName(ja)(1:1) .eq. 'C')ihsRad=2
         if(atomName(ja)(1:1) .eq. 'N')ihsRad=3
         if(atomName(ja)(1:1) .eq. 'O')ihsRad=4
         dAtMolWatMin2 = (hsRad(ihsRad)+hsRad(1))**2
         if(nbWpairListD2(jn) .gt. dAtMolWatMin2 .and.
     &     nbWPairListD2(jn) .lt. dAtMolWatMax2) then
         iw = ia-natomMOLnoH
         waterMolTag(iw)=ja      ! nearest solute MOLecule atom
c
         if(CONTROL)then
         d=sqrt((atomXYZMolWatnoH(ia*3-2)-atomXYZMolWatnoH(3*ja-2))**2
     &     +(atomXYZMolWatnoH(ia*3-1)-atomXYZMolWatnoH(3*ja-1))**2
     &     +(atomXYZMolWatnoH(ia*3)-atomXYZMolWatnoH(3*ja))**2)
c
         write(kanalp,*)'ia,iw,jn,ja:nbWpairListD2(jn):',
     &   ia,iw,jn,ja,sqrt(nbWpairListD2(jn)),d
         end if !C
c
        end if ! inSolvationSheLLCondition
         goto 201
        end if !ja .le. natomMOLnoH
        end do!jn
        end if!nnbWPairL(ia)
201     continue 
        end do !ia
c waterMolTag(ia) is set         
c
        if(CONTROL)then
        write(kanalp,*)'solvateMoL01:waterMolTag()=nearestSoluteMOLat:'
        do ia = natomMOLnoH+1,natMOLandWatnoH
        write(kanalp,*)'ia: waterMolTag(ia):',ia,waterMolTag(ia)
        end do !ia
        end if!C
c
        if(CONTROL)then
        natom=natomMOL
        nres=nresMOL
        write(kanalp,*)'solvateMoL01:Init:MOL+waterInSolvShel:'
        do ir=1,nres
        write(kanalp,*)'ir,startAt,stopAt(ir):',
     &  ir,startAtInRes(ir),stopAtInRes(ir)
        end do !ir
c
        write(kanalp,*)'nChain:',nChain
        do ich=1,nChain
        write(kanalp,*)'ich,startAtInCha(ich),stopAtInCha(ich):',
     &  ich,startAtInCha(ich),stopAtInCha(ich)
        end do !ich
c
        end if !C
c
c final collect WATer in solvation SHELL
c update newMOLEcule(oldMOL+watSheLL) atomXYZ(*) and PDBinfo
c update molecular topology data structure
c
        natom=natomMOL
        nres=nresMOL
        nWatMolAdd=0
        nWatMolAtAdd=0
c
        do iwM=1,nWatMolAdd01               
        if(waterMolTag(iwM) .gt. 0)then
        n3=3*natom
        wa1=natIn1watMol*(iwM-1)
        wa3=wa1*3        
        nWatMolAdd=nWatMolAdd+1
        nres=nres+1
        if(natom .gt. natomMAX)then
        write(kanalp,*)
     &  'solvateMoL01:ERROR:natom (mainMOL+watSheLL) >natomMAX:'
c
         write(kanalPStat,*)mError,
     &  'solvateMoL01:ERROR:natom (mainMOL+watSheLL) >natomMAX:' 
c
        stop
        end if
c
        if(nres .gt. nresMAX)then
        write(kanalp,*)
     &  'solvateMoL01:ERROR:(mainMOL+watSheLL) nres >nresMAX:',nresMAX
c
        write(kanalPStat,*)mError,
     &  'solvateMoL01:ERROR:(mainMOL+watSheLL) nres >nresMAX' 
        stop
        end if
c set atomXYZ(*) for WaterInSolvationSheLL
        do iwa = 1,natIn1watMol*3
        atomXYZ(n3+iwa)=atomXYZWat(wa3+iwa)
        end do !iwa 
c set waterAtomNAmes
        startAtInRes(nres) = natom+1
c
c pair12List,startPairL12,nPairL12 for waterAtoms
        do iwa = 1,natIn1watMol
         nPairL12(natom+iwa) = 2
         startPairL12(natom+iwa) = 
     &   startPairL12(natom+iwa-1)+nPairL12(natom+iwa-1) 
c moveAtomList(*) ! all WaterAtom are in movingList
         nmoveatom = nmoveatom + 1
         moveAtomList(nmoveatom)=natom + iwa        
         moveFlag(moveAtomList(nmoveatom)) = 1
         realAtomFlag(moveAtomList(nmoveatom)) = 1
c 
c addWaterXYZ to moveAtomXYZ(*)
         im3 = nmoveatom*3-3
         ia3 =  3*moveAtomList(nmoveatom)-3
         moveAtomXYZ(im3+1) = atomXYZ(ia3+1)
         moveAtomXYZ(im3+2) = atomXYZ(ia3+2)
         moveAtomXYZ(im3+3) = atomXYZ(ia3+3)
        end do !iwa
         pair12List(startPairL12(natom+1))=natom+2       ! O-H1, O-H2
         pair12List(startPairL12(natom+1)+1)=natom+3
         pair12List(startPairL12(natom+2))=natom+1       ! H1-O, H1-H2
         pair12List(startPairL12(natom+2)+1)=natom+3
         pair12List(startPairL12(natom+3))=natom+1       ! H2-H1, H2-O
         pair12List(startPairL12(natom+3)+1)=natom+2
c END pair12List,startPairL12,nPairL12 for waterAtoms    
c 
c add Water to global xyzPDBinfo.h :
        do iwa =1,natIn1watMol
        nWatMolAtAdd=nWatMolAtAdd+1
        natom=natom+1
	atomPWatTag(natom) = 2     ! 2=Water molecule
        resName(natom)=watResName
        resNameEx(natom)=resName(natom)
        resNumb(natom)=nres
        resNameRes(nres)=watResName
        resListSeq(nres)=watResName
        chNameRes(nres)='W'            !water chain name
        chName(natom)='W'
        atomName(natom)=atomNameWat(wa1+iwa)        
        atomNameEx(natom)=atomName(natom)
        atHvyNbInRes(nres)=1                ! one heawy atom in WaterMolecule
        atomNumb(natom) = natom
c TIP3P MODEL
        if(iwa .eq. 1)then
        ffAtomName(natom)='OW'
        atomQ(natom) = -0.8340         ! OW
        else 
        ffAtomName(natom)='HW'
        atomQ(natom) = 0.4170 
        end if
        atomQ0(natom) = atomQ(natom)
        atomBlockName(natom)='WS'
        end do!iwa
c
        stopAtInRes(nres) = natom
        startAtInRes(nres+1)=natom+1   
c
        if(CONTROL)then
        ir=nres
        write(kanalp,*)'iwM,ir,startAt,stopAt(ir):',
     &  iwM,ir,startAtInRes(ir),stopAtInRes(ir)
        end if !C 
c
        end if !waterMolTag(iwM) .gt. 0
        end do !iwM
c solvationSheLL WATer = newChain
c update nChain and startAtInCha(ich),stopAtInCha(ich)
        nChain = nChain + 1
        startAtInCha(nChain) = stopAtInCha(nChain-1) + 1
        stopAtInCha(nChain) = natom
        stopAtInMol(nMolec) = natom
c updated atomXYZ(*) consist of natomMOL - soluteAtoms
c                               + nWatMolAdd*natIn1watMol - waterMolAtoms
c
c define :nRestr1Seg and resEndRestr1(i2+1),resEndRestr1(i2+2):
c restr1ResHWinfo for solvationSheLLWaters
        nRestr1RHWSeg=1
        resEndRestr1RHW(2*nRestr1RHWSeg-1)=nresMOL+1
        resEndRestr1RHW(2*nRestr1RHWSeg) =
     &  resEndRestr1RHW(2*nRestr1RHWSeg-1)+nWatMolAdd - 1
c
        write(kanalp,*)'solvateMoL01: nRestr1RHWSeg:',nRestr1RHWSeg
        write(kanalp,*)' resEndRestr1RHW(*):',
     &  (resEndRestr1RHW(i),i=1,2*nRestr1RHWSeg)
c
        kanalXYZwSheLL = kanalwWatBrg 
        fileWatSheLL = "mol.watSheLLXYZ0.pdb"
        if(nMolNameLett .ge. 1)
     &  fileWatSheLL = 
     &  molNameArgLine(1:nMolNameLett)//".watSheLLXYZ0.pdb"
        open(unit=kanalXYZwSheLL,file=fileWatSheLL,form='formatted',
     &       status='unknown')
c
        write(kanalp,*)'nChain:',nChain
        do ich=1,nChain
        write(kanalp,*)'ich,startAtInCha(ich),stopAtInCha(ich):',
     &  ich,startAtInCha(ich),stopAtInCha(ich)
        end do !ich
c
        do ir=1,nres
        write(kanalp,*)'ir,startAt,stopAt(ir):',
     &  ir,startAtInRes(ir),stopAtInRes(ir)  
        end do !ir
c
        write(kanalXYZwSheLL,*)'PDB: initial molecule + waterSheLL:',
     &  ' dSHELL(A):',dSHELL
        write(kanalXYZwSheLL,*)'# nWatMolinSheLL:', nWatMolAdd
c
        do ich=1,nChain
        do ia = startAtInCha(ich),stopAtInCha(ich)
        ia3=ia*3-3
        write(kanalXYZwSheLL,7071)'ATOM  ',ia,atomName(ia),
     &        resName(ia),chName(ia),resNumb(ia),
     &        (atomXYZ(k+ia3),k=1,3), ffAtomNAme(ia),atomQ(ia)
         end do !ia
        write(kanalXYZwSheLL,'(a3)')'TER'
        end do !ich
        write(kanalXYZwSheLL,'(a3)')'END'    
c
        write(kanalp,*)'solvateMoL01: finish!:'
c         end if !C
        return 
c
7071    format(a6,i5,2x,a4,a4,a1,i4,4x,3f8.3,2x,a2,1x,f8.5) !orig PDB
c
        end
