c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c mdynSB                                                           *
c  Yury Vorobjev 2002                                              * 
c                2003                                              * 
c                2004                                              * 
c                2005                                              *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
 	subroutine inputMDSApar
c
        include 'xyzPDBsize.h'       
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h' 
        include 'xyzPDBcrd.h'
        include 'loopInfo.h'
        include 'mdRunPar.h'
        include 'simAnneal.h'
        include 'optionPar.h'
        include 'enForce.h'
        include 'nbondPairVCS.h'
        include 'filedat.h'
        include 'solvGSarray.h'
        include 'coulEnPar.h'
        include 'optimiz.h'
        include 'shake.h'
        include 'restrainInfo.h'
        include 'vdw12Par.h'
        include 'restrain1MHW.h'
        include 'restrain1RHW.h'
        include 'ligInfo.h'
        include 'restrDistA2.h'
        include 'compGeoEF.h'
        include 'hbond128.h'
        include 'replicaEx.h' 
        include 'rigBody01.h'
        include 'solvate01.h'
        include 'solvWBrg01.h'
        include 'ionShell.h'
        include 'solvGBorn.h'
        include 'dataSASdr.h'
        include 'ssbond.h'
        include 'engXYZtra.h'
        include 'modeXYZtra.h'
c
        integer kanalp
        logical CONTROL
        character*20 keyw
        integer keywl,keywt
        character*(charLenMAX) filenam
        character*(charLenMAX) valuec
        integer envVarCharMX_OPT
        character*(charLenMAX) HOMEBS_dirLoc
        character*(charLenMAX) resultFileDirLoc 
        integer valuei
        integer i,n,k
        integer nSAP,nMDREx
        real valuer
        real aSC,sc1,sc2
        character*80 line
        logical valuel,found
c
cx        kanalp = 6
        kanalp = kanalRunOut
        CONTROL = .false.
c
        write(kanalp,*)'inputMDSApar:  '
        write(kanalp,*)
c
c DIRECTORY and FILE NAMES
c
        call getenv("BIOPASEDHOME", HOMEBS_dir)
c
c extract VALID HOMEBS_dir string from "MDYN05HOME"
        envVarCharMX_OPT = charLenMAX  ! MAX character in envVariale
        ndirLett = 0
c extract rigth string for HOMEBS_dir
        call clean_space_right(HOMEBS_dir,envVarCharMX_OPT,
     &        HOMEBS_dirLoc,ndirLett)
c
        HOMEBS_dir(1:ndirLett) = HOMEBS_dirLoc(1:ndirLett) 
c
        write(kanalp,*)'HOMEBS_dir = ',HOMEBS_dir(1:ndirLett),
     &  ' ndirlett=',ndirlett
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c INIt forces scaling parameters
       fEngWF(1) = 1.0   ! vbond force weighting factor
       fEngWF(2) = 1.0   ! vangl
       fEngWF(3) = 1.0   ! impAng
       fEngWF(4) = 1.0   ! torsAng
       fEngWF(5) = 1.0   ! vdwR1
       fEngWF(6) = 1.0   ! coulR1
       fEngWF(7) = 1.0   ! coulR2
       fEngWF(8) = 1.0   ! retRFoces: restr1Eng+restr1MHWEng+restrDistA2Eng
       fEngWF(9) = 1.0   ! sovationGSmodel+compactGeoEn
       fEngWF(10)= 1.0   ! HBond128
       fEngWF(11)= 1.0
c 
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  DEFINE main  command file = MdynPar.inp :
c
cx        filenam = './MdynPar.inp'        ! in current job_dir
        filenam = mdynParFile            ! get via argLine
c
        write(kanalp,*)'MdynPar inputFile:',filenam
c
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c define absolute PATH to write all resultFiles *.pdb *.tra
       resultFileDir = ""
       nResFileDirLett = 0
       OPT_resFileDir = .false.
c
c define PATH from command file mdynParFile:
       keyw = 'resFileDir'
       keywl = 10
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found)then 
       OPT_resFileDir = .true.
       resultFileDir = valuec
c extract rigth string 
        call clean_space_right(resultFileDir,envVarCharMX_OPT,
     &        resultFileDirLoc,nResFileDirLett)
c
        resultFileDir(1:nResFileDirLett) = 
     &                  resultFileDirLoc(1:nResFileDirLett)
c
c add right slash for DIR name
        if(resultFileDir(nResFileDirLett:nResFileDirLett) .ne. 
     &    "/" ) then
        nResFileDirLett=nResFileDirLett+1
        resultFileDir(nResFileDirLett:nResFileDirLett) = "/"
       end if
c
       end if
        write(kanalp,*)'resultFileDir = ',
     &  resultFileDir(1:nResFileDirLett),
     &  ' nResFileDirLett=',nResFileDirLett
c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c DEFINE INPut data:
c atomXYZ().pdb file
cx       pdbfile = './molec.pdb' ! default=molec.pdb   !defualt in current job_dir
       pdbfile = xyzInpFile      ! get via argLine
c define Pathname from command file mdynParFile:
       keyw = 'molecFile'
       keywl = 9
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) pdbfile = valuec
c
       write(kanalp,*)'atomXYZ are taken : ',pdbfile
c
c resSeqFile        
       resSeqFile = './resSeq.inp'   ! default= ./resSeq.inp
c define NEW Pathname:
       keyw = 'resSeqFile'
       keywl = 10
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) resSeqFile = valuec
c
c       write(kanalp,*)'resSeqFile are taken : ', resSeqFile
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
c ffAtomTypeFile = HOMEBS_dir(1:nld)//'/dat/atmAAambff.dat' !defined in inputMDSApar.f
       ffAtomTypeFile = HOMEBS_dir(1:ndirLett)//'/dat/atmAAambff.dat'
c
       keyw = 'ffAtype'
       keywl = 7
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) ffAtomTypeFile = valuec            
c fullPathFilename
       write(kanalp,*)'ffAtomTypeFile :',ffAtomTypeFile
c       
c ffparametersFile:
       ffParFile = HOMEBS_dir(1:ndirLett)//'/dat/bsparBATV.dat'
c
       keyw = 'ffParFile'
       keywl = 9
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) ffParFile = valuec                 
c fullPathFilename
       write(kanalp,*)'ffParFile :', ffParFile
c
c solvAtmTypeFile: 
       solvAtmTypeFile = 
     &        HOMEBS_dir(1:ndirLett)//'/dat/solvGSPar_all_amb.dat'
c
       keyw = 'solvAtmTypeFile'
       keywl = 15
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) solvAtmTypeFile = valuec     
c
       write(kanalp,*)'solvAtmTypeFile : ',solvAtmTypeFile
c            
c solvAtmDatFile:
       solvAtmDatFile = HOMEBS_dir(1:ndirLett)//'/dat/solvGSPar.dat'
c
       keyw = 'solvAtmDatFile'
       keywl = 14
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) solvAtmDatFile = valuec                 
c
        write(kanalp,*)'solvAtmDatFile :',solvAtmDatFile
c
c pdbAmbZmFile: InPdb names --> amberZm names
       pdbAmbZmFile = 
     &      HOMEBS_dir(1:ndirLett)//'/dat/pdbAtName_ambZm.dat'
c
       keyw = 'pbdInAmbZmFile'
       keywl = 14
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) pdbAmbZmFile = valuec                 
c
        write(kanalp,*)'pdbAmbZmFile :',pdbAmbZmFile
c
c hInfoFile: info to add Hydrogens
       hInfoFile = 
     &      HOMEBS_dir(1:ndirLett)//'/dat/h_add.dat'
c
       keyw = 'hInfoFile'
       keywl = 9 
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) hInfoFile = valuec                 
c
        write(kanalp,*)'hInfoFile:', hInfoFile
c
c END fileName definition
c
c READ inPutInfo: OPT_ :
c
c * * * * * * * * * * * * * * *
c
        OPT_iHread=0
        OPT_Hread = .false.  ! atom H are missing in PDB ( LMOD01program)
        keyw = 'Hread'
        keywl = 5
        keywt = -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_Hread = .true.
        if(OPT_Hread)OPT_iHread=1
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        OPT_SSbondAuto=.true.     !atomatic corretion of CYS -->CYX in inPDBfile
	keyw = 'defSSbond'
	keywl = 9
	keywt = 2
	keywt = -1
	call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_SSbondAuto=.false.
c OPT_SSbondAuto=.false. ! inPDB should have right CYS/CYS/CYM residues names
c SSbonds will be constructed for the CYX residues
c ---------------------------------------------------------------------------
cXX        OPT_SSBonds=.true. 
cSS bomnds should be assigned by hand editing of PDBfile
        OPT_SSBonds=.false.
        if(OPT_SSBonds)then
        ssCYSname='CYX '     !CYS making SS
        ssAtname='SG  '
        dss2MAX_OPT=4.2**2   !default SS maxdistance to be included in SS list
        keyw = 'SSBonds'
        keywl = 7
        keywt = 2   ! -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
cx      if(found)OPT_SSBonds=.true. 
        if(found)then
        if(valuec(1:4) .ne. ssCYSname)ssCYSname=valuec(1:4)
        end if
        end if !OPT_SSBonds
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        fileWATBOX0=HOMEBS_dir(1:ndirLett)//'/dat/spc216.amb.dat'
        write(kanalp,*)'fileWATBOX0:',fileWATBOX0
c
	OPT_SolvateExWat=.false.
        iSolvateExWat=0
        iDefineMOLec=0
        dSOLVSheLL = 6.5 !4.5
        keyw = 'SolvateExWat'
        keywl = 12
        keywt = 1  
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_SolvateExWat = .true.
        iSolvateExWat=1
        dSOLVSheLL = valuer
        end if
c ------------------------------------------------------------------
        OPT_coIonShell=.false.
        keyw = 'coIonShell'
        keywl = 10
        keywt = -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_coIonShell = .true.
        end if
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        OPT_PDBresSeq = .true.   ! res Sequence are taken from the Inputmolec.pdb file
        OPT_resSeqInp = .false.  ! res sequence are taken from the Input file resSeq.inp
        keyw = 'resSeqInp'
        keywl = 9
        keywt = -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_resSeqInp = .true.
        OPT_PDBresSeq = .not.OPT_resSeqInp
c
c * * * * * * * * * * * * * * *
        OPT_AddLoop = .false.  ! all AtomXYZ exists in PDB 
        keyw = 'AddLoop'
        keywl = 7
        keywt = -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_AddLoop = .true.
c
c * * * * * * * * * * * * * * *
        OPT_fullProtMD = .true.
        keyw = 'fullProtMD'
        keywl = 10
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_fullProtMD=.true. 
c
         OPT_LoopMD = .not.OPT_fullProtMD
c * * * * * * * * * * * * * * *
        OPT_LoopMD = .false. ! F = fullProtMD
        keyw = 'MovingRes'
        keywl = 9
        keywt = -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_LoopMD = .true.
c
c needs file moveRes.inp to define residue List
c
c duplicate keyWords
        keyw = 'MoveRes'
        keywl = 7
        keywt = -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_LoopMD = .true.
c
        OPT_fullProtMD = (.not. OPT_LoopMD)  ! fullProtein MD
        OPT_fullProtMD = (.not. OPT_LoopMD)  ! fullProtein MD
        OPT_residMD = OPT_LoopMD
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        OPT_LigRes = .false. ! 
        ligResFile = './ligRes.inp'  !defualt file define LIGRES

        keyw = 'LigRes'   
        keywl = 6
        keywt = -1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_LigRes = .true.
c
c needs file ligRes.inp to define residue List
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        OPT_WRmolAtIntLig = .false.
        rAtInLig = 4.5   ! all molAtoms 4.5 from ligAtoms
        molAtXYZIntLigFile = './molAtXYZintLig.pdb'
c write file ./molAtXYZintLig.pdb - atoms in 5A vicinity from Lig atoms
c
         keyw = 'molAtIntLig'
        keywl = 11
        keywt = 1
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_WRmolAtIntLig = .true.
        rAtInLig = valuer
        end if
c
        write(kanalp,*)
     &  'inputMDSApar: OPT_WRmolAtIntLig:',OPT_WRmolAtIntLig
c
        if(rAtInLig .lt. 3.0) rAtInLig = 3.0
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *  
        OPT_doLigDock = 0      
        keyw = 'doLigDock'
        keywl = 9
        keywt = 0 
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_doLigDock = valuei
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        OPT_MDSA = .false.   !MD Simulated Annealing
        keyw = 'MDSA'
        keywl = 4
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_MDSA = .true.
c
c * * * * * * * * * * * * * * *
        OPT_SolvGS = .false.
        keyw = 'SolvGS'
        keywl = 6
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_SolvGS = .true.  
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c SAS calculation
        solMolRad_OPT = 2.8  ! default WatMol diameter
        dotden = 4.0         ! SAS for Wbrg and bornRad calculation
c
        OPT_SolvEFS = .false.
	OPT_SolvWbr = .false.
        fileEFSwbrXYZ = './watBrgSolvXYZ.pdb'
       if(nMolNameLett .ge. 1)then
       fileEFSwbrXYZ = molNameArgLine(1:nMolNameLett)//
     &                 '.watBrgSolvXYZ.pdb'
       end if
c
       if(OPT_resFileDir)then
       fileEFSwbrXYZ = resultFileDir(1:nResFileDirLett)//
     &                 'watBrgSolvXYZ.pdb'
       if(nMolNameLett .ge. 1)then
       fileEFSwbrXYZ = resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//
     &                 '.watBrgSolvXYZ.pdb' 
       end if
       end if
c
        keyw = 'SolvEFS'
        keywl = 7
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
	OPT_SolvEFS = .true.
	OPT_SolvWbr = .true.
	end if
c
        keyw = 'SolvWbrg'
        keywl = 8
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
	OPT_SolvEFS = .true.
        OPT_SolvWbr = OPT_SolvEFS
	end if
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        OPT_SolvSASHP = .false.
        nSAScall = 0
        keyw = 'SolvSAShp'
        keywl = 9
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_SolvSAShp = .true.
c*  * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        OPT_writeSASdotXYZ = .false.
c
        fileSASdotXYZ = 'molSASXYZ.pdb'
         if(nMolNameLett .ge. 1)then
        fileSASdotXYZ = molNameArgLine(1:nMolNameLett)//
     &                 '.molSASXYZ.pdb' 
        end if!
c
        if(OPT_resFileDir)then
        fileSASdotXYZ = resultFileDir(1:nResFileDirLett)//
     &                 'molSASXYZ.pdb'
        if(nMolNameLett .ge. 1)then
        fileSASdotXYZ = resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//
     &                 '.molSASXYZ.pdb'
        end if
        end if
c
        keyw = 'writeSASxyz' 
        keywl = 11
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_writeSASdotXYZ=.true.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        OPT_SolvGBorn = .false.
        keyw = 'SolvGBorn'
        keywl = 9
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_SolvGBorn = .true.
c *  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        OPT_zeroRot = .false.
        keyw = 'zeroRot'
        keywl = 7
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)OPT_zeroRot = .true.
c correct for selConsistency
        if(.not. OPT_fullProtMD) OPT_zeroRot = .false.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c LoopStartEnd FileName:
       if(OPT_LoopMD ) then
cx       loopfile = HOMEBS_dir(1:ndirLett)//'/job/loopMD.inp'
cx         loopFile = './moveRes.inp'
cx         loopFile = loopFile    ! in getArgLine
c
       keyw = 'moveResFile'
       keywl = 11
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) loopfile = valuec
       end if ! OPT_LoopMD
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c Hbond128
        nHBbondHxYMX = nHBbondHxYMAX
        nHBbondHxY = 0 
        OPT_HBond128 = .true. 
        fEngWF(10) = 1.0
        hB128WfScaleALL=1.5
        rcutHb128 = 3.5
c
        keyw = 'hBond128'
        keywl = 8 
        keywt = 1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_HBond128 = .true.
        fEngWF(10) = 1.0   
        hB128WfScaleALL=valuer
        rcutHb128 = 3.5               ! Hx ... Y
        end if !found
        hB128AgLpar(1) = -1.0         ! cosTet0
        hB128AgLpar(2) = 3.14159/6.0   ! dTeta= 30.0    
        hB128AgLpar(2) = cos(3.14159-hB128AgLpar(2))-hB128AgLpar(1)
        hB128AgLpar(2) = hB128AgLpar(2)**2  ! sigma^2
        hB128AgLpar(2) = 0.018
c Init hB128WfList(iHb)
        fEngWF(11) = fEngWF(10)  ! BB=1.0  BS=1.0  !default Optimal
c fEngWF(10)=SS, fEngWF(11)=(BS)hB128WfList()
        do i = 1,nHBbondHxYMAX
        hB128WfList(i) = fEngWF(10)*hB128WfScaleALL
        hB128TypeList(i) = 1
        end do  !i
c for eOp and initial MD run all hB128Angl, t.e. Back-Back,Side-Back,
c all Hb128Wfactors are equal to initial fEngWF(10)
c $hBond128=fEngWF(10)
c hB128WfList(i) = BB/BS in the initMDStart   
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c RESTRAINTS and COMpactization  FORCES
c
        OPT_harmAt1PosRst = .false.
        restr1AtConst = 0.0
        keyw = 'harmAt1PosRst'
        keywl = 13
        keywt = 1 
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_harmAt1PosRst = .true.  
        restr1AtConst = valuer
        end if !found 
c * * * * * * * * * * * * * * *
c REsrained Residues Segments file: defines: resEndRestr1(*)
cx       restr1File = './restrAt1.inp'
         restr1File = restrA1type
c
       if(OPT_harmAt1PosRst) then
c
       keyw = 'restr1File'
       keywl = 10
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) restr1File = valuec
       end if ! OPT_harmAt1PosRst
c * * * * * * * * ** * * * * * * * * * * * * * * * * * * * *
c harmonicWell Mol1 restraints
        OPT_HW1MolPosRst = .false.
        restr1MHWConst = 0.0
        sizeRestr1MHW =  4.0  ! 3.0 for smallLigand ! A  default
        keyw = 'harmMol1PosRst'
        keywl = 14
        keywt = 1 
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_HW1MolPosRst = .true.  
        restr1MHWConst = valuer
        end if !found 
c * * * * * * * * * * * * * * *
c REstrained Residues Segments file: defines: resEndRestr1(*)
c harmonicWell Mol1 
c
       restr1MHWFile = './restrMol1.inp'
       if(OPT_HW1MolPosRst) then
c
       keyw = 'harmRes1RsFile'
       keywl = 14
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) restr1MHWFile=valuec
       end if ! OPT_HW1MolPosRst 
c * * * * * * * * * * * * * * * * * * * *
c harmonicWell RES1 CMass restraints
        OPT_HW1ResPosRst = .false.
        restr1RHWConst = 0.00 ! 0.5  ! 
        sizeRestr1RHW =  2.8  ! RestrBox=2*sizeRestr1RHW for WATer
        keyw = 'harmRes1PosRst'
        keywl = 14
        keywt = 1 
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_HW1ResPosRst = .true.  
        restr1RHWConst = valuer
cx YV:2010 new SASdependent restrWatSheLLPot
cx        elseif (OPT_SolvateExWat) then
cx        OPT_HW1ResPosRst = .true.
cx        restr1RHWConst = 0.50      !default for water in solvSheLL
cx	sizeRestr1RHW = 1.70
        end if !found 
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c REstrained Residues Segments file: defines: resEndRestr1(*)
c harmonicWell Res1
       restr1RHWFile = './restrRes1.inp'
       if(OPT_HW1ResPosRst) then
c
       keyw = 'harmRes1RsFile'
       keywl = 14
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) restr1RHWFile=valuec
       end if ! OPT_HW1ResPosRst
c * * * * * * * * ** * * * * * * * * * * * * * * * * * * * * * * *
c a1-a2 DistantRestrA2
        OPT_restrDistA2 = .false.
        keyw = 'distRestrA2'
        keywl = 11
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_restrDistA2 = .true.
        end if !found
        write(kanalp,*)'inputMDSApar:OPT_restrDistA2:',OPT_restrDistA2
c
c distantRestrA2 file: defines:  a1  a2  dist12  hKonst
c
cx       restrA2File = "./distRestrA2.inp"
       restrA2File = restrA2type
c
       write(kanalp,*)'inputMDSApar:restrA2File:',restrA2File
       if(OPT_restrDistA2) then
c
       keyw = 'distRestrA2File'
       keywl = 15
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) restrA2File = valuec
       end if ! OPT_restrDistA2  
c
       write(kanalp,*)'inputMDSApar: restrA2File:',restrA2File
c
c COMPactizationForce 
c
c default parameters:
        compactGeoFpar(1) = 1.0 
        compactGeoFpar(2) = 15.0
        compactGeoFpar(3) = 1.0 
        compactGeoFpar(4) = 2.0
c
        OPT_CompactForce = .false.
        compactGeoEn=0.0 
        keyw = 'compactForce'
        keywl = 12
        keywt = 1 
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_CompactForce = .true.
        compactGeoFpar(1) = valuer
        end if !found
        write(kanalp,*)'inputMDSApar:OPT_CompactForce:',OPT_CompactForce
        if(OPT_CompactForce)then
        write(kanalp,*)'compactGeoF Par: aa,del1,bb,del2:',
     &   compactGeoFpar
        end if !
c
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c MD OPTIONs
c SHAKE OPT
        OPT_Shake = 0
        iShake = 0    !=0 noShake, =1 shakeH, =2 shakeAllBonds
       keyw = 'shake'    
       keywl = 5 
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) OPT_Shake = valuei
       iShake = OPT_Shake
c
       shTol = 0.0001      !shake accuracy
       shitMX = 100      !MAX numb shake iteration
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c OPT_RigidBodyMd
        OPT_RigidBodyMd = .false.
        rigBodyFile = './rigBody.inp'
        keyw = 'rigidBodyMd'
        keywl = 11
        keywt = -1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_RigidBodyMd = .true.
        end if !found
        write(kanalp,*)'inputMDSApar:OPT_RigidBodyMd:',OPT_RigidBodyMd
c
       if(OPT_RigidBodyMd) then
       keyw = 'rigBodyFile'
       keywl = 11
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) rigBodyFile = valuec
       end if ! OPT_RigidBodyMd
c
       write(kanalp,*)'inputMDSApar: restrA2File:',restrA2File
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c OPT_RigidBodyCADist
        OPT_RigidBodyCADist = .false.
        cAcARigBodyHK = 0.0
        keyw = 'rigBodyCAdist'
        keywl = 13
        keywt = 1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        OPT_RigidBodyCADist = .true.
        cAcARigBodyHK = valuer 
        end if !found
        write(kanalp,*)'inputMDSApar:OPT_RigidBodyCADist:',
     &  OPT_RigidBodyCADist,' cAcARigBodyHK:',cAcARigBodyHK
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c electrostaticEnOPTion:   
       coulVar_OPT = 2  ! neutralized background + eps(=r*cReps_OPT)
       cReps_OPT = 3.33 
c coulVar_OPT = 0,1 : standart Coulon, neutralized background (eps=1)
c
       keyw = 'coulEnVar'
       keywl = 9 
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) coulVar_OPT = valuei
c
       if(coulVar_OPT .eq. 2)then
       keyw = 'epsConstR'
       keywl = 9
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) cReps_OPT = valuer
       end if!coulVar_OPT .eq. 2
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c energyCalculation from inputPDB
       OPT_engCalc = .false.
c
       keyw = 'engCalc'
       keywl = 7
       keywt = -1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) OPT_engCalc = .true.
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c energyOPtimization PARAMETERS
       engEps = 0.001     !in kcal/mol energy Tolerance
       OPT_engOptim = .false.
c
       keyw = 'engOptim'
       keywl = 8
       keywt = -1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) OPT_engOptim = .true. 
c       
       nOptIter = 100    !max iteration 
       keyw = 'nOptStep'
       keywl = 8
       keywt = 0 
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) nOptIter = valuei
c       
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * *       
c MD RUN PARAMETERS:
c
c CUTOFF RADII, A
       rcutV = 8.0
       rbuffV = 1.00
       rbuffS = 1.00
       rbuffC = 1.00
c
       keyw = 'rcutV'
       keywl = 6
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) rcutV = valuer
       rcutS = rcutV    ! cutOff for solvation model
c
       rcutC = 14.0
       keyw = 'rcutC'
       keywl = 6
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) rcutC = valuer
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c SoftCore parameter rSC = rmvdd12*aSoftCore
       aSoftCoreMIN = 0.650  ! smallest SoftCore
       aSoftCore = aSoftCoreMIN
       keyw = 'aSoftCore'
       keywl = 9
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) then
        aSC = valuer
        sc1 = aSoftCoreMIN
        sc2 = 0.975
        aSoftCore=sc1*aSC + sc2*(1.0 - aSC)
        end if !f
c
c***********************************************************
c dielectric constant  parameters: default:
       epsMol_OPT = 6.00     ! eps inside protein
       epsSol_OPT = 80.0     ! eps in bulck solvent
       solMolRad_OPT = 2.8     ! water solvent molec radius
cx       AdcPar_OPT = 4.50  ! =1/2distFullScreen
c                          LazaridisModel DDDconst Aparam=8.0
cx       ndcPar_OPT = 6 ! 8  ! 4          ! LazaridisModel DDDconst n 
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c scaling (neutralization) of the Phosphate groups of NA
      	qPhosphSc_OPT = 1.00
        iqPhosphSc_OPT = 0
       keyw = 'qPhosphSc'
       keywl = 9 
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) then
        qPhosphSc_OPT = valuer
        iqPhosphSc_OPT = 1
        end if

c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  additional hydrophobic solvation
cx       atomHydPhSolv = 0.45  ! optimized additional to Lazaridis
         atomHydPhSolv = 0.0
c deltaGhydr = atomHydPhSolv*Vat**(2/3) = additional HPhobic(positive)
c
       ihydrPhSolv =  1  ! 0 !flag to do =1 modification of the atomHydPhSolv
       keyw = 'hydrPHsolv'
       keywl = 10
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) then
        atomHydPhSolv = valuer
        if(atomHydPhSolv .eq. 0.0)ihydrPhSolv = 0
        end if
c
        write(kanalp,*)'inputMDSApar: Mod:ihydrPhSolv:',ihydrPhSolv
        write(kanalp,*)'inputMDSApar: atomHydPhSolv:',atomHydPhSolv
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c       
c initial MDTemperature
       tempT0(1) = 50.0 
       tempT0(2) = tempT0(1)
       tempT0(3) = tempT0(1)   !! 27.02.2015
c
       keyw = 'initMDTemp'
       keywl = 10
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) then
        tempT0(1) = valuer
        tempT0(2) = tempT0(1)
        tempT0(3) = tempT0(1)
        endif  !! 27.02.2015
c
c target MDtemperature
        tempTg = 50.0 
        keyw = 'bathMDTemp'
        keywl = 10
        keywt = 1
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) tempTg = valuer
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ** **
c MD run 
c
c md integration time step
       deltat = 0.001  
       keyw = 'mdTimeStep'
       keywl = 10
       keywt = 1
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) deltat = valuer
c * * * * * * * * * * * * * * * * * * * * * * * * * * ** * * * * 
c mdStep
c       ntimeMX = 10000    ! 10 ps
        ntimeMX = 0
        OPT_doMD = .false. 
c       ntime0 = 0
        keyw = 'runMDnstep'
        keywl = 10
        keywt = 0
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) ntimeMX = valuei
c
c OPT_doMD * * * * * * * * ** * * * * * * * * * * * * * * * * * * **
       OPT_doMD = .false. 
       keyw = 'doMDyn'
       keywl = 6 
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) OPT_doMD = .true.
       if(OPT_doMD .and. ntimeMX .eq. 0)then
        ntimeMX = 100
       write(kanalp,*)'WARNING! ntimeMX is not defined! is set =100 '
       end if !
c
c OPT_mdRestart : read in XYZ and V from file
c
       OPT_mdRestart = .false.
       mdRestXYZVfile = './mdXYZVin.pdb'
       keyw = 'mdRestart'
       keywl = 9
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found)then
        OPT_mdRestart = .true.
       if(valuec(1:2) .ne. '  ')mdRestXYZVfile = valuec
       end if !
c
c * * * * * * * * * * ** * * * * *** * * * *** ** ****** * ** * * * * *
c twin range parList
c 0-R1 VDW(CoulR1)pairList + nearest(all atoms) GShel solvationModel
       ntimeR1 = 20        !update each ntimeR1 (smallest) timeStep
       keyw  = 'updateR1PL'
       keywl = 10
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) ntimeR1 = valuei 
c
       ntimeR2 = 2*ntimeR1    !update each ntimeR2 timeStep
       keyw  = 'updateR2PL'
       keywl = 10
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) ntimeR2 = valuei 
c
c multiple step MTS
       ntimeF1 = 1          ! R1 en/force update freq = satndart
       ntimeF2 = 2*ntimeF1    ! R2 en/force
       keyw = 'multStepF2'
       keywl = 10
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) ntimeF2 = valuei
c
       ntimeF3 = 4*ntimeF1    ! SLOW force (solvation)
       keyw = 'multStepF3'
       keywl = 10
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) ntimeF3 = valuei
c
c ensemble Type
       atype = 1         ! 0- free MD, 1- NTV
       keyw = 'NTV'
       keywl = 3 
       keywt = 0 
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) atype = valuei 
c
       tauTRF = 0.1      ! 0.4-0.5 ps ! default gromos standart
c
c traj files:   
       nwtra = 100          ! default write each 100 snap write ti trajFile          
       keyw = 'nwtra'
       keywl = 5 
       keywt = 0
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) nwtra = valuei
c
c OPT_WRxyzvtra : trajFile: WRite full_traj file
       OPT_WRxyzvtra = .false.
       xyzTraFile = './xyzMd.tra'             !trajFile:
c
       if(nMolNameLett .ge. 1)then
       xyzTraFile = molNameArgLine(1:nMolNameLett)//'.xyzMd.tra'
       end if
c
        if(OPT_resFileDir)then
        xyzTraFile = resultFileDir(1:nResFileDirLett)//
     &                 'xyzMd.tra'
        if(nMolNameLett .ge. 1)then
        xyzTraFile = resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//
     &                 '.xyzMd.tra'
        end if
        end if
c User defined fileName
       keyw = 'mdXYZtra'
       keywl = 8
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found)then
        OPT_WRxyzvtra = .true.
        if(valuec(1:2) .ne. '  ')xyzTraFile=valuec
       end if !
c
c energyTra File default:
       engTraFile = './engMd.tra'                  ! in current job_dir
       if(nMolNameLett .ge. 1)then
       engTraFile = molNameArgLine(1:nMolNameLett)//'.engMd.tra'
       end if 
c
        if(OPT_resFileDir)then
        engTraFile = resultFileDir(1:nResFileDirLett)//
     &                 'engMd.tra'
        if(nMolNameLett .ge. 1)then
        engTraFile = resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//
     &                 '.engMd.tra'
        end if
        end if
c
c User defined:
       keyw = 'engMdTra'
       keywl = 8
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) engTraFile = valuec
c
c mdResult XYZ snapshot files:
       pdbResFile = './molMdRes.pdb'               ! in current job_dir
       if(nMolNameLett .ge. 1)then
       pdbResFile = molNameArgLine(1:nMolNameLett)//'.molMdRes.pdb'
       end if
c
        if(OPT_resFileDir)then
        pdbResFile=resultFileDir(1:nResFileDirLett)//'molMdRes.pdb'
        if(nMolNameLett .ge. 1)then
        pdbResFile=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.molMdRes.pdb'
        end if
        end if 
c
       write(kanalp,*)'inputMDSApar.f:pdbResFile : ',pdbResFile
c
       keyw = 'pdbMdTra'
       keywl = 8
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) pdbResFile = valuec
c
       pdbFinFile = './molMdFin.pdb'               ! in current job_dir
       if(nMolNameLett .ge. 1)then
       pdbFinFile = molNameArgLine(1:nMolNameLett)//'.mdXYZVfin.pdb'
       end if
c
       if(OPT_resFileDir)then
        pdbFinFile=resultFileDir(1:nResFileDirLett)//'molMdFin.pdb'
        if(nMolNameLett .ge. 1)then
        pdbFinFile=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.molMdFin.pdb'
        end if
        end if 
c
       keyw = 'pdbMdFin'
       keywl = 8
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) pdbFinFile = valuec
c
       pdbOptFile = './molEnOpt.pdb'               ! in current job_dir
       if(nMolNameLett .ge. 1)then
       pdbOptFile = molNameArgLine(1:nMolNameLett)//'.molEnOpt.pdb' 
       end if
c
        if(OPT_resFileDir)then
        pdbOptFile=resultFileDir(1:nResFileDirLett)//'molEnOpt.pdb'
        if(nMolNameLett .ge. 1)then
        pdbOptFile=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.molEnOpt.pdb'
        end if
        end if
c
       keyw = 'pdbOptim'
       keywl = 8
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
       if(found) pdbOptFile = valuec
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c simulated Annealing
c
        if(OPT_MDSA) then
c define SA protocol
c       
cx       fileSAprotocol = './SAprotocol.inp' !default name 
       fileSAprotocol = saProtFile
c
       keyw = 'fileSAProt'
       keywl = 10
       keywt = 2
       call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found) fileSAprotocol = valuec 
c
        if(CONTROL)then
        write(kanalp,*)'SimAnnProtocol:file:',fileSAprotocol
        end if
c
c read fileSAProtocol
        open(unit=kanalInSAprot, file=fileSAprotocol, form='formatted',
     &       status= 'old') 
c
        nSAparMX = nSAparMAX
        nSAstep = 0
200     read(kanalInSAprot,'( a80 )', end=201 ) line
c
        write(kanalp,*)line
c
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 201
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 200
        end if 
c
c read number of SA steps
        read(line,'(i6)')nSAstep
c
201     continue
        if(nSAstep .gt. nSAstepMAX .or. nSAstep .le. 0 )then
        write(kanalp,*)'ERROR!!: nSAstep .gt. nSAstepMAX!!',nSAstep
        stop
        end if!nSAstep
c
c read protocol ntimeMX  tempTg(SA)
c
        nSAP = 0
300     read(kanalInSAprot,'( a80 )', end=301 ) line
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 301
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 300
        end if
c
        nSAP = nSAP+1
c                
        if(nSAP .gt. nSAstep) goto 301
c
        read(line, '(f10.0,1x,f8.1,1x,3(f6.1,1x))')
     &  (SAProtcol(nSAparMAX*nSAP-nSAparMAX+k),k=1,nSAparMAX)
c:nSAparMAX: nStep Temp vdwSC whbBB whbSS
c
c SAProtcol:  nStep  T  a(SC)  wf(hb128BB)  wf(hb128SB,SS)
c
        CONTROL = .true.
        if(CONTROL)then
        write(kanalp,*)'nSAP=',nSAP
        write(kanalp,*)line
        write(kanalp, '(f10.0,1x,f8.1,1x,3(f6.1,1x))')
     &  (SAProtcol(nSAparMAX*nSAP-nSAparMAX+k),k=1,nSAparMAX)
        end if
c
        goto 300
c
301     continue
c
        CONTROL = .true.
        if(CONTROL)then
        write(kanalp, *)'SimAnnProtocol: loaded'
        write(kanalp,*)'nSAstep : ',nSAstep
        write(kanalp,'(a40)')'ntimeMX, tempTg, vwwSC wfHbBB wfHbSS: '
c
        do k=1,nSAP
        write(kanalp,'(f10.0,1x,f8.1,1x,3(f6.1,1x))')
     &  (SAProtcol(nSAparMAX*k-nSAparMAX+i),i=1,nSAparMAX)
        end do!k
        end if
c
        close(unit=kanalInSAprot)      ! SImAnnProtocol.inp
c
        end if !MDSA
c -------------------------------------------------------------------------
c OPT_EssModeAnalys
        OPT_EssModeAnalys=.false.
        OPT_EssModeType=0        !NO essential mode analysis    
        keyw = 'EssModeAnalys'
        keywl = 13
        keywt = 2 
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found .and. OPT_MDSA)then
        
        OPT_EssModeAnalys=.true.
        OPT_EssModeType=1                          !default CA analysis
        if(OPT_fullProtMD)then
        OPT_zeroRot = .true.
        else
        OPT_zeroRot = .false.
        end if
        if(valuec(1:3) .eq. "PBB")OPT_EssModeType=2
        if(valuec(1:3) .eq. "PCA")OPT_EssModeType=1
        if(valuec(1:2) .eq. "CA") OPT_EssModeType=1
        if(valuec(2:3) .eq. "CA") OPT_EssModeType=1
        end if!
c
c --------------------------------------------------------------------------
c freeEnergy from mdTra
c
        OPT_freeEnergy = .false.
        keyw = 'FreeEnergy'
        keywl = 10
        keywt = 2
c
        call read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
        if(found)then
        if(OPT_EssModeAnalys) OPT_freeEnergy = .true.
        end if
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c SELF-CONSISTENCY:
        if(coulVar_OPT .eq. 3)coulVar_OPT = 2 !OPT 3  is out  30.08.2005
c
        if(OPT_SolvGBorn)then
        OPT_SolvGS = .false.
        OPT_SolvSASHP = .true.
        coulVar_OPT = 4
        end if !OPT_SolvGBorn
c
        if(OPT_SolvGS)then
        OPT_SolvGBorn = .false.
	if(.not.OPT_SolvateExWat)OPT_SolvSASHP = .false.
        end if !OPT_SolvGS        
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
        write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
        write(kanalp,*)'MD Statistics:'
        write(kanalp,*)
c
        write(kanalp,*)'MD initial Temp(K)      :', (tempT0(k),k=1,3) 
        if(atype .eq. 1 )then
        write(kanalp,*)'MD enesemble            :', 'NTV: const Temp'
        write(kanalp,*)'MD bathTemperature(K)   :', tempTg
        else 
        write(kanalp,*)'MD enesemble            :', 'NEV: const Energy'
        end if
c
        write(kanalp,*)'LoopModeler01: paramFile:','./LoopModeler01.inp'
        write(kanalp,*)'runMDnstep              :',ntimeMX
        write(kanalp,*)'SHAKE                   :',iShake
        write(kanalp,*)' = 0(no),1(H-),2(all bonds) '
        write(kanalp,*)'ffParFile               :',ffParFile
        write(kanalp,*)'solvAtmTypeFile         :',solvAtmTypeFile
        write(kanalp,*)'solvAtmDatFile          :',solvAtmDatFile
        write(kanalp,*)'OPT_LoopMD              :',OPT_LoopMD
        write(kanalp,*)'OPT_fullProtMD          :',OPT_fullProtMD
        write(kanalp,*)'coulVar_OPT             :',coulVar_OPT
        if(coulVar_OPT .eq. 2)
     &    write(kanalp,*)'cReps_OPT(epsConstR)    :',cReps_OPT
c        
        write(kanalp,*)'SoftCoreVDW(Coul) aSC   :',aSoftCore
        write(kanalp,*)'OPT_SolvateExWat        :',OPT_SolvateExWat
        if(OPT_SolvateExWat)
     &  write(kanalp,*)'ExplicitWaterShell      :',dSOLVSheLL
        write(kanalp,*)'OPT_SolvGS              :',OPT_SolvGS
	write(kanalp,*)'OPT_SolvSASHP           :',OPT_SolvSASHP
	write(kanalp,*)'OPT_SolvWbr             :',OPT_SolvEFS
        write(kanalp,*)'ihydrPhSolv             :',ihydrPhSolv
        write(kanalp,*)'hydrPHsolv+ (kcal/m)    :',atomHydPhSolv
        write(kanalp,*)'OPT_harmAt1PosRst       :',OPT_harmAt1PosRst
        write(kanalp,*)'restr1AtConst (kcal/A^2):',restr1AtConst
        write(kanalp,*)'OPT_HW1MolPosRst        :',OPT_HW1MolPosRst
        write(kanalp,*)'restr1MHWConst(kcal/A^2):',restr1MHWConst
        write(kanalp,*)'OPT_MDSA                :',OPT_MDSA
        if(OPT_MDSA)then
        write(kanalp,*)'SAprotocolFile          :',fileSAprotocol
        write(kanalp,*)'nSAstep                 :',nSAstep
        do k=1,nSAP
        write(kanalp,*)
     &       'ntimeMX, tempTg(K)                : '
     &       ,SAProtcol(2*k-1),SAProtcol(2*k)
        end do!k
        end if
c
        write(kanalp,*)'xyzTraFile              :',xyzTraFile
        write(kanalp,*)'engTraFile              :', engTraFile
        write(kanalp,*)'pdbResFile              :',pdbResFile
        write(kanalp,*)'runMDnstepMAX           :',ntimeMX
        write(kanalp,'(a26,f8.5)')'mdTimeStep(ps)          :',deltat
        write(kanalp,*)'updateR1PList after n   :',ntimeR1
        write(kanalp,*)'updateR2PList after n   :',ntimeR2
        write(kanalp,*)'md MTS  ntimeF1(stepDt) :',ntimeF1   
        write(kanalp,*)'md MTS  ntimeF2(*F1)    :',ntimeF2   
        write(kanalp,*)'md MTS  ntimeF3(*F2)    :',ntimeF3   
        write(kanalp,*)'pdbFile                 :', pdbfile 
        if(OPT_LoopMD ) then 
        write(kanalp,*)'LoopStartEndFile        :',loopfile
        end if
        if(OPT_EssModeAnalys)then
        valuec = "CA"
        if(OPT_EssModeType.eq.2)valuec = "Prot BackBone Atoms"
        write(kanalp,*)'Do EssModeAnalysis on:', valuec
        end if !
        write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
        write(kanalp,*) 
        write(kanalp,*)'inputMDSApar: Done!: ' 
c       
         return
	 end
c
c subroutine to read value of keyword from file
c 
         subroutine read_keyw(filenam,keyw,keywt,keywl,
     &                        found,valuel,valuei,valuer,valuec)
 
         character*(*) filenam
         character*(*) keyw
         character*(*) valuec
         integer keywt             !type -1,0,1,2 : logical,integer,real,char
         integer keywl
         logical valuel,found
         real valuer
         integer valuei
c
         include "charStringSiz.h" 
         include "output.h"
         include "kanalUse.h"
         include "statusMessg_mDyn.h"
c local         
         integer linetot
         integer startc,lenc
         integer kanalin
         character*20 word
cx         character*80 line
         character*(charLenMAX) line
         logical CONTROL
         integer kanalp
c
         CONTROL = .true. 
         kanalp = kanalRunOut
c
cx         if(CONTROL)then
cx         write(kanalp,*)'In read_keyw:',
cx     &    filenam,keyw,keywt,keywl
cx         end if
c
         valuel=.false.
         valuer=0.0
         valuei=0
         valuec = '     '
c  
         linetot = charLenMAX
         found=.false.
         startc=keywl+3
         lenc = linetot - startc
c
         kanalin = kanalInMdynPar
c
          open(unit=kanalin,file=filenam,form='formatted',
     &        status='old')
c
          rewind kanalin
c
100       continue
          read(kanalin,'(a256)',end=1001)line   !IMPORTANT WARNING !! charLenMAX=256 = format a256   ! 
          word=line(2:keywl+1)
c 
          if(line(1:1).eq.'$')then  !keyw read
          if(word.eq.keyw)then
          found=.true.

          if(keywt.eq.-1)valuel=.true.
          if(keywt.eq.0)then
          read(line(keywl+3:lenc),'(i9)')valuei
          end if
          if(keywt.eq.1)then
          read(line(keywl+3:lenc),'(f9.3)')valuer
          end if
          if(keywt.eq.2)then
          valuec(1:lenc) = line(keywl+3:linetot)
          end if
c
          end if
          end if
          if(found)then
          goto 1001
          else
          goto 100
          end if
 
1001      close(kanalin)
c
         if(CONTROL .and. found)then
         write(kanalp,*)'In read_keyw:',
     &   'file:', filenam,' keyw,keywt,keywl:',keyw,keywt,keywl,
     &   ' found:',found, ' valuel,valuei,valuer:',
     &   valuel,valuei,valuer,' valuec:', valuec
         end if
c
          return
          end
c
