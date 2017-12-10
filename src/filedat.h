cglobal variables 
c filenames  of data files
        character*(charLenMAX) HOMEBS_dir        ! HOMEdir name
        integer ndirLett
        common/filedat00/HOMEBS_dir,ndirLett
c
c  LIBrary files:
c
	character*(charLenMAX) ffAtomTypeFile    !convertor PDBname -> ffatomName
        character*(charLenMAX) ffParFile         ! ff parameters file
        character*(charLenMAX) resSeqFile        ! residueSeqFile
        character*(charLenMAX) pdbAmbZmFile      ! pdbAmbZmFile table to convert InPdb 
        character*(charLenMAX) hInfoFile         ! hInfoFile  to add H atoms
        character*(charLenMAX) fileWATBOX0       ! watBOX0pdb
c
	common/filedat01/ffAtomTypeFile,ffParFile,
     &             resSeqFile,pdbAmbZmFile,hInfoFile,fileWATBOX0
c
c myMoleculeXYZData files
c
        character*(charLenMAX) xyzInpFile        ! input XYZ.pdb for molecule
        common/filedat02/ xyzInpFile
c
        character*(charLenMAX) mdRestXYZVfile    ! mDyn RESTART file atomXYZ+velocities
        common/filedat03/mdRestXYZVfile  
c
c mDynProtocol files:
c
c mdynSB,simAnnProtocol,loopRes data files
        character*(charLenMAX) mdynParFile        ! mDynProtocolFile
        character*(charLenMAX) loopFile           ! movingResidueList
        character*(charLenMAX) saProtFile         ! simulatedAnnealing protocol
c
        common/filedat04/mdynParFile,loopFile,saProtFile
c
        character*(charLenMAX) restrA1type       ! defined A1 types of restraints
        character*(charLenMAX) restrA2type
        common/filedat05/restrA1type,restrA2type         
c
c RESULT files
c
        character*(charLenMAX) xyzTraFile       ! XYZV mdTrajectory  file
        character*(charLenMAX) engTraFile       ! energy mdTrajectory  file
        character*(charLenMAX) pdbResFile       ! XYZ.pdb mdTrajectory  files
        character*(charLenMAX) pdbFinFile       ! XYZVfinal.pdb mdResult file
        character*(charLenMAX) pdbOptFile       ! XYZenergyOpimization.pdb energy optimized XYZ.pdb
        common/filedat06/xyzTraFile,engTraFile,pdbResFile,pdbFinFile,
     &        pdbOptFile
c 
        character*(charLenMAX) molNameArgLine
        integer nMolNameLett
        common/filedat08/molNameArgLine,nMolNameLett
c
        character*(charLenMAX) fileWatSheLL     !XYZwater.pdb SolvationWaterMolecShell
        character*(charLenMAX) fileIonASXYZ     !XYZion.pdb   counterIons
        character*(charLenMAX) engTablFile1     !energy table file
        character*(charLenMAX) engTablFile2     !energy table file
        common/filedat09/fileWatSheLL,fileIonASXYZ,
     &                   engTablFile1,engTablFile2
c
	integer nRecPdb              !number of molResxxxx.pdb written
        common/filedat10/nRecPdb
c
        character*(charLenMAX) resultFileDir     !absolute name of Dir to write all result files
        integer nResFileDirLett
        common/filedat11/ resultFileDir, nResFileDirLett
c
