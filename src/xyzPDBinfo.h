c xyzPDBinfo.h     
c 2010: atom3ToLinkInREseq()
c
        character*(charLenMAX) pdbfile 
        character*8 atomNameEx(natomMAX)     !extendedName atm(1:4)+termTag(5:8)
        character*4 atomName(natomMAX)
        character*4 atomBlockName(natomMAX)  !blockName: P,R,B -NA
        character*4 resName(natomMAX)
        character*8 resNameEx(natomMAX)
        character*4 resNameRes(nresMAX) 
        character*6 resListSeq(nresMAX)
        character*1 chName(natomMAX)
        character*1 chNameRes(nresMAX)
        character*2 ffAtomName(natomMAX)    !ff atom code
        integer atHvyNbInRes(nresMAX)       !number of HeavyAtoms in Residue
        integer resNumb(natomMAX)
        integer atomNumb(natomMAX)
        integer natom,nres                   !current natom, nres
        integer natomMX,nresMX
cx       real    atomXYZ(3*natomMAX)        ! in xyzPDBcrd.h 
c
        real atomQ(natomMAX)                ! current atomQ
        real atomQ0(natomMAX)               ! original FF atomQ
        real atomRad(natomMAX)              ! atomRad for SAS
        real atomMass(natomMAX)
cx        real atomDconst(natomMAX)         !local DielConst in coulEnPar.h
        integer atomVDWtype(natomMAX)
c
        integer startAtInRes(nresMAX+1)
        integer stopAtInRes(nresMAX)
        character*9 resTypeInREseq(nresMAX) !BEG,ISO,INT,FIN + AA,NA,NAL residue type in sequence 
        character*12 atom3ToLinkInREseq(nresMAX) !3linkerAtomsRes(i)-->DUMM[res(i+1)
        integer resLoopFlag(nresMAX) 
c
	common/pdbInfo1/pdbfile
        common/pdbInfo2/atomName,atomNameEx,resNameEx,
     &         resName,resNameRes,chName,chNameRes,ffAtomName,
     &         resListSeq,atomBlockName
c
	common/pdbInfo3a/natomMX,natom,resNumb,atomNumb,
     &                  atomVDWtype,
     &                  atomMass,atomQ,atomQ0,atomRad
c
        common/pdbInfo3b/nresMX,nres,startAtInRes,stopAtInRes,
     &                   resTypeInREseq,atHvyNbInRes,
     &                   atom3ToLinkInREseq
c
        integer realAtomFlag(natomMAX)
        common/pdbInfo3c/realAtomFlag,resLoopFlag
c
        integer addHvyAtFlag(natomMAX)         ! 0/1 yes/no being added
        common/pdbInfo3d/ addHvyAtFlag
c
        integer nChain
        integer startAtInCha(nresMAX+1)
        integer stopAtInCha(nresMAX)
        integer startResInCha(nresMAX+1)
        integer stopResInCha(nresMAX)
        integer nMolec
        integer startAtInMol(nresMAX+1)
        integer stopAtInMol(nresMAX)
        integer startResInMol(nresMAX+1)
        integer stopResInMol(nresMAX)
        common/pdbInfo3e/nChain,startAtInCha,stopAtInCha,
     &                  startResInCha,stopResInCha,
     &                  nMolec,startAtInMol,stopAtInMol,
     &                  startResInMol,stopResInMol
c
c end
