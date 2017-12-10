c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c  bison                                                            *
c                                                                   *
c  Yury Vorobjev 2002                                               *
c                2005 SSbonds                                       *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c initialize  PART1: MoleculeTopology from AA sequence                     
c
	subroutine initMolecTopSeq01
c
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'filedat.h'
        include 'pair1234array.h'
        include 'loopInfo.h'
        include 'optionPar.h'
        include 'solvGSarray.h'
        include 'coulEnPar.h'
        include 'ssbond.h'
c
        integer nld,kanalp
        integer i,ia,n,ia3
        integer j,ja,ja3
        integer iHread
        real dss2
        logical CONTROL
        
c replica of array sizes
        natomMX = natomMAX
        nresMX = nresMAX
        npL12MX = npL12MAX
        npL13MX = npL13MAX
        npL14MX = npL14MAX
        npL123MX = npL123MAX
        npL1234MX = npL1234MAX
        nLoopMX = nLoopMAX
        nAtomInLoopMX=nAtomInLoopMAX
        nAnchorAtomMX =  nAnchorAtomMAX
        nImp1234MX = nImp1234MAX 
c
cx      kanalp = 6
        kanalp = kanalRunOut
c
        CONTROL = .true.
c
	write(kanalp,*) 'initMolecTopSeq01: START: ******************* '
        write(kanalp,*)
c
c AAcid Sequence:
c        resSeqFile = './resSeq.inp'
c
         write(kanalp,*)'PDBfile InPut: ',pdbfile
cx         write(kanalp,*)'resSeqFile :', resSeqFile
c
c define iHread
         iHread = 0
         if(OPT_Hread)iHread = 1
c
c Input residue Sequence data for protein from: resSeqFile
         if(OPT_resSeqInp)then
c
         write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
         write(kanalp,*)'There Is Not a complelete PDB file of the XYZ',
     &                  ' and RES sequence !!'
         write(kanalp,*)'Topology will be generated from resSeqFile:',
     &                   resSeqFile
         write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
c
	 call readResSeq(resSeqFile,nres,resListSeq,
     &                    nLoop,resLoopFlag,resEndLoop)
c
         end if !OPT_resSeqInp
c
c Input res sequence from XYZ.pdb file
         if(OPT_PDBresSeq)then    ! read resListSeq from molec.pdb IN file
c
         call readPDB05(pdbfile,natomMAX,nresMAX,iHread,
     &       natom,atomNumb,atomName,resName,chName,resNumb,
     &       atomXYZ,nres,resNameRes,chNameRes,
     &       atomNameEx,resNameEx,
     &       startAtInRes,stopAtInRes,realAtomFlag,
     &       resTypeInREseq,atHvyNbInRes,
     &       nChain,startAtInCha,stopAtInCha,
     &       startResInCha,stopResInCha,
     &       nMolec,startAtInMol,stopAtInMol,
     &       startResInMol,stopResInMol)
c
c copy to  addHxyz.h
           call readPDB05_nh
           write(kanalp,*)'readMolPDB: end of readPDB05_nh'
c
         write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
         write(kanalp,*)
     &   'WARNING!: resListSeq : are taken from INpdbFile:',pdbfile
         write(kanalp,*)
     &   'It is assumed that INpdbFile is COMPLETE(NO LOOPS withoutXYZ)'
         write(kanalp,*)'* * * * * * * * * * * * * * * * * * * * * * * '
c
c Correct RESnames and Atom Names from XYZ.pdb file to the amberZmatrix nameStyle
c call getAmbZmAtomNamef:
        write(kanalp,*)'initMolecTopSeq01: start: getAmbZmAtomNamef:'
c
        call getAmbZmAtomNamef(pdbAmbZmFile,
     &           natom,atomName,atomNameEx,resName,
     &           resNameEx,resNumb,nres,resNameRes)
 
c
        write(kanalp,*)'initMolecTopSeq01: end of getAmbZmAtomNamef:'
c
         do i = 1,nres
         resListSeq(i) = resNameRes(i)(1:4) 
         resLoopFlag(i) = 0   
         end do !nres
c
         end if ! OPT_PDBresSeq
c
c read Zmatrx of RESidues and assemble zMatrix for the AASequence
c define atomPDB info put it into Global xyzPDBinfo.h
c calculate pair12List and put it in Global array
c
         call  getZMatrixREseq(nres,resListSeq,
     &          resTypeInREseq,atom3ToLinkInREseq,atomBlockName,
cx     &              resTypeInREseq,atomBlockName,
     &              natom,atomName,ffatomName,
     &              resName,resNumb,atomNameEx,
     &              startAtInRes,stopAtInRes,
     &              nChain,startAtInCha,stopAtInCha,
     &              startResInCha,stopResInCha,
     &              nMolec,startAtInMol,stopAtInMol,
     &              startResInMol,stopResInMol,
     &              pair12List,startPairL12,
     &              nPairL12,atomQ)
c
c subroutine readMolPDB:
c 1)  read (good or  bad incomplete) PDB file of ProteinCore
c     addXYZ of missing HeavyAtoms(sideChains) in Residues
c     addXYZ of missing Hydrogens
c 2)  put XYZ in the right ORDER as it is in the resSeqZmatrix
c 3)  defines xyzPDBinfo.h + xyzPDBcrd.h
c 4)  defines:
c        realAtomFlag(ia) for ProteiCore atoms with XYZ
c
         call readMolPDB
c
         call getSSbonds 

c
c define ALL Loop atoms 
c define atomInLoopList() - list of atoms in all LOOPs to be added
c                         - it works if nLoop > 0 and  OPT_AddLoop=true
c
c read loopMD.inp (moveRes.inp)  file to define moving atom List
        nLoop = 0      !default i.e. noLoop, fullProtMD
        if(OPT_LoopMD)then
        call readLoopInf(loopFile,nLoopMX,nLoop,resEndLoop)
c all loops are moving
        nMoveLoop = nLoop
        do i = 1,nLoop
        moveLoopList(i) = i
        end do !i 
        end if !OPT_LoopMD
c
        addAnchorAtF = 0     ! no anchorAtoms
c
        call inLoopAtomList(addAnchorAtF,
     &           nLoopMX,nLoop,resEndLoop,
     &           natom,atomName,resName,chName,resNumb,
     &           nres,resNameRes,chNameRes,
     &           atomNameEx,startAtInRes,
     &           nAtomInLoop,nAtomInLoopMX,atomInLoopList,
     &           startAtInLoop,stopAtInLoop,loopType,
     &           nAnchorAtom,nAnchorAtomMX,anchorAtomList,
     &           atomXYZ,refPosAnchorAtom)
c
c store initial ff-atomic Q 
         do i=1,natom
         atomQ0(i) = atomQ(i)     ! atomQ0() original FF atomQ
         end do !
c neutralize Q Phosphote groups
         if (iqPhosphSc_OPT .eq. 1)then
         write(kanalp,*)
         write(kanalp,*)'SCALE Q of NA Phosphate groups: qPhosphSc=',
     &                   qPhosphSc_OPT
         write(kanalp,*)
         do i=1,natom
         if(atomBlockName(i)(1:1) .eq. 'P')
     &         atomQ(i) = atomQ0(i)*qPhosphSc_OPT
         write(kanalp,7071)'ATOM',i,atomName(i),resName(i),chName(i),
     &   resNumb(i),
     &   atomXYZ(3*i-2),atomXYZ(3*i-1),atomXYZ(3*i),atomBlockName(i),
     &   atomQ(i),atomQ0(i)
         end do!i
         end if !iqPhosphSc_OPT
c
c NEXT PART2 of initMOlecTopology :
c define:  atomMovList(ia) and realAtomFlag(ia)
c and initMolecTopSeq02 
c for a GIVEN loopList and OPT_FullProt
c all 12,13,14 etc. list will be calculated for a given atomMovList()
c     and realAtomFlag(ia)
c assuming that all other atoms are fixed
c
	write(kanalp,*) 'initMolecTopSeq01: FINISH: successfully !!! '
        write(kanalp,*)
c
7071    format(a6,i5,2x,a4,a4,a1,i4,4x,3f8.3,1x,a4,1x,2f8.4) !orig PDB
c
	return
        end
