c global variables
c Info: atoms in ligand
        character*(charLenMAX) ligResFile     ! defines Lig residues
c file name ligRes.inp   : default
        character*(charLenMAX) molAtXYZIntLigFile  ! PDB file for protAtom local to Lig Atoms
        real rAtInLig
c molAtXYZintLig.pdb : default
        integer nLigdMAX
        parameter (nLigdMAX = 1)
        integer resStEndLig(2*nLigdMAX)
        integer nAtomInLigMAX
        parameter (nAtomInLigMAX=natomMAX/4) 
        integer atomInLigList(nAtomInLigMAX)  ! atomInLigList(iaL) = global iaNumb=iatGL
c                                               ! iaL = 1,..,nAtomInLig  
        integer atInLigFlag(natomMAX)
        integer nDockVarMAX
        parameter (nDockVarMAX = 10)          ! number of best dockVariants
        integer nDockVar
c
        real atomLigXYZ(3*nAtomInLigMAX)
        real atomLigXYZ0c(3*nAtomInLigMAX)     ! ligXYZ relative ligGeoCentr from PDB
        real atomLigXYZcGeo0(3)                ! initial from PDB GeoCentr
        real atomLigXYZcGeo(3)                 ! current geoCenter
c
        real atomLigXYZdockVar(nDockVarMAX*3*nAtomInLigMAX)
	real engLigdockVar(nDockVarMAX)
c
        integer nLigd,nAtomInLig               ! nLigd - number of Ligands
        integer nLigdMX
        integer nAtomInLigMX
c
        common/ligInfo01/ ligResFile
        common/liginfo02/rAtInLig,molAtXYZIntLigFile
c
	common/ligInfo03/ nLigd,nAtomInLig,resStEndLig,
     &                   atomInLigList,atInLigFlag
c
        common/ligXYZ/ atomLigXYZ0c,atomLigXYZ,
     &                 atomLigXYZdockVar
c
        common/ligXYZcG/ atomLigXYZcGeo0,atomLigXYZcGeo
c
c ligandAt-moleculeAt pairList
c
        integer makeVdWLig
        integer makeCLLig
        integer makeSLLig
c
        integer ncallNBPLLig
        integer nnbpLVLigMAX
        parameter (nnbpLVLigMAX=nAtomMAX*200)  !~ 8 A cutOff
        integer nbpairListVLig(nnbpLVLigMAX)
        integer nnbPairLVLig(nAtomMAX)
        integer startnbPairLVLig(nAtomMAX)
        integer nnbpLVLigT                     ! current tot numb of neigb in the List
        integer nnbpLVLigMX
        common/ligNbondVL/nbpairListVLig,startnbPairLVLig,
     &                   nnbPairLVLig,nnbpLVligT,ncallNBPLLig
c
        integer nnbPairLCLig(nAtomMAX)    
        integer startnbPairLCLig(nAtomMAX)
        integer nnbpLCLigMAX
        parameter (nnbpLCLigMAX=nAtomMAX*400)  ! ~ 14 A cutOff
        integer nbpairListCLig(nnbpLCLigMAX)
        integer nnbpLCligT                     ! total number of neigbours in the List
        integer nnbpLCligMX
        common/ligNbondCL/nbpairListCLig,startnbPairLCLig,
     &                   nnbPairLCLig
c
        integer nnbPairLSLig(nAtomMAX)
        integer startnbPairLSLig(nAtomMAX)
        integer nnbpLSLigMAX
        parameter (nnbpLSLigMAX=nAtomMAX*300)
        integer nbpairListSLig(nnbpLSLigMAX)
        integer nnbpLSLigT                     ! total number of neigbours in the List
        integer nnbpLSLigMX
        common/ligNbondL/nbpairListSLig,startnbPairLSLig,
     &                   nnbPairLSLig
c
c Lig-molecule energies:
        real eVbondDefLg,eVangDefLg
        real eImpDefLg,eTorsDefLg
        real engVDWR1Lg,engCOULR1Lg
        real engCOULR2Lg
        real restr1MHWEngLg,restr1EngLg
        real molSolEnLg,eGeoDefLg
        real eGeoRstLg,engCOULLg
        real hBHxYeng128Lg
        real engPOTENTLg
c
        common/ligMolEng/ eVbondDefLg,eVangDefLg,eImpDefLg,eTorsDefLg,
     &                    engVDWR1Lg,engCOULR1Lg,engCOULR2Lg,
     &                    restr1MHWEngLg,restr1EngLg,molSolEnLg,
     &                    eGeoDefLg,eGeoRstLg,engCOULLg,engPOTENTLg,
     &                    hBHxYeng128Lg
c end  
