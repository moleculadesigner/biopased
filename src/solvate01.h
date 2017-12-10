c Yury Vorobje 2005
c solvationBox + solvation shell water
c
        integer nWatInBox0MAX
        parameter (nWatInBox0MAX=216)
        integer natIn1watMol
        parameter (natIn1watMol=3)
        character*4 watAtName(natIn1watMol)
c
        integer nWaterMolMAX
        parameter (nWaterMolMAX = 15000)
        integer nAtomWatMAX
        parameter (nAtomWatMAX = nWaterMolMAX*natIn1watMol)
c
        integer nWatMolAdd
        integer nWatMolAtAdd
        real atomXYZWat(3*nAtomWatMAX)
        integer resNumbWat(nAtomWatMAX)
        character*4 atomNameWat(nAtomWatMAX)
        character*4 resNameWat(nAtomWatMAX)
        integer waterMolTag(nWaterMolMAX)
c
c soluteMOL+Water Heavy atoms: noH for MOLECatoms and WATers
cx        integer natomMOLnoH
cx        integer natMOLandWatnoH
cx        real atomXYZMolWatnoH(3*nAtomWatMAX)
cx        character*4 atomNamenoH(nAtomWatMAX)
cx        character*4 resNamenoH(nAtomWatMAX)
cx        integer resNumbnoH(nAtomWatMAX)
c
c neihbourwPairList
cx        integer ncallwPL
cx        real rCUTwPL
c local WpairList dataStructure
cx        integer nnbwPLMAX
cx        parameter (nnbwPLMAX = nAtomWatMAX*200)    !
cx        integer nnbwPairL(natomMAX)
cx        integer startnbwPairL(natomMAX)
cx        integer nbWpairList(nnbwPLMAX)
cx        real nbWpairListD2(nnbwPLMAX)
cx        integer nnbWpLT,nnbWpLMX
c
c WatInBOX0
c
        integer nWatAtInBox0Max
        parameter (nWatAtInBox0Max=natIn1watMol*nWatInBox0MAX)
        integer nWatAtInBox0,nWatMoLInBox0
        real watXYZinBox0(3*nWatAtinBox0MAX)
        character*4 watAtNameInBox0(nWatAtinBox0MAX)
        character*4 watResNameInBox0(nWatAtinBox0MAX)
        integer     watResNumbInBox0(nWatAtinBox0MAX)
        character*1 watchNameInBox0(nWatAtinBox0MAX)
        integer watAtNumbInBox0(nWatAtinBox0MAX)
c 
        real dSOLVSheLL
        common/solvbox01/dSOLVSheLL
c
	integer natomSMOL,nresSMOL    !nAt in solute molecule
	common/solvate02/natomSMOL,nresSMOL,
     &        nWatMolAdd,nWatMolAtAdd 
c
