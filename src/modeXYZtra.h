c essential mode analysis: Mode: EssentialMode:FreeEnregy
c
c average sturucture
        real atomXYZtraAve(3*natomMAX)
        real atomXYZRmsdAve(natomMAX)
        real atomXYZRmsdTra(nStepTraMAX)    ! averageRMSD for snaps
        common/modeXYZ00/atomXYZtraAve,atomXYZRmsdAve,
     &        atomXYZRmsdTra
c
c matrixPositionalFluctuations
c
        integer nAtomTypeMPFMAX       ! atomTypes to be included in MatrixPosFluct
        parameter (nAtomTypeMPFMAX=4) ! CA,N,C,O for proteins
        integer nAtomTypeMPF       ! default =1 CA atoms; =4 [CA,N,C,O]
        integer OPT_EssModeType    !=1 CA of prot backbone for MatrixPositionaFluct analysis
c                                  !=2 CA N C O all ProtBB atoms
        character*4 atomTypeMPFList(nAtomTypeMPFMAX) ! =[CA,N,C,O for proteins]
c
c proteinBackBone atoms List
        integer atomNumbProtBBList(natomMAX)   !list of protein backbone atoms
        integer nAtomProtBB
        integer nProtBBatRes
        common/protBB00/ OPT_EssModeType,nAtomTypeMPF,
     &                 nProtBBatRes,nAtomProtBB,atomNumbProtBBList
c
        integer nAtomMPFMAX       !max number of atom for MatrixPositionaFluct (MPF)
cx      parameter(natomMPFMAX=natomMAX/10)  ! ~ 250 resProtein for [CA,N,C,O] 
cx IMportant RAMemory consuming parameter RAM ~ 10*10*(natomMPFMAX*natomMPFMAX) bite
        parameter(natomMPFMAX=550)  ! ~ 125 resProtein for [CA,N,C,O]
c
        integer atomNumbMPFList(natomMPFMAX)  !atom numbers list in MatrixPosFluct analysis
        integer nAtomMPFMX         ! = natomMPFMAX
        integer nAtomMPFList       ! current number of atoms included in MPF analysis
        common/modeXYZ01/atomTypeMPFList,atomNumbMPFList
c
        real atomXiXjAve(3*natomMPFMAX,3*natomMPFMAX)
        common/modeXYZ02/atomXiXjAve
c        
        integer nrotJac                !nRotationJacobi diagonalization
        real lambdaMPF(3*natomMPFMAX)
        real essModeRMSDperAtom(3*natomMPFMAX)
        real essModeFreq(3*natomMPFMAX)
        real eigVectMPF(3*natomMPFMAX,3*natomMPFMAX)
        real bJacobi(3*natomMPFMAX),zJacobi(3*natomMPFMAX) ! tempArray for jacobiGen
        real essModeEntropy
	common/modeXYZ03/lambdaMPF,essModeRMSDperAtom,
     &                   essModeFreq,essModeEntropy
        common/modeXYZ04/eigVectMPF
c
        integer nModeXYZmAX
        parameter (nModeXYZmAX=10)  ! max number of essModes reconsructed
        integer nModeSnapMAX        ! max Snap permode:2,4
        parameter (nModeSnapMAX=4)  ! xyz+/xyz- eigVector
        integer nModeXYZ
        integer nModeXYZSnap
        real atEssModeXYZ(nModeXYZmAX*nModeSnapMAX*3*nAtomMAX) 
        common/modeXYZ05/nModeXYZ,nModeXYZSnap,atEssModeXYZ
cEND
