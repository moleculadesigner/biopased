c data structure for implicit waterBridges solvation model
c
	integer nWBrgMAX
        parameter (nWBrgMAX = 1000)
c
        real dipMWEFS_OPT      ! 0.477
        real eBindMIN_OPT      ! 0.6
        real wBrgEpotQMAX_OPT   ! -2.50 kcal
        real epsEFSMod,eLectScaleEFSMod
        real RcutEFSMAX 
        common/solWBrg00/dipMWEFS_OPT,epsEFSMod,eLectScaleEFSMod,
     &                   eBindMIN_OPT,RcutEFSMAX,wBrgEpotQMAX_OPT
c
        integer nWBrgNow
        real eFieldMX(nWBrgMAX)
        integer ideFMX(nWBrgMAX)
        common/solWBrg01/nWBrgNow,ideFMX,eFieldMX
        real wBrgXYZ(3*nWBrgMAX)
        real wBrgWmoLXYZ(3*3*nWBrgMAX)
        real wBrgAtomQ(3*nWBrgMAX)
        real wBrgDip(3*nWBrgMAX)   ! dipMom vector of wBridge waterMolecule
        real wBrgEpot(nWBrgMAX)    !pot Eng of wBrg dipol in the molEfield in build up
c                                   procedure
        real wBrgEpotQ(nWBrgMAX)  !pot Eng of wBrg q-q electEng
        real scaleWBrgEng(nWBrgMAX) !scalingCoeff
c
        common/solWBrg02/wBrgXYZ,wBrgDip,wBrgEpot,wBrgEpotQ,
     &                   scaleWBrgEng 
c
        common/solWBrg03/wBrgWmoLXYZ,wBrgAtomQ
c
        character*4 wBrgAtomName(3*nWBrgMAX),wBrgResName(nWBrgMAX)
        common/solWBrg04/wBrgResName,wBrgAtomName
c
c RESult:
        character*(charLenMAX) fileEFSwbrXYZ     !XYZwaterBrgidge.pdb file = watBrgSolvXYZ.pdb
        common/solWBrg05/fileEFSwbrXYZ  
c
