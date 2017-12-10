c global variables
c
c nonBondedPairList for VdW, Coul, SolvationGS
c for given nmoveatom of the movingAtomList()
c movingAtom List
c
c        integer nmoveatom
c        integer moveAtomList(natomMAX)
c	common/moveAtom/nmoveatom,moveAtomList
c
        integer ncallNBPL     ! call counter for nonbondListVCS

c nonbonded neighbourList for GS solvation model
        integer makeSL                 ! 0/1 update flag
        real rcutS,rbuffS
        integer nnbpLSMAX
        integer nnbpLSMX
        integer nnbpLST
        parameter (nnbpLSMAX = natomMAX*300)  ! ~ 8 rcutOff
        integer nbpairListS(nnbpLSMAX)
        integer startnbPairLS(natomMAX)
        integer nnbPairLS(natomMAX)
c
	common/pairLS/nnbpLSMX,nbpairListS,
     &                startnbPairLS,nnbPairLS,
     &                rcutS,rbuffS,makeSL,nnbpLST
c 
c nonbonded neighbourList for VdW interections 
c 12,13,14 - pairs are excluded,
        integer makeVdW                ! 0/1 update flag
        real rcutV,rbuffV
        integer nnbpLVMAX
        integer nnbpLVMX
        integer nnbpLVT
        parameter (nnbpLVMAX = natomMAX*300)  ! ~ 8 rcutOff
        integer nbpairListV(nnbpLVMAX)
        integer startnbPairLV(natomMAX)
        integer nnbPairLV(natomMAX)
c
	common/pairLVDW/nnbpLVMX,nbpairListV,
     &          startnbPairLV,nnbPairLV,rcutV,rbuffV,
     &          makeVdW,nnbpLVT,ncallNBPL
c 
c nonbonded neighbourList for Coul interections
        integer makeCL                 !1/0 update flag
        real rcutC,rbuffC
        integer nnbpLCMAX
        integer nnbpLCMX
        integer nnbpLCT
        parameter (nnbpLCMAX = natomMAX*400)  ! ~ 14 rcutOff
        integer nbpairListC(nnbpLCMAX)
        integer startnbPairLC(natomMAX)
        integer nnbPairLC(natomMAX)
c
	common/pairCL/nnbpLCMX,nbpairListC,startnbPairLC,
     &                nnbPairLC,rcutC,rbuffC,makeCL,nnbpLCT
c
cend
