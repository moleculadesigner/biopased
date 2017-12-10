c* initialize: all nonBondList    for given:  atomXYZs()       
c              and given flags: makeVdWs,makeCLs,makeSLs                               
c
	subroutine initNonBondList(atomXYZs,makeVdWs,makeCLs,makeSLs)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'movingAtom.h'
c
        real atomXYZs(*)
        integer makeVdWs,makeCLs,makeSLs
c
	if(makeVDWs .eq. 0 ) return
c default
c        rcutV = 8.0
c        rcutC = 14.0
c        rcutS = 6.0
c        rbuffV = 0.75
c        rbuffS = 1.25
c        rbuffC = 1.25
c
c neighbour list
c
        call nonbondListVCS(ncallNBPL,rcutV,rcutC,rcutS,
     &           atomXYZs,atomQ,
     &           rbuffV,rbuffC,rbuffS,
     &           makeVdWs,makeCLs,makeSLs,
     &           natom,atomNumb,atomName,resName,chName,resNumb,
     &           nres,resNameRes,atomBlockName,
     &           atomNameEx,startAtInRes,
     &           nmoveatom,moveAtomList,moveFlag,
     &           pair12List,startPairL12,nPairL12,
     &           pair13List,startPairL13,nPairL13,
     &           pair14List,startPairL14,nPairL14,
     &           nbpairListV,startnbPairLV,nnbPairLV,nnbpLVT,nnbpLVMX,
     &           nbpairListC,startnbPairLC,nnbPairLC,nnbpLCT,nnbpLCMX,
     &           nbpairListS,startnbPairLS,nnbPairLS,nnbpLST,nnbpLSMX)
c
c
         return
         end
