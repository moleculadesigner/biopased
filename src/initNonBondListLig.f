c* initialize: all nonBondListLig    for given:  atomXYZs()       
c              and given flags: makeVdWs,makeCLs,makeSLs          
c USE: nonBondListLig  to calculate Ligand Eng on the Molecule
c      if MOving RESlist > LIGresList                     
c
	subroutine initNonBondListLig(atomXYZs,makeVdWs,makeCLs,makeSLs)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'ligInfo.h'
        include "output.h" 
c
        logical CONTROL
        real atomXYZs(*)
        integer makeVdWs,makeCLs,makeSLs
        integer kanalp
c
	kanalp = kanalRunOut
        CONTROL = .false.
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
         nnbpLVLigMX = nnbpLVLigMAX
         nnbpLCLigMX = nnbpLCLigMAX
         nnbpLSLigMX = nnbpLSLigMAX
c
         if(CONTROL)then
         write(kanalp,*)'initNonBondListLig:ST: ncallNBPLLig:',
     &             ncallNBPLLig
         end if ! CL
c
         call nonbondListVCS(ncallNBPLLig,rcutV,rcutC,rcutS,
     &           atomXYZs,atomQ,
     &           rbuffV,rbuffC,rbuffS,
     &           makeVdWs,makeCLs,makeSLs,
     &           natom,atomNumb,atomName,resName,chName,resNumb,
     &           nres,resNameRes,atomBlockName,
     &           atomNameEx,startAtInRes,
     &           nAtomInLig,atomInLigList,atInLigFlag,
     &           pair12List,startPairL12,nPairL12,
     &           pair13List,startPairL13,nPairL13,
     &           pair14List,startPairL14,nPairL14,
     &           nbpairListVLig,startnbPairLVLig,nnbPairLVLig,
     &           nnbpLVLigT,nnbpLVLigMX,
     &           nbpairListCLig,startnbPairLCLig,nnbPairLCLig,
     &           nnbpLCLigT,nnbpLCLigMX,
     &           nbpairListSLig,startnbPairLSLig,nnbPairLSLig,
     &           nnbpLSLigT,nnbpLSLigMX)
c
          if(CONTROL)then
          write(kanalp,*)'initNonBondListLig:Fin: ncallNBPLLig:',
     &    ncallNBPLLig,' nnbpLVLigT,nnbpLCLigT,nnbpLSLigT:',
     &    nnbpLVLigT,nnbpLCLigT,nnbpLSLigT
c
          write(kanalp,*)'initNonBondListLig : FINISH!'
          end if !CL
c
         return
         end
