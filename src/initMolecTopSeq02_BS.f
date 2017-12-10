c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2003                                               *
c             PWatEx:2009 separated Temp of Solute and Solvent molecules
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c initialize PART2 of: MoleculeTopology from AA sequence                     
c define molecTopology:  
c 1) define: moveAtomList(ia)
c 2) for a given realAtomFlag(ia), moveAtomList(ia) and pair12List()
c calculate:  bondPair12(), neighbours 13, 14; triplets 123, quartets 1234
c
	subroutine initMolecTopSeq02
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
        include 'movingAtom.h'
        include 'shake.h'
	include 'enForce_PW.h'
c
clocal
        integer ivar1234,ia,iam
	integer scPWFDG(3)
        integer kanalp
        logical CONTROL
c        
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
c initialize:
cx        kanalp = 6
        kanalp = kanalRunOut
        CONTROL = .true.
c
        ivar1234 = 0
        if(OPT_LoopMD)ivar1234 = 1
c
c define molecTopology:  for given realAtomFlag(ia), moveAtomList(ia)
c and pair12List()
c calculate:  bondPair12(), neighbours 13, 14; triplets 123, quartets 1234
        call vbondListAAseq(ivar1234,atomXYZ,realAtomFlag,
     &           natom,atomNumb,atomName,resName,chName,resNumb,
     &           nres,resNameRes,atomBlockName,
     &           atomNameEx,startAtInRes,
     &           nmoveatom,moveAtomList,moveFlag,
     &           pair12List,startPairL12,nPairL12,npL12MX,
     &           pair13List,startPairL13,nPairL13,npL13MX,
     &           pair14List,startPairL14,nPairL14,npL14MX,
     &           pair14VList,startPairL14V,nPairL14V,
     &           bond12List,nbond12,nbond12noH,
     &           trip123List,nTrip123,npL123MX,
     &           quar1234List,nQuar1234,npL1234MX,
     &           quarImp1234L,nImp1234,nImp1234MX)

c
c define SHake related variables:
        if(iShake .eq. 0) nbondShake = 0
        if(iShake .eq. 1) nbondShake = nbond12-nbond12noH
        if(iShake .ge. 2) nbondShake = nbond12
c number degrees of freedom of whole system
        nAtFreeDg(1) = 3*nmoveatom - nbondShake 
c nAtFreeDg(*): 1,2,3 = totalSystem, ProteinAtoms, WaterSolvent
        nAtFreeDg(2)=0
	nAtFreeDg(3)=0
        scPWFDG(1)=3
	nMovingAtom(1)=0
	nMovingAtom(2)=0
	nMovingAtom(3)=0
        do iam = 1,nmoveatom
	ia = moveAtomList(iam)
	scPWFDG(2)=3
	scPWFDG(3)=2
        if(atomPWatTag(iam).eq.1)then
	if(iShake .ge. 1)scPWFDG(2) = 2
	nAtFreeDg(2)=nAtFreeDg(2) + scPWFDG(2)
	nMovingAtom(2)=nMovingAtom(2) + 1
	end if!
	if(atomPWatTag(iam).eq.2)then
	nAtFreeDg(3)=nAtFreeDg(3) + scPWFDG(3)
	nMovingAtom(3)=nMovingAtom(3) + 1
	end if!
	end do!ia
c 
	nAtFreeDg(1) = nAtFreeDg(2)+nAtFreeDg(3)
	nMovingAtom(1)=nMovingAtom(2)+nMovingAtom(3)
c
        write(kanalp,*)
     &	'initMolecTopSeq02: nAtFreedomDegree(tot,solute,solvent):',
     &   nAtFreeDg
        write(kanalp,*)
     &  'initMolecTopSeq02: nMovingAtom(tot,solute,solvent):',
     &   nMovingAtom
c
	return
        end
