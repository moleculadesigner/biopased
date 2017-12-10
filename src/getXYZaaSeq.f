c generate XYZ of AminoAcid sequence of nresLoop Residues 
c          with given Fi,Psi of residues
c
        subroutine getXYZaaSeqFiPsi(nresLoop,resNameLoop,          
     &              FiPsiLoop,natLoop,atLoopXYZ,
     &              atQLoop,atNameLoop,ffatNameLoop,
     &              startAtResLoop,atLoopCMatrx)         
c
c IN:
c nresLoop - number Res in Loop
c resNameLoop(ir), ir=1,nresLoop
c FiPsiLoop(2*(nresLoop+2)) - Fi,Psi(ir)
c     FiPsi(1)=0.  FiPsi(2) = Psi(-1 res) defines DUMM atoms
c     FiPsi(3)=Fi(ir=1), FiPsi(4)=Psi(ir=1) 
c     FiPsi(2*(nresLoop+1)+1,2) = Fi,Psi of terminal GLY (additionalREs)
c
c OUT:
c  natLoop - total NAtoms
c atLoopXYZ(*),atNameLoop()-atNames
c resNameLoop(*) resNames for atoms
c ffatNameLoop(*) ff atom names
c startAtResLoop(*) - starting atoms of residue ir
c atLoopCMatrx(*) - connectivityMatrix
c 
cx        implicit none
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'zmatrix.h'
        include 'zmatrixtmp.h'
        include 'filedat.h'
c
        integer nresLoop
        character*4 resNameLoop(*)
        character*4 atNameLoop(*)
        real FiPsiLoop(*)
        integer natLoop
        real atLoopXYZ(*)
        character*2 ffatNameLoop(*)
        integer startAtResLoop(*)
        real atQLoop(*)
        integer atLoopCMatrx(*)
c local
        integer nres
        integer nAtPairCycMX
        integer ir,ia,ia3,ir2
        integer nz,k,i,i2,i5
        integer ia4,nL3,nL4,ns
        integer zmflag
        integer kanalp
        logical CONTROL
        integer unitResLib,unitResLibB,unitResLibE
        integer unitRes
        character*80 line
        integer natinres
        integer dumFlag
c
        kanalp = kanalRunOut 
        CONTROL = .true.
c
       zmfile =    HOMEBS_dir(1:ndirLett)//'/dat/bs_int_all94.dat'
       zmfileBEG = HOMEBS_dir(1:ndirLett)//'/dat/bs_one_all94.dat'
       zmfileFIN = HOMEBS_dir(1:ndirLett)//'/dat/bs_fin_all94.dat'
c
       unitResLib = 21
       unitResLibB = 22
       unitResLibE = 23
c
       open(unit=unitResLib,file=zmfile,form='formatted',
     &      status='old')
c
       open(unit=unitResLibB,file=zmfileBEG,form='formatted',
     &      status='old')
c
       open(unit=unitResLibE,file=zmfileFIN,form='formatted',
     &      status='old')
c
       if(CONTROL)then
       write(kanalp,*)'getXYZaaSeqFiPsi : read Files: '
       write(kanalp,*)'ZMfile:',zmfile
       write(kanalp,*)'ZMfileBEG:',zmfileBEG
       write(kanalp,*)'ZMfileFIN:',zmfileFIN
       end if !C
c
c add GLY at the end of original aaLoop
       nres = nresLoop+1
       resNameLoop(nres) = 'GLY '
c
c read RESidues and assembe zMatrix for the AASequence
c
       natLoop= 0
       dumFlag = 0
       natPairCyc = 0
       natPairCycMX = natPairCycMAX
       do ir = 1,nres
c
       unitRes = unitResLib   ! INTernal resType
c
       startAtResLoop(ir) = natLoop+1
           
       call readResLibSB(unitRes,resNameLoop(ir),
     &               dumFlag,natLoop,natinres,
     &               molZmatrx1Tmp,molZmatrx2Tmp,molZmatrx3Tmp,
     &               atQLoop,ffatNameLoop,
     &              natPairCycMX,natPairCyc,atPairCycList)
c
c get Fi,Psi into molZmatrx2Tmp() of residue
c       if(ir .eq. 1) then
c       molZmatrx3Tmp(12) = FiPsiLoop(2)  ! at 4 DUMM - standart
c       end if 
c
       ir2 = 2*ir
       do ia = startAtResLoop(ir),startAtResLoop(ir)+natinres
c loop over atoms in res = ir
       ia3=3*ia-3
       if(molZmatrx1Tmp(ia3+1)(1:3) .eq. 'C  ')
     &    molZmatrx3Tmp(ia3+3) = FiPsiLoop(ir2+1)   ! Fi of the RES i
c
       if(molZmatrx1Tmp(ia3+1)(1:3) .eq. 'HA ')
cX     &    molZmatrx3Tmp(ia3+3) = FiPsiLoop(ir2+1) - 120.0  !Daa; Fi of the RES i
     &    molZmatrx3Tmp(ia3+3) = FiPsiLoop(ir2+1) + 120.0   !Laa; Fi of the RES i
c
       if(molZmatrx1Tmp(ia3+1)(1:3) .eq. 'CB ')
cX     &    molZmatrx3Tmp(ia3+3) = FiPsiLoop(ir2+1) + 120.0 !Daa;Fi of the RES i
     &    molZmatrx3Tmp(ia3+3) = FiPsiLoop(ir2+1) - 120.0 !Laa;Fi the RES i
c      
       if(molZmatrx1Tmp(ia3+1)(1:3) .eq. 'O  ')
     &    molZmatrx3Tmp(ia3+3) = FiPsiLoop(ir2+2)+180.0   ! Psi of the RES i
c
       if(molZmatrx1Tmp(ia3+1)(1:3) .eq. 'N  ')
     &    molZmatrx3Tmp(ia3+3) = FiPsiLoop(ir2)   ! Psi of the RES i-1
c
        end do !ia
c
        dumFlag = 1   !delete DUMM atoms
c
       end do! ir
c
c remove all atoms of the C-end CLY, keep N,CA,C and make them DUMM
       natLoop = startAtResLoop(nres)-1 
       do ia = startAtResLoop(nres),startAtResLoop(nres)+natinres
       ia3 = 3*ia-3
       ia4 = 4*ia-4
c
       if(molZmatrx1Tmp(ia*3-2)(1:3) .eq. 'N  ' .or.
     &    molZmatrx1Tmp(ia*3-2)(1:3) .eq. 'CA ' .or.
     &    molZmatrx1Tmp(ia*3-2)(1:3) .eq. 'C  ' ) then
    
       ns = 0 
       if(molZmatrx1Tmp(ia*3-2)(1:3) .eq. 'C  ' ) ns=1
c
       natLoop = natLoop + 1
       nL4 = natLoop*4-4
       nL3 = natLoop*3-3
c
       molZmatrx1Tmp(nL3+1)='DUMM'
       molZmatrx1Tmp(nL3+2)='DU  '
       molZmatrx1Tmp(nL3+3)=molZmatrx1Tmp(ia3+3)
       do k = 1,3
       molZmatrx3Tmp(nL3+k) = molZmatrx3Tmp(ia3+k) 
       end do!k
       molZmatrx2Tmp(nL4+1) = natLoop
       molZmatrx2Tmp(nL4+2) = molZmatrx2Tmp(ia4+2)-ns
       do k = 3,4
       molZmatrx2Tmp(nL4+k) = molZmatrx2Tmp(ia4+k)
       end do !k
c       
       end if  
       end do !ia
c
       close(unitResLib) 
       close(unitResLibE) 
       close(unitResLibB) 
c
        call  zmcoord01(natLoop,
     &    molZmatrx1Tmp,molZmatrx2Tmp,molZmatrx3Tmp,
     &    atLoopXYZ,atLoopCMatrx,atNameLoop)
c
c        
       return
       end
