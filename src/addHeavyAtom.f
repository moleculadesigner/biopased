c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c Yurii Vorobjev  2003                                            *
c                                                                 *
c                                                                 *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
c add XYZ of Heavy-atoms to PDB atoms of the PDB file 
c
c PROTEIN sideChain atoms  ONLY!
c         atoms N,CA,C,O ussumed to be existed
c
c  USE: 
c        xyzPDBinfo.h
c        xyzPDBcrd.h
c        zmatrix.h
c        residue MolTopo ZMatrix Lib
c  RESULT:
c        xyzPDBinfo.h
c        xyzPDBcrd.h
c  include XYZ and PDB info for all Heavy Atoms of MOLECULE
c	
        subroutine addHeavyAtom
c
c..................................................................
       include 'xyzPDBsize.h'
c
       include "output.h"
       include "kanalUse.h"
       include "statusMessg_mDyn.h"
c
       include 'xyzPDBinfo.h'
       include 'zmatrix.h'
       include 'zmatrixtmp.h' 
       include 'addHxyz.h'
       include 'optionPar.h' 
       include 'filedat.h'
c local copy of pdbinfo.h
        character*8 atomNameEx0(natomMAX)     !extendedName atm+termTag
        character*4 atomName0(natomMAX)
        character*4 atomBlockName0(natomMAX)  !blockName: P,R,B -NA
        character*4 resName0(natomMAX)
        character*8 resNameEx0(natomMAX)
        character*2 ffAtomName0(natomMAX)
        character*1 chName0(natomMAX) 
        real atomQ00(natomMAX)
        real atomXYZ0(3,natomMAX)
        integer addHvyAtFlag0(natomMAX)
        integer realAtomFlag0(natomMAX)
c getXYZaaSeqFiPsi() parameters
        integer irMAX
        parameter (irMAX = 1)
        integer natresMAXloc
        parameter (natresMAXloc = 100)
c
        integer nresLoop
        integer natLoop
        real FiPsiLoop(2*(irMAX+2))
        character*4 resNameLoop(irMAX+1)
        character*4 atNameLoop(natresMAXloc)
        real atLoopXYZ(3*natresMAXloc)
        character*2 ffatNameLoop(natresMAXloc)
        real atQLoop(natresMAXloc)
        integer startAtResLoop(nresMAX)
        integer atLoopCMatrx(5*natresMAXloc)
c  
        integer nHvyAddRes(nresMAX)
c
        integer nrefAtMAX
        parameter (nrefAtMAX = 4)      
        character*2 refAtName(nrefAtMAX)
        integer refAtNumb0(nrefAtMAX)
        integer refAtNumb1(nrefAtMAX)
        real refAtXYZ0(3*nrefAtMAX)
        real refAtXYZ1(3*nrefAtMAX)
        real atHvyXYZ0(3*natresMAXloc)
        real atHvyXYZ1(3*natresMAXloc)
        character*4 atNameHvy(natresMAXloc)
        integer nHvyAt
        integer nHvyAddTot
        integer nrefAt0,nrefAt1
c
        logical addHvyDO
        integer kanalPDB
        character*(charLenMAX) fileHvyAt
c
        integer natom0
        logical CONTROL1
        logical CONTROL0
        integer resPosit
        integer inCha,Nend,Cend
        integer nDUM
        logical delDUM
        integer i,j,k,ia
        integer ir,ich
        integer i1,i3,ia3
        integer ias,iaf,ia0
        integer nrf,nrf3
	integer kanalp
        data refAtName/'N ','CA','C ','O '/
c
c initialize
        CONTROL1= .false.
        CONTROL0 = .true.
cx        kanalp = 6
        kanalp = kanalRunOut 
c
        resPosit = 0  ! insideChain=0, Nend=1, Cend=2
        inCha=0
        Nend=1
        Cend=2
c local copy of initial xyzPDBinfo.h
        natom0 = natom_nh
c
        do i = 1,natom0
        addHvyAtFlag0(i) = addHvyAtFlag(i)
        atomNameEx0(i) = atomNameEx_nh(i)
        atomName0(i) = atomName_nh(i) 
        resName0(i) = resName_nh(i)
        resNameEx0(i) = resNameEx_nh(i)
        realAtomFlag0(i) = realAtomFlag_nh(i)
        chName0(i) = chName_nh(i)
c
        atomXYZ0(1,i)=atomXYZ_nh(1,i)
        atomXYZ0(2,i)=atomXYZ_nh(2,i)
        atomXYZ0(3,i)=atomXYZ_nh(3,i)
        end do !i
c      
c ADD HEAVY atoms start:!
        nHvyAddTot = 0
cx        do ir = 1,nres
         do ich = 1,nChain
         do ir = startResInCha(ich),stopResInCha(ich)
c resPosition
         resPosit=inCha
         if(ir .eq. startResInCha(ich))resPosit=Nend
         if(ir .eq. stopResInCha(ich))resPosit=Cend
c check residue ir is complete?
        nHvyAddRes(ir) = 0
        addHvyDO = .false.
        if(atHvyNbInRes(ir) .ne. atHvyNbInResZm(ir))then  !if1
c add SideChainAtoms
         addHvyDO = .true.
        nHvyAddRes(ir) = atHvyNbInResZm(ir)-atHvyNbInRes(ir)
c         
        if(CONTROL0)then
        write(kanalp,*)'WARNING! addHeavyAtom: Miss heavy At in resN:'
     &  ,ir,' atPDB, atZMlib:', atHvyNbInRes(ir), atHvyNbInResZm(ir)
        end if !C
c
        end if !if1
c
c 1) define INdata for getXYZaaSeqFiPsi()
c    calculate atXYZ from Zm template
c
        if(addHvyDO)then   !if2
        nresLoop = 1
        resNameLoop(1) = resNameRes(ir)
c
        FiPsiLoop(1) = 0.0      ! DUMM
        FiPsiLoop(2) = 180.0    ! Psi (-1) res defines DUMM at
        FiPsiLoop(3) = 180.0    ! Fi(ir=1)
        FiPsiLoop(4) = 180.0    ! Psi(ir=1)
        FiPsiLoop(5) = 180.0    ! term GLY   
        FiPsiLoop(6) = 180.0    ! term GLY
c
        call getXYZaaSeqFiPsi(nresLoop,resNameLoop,
     &              FiPsiLoop,natLoop,atLoopXYZ,
     &              atQLoop,atNameLoop,ffatNameLoop,
     &              startAtResLoop,atLoopCMatrx)
c
c define reference atoms
        nrefAt0 = 0
        nrefAt1 = 0
        do j = 1,nrefAtMAX
c
        do ia = 1,natLoop
        if(atNameLoop(ia)(1:2) .eq. refAtName(j))then
        nrefAt0 = nrefAt0+1
        refAtNumb0(nrefAt0) = ia
        i3 = 3*nrefAt0-3
        ia3 = 3*ia-3
        refAtXYZ0(i3+1)=atLoopXYZ(ia3+1)
        refAtXYZ0(i3+2)=atLoopXYZ(ia3+2)
        refAtXYZ0(i3+3)=atLoopXYZ(ia3+3)
        end if
        end do!ia
c
c define refAtoms on Protein
        do ia = startAtInRes_nh(ir),stopAtInRes_nh(ir)
c
        if(atomName0(ia)(1:2) .eq. refAtName(j))then
        nrefAt1 = nrefAt1+1
        refAtNumb1(nrefAt1) = ia
        i3 = 3*nrefAt1-3
        ia3 = 3*ia-3
        refAtXYZ1(i3+1)=atomXYZ0(1,ia)
        refAtXYZ1(i3+2)=atomXYZ0(2,ia)
        refAtXYZ1(i3+3)=atomXYZ0(3,ia)
        end if ! atomName_nh(ia) .eq.
c
        if(CONTROL1)then
        write(kanalp,*)
     &  'addHeavyAtom: j,refAtName(j),ia,atomName0(ia):',
     &   j,refAtName(j),ia,atomName0(ia)
        write(kanalp,*)'nrefAt1 : ',nrefAt1
        end if
c
        end do!ia
c
        end do!j 
c
        if(CONTROL0)then
        write(kanalp,*)
     &  'addHeavyAtom: res,startAtInRes_nh(ir),stopAtInRes_nh(ir):',
     &   ir, startAtInRes_nh(ir),stopAtInRes_nh(ir)
        end if !C
c
        if(nrefAt1 .ne. nrefAtMAX)then
        write(kanalp,*)'WARNING! Missing BackBone Atoms ! ',
     &  ' nrefAt0,nrefAt1,nrefAtMAX:',nrefAt0,nrefAt1,nrefAtMAX
        write(kanalp,*)' refAtNumb1():',(refAtNumb1(k),k=1,nrefAt1)
        write(kanalp,'(a12,6a4)')
     &  ' refAtNames:', (atomName_nh(refAtNumb1(k)),k=1,nrefAt1)
        write(kanalp,'(a25,6a4)')
     &  'refAtNames(EXpect):', 
     &   (atomName_nh(refAtNumb1(k)),k=1,nrefAtMAX)
        write(kanalp,*)'FATAL ERROR: program ABorted!!'
        write(kanalp,*)
     &  'WARNING!check atoms with EQUAL names in the RES in INit.pdb!'
c
          write(kanalPStat,*)mError,
     &    ' (-c file) has ATOMs missing in the molTopLIB:bs*.dat files',
     &    ' or (-c file) has ATOMs with EQUAL names in the same RES '

         write(kanalPStat,*)mError,
     &    ' see (-o file runOut ) for details '
c
        stop
        end if
c
        if(CONTROL0)then
        write(kanalp,*)'addHeavyAtom: refAt0:',refAtNumb0 
        write(kanalp,*)'addHeavyAtom: refAtXYZ0:'
        write(kanalp,'(4(3f7.2),1x)')refAtXYZ0
c
        write(kanalp,*)'addHeavyAtom: refAt1:',refAtNumb1 
        write(kanalp,*)'addHeavyAtom: refAtXYZ1:'
        write(kanalp,'(4(3f7.2),1x)')refAtXYZ1
c
        end if!C
c
c define subset of Heavy Atoms of the atLoopXYZ()
c in local XYZ0 system
        nHvyAt = 0
        nDUM = 0
        do ia = 1,natLoop
c
        if(resPosit .eq. inCha .or. resPosit .eq. Nend)then
        if(atNameLoop(ia)(1:2) .ne. 'DU' 
     &     .and. atNameLoop(ia)(1:1) .ne. 'H')then
        nHvyAt = nHvyAt + 1
        i3 = 3*nHvyAt - 3
        ia3 = 3*ia-3
        atHvyXYZ0(i3+1) = atLoopXYZ(ia3+1) 
        atHvyXYZ0(i3+2) = atLoopXYZ(ia3+2) 
        atHvyXYZ0(i3+3) = atLoopXYZ(ia3+3) 
        atNameHvy(nHvyAt) = atNameLoop(ia)
        end if !hvyAt
        end if !resPos .eq. inCha
c
        if(resPosit .eq. Cend)then
        delDUM = .true. 
c
        if(atNameLoop(ia)(1:2) .eq. 'DU')then
        nDUM=nDUM+1
        if(nDUM .eq. 4)delDUM = .false.  ! it is OXT terminal
        end if !atNameLoop(ia)(1:2) .eq. 'DU'
c
        if( .not. delDUM  .or.  (atNameLoop(ia)(1:2) .ne. 'DU'
     &     .and. atNameLoop(ia)(1:1) .ne. 'H') )then  ! hvyAt
        nHvyAt = nHvyAt + 1
        i3 = 3*nHvyAt - 3
        ia3 = 3*ia-3
        atHvyXYZ0(i3+1) = atLoopXYZ(ia3+1) 
        atHvyXYZ0(i3+2) = atLoopXYZ(ia3+2) 
        atHvyXYZ0(i3+3) = atLoopXYZ(ia3+3) 
        atNameHvy(nHvyAt) = atNameLoop(ia)
        if(.not. delDUM) atNameHvy(nHvyAt) = "OXT"
c
        end if !hvyAt
        end if !resPos .eq. Cend 
c
        end do !ia
c
c put the atHvyXYZ0() on the refAtXYZ1() positions
        call rotMolBy3P(refAtXYZ0(1),refAtXYZ0(4),
     &                  refAtXYZ0(7),refAtXYZ0(4),
     &                  nHvyAt,atHvyXYZ0,
     &                  refAtXYZ1(1),refAtXYZ1(4),
     &                  refAtXYZ1(7),refAtXYZ1(4),
     &                  atHvyXYZ1)
c
        end if !if2
c
c insert into _nh arrays new Heavy atom
        ias = startAtInRes_nh(ir)+nHvyAddTot
        iaf = stopAtInRes_nh(ir)+nHvyAddTot+nHvyAddRes(ir)
c      
        do ia = ias,iaf
        i1 = ia-ias+1
        i3 = 3*i1-3
        ia3=3*ia-3
        ia0 = ia-nHvyAddTot
c
        if(addHvyDO)then  !if3
        atomName_nh(ia) = atNameHvy(i1)
        atomNameEx_nh(ia) = 
     &   atomName_nh(ia)//resNameEx0(startAtInRes_nh(ir))(5:8)
        resNameEx_nh(ia) = 
     &   resNameRes_nh(ir)//resNameEx0(startAtInRes_nh(ir))(5:8)
        resName_nh(ia) = resNameRes_nh(ir)
c
        nrf = 0
        if(atNameHvy(i1)(1:3) .eq. 'N  ' )nrf=1
        if(atNameHvy(i1)(1:3) .eq. 'CA ' )nrf=2
        if(atNameHvy(i1)(1:3) .eq. 'C  ' )nrf=3
        if(atNameHvy(i1)(1:3) .eq. 'O  ' )nrf=4
         if(nrf .gt. 0)then !mainChain
c save original backbone atoms
        nrf3=3*nrf-3
        atomXYZ_nh(1,ia) = refAtXYZ1(nrf3+1)
        atomXYZ_nh(2,ia) = refAtXYZ1(nrf3+2)
        atomXYZ_nh(3,ia) = refAtXYZ1(nrf3+3)
         else
c add calculated sidechain atoms
        atomXYZ_nh(1,ia) = atHvyXYZ1(i3+1)
        atomXYZ_nh(2,ia) = atHvyXYZ1(i3+2)
        atomXYZ_nh(3,ia) = atHvyXYZ1(i3+3)
        end if !mainChain
c
        else      ! NoAddHvy, keep OldAtoms
        atomName_nh(ia) = atomName0(ia0)
        atomXYZ_nh(1,ia) = atomXYZ0(1,ia0)        
        atomXYZ_nh(2,ia) = atomXYZ0(2,ia0)        
        atomXYZ_nh(3,ia) = atomXYZ0(3,ia0)        
        atomNameEx_nh(ia) = atomNameEx0(ia0)
        resNameEx_nh(ia) = resNameEx0(ia0)
        end if !if3
c
        resName_nh(ia) = resNameRes_nh(ir)
        resNumb_nh(ia) = ir
        chName_nh(ia) = chNameRes_nh(ir)
        realAtomFlag_nh(ia) = 1
c
        end do!ia
c
c correct start/stop
        startAtInRes_nh(ir) = startAtInRes_nh(ir)+nHvyAddTot
        nHvyAddTot = nHvyAddTot + nHvyAddRes(ir)
        stopAtInRes_nh(ir) = stopAtInRes_nh(ir) + nHvyAddTot
c
         end do !ir loop
         end do !ich
c
        natom_nh = natom0 + nHvyAddTot
c
        write(kanalp,*)' Finish addHeavyAtom !'
c
        kanalPDB = 11
        fileHvyAt = 'molAddHvyAt.pdb'
c
        if(nMolNameLett .ge. 1)then
        fileHvyAt=
     &  molNameArgLine(1:nMolNameLett)//'.molAddHvyAt.pdb'
        end if
c
        if(OPT_resFileDir)then
        fileHvyAt=resultFileDir(1:nResFileDirLett)//'molAddHvyAt.pdb'
        if(nMolNameLett .ge. 1)then
        fileHvyAt=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.molAddHvyAt.pdb'
        end if
        end if
c
        open(unit=kanalPDB,file=fileHvyAt,form='formatted',
     &       status='unknown')
        write(kanalPDB,'(a36)')'REMARK: molec.pdb + added HeavyAtoms'
        do ia = 1,natom_nh
        write(kanalPDB,7071)'ATOM  ',ia,atomName_nh(ia),resName_nh(ia),
     &  chName_nh(ia),resNumb_nh(ia),(atomXYZ_nh(k,ia),k=1,3)
        end do!ia
c
        write(kanalPDB,'(a6)')'TERCHA'
        write(kanalPDB,'(a6)')'ENDMOL' 
c
        close(kanalPDB)
c
7071    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB  
c
	return
	end
