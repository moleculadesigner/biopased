c READ in PDBXYZ 
c PDB file must be organized according PDB format/rule
c molecule/chain atoms terninated by TER line
c SPECIAL:
c water HOH residue or ION residue doNOT require TER line
c
c atmXYZ() real*4
c      read BAD PDB, some residues (atoms) are missing, no XYZ
c
        subroutine readPDB05(pdbfile,natomMAX,nresMAX,iHread,
     &       natom,atomNumb,atomName,resName,chName,resNumb,
     &       atmXYZ,nres,resNameRes,chNameRes,
     &       atomNameEx,resNameEx,
     &       startAtInRes,stopAtInRes,realAtomFlag,
     &       resTypeREseq,atHvyNbInRes,
     &       nChain,startAtInCha,stopAtInCha,
     &       startResInCha,stopResInCha,
     &       nMolec,startAtInMol,stopAtInMol,
     &       startResInMol,stopResInMol)
c
c iHread - flag 0/1 NO/YES to read Hydrogen from initial XYZ.pdb
c
c atomNumb(ia) - original atom N in PDBfile for ia in the sequential list
c                1,2,..,
c natom - MAX number of atoms, ATOM number as it is in the PDB file
c nres  - MAX number of RESiduesNUMB line 
c
        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        character*(*) pdbfile
        integer natomMAX,natom
        integer iHread
        character*4 atomName(*)
        character*4 resName(*)
        character*1 chName(*)
        real*4 atmXYZ(*)
        integer atomNumb(*),resNumb(*)
        integer nres,nresMAX
        character*4 resNameRes(*)
        character*1 chNameRes(*) 
        character*8 atomNameEx(*)
        character*8 resNameEx(*)
        integer startAtInRes(*)
        integer stopAtInRes(*)
        integer realAtomFlag(*)
        character*9 resTypeREseq(*)
        integer atHvyNbInRes(*)
        integer nChain
        integer startAtInCha(*)
        integer stopAtInCha(*)
        integer nMolec
        integer startAtInMol(*)
        integer stopAtInMol(*)
        integer startResInCha(*)
        integer stopResInCha(*)
        integer startResInMol(*)
        integer stopResInMol(*)
c local variables
        integer nAtTagMAX
        parameter (nAtTagMAX = 2)
        character*1 atomTagList(nAtTagMAX)
        integer itg
        character*6 lastHead
c
        character*6 head
        character*6 TERCHA,ENDMOL
        character*6 resnGen,resnGenM1
        character*1 pdbVarA,pdbVarS
        character*1 pdbVarB,pdbVar
        character*4 atNameLn
        character*6 STARTM,ENDMDL
        logical acceptA,acceptS
        logical watIon,lastWatIon
c
        integer i,j,k,k3  
        integer nres1
        character*80 line
        integer kanalXYZ,kanalp
        logical readAtom,newChain
        logical CONTROL,CONTROL0
c
        kanalXYZ=kanalInPdb
        kanalp=kanalRunOut
c
        CONTROL0 = .true.
        CONTROL  = .false. 
c
        write(kanalp,*)'In readPDB05:pdbfile:',pdbfile 
c init
        TERCHA = "TERCHA"
        ENDMOL = "ENDMOL"
c
        STARTM = "MODEL "
        ENDMDL = "ENDMDL"         ! end of MODEL - multiple conformation pdb file
c
        pdbVarA = 'A'
        pdbVarB = 'B'
        pdbVarS = ' '
        resnGen = '      '
        acceptA = .false.
        acceptS = .false.
        atNameLn = '    '
        resnGenM1 = resnGen  !yv
c
        atomTagList(1) = pdbVarA
        atomTagList(2) = pdbVarB       
        lastHead = '      '
        lastWatIon = .false.
c init PDB
        do k = 1,natomMAX
        resName(k) = '    '
        atomName(k) = '    '
        atomNameEx(k) = '        '
        chName(k) = ' '
        resNumb(k) = 0
        atomNumb(k) = 0    !if =0, atomXYZ does not exist
        atmXYZ(3*k-2) = 0.0
        atmXYZ(3*k-1) = 0.0
        atmXYZ(3*k) = 0.0
        realAtomFlag(k) = 0
        end do!k
        do k=1,nresMAX
        startAtInRes(k)=0
        stopAtInRes(k)=0
        startAtInCha(k)=0
        stopAtInCha(k)=0
        startAtInMol(k)=0 
        stopAtInMol(k)=0
        atHvyNbInRes(k)=0
        end do
c initialization
        nres = 0
        nres1 = 0
        nChain = 0
        nMolec = 0
c
       open(unit=kanalXYZ,file=pdbfile,form='formatted',
     &      status='old' )

c READ PDB file and assign to variables
           k=0
c read XYZ
          rewind kanalXYZ

400       read(kanalXYZ,'(a80)',end=401)line
c
          write(kanalp,*)line
c
           readAtom=.false.
           watIon = .false.
           newChain = .false.
c
           if(line(1:4).eq.'ATOM' )then   !ATOM
c
           readAtom=.true. 
c
           pdbVar = line(17:17)
           resnGen = line(22:27)    !ChaTag+ResNumb
c
 
          if(iHread .eq. 0 .and. line(14:14) .eq. 'H')readAtom=.false.
c
           if(resnGen .eq. resnGenM1)then
c the same Cha+ResNumb
            if(pdbVar .eq. pdbVarB )readAtom=.false. 
c
           if(atNameLn .ne. line(13:16))then
           acceptA = .false.
           acceptS = .false.
           else 
           if(acceptA .or. acceptS)readAtom=.false.
           end if !newAtName
c
           if((.not. acceptA) .and. (.not. acceptS))atNameLn=line(13:16)
c
           else
           readAtom=.true. 
           end if !
c
           if(CONTROL)then
           write(kanalp,'(a68,a12,a2,L1)')
     &     line(1:68),' OrigPDBline','R=',readAtom
           end if
           end if !ATOM
c
           if(readAtom)then
           k=k+1
           k3=3*k-3
c
           resnGen = line(22:27)
c
           if(pdbVar .eq. pdbVarA)acceptA = .true.
           if(pdbVar .eq. pdbVarS)acceptS = .true.
c
c control
           if(k.gt.natomMAX)then
           write(kanalp,*)'ERROR!: Too much atoms in PDB file'
           write(kanalp,*)'       Allowed atomMAX = ',natomMAX
           write(kanalp,*)'       increase parameter (atomMAX)...'
c
           write(kanalPStat,*)mError,
     &    ' (-c file) has to much atoms, increase parameter (atomMAX)',
     &    ' in xyzPDBsize.h file '
c
           stop
           end if
  
           chName(k) = ' ' 
           read(line,7071)head,
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atmXYZ(k3+j),j=1,3)
           if(chName(k) .eq. ' ')chName(k)='A'
           realAtomFlag(k) = 1
ccontrol
           if(CONTROL)then
           write(kanalp,7071)head,
     &     atomNumb(k),atomName(k),resName(k),chName(k),resNumb(k),
     &     (atmXYZ(k3+j),j=1,3)
           end if !CONTROL

           if(k .eq. 1)then
           nres1 = resNumb(k)             ! number of res1 in PDB
cx          nres = resNumb(k) - nres1+1    ! nres = 1 for katom=1
           nres = 1
           resNameRes(nres) = resName(k)
           chNameRes(nres) = chName(k)
           startAtInRes(nres) = 1
c
           if(atomName(k)(1:1) .ne. 'H')then
           atHvyNbInRes(nres) = atHvyNbInRes(nres) + 1
           end if! noH
c
           resTypeREseq(nres) = 'BEG'
           startAtInCha(1) = 1
           startAtInMol(1) = 1
           startResInCha(1) = 1
           startResInMol(1) = 1
           end if !k=1
c
           if( k .gt. 1 ) then  

           if(resnGen .ne. resnGenM1)then  !next REs starts
           nres = nres + 1
           resNameRes(nres) = resName(k)
           chNameRes(nres) = chName(k)
           startAtInRes(nres) = k
           if(nres .gt. 1)stopAtInRes(nres-1) = k-1
c
           if(resTypeREseq(nres-1)(1:3) .eq. 'FIN')then
           resTypeREseq(nres) = 'BEG'
c
           else
           resTypeREseq(nres) = 'INT'
           end if !resTypeREseq(nres-1) .eq. 'FIN'
c
           end if !nextRes
c
           if(CONTROL)then
           write(kanalp,*)'pdbRead05:resnGen resnGenM1,nres:',
     &     resnGen,resnGenM1,nres
           end if! CONTROL
c
           if(atomName(k)(1:1) .ne. 'H')then
           atHvyNbInRes(nres) = atHvyNbInRes(nres) + 1
           end if! noH
c
           end if !k > 1
c
           resNumb(k) = nres
c          
           lastHead = 'ATOM  '
           resnGenM1 = resnGen  !eyv
c
           if((resName(k)(1:3) .eq. 'HOH'
     &       .and. atomName(k) .eq. 'O   ') .or.
     &          (resName(k)(1:3) .eq. 'ION'))then 
           watIon = .true.
           lastWatIon = .true.
           end if !'HOH'
c
           end if !readatom 
cedyv15
cx           write(kanalp,*)
cx     &     'pdbRead05:lastHead,lastWatIon,watIon,readAtom,newChain:',
cx     &      lastHead,lastWatIon,watIon,readAtom,newChain
c
         if((line(1:3).eq. TERCHA(1:3) .and. (.not. lastWatIon))
     &        .or. watIon ) then !CHAIN/MOLec Termination
          newChain = .true.
          lastHead = 'TERCHA'
          end if
       if(line(1:3).eq. TERCHA(1:3) .and. (.not. lastWatIon))then
          lastWatIon = .false. 
       end if !TERCHA
c
         if(newChain)then
c TERCHA = "TERCHA"
           nChain = nChain+1
           resTypeREseq(nres) = 'FIN'
c
           if(nres .gt. 1)then
           if(resTypeREseq(nres-1)(1:3) .eq. 'FIN')
     &     resTypeREseq(nres) = 'BEG'
c
           if(nChain .gt. 1 .and.
     &     (stopResInCha(nChain-1) .eq. startResInCha(nChain-1)))
     &     resTypeREseq(nres) = 'BEG' 
c
           else
           resTypeREseq(nres) = 'BEG'
           end if !nres .gt. 1
c
           stopAtInCha(nChain) = k
           startAtInCha(nChain+1) = k+1
           stopResInCha(nChain) = nres
           startResInCha(nChain+1) = nres+1
c
           if(nChain .gt. 1)then
           if(stopResInCha(nChain) .eq. startResInCha(nChain))
     &     resTypeREseq(nres) = 'BEG'
           end if ! nChain .gt. 1
c
           newChain = .false.
c
           end if ! newChain 
c          
           if((line(1:3) .eq. ENDMOL(1:3)) .or.
     &        (line(1:6) .eq. ENDMOL(1:6)) .or.
     &        (line(1:6) .eq. ENDMDL ))then !Molec Termination
c ENDMOL = "ENDMOL"
           nMolec = nMolec+1
           stopAtInMol(nMolec) = k
           startAtInMol(nMolec+1) = k+1
           stopResInMol(nMolec) = nres
           startResInMol(nMolec+1) = nres+1
c
c correct nChain number if no TER line before END
c
           if(lastHead .eq. 'ATOM  '.and.(.not. lastWatIon))then
           nChain = nChain + 1
           stopAtInCha(nChain) = k
           startAtInCha(nChain+1) = k+1
           stopResInCha(nChain) = nres
           startResInCha(nChain+1) = nres+1            
           end if ! lastHead .eq. 'ATOM  '
c
           end if ! 'ENDMOL'
c
c end of MODEL
           if((line(1:3) .eq. ENDMOL(1:3)) .or.
     &        (line(1:6) .eq. ENDMDL ))then
            goto 401
            end if !endofMODELXYZ
c
           goto 400

c end of file is reached:
401        natom=k
           startAtInRes(nres+1) = natom+1
           stopAtInRes(nres) = natom
c
           startAtInCha(nChain+1) = natom+1
           stopAtInCha(nChain) = natom
           stopAtInMol(nMolec+1) = natom+1
           stopAtInMol(nMolec) = natom
c
           if(natom .eq. 0)then
           write(kanalp,*)'ERROR!! molec.pdb does not has ATOM records'
c
           write(kanalPStat,*)mError,
     &    ' (-c file) has no  ATOM records .!'
           stop
           end if
c correct startAtInCha(*)
           startAtInCha(1)=1
           do i=1,nChain
           stopAtInCha(i)=stopAtInRes(stopResInCha(i))
           startAtInCha(i+1)=stopAtInCha(i)+1
           end do !ich
c
           do k=1,natom
c correct PDB atomNames, delete pdbVarTag
           do itg = 1,nAtTagMAX
            if(atomName(k)(4:4) .eq. atomTagList(itg))then
            atomName(k)(4:4) = ' '
            end if
            end do !itg
c
           atomNameEx(k)=atomName(k)//'    ' 
           end do!k
c
          if(nres .gt. nresMAX)then
          write(kanalp,*)'ERROR!:pdbRead: nresMAX is low !'
c
           write(kanalPStat,*)mError,
     &    ' (-c file) nresMAX is low! , increase param nresMAX !'
c
          stop
          end if
c
          if(CONTROL)then
          write(kanalp,*)'pdbRead05:  nres1,nres:',nres1,nres
          do i=1,nres
          write(kanalp,*)
     &     'nres:',i,' startAtInRes(i):',startAtInRes(i),
     &     ' stopAtInRes(i):',stopAtInRes(i),
     &     ' resTypeREseq:', resTypeREseq(i)
          end do
c
          write(kanalp,*)'pdbRead05:  startAtInCha(i): '
          do i=1,nChain
          write(kanalp,*)
     &     'nCha:',i,' startAtInCha(i):',startAtInCha(i),
     &     ' stopAtInCha(i):',stopAtInCha(i),
     &     ' start/stopRes:', startResInCha(i),stopResInCha(i)
          end do !i
c
          write(kanalp,*)'pdbRead05: startAtInMol(i): '
          do i=1,nMolec
          write(kanalp,*)
     &     'nMol:',i,' startAtInMol(i):',startAtInMol(i),
     &     ' stopAtInMol(i):',stopAtInMol(i),
     &     ' start/stopRes:', startResInMol(i),stopResInMol(i)
          end do !i
          end if !Control
c
c make extended atom Names
           do i = 1,nres
c
           do k=startAtInRes(i),stopAtInRes(i)
           if(resTypeREseq(i)(1:3) .eq. 'BEG')then
           atomNameEx(k)=atomName(k)//'NT  '
           end if
           if(resTypeREseq(i)(1:3) .eq. 'FIN')then
           atomNameEx(k)=atomName(k)//'CT  '
           end if
           end do !k
c
           end do !i nres
c resNameEx()
           do k = 1,natom
           resNameEx(k)(1:4) = resName(k)
           resNameEx(k)(5:8) = atomNameEx(k)(5:8)
           end do !k 
c contorl
           if(CONTROL0)then 
           write(kanalp,*)'PDBread05:pdbXYZInPut: Natom:',natom
           write(kanalp,*)'pdbLine,atomNumb,atomName,resName,etc.'
           head = 'ATOM  '
           do k=1,natom
           k3=3*k-3
c
           write(kanalp,7072)
     &     k, head,atomNumb(k),atomName(k),resName(k),chName(k),
     &     resNumb(k), (atmXYZ(k3+j),j=1,3),atomNameEx(k),resNameEx(k)
           end do !k
c
         write(kanalp,*)
     &  'pdbRead05:startAtInRes(k) stopAtInRes(k) resType atHvyNbInRes'
           do k=1,nres
           if(startAtInRes(k) .gt. 0)then
           write(kanalp,'(i6,a6,a4,i5,1x,i5,1x,a3,i5)')
     &     k,resNameRes(k),chNameRes(k),startAtInRes(k),stopAtInRes(k),
     &     resTypeREseq(k),atHvyNbInRes(k)
           else
           write(kanalp,'(i6,a6,a4,i5,1x,i5,1xa3)')
     &     k,'      ','    ',startAtInRes(k),stopAtInRes(k),
     &     resTypeREseq(k)
           end if!startAtInRes(k) .gt. 0 !nonLoop 
           end do!k
c
          write(kanalp,*)'pdbRead05:  startAtInCha(i): '
          do i=1,nChain
          write(kanalp,*)
     &     'nCha:',i,' startAtInCha(i):',startAtInCha(i),
     &     ' stopAtInCha(i):',stopAtInCha(i),
     &     ' start/stopRes:', startResInCha(i),stopResInCha(i)
          end do !nCha
c
          write(kanalp,*)'pdbRead05: startAtInMol(i): '
          do i=1,nMolec
          write(kanalp,*)
     &     'nMol:',i,' startAtInMol(i):',startAtInMol(i),
     &     ' stopAtInMol(i):',stopAtInMol(i),
     &     ' start/stopRes:', startResInMol(i),stopResInMol(i)
          end do !nMol
c
           write(kanalp,*)'END of pdbRead05'
           end if !control
c
          close(kanalXYZ)
c
7071    format(a6,i5,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB
7072    format(i5,1x,a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,2x,a8,1x,a8) !(i5,1x)+ PDB
c
            return
            end
