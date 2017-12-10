c
c Yurii Vorobjev,  2003, 2005
c
c read initial bad(incomplete) MolPDB       
c                             (to insert Loops )      
c define realAtomFlag(ia)
c add missing Protein Side Chain heavyAtoms
c add H
c sort out atomXYZ(*) in the right order = zMatrix from resSeq
c
	subroutine readMolPDB  
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
        include 'optionPar.h'
        include 'addHxyz.h'
c local 
        character*(charLenMAX) fileAllAtPDB
        integer kanalPDB
        integer iHread
        integer i,ia,ir,ja,k,imiss
        integer kanalp
        logical find
        logical CONTROL
        logical CONTROL0 
c replica of array sizes
        natomMX = natomMAX
        nresMX = nresMAX
c
        kanalp = kanalRunOut
        CONTROL0 = .true.
        CONTROL = .false.
c        
cx       pdbfile = './molec.pdb'
c
         write(kanalp,*)'readMolPDB: PDBfile: ',pdbfile
c
         iHread=0
         if(OPT_Hread)iHread=1   ! readPDB with Hydrogens
c
         write(kanalp,*)'in readMolPDB:OPT_Hread:',OPT_Hread
c
        if(iHread .eq. 0)then
c 1) Correct xyz.pdb RESnames and Atom Names to amberZmatrix style
c
        write(kanalp,*)'readMolPDB: start  getAmbZmAtomNamef:'
c
        call getAmbZmAtomNamef(pdbAmbZmFile,
     &           natom_nh,atomName_nh,atomNameEx_nh,
     &           resName_nh,resNameEx_nh,
     &           resNumb_nh,nres_nh,resNameRes_nh) 
c
        write(kanalp,*)'readMolPDB: end of getAmbZmAtomNamef:'
c
c 2) Add XYZ of missing HeavyAtoms(relative MolTopo) to the PDB
        call addHeavyAtom
c
c 3)add Hydrogens 
        write(kanalp,*)'readMolPDB: start of assign_htg:'
c
        do i = 1,natom_nh
         head_nh(i) = 'ATOM  '
        end do !i
c
c3.1) read hInfoFile and create dataStructure assign_htg.h
        call assign_htg(hInfoFile,
     &         natom_nh,atomName_nh,resNameEx_nh,
     &         resNumb_nh,chName_nh)
c
c3.2) using dataStructure assign_htg.h: calculate h-atoms XYZ
c 
        write(kanalp,*)'readMolPDB: start  add_hatom:'
c
         call add_hatom
     &    (natom_nh,head_nh,
     &     atomName_nh,resName_nh,resNumb_nh,chName_nh,atomXYZ_nh,
     &     nres_nh,startAtInRes_nh,stopAtInRes_nh,resNameRes_nh,
     &     natom_h,head_h,
     &     atomName_h,resName_h,resNumb_h,chName_h,atomXYZ_h,
     &     nres_h,startAtInRes_h,stopAtInRes_h,resNameRes_h)
c      
        end if !iHread .eq. 0 
c
c assign Hydr to current pdb_info
c assign atomXYZ(*)  according to the MolTop (from Zmatrix) order
c it is assumed that RESIDUES have the same (right) order in PDB and MolTop
c
         imiss=0
         do ia = 1,natom
         atomNumb(ia) = ia
         chName(ia) = 'A'
c
         realAtomFlag(ia)=0
         ir=resNumb(ia)
         resNameRes(ir) = resName(ia)
c
         find = .false.
c
         if(iHread .eq. 0)then    ! process atomXYZ_h() with added Hydrogens
         if(startAtInRes_h(ir) .ge. 1)then
c
         do ja=startAtInRes_h(ir),stopAtInRes_h(ir)
         if(atomName(ia) .eq. atomName_h(ja))then
         atomXYZ(3*ia-2) = atomXYZ_h(1,ja)  
         atomXYZ(3*ia-1) = atomXYZ_h(2,ja)  
         atomXYZ(3*ia)   = atomXYZ_h(3,ja)
c
         if(CONTROL)then
         write(kanalp,*)'find: ia,ir,ja:',ia,ir,ja,
     &    ' atomXYZ:',atomXYZ_h(1,ja),atomXYZ_h(2,ja),atomXYZ_h(3,ja)
         end if!C
c
         find = .true.
         realAtomFlag(ia)=1
         goto 201 
         end if ! atomName(ia) .eq. (ja)
         end do !ja
         end if !startAtInRes_h(ir) .ge. 1
c
         if(.not. find)then
         imiss = imiss+1
c
         write(kanalp,'(a14,i3,1x,a1,i4,1x,a4,a4,
     &                  i4,1x,i4,1x,i4,1x,a26)')
     &   'WARNING! atom:',imiss,': ',
     &     ia,atomName(ia),resName(ia),ir,
     &     startAtInRes_h(ir),stopAtInRes_h(ir),
     &   ' is MISSING in the PDB !!'
c
         end if !.not. find
c
         end if ! (iHread .eq. 0)
c
         if(iHread .eq. 1)then  ! readPDB with Hydrogens
         if(startAtInRes_nh(ir) .ge. 1)then
c
         do ja=startAtInRes_nh(ir),stopAtInRes_nh(ir)
c
         if(CONTROL)then
         write(kanalp,*)'comp: ia,ir,ja:',ia,ir,ja,
     &  ' aNamePDB:',atomName_nh(ja),' aNameZm:',atomName(ia),
     &  ' atomXYZpdd:',
     &   atomXYZ_nh(1,ja),atomXYZ_nh(2,ja),atomXYZ_nh(3,ja)
         end if!C
c
         if(atomName(ia) .eq. atomName_nh(ja))then
         atomXYZ(3*ia-2) = atomXYZ_nh(1,ja)  
         atomXYZ(3*ia-1) = atomXYZ_nh(2,ja)  
         atomXYZ(3*ia)   = atomXYZ_nh(3,ja)
c
         if(CONTROL)then
         write(kanalp,*)'find: ia,ir,ja:',ia,ir,ja,
     &  ' aNamePDB:',atomName_nh(ja),' aNameZm:',atomName(ia),
     &  ' atomXYZpdd:',
     &   atomXYZ_nh(1,ja),atomXYZ_nh(2,ja),atomXYZ_nh(3,ja)
         end if!C
c
         find = .true.
         realAtomFlag(ia)=1
         goto 201 
         end if ! atomName(ia) .eq. (ja)
         end do !ja
         end if !startAtInRes_nh(ir) .ge. 1
c
cx       startAtInRes(ir) = startAtInRes_h(ir) ! true startAtInRes() is 
cx                            defined by call  getZMatrixAAseq( )
c         
         if(.not. find)then
         imiss = imiss+1
c
         write(kanalp,'(a14,i3,1x,a1,i4,1x,a4,a4,
     &                  i4,1x,i4,1x,i4,1x,a26)')
     &   'WARNING! atom:',imiss,': ',
     &     ia,atomName(ia),resName(ia),ir,
     &     startAtInRes_nh(ir),stopAtInRes_nh(ir),
     &   ' is MISSING in the PDB !!'
c
         end if !.not. find
         end if ! (iHread .eq. 1) 
c
201     continue
c
         end do !ia
c
         if( .not. OPT_AddLoop .and. imiss .gt. 1) then
         write(kanalp,*)
     &   'ERROR!: initialPDB has MISSING atom ',
     &   ' (relative to the molTOPO)'
         stop
c
         end if ! .not. OPT_AddLoop
c
         if(CONTROL0) then
         write(kanalp,'(a40)')
     &   'readMolPDB: reOrdered TOPO: xyzPDBinfo.h'
         write(kanalp,*) 
     &    '      Nat  aNam rNam ch resN   X  Y  Z   realAtFlag ff  Qat'
c
         do ia=1,natom
         write(kanalp,7071)
     &     'ATOM  ',atomNumb(ia),atomName(ia),resName(ia),chName(ia),
     &     resNumb(ia), (atomXYZ(ia*3-3+k),k=1,3),realAtomFlag(ia),
     &     ffAtomName(ia),atomQ(ia)
         end do !ia
c
         end if !Contr0
c  write prepared allAtomXYZ.pdb
         kanalPDB = 11
         fileAllAtPDB='molAllAtXYZ.pdb'
c
         if(nMolNameLett .ge. 1)then
        fileAllAtPDB=
     &  molNameArgLine(1:nMolNameLett)//'.molAllAtXYZ.pdb'
        end if
c
        if(OPT_resFileDir)then
        fileAllAtPDB=resultFileDir(1:nResFileDirLett)//'molAllAtXYZ.pdb'
        if(nMolNameLett .ge. 1)then
        fileAllAtPDB=resultFileDir(1:nResFileDirLett)//
     &  molNameArgLine(1:nMolNameLett)//'.molAllAtXYZ.pdb'
        end if
        end if
c
         open(unit=kanalPDB,file=fileAllAtPDB,form='formatted',
     &         status='unknown')
        write(kanalPDB,'(a40)')
     &        'REMARK: molec.pdb+restoredHvy+Hydr atoms'
c
         do ia=1,natom
         write(kanalPDB,7071)
     &     'ATOM  ',atomNumb(ia),atomName(ia),resName(ia),chName(ia),
     &     resNumb(ia), (atomXYZ(ia*3-3+k),k=1,3),realAtomFlag(ia),
     &     ffAtomName(ia),atomQ(ia)
         end do !ia
c
        write(kanalPDB,'(a6)')'TERCHA'
        write(kanalPDB,'(a6)')'ENDMOL'
        close(kanalPDB)
c
c
7071    format(a6,i5,2x,a4,a4,a1,i4,4x,3f8.3,i2,1x,a2,1x,f8.5) !PDB+(i2,1x,a2)
 	 return
         end
c
