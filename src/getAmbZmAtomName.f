C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     getAmbZmAtomName()                                                     *
C* corrects the InPut XYZ.pdb names to the amber Zmatrix atom names           *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getAmbZmAtomNamef(pdbAmbZmFile,
     &              natom,atmName,atmNameEx,atmResName,
     &              atmResNameEx,atmResNumb, nres,resNameRes)

C* INput: 
C* PdbAmbZmFile - fullPath/fileName table of CONVERSION PDBname --> AmbZmAtomName
C* natom,
C* nres
C* atmName(ia) - atomNames 
C* atmResName(ia) - resName for atom ia
C* RESULT:
C*    atmName(ia) - corrected     
C*    atmResName(ia) - corrected
C*    resNameRes(ir) - corrected                 
c
	implicit none
	character*(*) pdbAmbZmFile
	integer natom,nres
	character*4 atmName(*)
	character*8 atmNameEx(*)
        character*4 atmResName(*)
        character*8 atmResNameEx(*)
        character*4 resNameRes(*)
        integer atmResNumb(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
	include "ssbond.h"
c local
        integer ntbmax                  !convertionTablsiz.h
        parameter (ntbmax=2047)         !MaxNumb of table lines
        character*4 resnamInR(ntbmax)
        character*4 resnamOutR(ntbmax)
        character*4 resnamInA(ntbmax)
        character*4 atmNameInA(ntbmax)
        character*4 atmNameInAEx(ntbmax)
        character*4 atmNameOutA(ntbmax)
c
        character*80 line
        character*4 resNa4,atmNa4
        character*4 atmNaEx4
        character*4 resNa4Out,atmNa4Out
        logical findR,findA
        integer i,ia,ir
        integer ntbrecR,ntbrecA
        integer fcontr
        integer resNu
        integer kanalp
        integer linematch 
        logical find
	logical replace
        logical CONTROL
c initialize
cx        kanalp = 6
        kanalp = kanalRunOut
        CONTROL = .true. 
c
        if(CONTROL)then
          write(kanalp,*)'In getAmbZmAtomName:'
          write(kanalp,*)'pdbAmbZmFile:',pdbAmbZmFile
        end if
c read in pdbAtName_ambZm.dat
        open(unit=11, file=pdbAmbZmFile, form='formatted',
     &       status= 'old')      
c
        ntbrecR = 0
        ntbrecA = 0
c
100     read(11,'( a80 )', end=101 ) line
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        if(line(1:4) .eq. 'RESI')then
        ntbrecR=ntbrecR+1
c
        if(ntbrecR .gt. ntbmax) then
        write(kanalp, *)'ERROR!in getAmbZmAtomName:(RESI): ',
     &  ' ntbmax (param) is SMALL '
        stop
        end if
c 
        read(line, '(10x,a4,1x,a4)')
     &  resnamInR(ntbrecR),resnamOutR(ntbrecR)
c
        if(CONTROL)then
        write(kanalp,'(10x,a4,1x,a4,a10,i4)')
     &  resnamInR(ntbrecR),resnamOutR(ntbrecR),' ntbrecR:',ntbrecR
        end if
c
        end if ! RESI
c
        if(line(1:4) .eq. 'ATOM')then
        ntbrecA=ntbrecA+1
c
        if(ntbrecA .gt. ntbmax) then
        write(kanalp, *)'ERROR!in getAmbZmAtomName:(ATOM): ',
     &  ' ntbmax (param) is SMALL '
        stop
        end if
c 
        read(line, '(10x,a4,a4,a4,1x,a4)')
     &  atmNameInA(ntbrecA),atmNameInAEx(ntbrecA),
     &  resnamInA(ntbrecA),atmNameOutA(ntbrecA)
c
        if(CONTROL)then
        write(kanalp,'(10x,a4,4x,a4,1x,a4,a10,i4)')
     &  atmNameInA(ntbrecA),resnamInA(ntbrecA),atmNameOutA(ntbrecA),
     &  'ntbrecA:', ntbrecA
        end if
c
        end if ! ATOM
c
        goto 100
 
101     continue
        if(CONTROL)then
        write(kanalp,*)'endFile***********************'
        end if !C
        close(unit=11)
c
c correct RESIDUE NAMES and Heawy atom NAMES:
c
        do i =  1,natom
c
c find RESIDUE:      
c init
          resNa4 = atmResName(i)
          atmNa4 = atmName(i)
          atmNaEx4 = atmNameEx(i)(5:8)
          resNa4Out = resNa4
          atmNa4Out = atmNa4
          resNu = atmResNumb(i)
          findR = .false.
          findA = .false.
c
          do ir = 1,ntbrecR
          if(resNa4 .eq. resnamInR(ir))then
          findR = .true.
          resNa4Out = resnamOutR(ir)
          end if          
          end do!ir   
c
          do ia = 1,ntbrecA
c
          if(atmNa4 .eq. atmNameInA(ia) .and.
     &       atmNaEx4 .eq. atmNameInAEx(ia) .and.
     &       resNa4 .eq. resnamInA(ia) ) then
             atmNa4Out = atmNameOutA(ia)
             findA = .true.
          end if
c
          end do !ia
c change arrays PDB info
          if(findR)then
	  replace = .true.
c CYS residues are replaced if OPT_SSbondAuto=.true.
          if((atmResName(i)(1:3) .eq. "CYS" .or. 
     &        atmResName(i)(1:3) .eq. "CYX" .or.
     &        atmResName(i)(1:3) .eq. "CYM") .and.
     &                .not.OPT_SSbondAuto)replace =.false.
c
          if(replace)then
          atmResName(i) = resNa4Out
          atmResNameEx(i)(1:4) = resNa4Out
          resNameRes(resNu) = resNa4Out 
	  end if!replace
          end if !findR
c
          if(findA)then
          atmName(i) = atmNa4Out
          atmNameEx(i)(1:4) = atmNa4Out 
          end if
c
        if(CONTROL .and. (findR .or. findA ))then
        write(kanalp,*)'Result: of  getAmbZmAtomName:'
        write(kanalp,'(a6,i5,2x,a8,a4,a8,a1,i4)')
     &  'ATOM  ',i,atmNameEx(i),atmResName(i),atmResNameEx(i),
     &  ' ',atmResNumb(i)
        end if!C
        end do!i
c
        if(CONTROL)then
        write(kanalp,*)'End of  getAmbZmAtomName:'
        end if
c
	return
        end
