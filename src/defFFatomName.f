C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     defFFatomName()                                                        *
C* Defines  the ForceField atom names from PDB atom name                      *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine defFFatomName (ffAtomTypeFile,
     &              natom,atmNameEx,atmResName,chname,
     &              ffAtomName,atomCharge)
 
C* ffAtmTypeFile - fullPath/fileName table of CONVERSION PDBname --> ffAtomName
C* natom
C* atmNameEx(ia) - atomNames Extended
C* atmResName(ia) - resName for atom ia
C* chname(ia) - chainName
C* RESULT:
C*    ffAtomName(ia) - ForceFieled specified atom name
C*    atomCharge(ia) - ForceFieled specified atomic charge
c
	implicit none
	character*(*) ffAtomTypeFile
	integer natom
	character*8 atmNameEx(*)
        character*4 atmResName(*)
        character*1 chname(*)
        character*2 ffAtomName(*)
        real atomCharge(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer ntbmax                  !convertionTablsiz.h
        parameter (ntbmax=2047)         !MaxNumb of table lines 
        character*4 atnam_tb(ntbmax)      
        character*4 term_tb(ntbmax)
        character*4 res_tb(ntbmax)
        character*4 rnum_tb(ntbmax)
        character*1 chnam_tb(ntbmax)
        character*2 atomFFname_tb(ntbmax)
        real    atomFFq_tb(ntbmax)
c hashTable local
        integer nhashmx
        parameter (nhashmx = ntbmax )
        integer ntabl_ihash(nhashmx)
        integer hashtblink(nhashmx)
c
        character*80 line
        character*4 rnum4,atName4
        character*8 resn8
        integer i
        integer ntbrec,fcontr
        integer kanalp
        integer linematch 
        logical find
        logical CONTROL,CONTROL0
        logical OPT_hash

c initialize
cx        kanalp = 6
        kanalp = kanalRunOut
c
        CONTROL = .false.
        CONTROL0 = .true.
        OPT_hash = .true.

c read in pdbAtom_ambff.tab 
        open(unit=11, file=ffAtomTypeFile, form='formatted',
     &       status= 'old')      
c
        ntbrec = 0
        if(CONTROL0)then
          write(kanalp,*)'In defFFatomName:'
          write(kanalp,*)'ffAtomTypeFile : ',ffAtomTypeFile
        end if
c
100     read(11,'( a80 )', end=101 ) line
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        ntbrec=ntbrec+1
        if(ntbrec .gt. ntbmax) then
        write(kanalp, *)'ERROR!in ffAtomTypeFile: ',
     &  ' ntbmax (param) is SMALL '
        stop
        end if
 
        read(line, '(a4,2x,a4,a4,1x,a1,2x,a2,2x,f8.5)')
     &  atnam_tb(ntbrec),term_tb(ntbrec),res_tb(ntbrec),
     &  chnam_tb(ntbrec),atomFFname_tb(ntbrec),atomFFq_tb(ntbrec) 
c
        rnum_tb(ntbrec)='    '
 
        if(CONTROL)then
        write(kanalp,'(a4,2x,a4,a4,1x,a1,2x,a2,2x,f8.5,1x,i6)')
     &  atnam_tb(ntbrec),term_tb(ntbrec),res_tb(ntbrec),
     &  chnam_tb(ntbrec),atomFFname_tb(ntbrec),atomFFq_tb(ntbrec),
     &  ntbrec
        end if
 
        goto 100
 
101     continue
        close(unit=11)
c end of readIn ffAtomTypeFile
c
c calculate hashTable for the dictionaryRecords
         fcontr = 0      !test print
c 
         call hash_dict_tabl
     &    (ntbrec,atnam_tb,term_tb,res_tb,rnum_tb,chnam_tb,
     &     ntabl_ihash,hashtblink,nhashmx,fcontr)
c
c define ffAtomName
          if(CONTROL0)then
          write(kanalp,*)'defFFatomName: ambNam  Q'
          write(kanalp,*)'atmNameEx resnEx ffAtomName atomCharge'
          end if
c
          fcontr = 0  !test print
          do i=1,natom
          atName4 = atmNameEx(i)(1:4)
          resn8 = atmResName(i)//atmNameEx(i)(5:8)    
          rnum4 = '7777'   ! no right resNumb
c
          if(.not.OPT_hash)then
          call find_match_word(atName4,resn8,rnum4,
     &    chname(i),ntbrec,atnam_tb,term_tb,res_tb,rnum_tb,
     &    chnam_tb,linematch,fcontr,find)
 
          else
 
          call find_match_hashtb2(atName4,resn8,rnum4,
     &        chname(i),ntbrec,atnam_tb,term_tb,res_tb,rnum_tb,
     &        chnam_tb,ntabl_ihash,hashtblink,nhashmx,
     &        linematch,fcontr,find)
          end if
c
          if(find)then
          ffAtomName(i) = atomFFname_tb(linematch) 
          atomCharge(i) = atomFFq_tb(linematch)
          end if 
c
	if(CONTROL0)then
        write(kanalp,'(i5,1x,a4,a8,2x,a2,2x,f8.5,a6,i6)')i,
     &  atmNameEx(i)(1:4),resn8,ffAtomName(i), atomCharge(i),
     &  ' tbline: ',linematch
        end if
          end do !i
c
         write(kanalp,*)'End of defFFatomName:'
c
	return
        end
