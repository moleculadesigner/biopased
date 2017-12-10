C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C* getBondDefPar()                                                            *
C* Defines  the ForceField bondDef parameters                                 *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getBondDefPar(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              bond12List,nbond12,bond12ParL)
c
c InPut:
c       ffParFile - ffParameters file 
c       ffAtomName(ia) - FFatomName to search table
c RESULT: bond12ParL(2*ib-1),(2*ib) = (Kb,b0) for gromos vbond function for 
c                         bond list in bond12List(ib),ib=1,..,nbond12 
        implicit none
        character*(*) ffParFile     
        integer natom
        character*8 atomNameEx(*)
        character*4 ResName(*)
        character*1 chName(*)
        character*2 ffAtomName(*)
        integer bond12List(*)
        integer nbond12
        real bond12ParL(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer ntbmax                  !convertionTablsiz.h
        parameter (ntbmax=2047)         !MaxNumb of table lines
c
        integer tlen
        parameter (tlen = 5)        
        character*(tlen) bondPairName_tb(ntbmax)
        real    kb_tb(ntbmax)
        real    b0_tb(ntbmax)
        character*20 line_tb(ntbmax) 
c hashTable local
        integer nhashmx
        parameter (nhashmx = ntbmax )
        integer ntabl_ihash(nhashmx)
        integer hashtblink(nhashmx)    
        integer nhash
c
        character*80 line
        character*(tlen) aabb,bbaa    
        character*1 s1
        character*2 aa,bb
        character*8 recBlock
        character*4 endBlock 
        integer i,i2,len,nrb
        integer ntbrec,contr
        integer nstr1
        integer kanalp
        integer linematch
        logical find
        logical myblock
        logical CONTROL,CONTROL1
        logical OPT_hash
c
c initialize
c        kanalp = 6
        kanalp = kanalRunOut
c
        CONTROL = .true. 
        CONTROL1 = .false.
        OPT_hash = .true.

c read in ffParfile
        open(unit=11, file=ffParFile, form='formatted',
     &       status= 'old')
c
        ntbrec = 0
        if(CONTROL)then
          write(kanalp,*)'In defBondDefPar:'
          write(kanalp,*)'ffParFile : ',ffParFile
        end if
c
        recBlock = '$BONDPAR'
        nrb = 8
        endBlock = '$END'
        myblock = .false.
c
100     read(11,'( a80 )', end=101 ) line
c
        if(CONTROL1)then
        write(kanalp,*)line
        end if
        if(myblock .and. (line(1:4) .eq. endBlock) ) goto 101
        if(line(1:nrb) .eq. recBlock )then
         myblock = .true.  
         goto 100
         end if 
c
        if(.not.myblock)goto 100
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        if(myblock)then
        ntbrec=ntbrec+1
        if(ntbrec .gt. ntbmax) then
        write(kanalp, *)'ERROR!in defBondPar: ',
     &  ' ntbmax (param) is SMALL '
        stop
        end if
c 
        read(line, '(a5,2x,f8.2,f8.5)')
     &  bondPairName_tb(ntbrec),kb_tb(ntbrec),b0_tb(ntbrec)
c
        if(CONTROL1)then
	write(kanalp,'(a5,2x,f8.2,f8.5,1x,i6)')
     &  bondPairName_tb(ntbrec),kb_tb(ntbrec),b0_tb(ntbrec),ntbrec
        end if
c
        end if !myblock
c
        goto 100
c
101     continue
        close(unit=11)
c end of readIn ffParFile
c
c calculate hashTable for dictionary table
       len = tlen
       nstr1 = len
       contr = 0     !control print
c init hashtable
       do i = 1,nhashmx
       hashtblink(i) = 0
       ntabl_ihash(i) = 0
       end do!i
c 
       do i = 1,ntbrec
       line_tb(i) = bondPairName_tb(i)   ! string length
c calculate hash number from string strtabl1(1:nstr1) simbols
       call hashCalculator(bondPairName_tb(i),nstr1,nhashmx,nhash,contr)
 
c hashtblink() - linkedList to resolve degenerate nhash
c ntabl_ihash() - head() of linked list for given nhash
 
       hashtblink(i) = ntabl_ihash(nhash)
       ntabl_ihash(nhash) = i
 
       if(CONTROL1)then
       write(kanalp,'(a27,i6,i6,a5,a5)')
     & 'hash_dict_tabl:itab,nhash:',i,nhash,' str:',
     &      bondPairName_tb(i)
       end if

       end do !i
c
c find match for bondPairs
        s1 = "-"
        if(CONTROL1)then
        write(kanalp,*)'getBondDefPar: result: ',
     &  ' Kb(gromos) b0   kb(datTable)   b0'
        end if
c
        do i = 1, nbond12
        i2 = i*2
        aa = ffAtomName(bond12List(i2-1))
        bb = ffAtomName(bond12List(i2))
        aabb = aa//s1//bb
        bbaa = bb//s1//aa
c
        contr = 0  !
c
        call find_match_hashtb0(aabb,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
        if ( .not. find )then
        call find_match_hashtb0(bbaa,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
        end if
c
c default value
        bond12ParL(i2-1) = 0.0   ! kb
        bond12ParL(i2) =   1.30  ! b0
c
        if(find)then
c read Kharm/2 from ffdata table = kb_tb(line) 
        call getInternalKb
     &  (kb_tb(linematch),b0_tb(linematch),bond12ParL(i2-1))
c 
	bond12ParL(i2) = b0_tb(linematch)
c
        if(CONTROL)then
        write(kanalp,'(3i6,1x,a5,1x,a5,1x,2(f8.2,f8.5))')
     &   i, bond12List(i2-1),bond12List(i2), 
     &  aabb,bbaa, bond12ParL(i2-1),bond12ParL(i2)
     &  ,kb_tb(linematch),b0_tb(linematch)
        end if !control

        else
        if(CONTROL)then
        write(kanalp,'(3i6,1x,a5,1x,a5,1x,a40)')
     &   i, bond12List(i2-1),bond12List(i2),
     &  aabb,bbaa, 'WARNING! bondPar is not found: Kb = 0 !'
c
        write(kanalPStat,*)mError,
     &  'WARNING! bondPar is not found: see runOUT file '
c
        end if!Control
c
        end if !find

        end do!i
c        
        return
	end
c
	subroutine getinternalKb(kharm05,b0,kb)
c convert Kharm from FFdata to kb convenient for energy calculation
c Ebond = 0.25kb*[b**2 - b0**2]**2   - gromos function
	implicit none
	real kharm05,kb,b0
c
        kb = (2.0*kharm05)/(2.0*b0**2)
c
	return
	end
