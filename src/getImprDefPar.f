C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     getImpDefPar()                                                        *
C* Defines  the ForceField impanglDefParameters                                  *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getImpDefPar(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              quarImp1234L,nImp1234,impAng1234Par)
c
c InPut:
c       ffParFile - ffParameters file 
c       ffAtomName(ia) - FFatomName to search table
c RESULT: impAng1234Par(Kb,a0) for list of anglTriplets in the quarImp1234L(i),i=1,..,nImp1234
c
        implicit none
        character*(*) ffParFile     
        integer natom
        character*8 atomNameEx(*)
        character*4 ResName(*)
        character*1 chName(*)
        character*2 ffAtomName(*)
        integer quarImp1234L(*)
        integer nImp1234
        real impAng1234Par(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer ntbmax                  !convertionTablsiz.h
        parameter (ntbmax=2047)         !MaxNumb of table lines
        integer tlen
        parameter (tlen = 11)
        character*(tlen) impQuarName_tb(ntbmax)
        character*(tlen) bbw(6)
        real    kb_tb(ntbmax)
        real    a0_tb(ntbmax)
        character*20 line_tb(ntbmax) 
c hashTable local
        integer nhashmx
        parameter (nhashmx = ntbmax )
        integer ntabl_ihash(nhashmx)
        integer hashtblink(nhashmx)    
        integer nhash
c
        character*80 line
        character*(tlen) aabbcc,ccbbaa    
        character*1 s1
        character*2 aa,bb,cc,dd
        character*9 recBlock
        character*4 endBlock 
        integer i,i2,i4,len
        integer iw,nw
        integer ntbrec,contr
        integer nstr1
        integer nlrcb
        integer kanalp
        integer linematch
        logical find
        logical myblock
        logical CONTROL,CONTROL1
        logical OPT_hash
c
c initialize
        kanalp = kanalRunOut 
        CONTROL = .true.
        CONTROL1 = .false. 
        OPT_hash = .true.
c read in ffParfile
        open(unit=11, file=ffParFile, form='formatted',
     &       status= 'old')
c
        ntbrec = 0
        if(CONTROL)then
          write(kanalp,*)'In getImprDefPar:'
          write(kanalp,*)'ffParFile : ',ffParFile
        end if
c
        recBlock = '$IMPROPER'
        nlrcb=9
        endBlock = '$END'
        myblock = .false.
        len = tlen
c
100     read(11,'( a80 )', end=101 ) line
c
        if(CONTROL1)then
        write(kanalp,*)line
        end if
        if(myblock .and. (line(1:4) .eq. endBlock) ) goto 101
        if(line(1:nlrcb) .eq. recBlock )then
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
        read(line, '(a11,3x,f9.3,f9.3)')
     &  impQuarName_tb(ntbrec),kb_tb(ntbrec),a0_tb(ntbrec)
c
        if(CONTROL1)then
	write(kanalp,'(a11,3x,f9.3,f9.3,1x,i6)')
     &  impQuarName_tb(ntbrec),kb_tb(ntbrec),a0_tb(ntbrec),ntbrec
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
       nstr1 = len
       contr = 0     !control print
c init hashtable
       do i = 1,nhashmx
       hashtblink(i) = 0
       ntabl_ihash(i) = 0
       end do!i
c 
       do i = 1,ntbrec
       line_tb(i) = impQuarName_tb(i)   ! string length
c calculate hash number from string strtabl1(1:nstr1) simbols
       call hashCalculator(impQuarName_tb(i),nstr1,nhashmx,nhash,contr)
 
c hashtblink() - linkedList to resolve degenerate nhash
c ntabl_ihash() - head() of linked list for given nhash
 
       hashtblink(i) = ntabl_ihash(nhash)
       ntabl_ihash(nhash) = i
 
       if(CONTROL1)then
       write(kanalp,'(a27,i6,i6,a11,a11)')
     & 'hash_dict_tabl:itab,nhash:',i,nhash,' str:',
     &  impQuarName_tb(i)
       end if

       end do !i
c
c find match for bondPairs
        s1 = "-"
        if(CONTROL)then
        write(kanalp,*)'getImpAngPar:    Kdef, A0'
        end if
c
        do i = 1, nImp1234
        i4 = i*4
        i2 = 2*i
c
c default values Planar:
	impAng1234Par(i2-1) = 40.0 ! kcal/rad^2            
	impAng1234Par(i2) = 0.0               
c
        aa = ffAtomName(quarImp1234L(i4-3))
        bb = ffAtomName(quarImp1234L(i4-2))
        cc = ffAtomName(quarImp1234L(i4-1))
        dd = ffAtomName(quarImp1234L(i4))
c*
c        aabbccdd = aa//s1//bb//s1//cc//s1//dd
c make all permutations of aa,cc,dd
        nw = 6
        bbw(1) = aa//s1//bb//s1//cc//s1//dd 
        bbw(2) = aa//s1//bb//s1//dd//s1//cc
        bbw(3) = cc//s1//bb//s1//dd//s1//aa
        bbw(4) = cc//s1//bb//s1//aa//s1//dd
        bbw(5) = dd//s1//bb//s1//aa//s1//cc
        bbw(6) = dd//s1//bb//s1//cc//s1//aa
c
        contr = 0  !
c
        do iw = 1,nw
        call find_match_hashtb0(bbw(iw),len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
        if(find)goto 201
        end do !iw
c
201     if(find)then
	impAng1234Par(i2-1) = kb_tb(linematch)
	impAng1234Par(i2) = a0_tb(linematch)
c
        if(CONTROL1)then
        write(kanalp,'(i6,1x,4i6,1x,a11,1x,2(f9.3,f9.3))')
     &  i, quarImp1234L(i4-3),quarImp1234L(i4-2),
     &     quarImp1234L(i4-1),quarImp1234L(i4), 
     &  bbw(1), impAng1234Par(i2-1),impAng1234Par(i2)
     &  ,kb_tb(linematch),a0_tb(linematch)
        end if !control

        else
        if(CONTROL)then
        write(kanalp,'(i6,1x,4i6,1x,a11,1x,a49))')
     &  i, quarImp1234L(i4-3),quarImp1234L(i4-2),
     &     quarImp1234L(i4-1),quarImp1234L(i4),
     &  bbw(1), 'WARNING!:Imp are not found: deflt: KIm=40,planar'
c
cx        write(kanalPStat,*)mError,
cx     &  ' WARNING!:Imp are not found: deflt are taken: KIm=40,planar'
c
        end if !control

        end if !find
        end do!i
c        
        return
	end
