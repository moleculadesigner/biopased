C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     getVangDefPar()                                                        *
C* Defines  the ForceField anglDefParameters                                  *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getVangDefPar(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              trip123List,nTrip123,ang123ParL)
c
c InPut:
c       ffParFile - ffParameters file 
c       ffAtomName(ia) - FFatomName to search table
c RESULT: ang123ParL(Kb,a0) for list of anglTriplets in the trip123List(i),i=1,..,nTrip123
c
        implicit none
        character*(*) ffParFile     
        integer natom
        character*8 atomNameEx(*)
        character*4 ResName(*)
        character*1 chName(*)
        character*2 ffAtomName(*)
        integer trip123List(*)
        integer nTrip123
        real ang123ParL(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer ntbmax                  !convertionTablsiz.h
        parameter (ntbmax=2047)         !MaxNumb of table lines
        integer tlen
        parameter (tlen = 8)
        character*(tlen) angTriplName_tb(ntbmax)
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
        character*2 aa,bb,cc
        character*8 recBlock
        character*4 endBlock 
        integer i,i2,i3,len
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
          write(kanalp,*)'In defBondDefPar:'
          write(kanalp,*)'ffParFile : ',ffParFile
        end if
c
        recBlock = '$ANGLPAR'
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
        if(line(1:8) .eq. recBlock )then
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
c
        write(kanalPStat,*)mError,
     &    ' defBondPar: ntbmax (param) is SMALL '
        stop
        end if
c 
        read(line, '(a8,3x,f9.3,f9.3)')
     &  angTriplName_tb(ntbrec),kb_tb(ntbrec),a0_tb(ntbrec)
c
        if(CONTROL1)then
	write(kanalp,'(a8,3x,f9.3,f9.3,1x,i6)')
     &  angTriplName_tb(ntbrec),kb_tb(ntbrec),a0_tb(ntbrec),ntbrec
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
       line_tb(i) = angTriplName_tb(i)   ! string length
c calculate hash number from string strtabl1(1:nstr1) simbols
       call hashCalculator(angTriplName_tb(i),nstr1,nhashmx,nhash,contr)
 
c hashtblink() - linkedList to resolve degenerate nhash
c ntabl_ihash() - head() of linked list for given nhash
 
       hashtblink(i) = ntabl_ihash(nhash)
       ntabl_ihash(nhash) = i
 
       if(CONTROL1)then
       write(kanalp,'(a27,i6,i6,a8,a8)')
     & 'hash_dict_tabl:itab,nhash:',i,nhash,' str:',
     &  angTriplName_tb(i)
       end if

       end do !i
c
c find match for bondPairs
        s1 = "-"
        if(CONTROL1)then
        write(kanalp,*)'getVangDefPar: result: '
        write(kanalp,*)'                        ',
     &  ' Kgrom,cosA0: Kang(ffdat),A0'
        end if
c
        do i = 1, nTrip123
        i3 = i*3
        i2 = 2*i
        aa = ffAtomName(trip123List(i3-2))
        bb = ffAtomName(trip123List(i3-1))
        cc = ffAtomName(trip123List(i3))
        aabbcc = aa//s1//bb//s1//cc
        ccbbaa = cc//s1//bb//s1//aa
c
        contr = 0  !
c
        call find_match_hashtb0(aabbcc,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
        if ( .not. find )then
        call find_match_hashtb0(ccbbaa,len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,linematch,contr,find)
        end if
c
c defualt values:
	ang123ParL(i2-1) = 0.0    ! kb_tb(linematch)
	ang123ParL(i2) =  -0.334  ! default cos(a0), a0=109.5 
c
        if(find)then
c read  Kang(harm)/2  from fftable
        call getInternalAngPar(
     &       kb_tb(linematch),a0_tb(linematch),
     &       ang123ParL(i2-1),ang123ParL(i2))
c
        if(CONTROL1)then
        write(kanalp,'(4i6,1x,a8,1x,a8,1x,2(f9.3,f9.3))')
     &   i, trip123List(i3-2),trip123List(i3-1),trip123List(i3), 
     &  aabbcc,ccbbaa, ang123ParL(i2-1),ang123ParL(i2)
     &  ,kb_tb(linematch),a0_tb(linematch)
        end if !control

        else
        if(CONTROL)then
           write(kanalp,'(4i6,1x,a8,1x,a8,1x,a45)')
     &   i, trip123List(i3-2),trip123List(i3-1),trip123List(i3),
     &  aabbcc,ccbbaa,
     &  'WARNING! AnglDefPar are not found: KangDef=0 '
        end if !control
        end if !find

        end do!i
        

        return
	end
c
	subroutine getInternalAngPar(Kh05,a0,kcs,cs0)
c Kh05 = Kharm/2
c convert 0.5kh(a0-a)**2 - classical function
c to gromos function Eangl= 0.5*kcs*[cosA - cosA0]**2
c
	implicit none
	real Kh05,a0,kcs,cs0
        real kt,a0r,a1,ar
        kt=0.596
        ar=3.1415927/180.0
        a0r=a0*ar
        cs0 = cos(a0r)
        a1 = sqrt(kt/(2.0*Kh05))
c
        if(kh05 .le. 0.0)then
        kcs=0.0
        else
        kcs = 2.0*kt/( (cos(a0r+a1) - cs0)**2+
     &        (cos(a0r-a1) - cs0)**2)
        end if
c
        return
        end
