C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     getTorsPar()                                                           *
C* Defines  the ForceField torsionParameters                                  *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getTorsPar(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              quar1234L,nQuar1234,quar1234Par,quar1234nPar)        
c
c InPut:
c       ffParFile - ffParameters file 
c       natom,atomNameEx,ResName,chName : PDB info
c       ffAtomName(ia) - FFatomName to search table
c       the quar1234L(i),i=1,..,nQuar1234 : the QuartetList
c RESULT: quar1234Par(16*nQuar1234) - torsionFF parameters for list of quartets 
c         pass,Vt/2,delta,nFi - (printed) for each torsHarmonics,
c         pass,Vt/2/pass,cos(delta),nFi - finally in array
c        4- torsionHarmanics is possible.
c        quar1234nPar(iQuart) - number of torsHarmonics for the torsAngl
c
        implicit none
        character*(*) ffParFile     
        integer natom
        character*8 atomNameEx(*)
        character*4 ResName(*)
        character*1 chName(*)
        character*2 ffAtomName(*)
        integer quar1234L(*)
        integer nQuar1234
        real quar1234Par(*)
        integer quar1234nPar(*)
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
        character*(tlen) quarName_tb(ntbmax)
        integer pass_tb(ntbmax)
        real    Vt_tb(ntbmax)
        real    delta_tb(ntbmax)
        real    nFi_tb(ntbmax)
        character*20 line_tb(ntbmax) 
c hashTable local
        integer nhashmx
        parameter (nhashmx = ntbmax )
        integer ntabl_ihash(nhashmx)
        integer hashtblink(nhashmx)    
        integer nhash
c
        character*80 line
        integer nwordMAX
        parameter (nwordMAX = 8)
        character*(tlen) abcdw(nwordMAX)    
        integer nlinematchMAX,nlinematchMX
        parameter (nlinematchMAX = 4)
        integer linematch(nlinematchMAX)
        integer nlinematch
        integer nTorsHarmMAX,nTorsParAng
c
        character*1 s1
        character*2 aa,bb,cc,dd,xx
        character*8 recBlock
        character*4 endBlock 
        integer i,i4,i16,k,k4,t
        integer len
        integer iw,nw
        integer ntbrec,contr
        integer nstr1
        integer nlrcb
        integer kanalp
        logical find
        logical myblock
        logical CONTROL,CONTROL1,CONTROL0
        logical OPT_hash
c
c initialize
        kanalp = kanalRunOut
        CONTROL0 = .true.    ! alltime
        CONTROL =  .false.    ! test ! false
        CONTROL1 = .false.   ! test ! false
        OPT_hash = .true.
        xx = 'X '
        nlinematchMX = nlinematchMAX
        nTorsParAng = 4
        nTorsHarmMAX = 4
c read in ffParfile
        open(unit=11, file=ffParFile, form='formatted',
     &       status= 'old')
c
        ntbrec = 0
        if(CONTROL0)then
          write(kanalp,*)'In getTorsPar:'
          write(kanalp,*)'ffParFile : ',ffParFile
          write(kanalp,*)
     &  'nTang nTlib   a-b-c-d      i-j-k-l       nPass Vt/2  delt nFi'        
        end if
c
        recBlock = '$TORSPAR'
        nlrcb=8
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
        write(kanalp, *)'ERROR!in defTorsPar: ',
     &  ' ntbmax (param) is SMALL '
c
        write(kanalPStat,*)mError,
     &    ' defTorsPar: ntbmax (param) is SMALL ' 
c
        stop
        end if
c 
        read(line, '(a11,2x,i3  ,2x,f6.3,6x,f8.1,8x,f8.1)')
     &  quarName_tb(ntbrec),pass_tb(ntbrec),Vt_tb(ntbrec),
     &  delta_tb(ntbrec),nFi_tb(ntbrec)
c
        if(CONTROL1)then
	write(kanalp,'(a11,2x,i3  ,2x,f6.3,6x,f8.1,8x,f8.1)')
     &  quarName_tb(ntbrec),pass_tb(ntbrec),Vt_tb(ntbrec),
     &  delta_tb(ntbrec),nFi_tb(ntbrec)
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
       line_tb(i) = quarName_tb(i)   ! string length
c calculate hash number from string strtabl1(1:nstr1) simbols
       call hashCalculator(quarName_tb(i),nstr1,nhashmx,nhash,contr)
 
c hashtblink() - linkedList to resolve degenerate nhash
c ntabl_ihash() - head() of linked list for given nhash
 
       hashtblink(i) = ntabl_ihash(nhash)
       ntabl_ihash(nhash) = i
 
       if(CONTROL1)then
       write(kanalp,'(a27,i6,i6,a11,a11)')
     & 'hash_dict_tabl:itab,nhash:',i,nhash,' str:',
     &  quarName_tb(i)
       end if

       end do !i
c
c find match for quartet  
        s1 = "-"
        nw = nwordMAX
c
        do i = 1, nQuar1234
c
c init
        i16 = nTorsHarmMAX*nTorsParAng*(i-1)
        i4 = 4*i-4
c  nTorsHarmMAX=4, nTorsParAng=4 defined in the pair1234array.h
        quar1234nPar(i) = 0                   
c number of torsHarm for this torsAng
        do k = 1,nTorsHarmMAX*nTorsParAng
        quar1234Par(i16 + k) = 0.0
        end do !k
c
        aa = ffAtomName(quar1234L(i4+1))
        bb = ffAtomName(quar1234L(i4+2))
        cc = ffAtomName(quar1234L(i4+3))
        dd = ffAtomName(quar1234L(i4+4))
c
c make all possible words of the torsLib from aa,bb,cc,dd
        abcdw(1) = aa//s1//bb//s1//cc//s1//dd 
        abcdw(2) = dd//s1//cc//s1//bb//s1//aa
        abcdw(3) = xx//s1//bb//s1//cc//s1//xx
        abcdw(4) = xx//s1//cc//s1//bb//s1//xx
        abcdw(5) = xx//s1//bb//s1//cc//s1//dd
        abcdw(6) = xx//s1//cc//s1//bb//s1//aa
        abcdw(7) = xx//s1//xx//s1//cc//s1//dd
        abcdw(8) = xx//s1//xx//s1//bb//s1//aa
c
        if(CONTROL)then
        write(kanalp,'(a11,i6,1x,8(a11,1x))')'quartet:i:',i,abcdw
        end if
c
        contr = 0  !
cx        contr = 1  ! test
c
        do iw = 1,nw
        if(CONTROL)then
        write(kanalp,*)'find_match: iw,abcdw(iw):',iw,abcdw(iw)
        end if
c
        nlinematch = 0
        call find_match_hashtb01(abcdw(iw),len,
     &         ntbrec,line_tb,
     &         ntabl_ihash,hashtblink,nhashmx,
     &         nlinematch,nlinematchMX,linematch,contr,find)
c
        if(find)goto 201
        end do !iw
c
c assign torsPar from the lines
201     if(find)then
        quar1234nPar(i) = nlinematch             !number of torsHarmonics
c
        if(CONTROL)then
        write(kanalp,*)
     &  'getTorsPar:TorsHarm=nlinematch:',quar1234nPar(i)
        end if !C
c
        do k = 1,nlinematch
        k4 = 4*(k-1)
	quar1234Par(i16+k4+1) = pass_tb(linematch(k))
	quar1234Par(i16+k4+2) = Vt_tb(linematch(k))
	quar1234Par(i16+k4+3) = delta_tb(linematch(k))
	quar1234Par(i16+k4+4) = nFi_tb(linematch(k))
c
        if(CONTROL)then
        write(kanalp,
     &  '(i4,1x,i4,1x,a11,1x,
     &    4(i4,1x),1x,f3.1,1x,f6.3,2x,f5.1,1x,f5.1)')         
     &  i, linematch(k),
     &  abcdw(1), (quar1234L(i4+t),t=1,4),
     &            (quar1234Par(i16+k4+t),t=1,4)
        end if !control
c
c process for convinience to USE in sub: torsanglenf
	quar1234Par(i16+k4+2) = Vt_tb(linematch(k))
     &                          /pass_tb(linematch(k))
	quar1234Par(i16+k4+3) = cos(0.01745*delta_tb(linematch(k)))
	quar1234Par(i16+k4+4) = abs(nFi_tb(linematch(k)))
c        
        if(CONTROL1)then
        write(kanalp,
     &  '(i4,a1,1x,i3,1x,a11,1x,
     &    4(i4,1x),1x,f3.1,1x,f6.3,2x,f5.1,1x,f5.1)')         
     &  i,'*',linematch(k),
     &  abcdw(1), (quar1234L(i4+t),t=1,4),
     &            (quar1234Par(i16+k4+t),t=1,4)
        end if !control
c
        end do !k
c
        else

        if(CONTROL0)then
        write(kanalp,'(i4,1x,4x,1x,a11,1x,4(i4,1x),a40)')         
     &  i, abcdw(1),(quar1234L(i4+t),t=1,4),
     &  ' WARNING!: par not found! eTors=0 '
        end if !CONTROL0
c
        end if !find
c
        end do!i
c        
        return
	end
