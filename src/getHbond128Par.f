c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  getHbond128Par()                                                          *
c*  Defines  the ForceField  parameters for Hb Hx-Y vdw12-8                   *
c*                                                                            *
c*     Yury Vorobjev 2004                                                     *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getHbond128Par(ffParFile,nHXhb128MX,nYYhb128MX,
     &             nHXhb128,nYYhb128,hb128HxList,hb128YYList,
     &             nHbond128Types,hB128pariHjY,hB128parRE,
     &             hb128parAB)
c
c InPut:
c       ffParFile - ffParameters file 
c
c RESULT: nHXhb128,nYYhb128 - numbers of HXdonors, acceptors YY 
c         hb128HxList() - list of atNameFF for HxDonors
c         hb128YYList() - list of atNameFF for YY acceptors
c         hB128par(2*nHXhb128*nYYhb128) - rMin, eMin vdw12-8
c
        implicit none
        character*(*) ffParFile     
        integer nHXhb128MX,nYYhb128MX
        integer nHXhb128,nYYhb128
        character*2 hb128HxList(*)
        character*2 hb128YYList(*)
        integer nHbond128Types
        integer hB128pariHjY(*)
        real hB128parRE(*)
        real hb128parAB(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
clocal
        integer getPointHB128ParTable
        integer ntbrec 
        character*80 line
        character*9 recBlock
        character*4 endBlock 
        integer nFormatMX 
        integer k,kShift 
        integer iHx,jYY        
        integer nrb
        real rm4,rm6,rm8
        logical HBPOT12_8,HBPOT12_6
        integer kanalp
        logical myblock
        logical CONTROL,CONTROL1
c
c initialize
        kanalp = kanalRunOut
        CONTROL  = .true.  
        CONTROL1 = .false.
c
        HBPOT12_8 = .false.
        HBPOT12_6 = .true.
c
c read in ffParfile
       if(CONTROL)then
       write(kanalp,*)'getHbond128Par :01: readFile:',ffParFile
       end if
c
        open(unit=11, file=ffParFile, form='formatted',
     &       status= 'old')
c
c READ ATOMMASS:
        recBlock = '$HBOND128'
        nrb = 9
        endBlock = '$END'
        myblock = .false.
        ntbrec = 0
        nFormatMX=10
c
100     read(11,'( a80 )', end=101 ) line
c
        if(CONTROL1)then
        write(kanalp,*)line
        write(kanalp,*) 'myblock :', myblock 
        end if
c
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
c 
        if(line(1:6) .eq. 'ATHX_N')then
        read(line, '(6x, i4)') nHXhb128
        if(nHXhb128 .gt. nHXhb128MX)then
        write(kanalp,*)'getHbond128Par:ERROR: nHXhb128 too large',
     &  ' increase nHXhb128MAX in hbond128.h ' 
        STOP
        end if !nHXhb128 .gt. nHXhb128MX
        end if ! ATHX_N
c
        if(line(1:6) .eq. 'ATYY_N')then
        read(line, '(6x, i4)') nYYhb128
        if(nYYhb128 .gt. nYYhb128MX)then
        write(kanalp,*)'getHbond128Par:ERROR: nYYhb128 too large',
     &  ' increase nYYhb128MAX in hbond128.h '
        STOP
        end if !nYYhb128 .gt. nYYhb128MX
        end if ! ATYY_N
c 
        if(CONTROL)then
        write(kanalp,*)
     &  'getHbond128Par:nHXhb128,nYYhb128:',nHXhb128,nYYhb128 
        end if
c
        if(line(1:6) .eq. 'ATHX_L')then
        read(line, '(6x, 15(2x,a2))')(hb128HxList(k),k=1,nHXhb128)
        end if ! ATHX_N
c
        if(line(1:6) .eq. 'ATYY_L')then
        read(line, '(6x, 15(2x,a2))')(hb128YYList(k),k=1,nYYhb128)
        end if ! ATYY_N
c
        if(line(1:6) .eq. 'HBPAR ')then
        read(line, '(6x,2i4)')iHx,jYY   
c       kShift= 2*(nYYhb128*(iHx-1) + jYY-1)
        kShift =  getPointHB128ParTable(iHx,jYY,nYYhb128)
c
        if(CONTROL1)then
        write(kanalp,*)'getHbond128Par: kShift=',kShift
        end if
c
        hB128pariHjY(kShift+1)=iHx
        hB128pariHjY(kShift+2)=jYY
        read(line, '(14x,2f10.6)')
     &  hB128parRE(kShift+1),hB128parRE(kShift+2) 
c calculate A12, B8
c
         if(HBPOT12_8)then
         rm4 =  hB128parRE(kShift+1)**4
         rm8 = rm4**2
         hb128parAB(kShift+1) = hB128parRE(kShift+2)*rm8*rm4
         hb128parAB(kShift+2) = hB128parRE(kShift+2)*rm8*1.50
         end if
c
         if(HBPOT12_6)then
         rm6 = hB128parRE(kShift+1)**6
         hb128parAB(kShift+1) = hB128parRE(kShift+2)*rm6*rm6
         hb128parAB(kShift+2) = hB128parRE(kShift+2)*rm6*2.00
         end if
c
         end if !HBPAR
c
        end if !myblock
c
        goto 100
c
101     continue
c
c end of readIn ffParFile
c
        nHbond128Types = nHXhb128*nYYhb128 
        if(CONTROL)then
        write(kanalp,*)'getHbond128Par: '
        write(kanalp,'(a6, 15(2x,a2))')
     &  'ATHX_N',(hb128HxList(k),k=1,nHXhb128) 
        write(kanalp,'(a6, 15(2x,a2))')
     &  'ATYY_L',(hb128YYList(k),k=1,nYYhb128) 
        write(kanalp,*)'Y-HX   Hn  Ym Rmin(A)   Em(kcal/m) : vdw 12-8'
        do k=1,nHbond128Types
        write(kanalp,'(a6,2i4,1x,2f10.6)')
     & 'HBPAR ',hB128pariHjY(2*k-1),hB128pariHjY(2*k),
     &          hB128parRE(2*k-1),hB128parRE(2*k)
        end do!k 
        end if!CONTROL
c
        if(CONTROL1) STOP
c
        return
	end
c
	integer function getPointHB128ParTable(iHx,jYY,nYYhb128)
        implicit none
        integer iHx,jYY,nYYhb128
c
        getPointHB128ParTable = 2*(nYYhb128*(iHx-1) + jYY-1)
c
        return
        end
