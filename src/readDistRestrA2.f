c irestDistA2 read restrDistInfo 
c reads inf from file and calculates data for restrDistA2.h
c
	subroutine readRestrDistA2
     &             (restFile,nRestrDistA2MX,nRestrDistA2,
     &                       restrDistA2List,restrDistA2DistHK)
c
c InPut:
c restrfile - "./distRestrA2.inp"  file name with input information
c #restr residues (define restricted MOLECule)
c (a6,i6,i6,1x,f6.2,1x,f6.2)
c RESTA2   AtNumb1  AtNumb2  Dist12   hK12 (kcal/A**2)
c RESTA2   AtNumb1  AtNumb2  Dist12   hK12
c end
c -------------------
c nRestrDistA2MX - MAX number of atomPairs           
c nRestrDistA2   - number of atomPairs                       
c OUT:
c restrDistA2List(2*nRestrDistA2) - pairs of atom List                               
c restrDistA2DistK(2*nRestrDistA2) - Dist12,HKonst for atom pair
c
        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        character*(*) restFile
        integer nRestrDistA2MX,nRestrDistA2
        integer restrDistA2List(*)
        real restrDistA2DistHK(*)
c local
        character*80 line
        character*6 keyw,keywf
        integer nlkey
        integer i,i2,nline
        integer kanalInp
        integer kanalp
        logical CONTROL
c      
        kanalp = kanalRunOut
        CONTROL = .true.
        keyw = 'RESTA2' 
        nlkey = 6
c
        write(kanalp,*)'readDistRestrA2:nRestrDistA2MX:',
     &  nRestrDistA2MX
        write(kanalp,*)'readDistRestrA2: restFile:',restFile
c     
        do i=1,nRestrDistA2MX
        i2=2*i-2
        restrDistA2List(i2+1)=0
        restrDistA2List(i2+2)=0
        restrDistA2DistHK(i2+1)=0.0
        restrDistA2DistHK(i2+2)=0.0
        end do!i
c
        nRestrDistA2=0
c 
c read in file  restDistA2.inp
c
        kanalInp = 27
        open(unit=kanalInp, file=restFile, form='formatted',
     &       status= 'old') 
c
	nline = 0
 100    read(kanalInp,'( a80 )', end=101 ) line
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        nline=nline+1
	if(line(1:nlkey) .eq. keyw)then
        i2 = 2*nRestrDistA2
        nRestrDistA2 = nRestrDistA2 + 1 
c
        if(nRestrDistA2 .gt. nRestrDistA2MX)then
        write(kanalp,*)'readDistRestrA2:ERROR!nRestrDistA2 too large:',
     &  nRestrDistA2
        stop
        end if
c
        read(line, 707, end=101)keywf,
     &  restrDistA2List(i2+1),restrDistA2List(i2+2),
     &  restrDistA2DistHK(i2+1),restrDistA2DistHK(i2+2)
        end if !RESTA2
        goto 100

101     continue
c
        close(unit=kanalInp)
c
        if(CONTROL)then
        write(kanalp,*)'readDistRestrA2: nRestrDistA2:',nRestrDistA2
        do i=1,nRestrDistA2 
        i2=i*2-2
        write(kanalp,707)keyw,
     &  restrDistA2List(i2+1),restrDistA2List(i2+2),
     &  restrDistA2DistHK(i2+1),restrDistA2DistHK(i2+2)
        end do!
        write(kanalp,*)'readDistRestrA2 : finish OK'
        end if!control
c
707     format(a6,i6,i6,1x,f6.2,1x,f6.2)
c
	return
        end
