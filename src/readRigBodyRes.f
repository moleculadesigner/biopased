c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2004                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c read rigidBodyRes
	subroutine readRigBodyRes 
     &             (rgbFile,nRigBodyMX,nRigBody,rigBodyStEndRes)
c
c InPut:
c rgbFile - file name with start/end Res of rigidBody segments
c #rigibBodyRes
c RIGB01  11  16
c RIGB02  47  59
c end
c -------------------
c nRigBodyMX - MAX number of rigidBodySegment of residues
c
c OUT:
c nRigBody - number of rigidBody Segment
c rigBodyStEndRes(2*nRigBody) - nkN,nkC: start/end ResNumb of the rigBodySegments
c
        implicit none
        character*(*) rgbFile
        integer nRigBodyMX,nRigBody    
        integer rigBodyStEndRes(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        character*80 line
        character*6 keyw
        integer i,i2,nline
        integer nlkey
        integer kanalp
        integer kanalRigB
        logical fext
        logical CONTROL
c      
        kanalp=kanalRunOut
        CONTROL = .true.
        keyw = 'RIGB' 
        nlkey = 4
        kanalRigB = 22
c     
        do i=1,nRigBodyMX
        i2=2*i-1
        rigBodyStEndRes(i2)=0
        rigBodyStEndRes(i2+1)=0
        end do
c
        nRigBody = 0
c 
        write(kanalp,*)'readRigBodyRes:01: rgbFile:',rgbFile
c
c read in file  rigBodyRes.inp    
        inquire(file=rgbFile, exist=fext)
        if(.not. fext) then
        write(kanalp,*) 'ERROR:readRigBodyResFile:',rgbFile,
     &  ' does not exist'
        stop
        end if
c
        open(unit=kanalRigB, file=rgbFile, form='formatted',
     &       status= 'old') 
c
	nline = 0
 100    read(kanalRigB,'( a80 )', end=101 ) line
c
        write(kanalp,*)line
c
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        nline=nline+1
	if(line(1:nlkey) .eq. keyw)then
        i2 = 2*nRigBody
        nRigBody = nRigBody + 1
cx        read(line, '(6x,1x,i4,1x,i4)')
         read(line, '(6x,i4,i4)') 
     &   rigBodyStEndRes(i2+1),rigBodyStEndRes(i2+2)
        end if !RIGB 
        goto 100

101     continue
        close(unit=kanalRigB)
c
        if(CONTROL)then
        write(kanalp,*)'readRigBodyRes:nRigBody:',nRigBody      
        do i=1,nRigBody
        i2=i*2-2
        write(kanalp,*)'rigBodyStEndRes():',
     &  rigBodyStEndRes(i2+1),rigBodyStEndRes(i2+2)  
        end do!
        end if!control
c
	return
        end
