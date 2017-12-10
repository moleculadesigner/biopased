c read restr1MolInfo
c reads inf from file and calculates data for restrain1HW.h
c
	subroutine readRestr1MolHW
     &             (restFile,nRestr1SegMX,nRestr1Seg,resEndRestr1)
c
c InPut:
c restrfile - file name with input loopEndRes information
c #restr residues (define restricted MOLECule)
c (a6,i4,i4)
c RESTML   1  16
c RESTML  47  59
c end
c -------------------
c nRestr1SegSegMX - MAX number of Segment of residues
c nrestSeg  - number of restrained Segment of residues
c OUT:
c resEndReRestr1(2*nRestr1Seg) - nkN,nkC - pairs of numbers from the res N and C ENDs
c
        implicit none
        character*(*) restFile
        integer nRestr1Seg,nRestr1SegMX
        integer resEndRestr1(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        character*80 line
        character*6 keyw
        integer nlkey
        integer i,i2,nline
        integer kanalp
        integer kanalIn
        logical CONTROL
c      
        kanalp = kanalRunOut
        CONTROL = .true.
        keyw = 'RESTML' 
        nlkey = 6
c     
        do i=1,nRestr1SegMX
        i2=2*i-2
        resEndRestr1(i2+1)=0
        resEndRestr1(i2+2)=0
        end do
c
        nRestr1Seg = 0
c 
        write(kanalp,*)'readRestr1MolHW: restFile:',restFile
c
        kanalIn =  kanalInRestrALLType
c* read in file         
        open(unit=kanalIn, file=restFile, form='formatted',
     &       status= 'old') 
c
	nline = 0
 100    read(kanalIn,'( a80 )', end=101 ) line
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        nline=nline+1
	if(line(1:nlkey) .eq. keyw)then
        i2 = 2*nRestr1Seg
        nRestr1Seg = nRestr1Seg + 1
        read(line, '(6x,2i4)')resEndRestr1(i2+1),resEndRestr1(i2+2)
        end if !REST 
        goto 100

101     continue
        close(unit=kanalIn)
c
        if(CONTROL)then
        write(kanalp,*)'readRestr1MHW:nRestr1Seg:',nRestr1Seg      
        do i=1,nRestr1Seg
        i2=i*2-2
        write(kanalp,*)'resEndRestr1():',
     &  resEndRestr1(i2+1),resEndRestr1(i2+2)  
        end do!
        end if!control
c
	return
        end
