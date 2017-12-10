c read LoopInfo or list of RES segments to move by MD/eOptim
c 
	subroutine readLoopInf(loopFile,nLoopMX,nLoop,resEndLoop)
c
c InPut:
c loopfile - file name with input loopEndRes information
c #inLOOP residues
c MOVRES  11  16
c MOVRES  47  59
c end
c -------------------
c nLoopMX - MAX number of Loops (resSegments)
c nLoop    - numbers of Loops
c OUT:
c resEndLoop(2*nLoop) - nkN,nkC - pairs of numbers from the res N and C ENDs
c
        implicit none
        character*(*) loopFile
        integer nLoop,nLoopMX
        integer resEndLoop(*)
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
        logical CONTROL
        logical fext
      
        kanalp=kanalRunOut 
        CONTROL = .true.
        keyw = 'MOVRES'
        nlkey = 6
c     
        do i=1,nLoopMX
        i2=2*i-2
        resEndLoop(i2+1)=0
        resEndLoop(i2+2)=0
        end do
c
        nLoop=0
c
c* read in file  loop.inp
        inquire(file=loopFile, exist=fext)
        if(.not. fext) then
        write(kanalp,*)
     &    'ERROR:UseINp movingRes inp file:',loopFile,' does not exist'
c
         write(kanalPStat,*)mError,
     &   ' movingRes inp file  ', loopFile,' does not exist'
 
        stop
        end if       
c
        open(unit=kanalInMoveRes, file=loopFile, form='formatted',
     &       status= 'old') 
c
	nline = 0
 100    read(kanalInMoveRes,'( a80 )', end=101 ) line
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        nline=nline+1
	if(line(1:nlkey) .eq. keyw)then
        i2 = 2*nLoop
        nLoop = nLoop + 1
        read(line, '(6x,2i4)')resEndLoop(i2+1),resEndLoop(i2+2)
        end if !LOOP 
        goto 100
c
101     continue
        close(unit=kanalInMoveRes)
c
        if(CONTROL)then
        write(kanalp,*)'readLoopInf:nLoop:',nLoop      
        do i=1,nLOOP
        i2=i*2-2
        write(kanalp,*)'resEndLoop():',
     &  resEndLoop(i2+1),resEndLoop(i2+2)  
        end do! i
        end if!control
c
	return
        end
