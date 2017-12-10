c read restr1ResHWinfo
c reads info from file and calculates data for restrain1RHW.h
c
	subroutine readRestr1ResHW
     &             (init,restFile,nRestr1SegMX,nRestr1Seg,resEndRestr1)
c
c InPut:
c init =0/1 if =0 initialization to zero nRestr1Seg,resEndRestr1(*)
c restrFile - file name with input resr1ResCM information
c #restr residues (define restricted REsidues for REsid CMass)
c (a6,i4,i4)
c RESTRS   1  16
c RESTRS  47  59
c end
c -------------------
c nRestr1SegSegMX - MAX number of Segment of residues
c nrestSeg  - number of restrained Segment of residues
c OUT:
c resEndReRestr1(2*nRestr1Seg) - nkN,nkC - pairs of start/stop 
c                                res numbers from the res N and C ENDs
c
        implicit none
        character*(*) restFile
        integer init
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
        logical fext,CONTROL
c      
        kanalp=kanalRunOut 
        CONTROL = .true.
        keyw = 'RESTRS' 
        nlkey = 6
c     
        if(init .eq. 0)then
        nRestr1Seg = 0
        do i=1,nRestr1SegMX
        i2=2*i-2
        resEndRestr1(i2+1)=0
        resEndRestr1(i2+2)=0
        end do
        return
        end if !init
c
        write(kanalp,*)'readRestr1ResHW: restFile:',restFile
c
        kanalIn =  kanalInRestrALLType 
c
c* read in file  restr1A.inp
        inquire(file=restFile, exist=fext)
        if(.not. fext) then
        write(kanalp,*)
     &  'WARNING!:file (restrHWRes1.inp)',restFile,' does not exist '
     &  ,'ALL waterMol in hShell are restrained in 3.5A box by default'
c
cx         write(kanalPStat,'(a30,a20,a15)')
cx     &  'WARNING!:file (restrHWRes1.inp)',restFile,' does not exist'
cx         write(kanalPStat,*)
cx     &  ' ALL waterMolec in hydrShell are restrained by default'
c*
cx         nRestr1Seg = 0
cx         resEndRestr1(1)=0
cx         resEndRestr1(2)=0
         end if! .not.fext
c*       
         if(fext)then
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
        read(line, '(6x,2i4)')resEndRestr1(i2+1),resEndRestr1(i2+2)
        if(resEndRestr1(i2+1) .gt. 0 
     &  .and. resEndRestr1(i2+2) .gt. 0)nRestr1Seg = nRestr1Seg + 1
        end if !REST 
        goto 100

101     continue
        close(unit=kanalIn)
c
        end if !fext
c
        if(CONTROL)then
        write(kanalp,*)'readRestr1REsHW:nRestr1Seg:',nRestr1Seg      
        do i=1,nRestr1Seg
        i2=i*2-2
        write(kanalp,*)'resEndRestr1():',
     &  resEndRestr1(i2+1),resEndRestr1(i2+2)  
        end do!
        end if!control
c
	return
        end
