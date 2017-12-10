c read restrAt1Inf
c modified verson: 2009
c each res-res segment has its own kharm constant
c
	subroutine readRestr1AtInf09
     &             (restFile,nRestr1SegMX,nRestr1Seg,resEndRestr1,
     &              resResSegKhConst,resRest1AtLine)
c
c InPut:
c restrfile - file name with input loopEndRes information
c #restr residues
c (6x,2i4,1x,f8.6,1x,a40)
c 123456789012345678901234567890...
c         r1- r2  Kharm   atomList
c xxxxxxIIIIiiiixffffffffxAAAAAAAAAA
c RESTAT   1   4  0.001   NAB         :NuclAcidBaseAtoms are Restrained
c RESTAT   1 219  1.0     PBB CA      :ProtBackBone CA
c RESTAT   1  16  0.5     ALL         : ALL atoms
c RESTAT  47  59  0.05               : ALL atoms
c end
c -------------------
c nRestr1SegSegMX - MAX number of Segment of residues
c nrestSeg  - number of restrained Segment of residues
c OUT:
c resEndReRestr1(2*nRestr1Seg) - nkN,nkC - pairs of numbers from the res N and C ENDs
c resResSegKhConst(nRestr1Seg) - KharmConstatn for atoms of the res-res Segment
c resRest1AtLine(iRestrSeg) : defines Names of Restricted atoms
c                           : ALL; PBB; PBB CA C O N
c                           : ALL - all atoms of residue are restrained
c                           : PBB :  backBone (all) are restrained
c                           : PBB CA : selected names in BackB restrained
c                           : PBB CA C O N  :selected names restrained
c
cx        implicit none
        include 'xyzPDBsize.h'
        character*(*) restFile
        integer nRestr1Seg,nRestr1SegMX
        integer resEndRestr1(*)
	real resResSegKhConst(*)
        character*40 resRest1AtLine(*)
c
        include "xyzPDBinfo.h"
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        character*80 line
        character*6 keyw
        integer i,i2,nline
        integer nlkey
        integer kanalp
        integer kanalIn
        integer nLtInLineMX
        logical fext
        logical CONTROL
c      
        nLtInLineMX = 40  !max letters in line to analysed
        kanalp=kanalRunOut
        CONTROL = .true.
        keyw = 'RESTAT' 
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
        write(kanalp,*)'readRestrAtInfo: restFile:',restFile
c
        kanalIn =  kanalInRestrALLType
c* read in file  restr1A.inp  
        inquire(file=restFile, exist=fext)
        if(.not. fext) then
        write(kanalp,'(a30,a20,a15)')
     &  'WARNING!:file (restrAt1.inp)',restFile,' does not exist'
c
         write(kanalPStat,'(a30,a20,a15)')
     &   'WARNING!: file (restrAt1.inp) ', restFile,' does not exist'
         write(kanalPStat,*)
     &   ' ALL ProtBackBone atoms are restrained by default' 
cx       stop
c default:
         nRestr1Seg = 1
         resEndRestr1(1) = 1
         resEndRestr1(2) = nres
         resRest1AtLine(nRestr1Seg) = "PBB "
	 resResSegKhConst(1) = 0.0001
        end if
c
        if(fext)then
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
cx        read(line, '(6x,2i4)')resEndRestr1(i2+1),resEndRestr1(i2+2)
cx        read(line, '(6x,2i4,1x,f8.6,1x,a40)')             !06.02.2015
        read(line, '(6x,2i4,1x,f8.6,1x,a12)')               !06.02.2015
     &	    resEndRestr1(i2+1),resEndRestr1(i2+2),
     &      resResSegKhConst(nRestr1Seg),resRest1AtLine(nRestr1Seg)
        end if !REST 
        goto 100

101     continue
        close(unit=kanalIn)
c
        end if !fext
c
        if(CONTROL)then
        write(kanalp,'(a31,1x,i2)')
     &  'readRestr1r1Inf:nRestr1Seg:',nRestr1Seg      
        do i=1,nRestr1Seg
        i2=i*2-2
        write(kanalp,'(a20,i4,i4,1x,f8.6,1x,a20)')'resEndRestr1():',
     &  resEndRestr1(i2+1),resEndRestr1(i2+2),
     &  resResSegKhConst(i),resRest1AtLine(i)  
        end do!
        end if!control
c
	return
        end
