c read LigInfo or list of RES segment of Ligand
c 
	subroutine readLigInf(ligFile,nLigMX,nLig,resEndLig)
c
c InPut:
c ligfile - file name with input ligEndRes information
c #inLIG residues
c LIGRES  11  16
c LIGRES  47  59
c end
c -------------------
c nLigMX - MAX number of Ligs (resSegments)
c nLig    - numbers of Ligands; Lig is defined as a set of consequtive REsidues
c
c OUT:
c resEndLig(2*nLig) - nkN,nkC - pairs of numbers from the res N and C ENDs
c
        implicit none
        character*(*) ligFile
        integer nLig,nLigMX
        integer resEndLig(*)
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
        integer kanalp,kanal11
        logical CONTROL
      
        kanalp=kanalRunOut 
        kanal11=kanalInRestrALLType
        CONTROL = .true.
        keyw = 'LIGRES'
        nlkey = 6
c     
        do i=1,nLigMX
        i2=2*i-2
        resEndLig(i2+1)=0
        resEndLig(i2+2)=0
        end do
c
        nLig = 0
c
c* read in file  ligRes.inp      
        open(unit=kanal11, file=ligFile, form='formatted',
     &       status= 'old') 
c
	nline = 0
 100    read(kanal11,'( a80 )', end=101 ) line
        if(line(1:3) .eq. 'end' .or.
     &      line(1:3) .eq. 'END')goto 101
c skip if comment characters ! or #
        if( line(1:1) .eq. '!' .or. line(1:1) .eq. '#')then
        goto 100
        end if
c
        nline=nline+1
	if(line(1:nlkey) .eq. keyw)then
        i2 = 2*nLig
        nLig = nLig + 1
c
        if(nLig .gt. nLigMX)then
        write(kanalp,*)'ERROR! readLigInf: nLigMAX is LOW!'
        stop
        end if
c
        read(line, '(6x,2i4)')resEndLig(i2+1),resEndLig(i2+2)
        end if !keyw 
        goto 100
c
101     continue
        close(unit=kanal11)
c
        if(CONTROL)then
        write(kanalp,*)'readLigInf:nLig:',nLig,'ligFile:',ligFile
        do i=1,nLig 
        i2=i*2-2
        write(kanalp,*)'resEndLig():',
     &  resEndLig(i2+1),resEndLig(i2+2)  
        end do!
        end if!control
c
	return
        end
