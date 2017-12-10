c init PWatAtomTag 
c
	subroutine initPWatAtomTag
c
c READ in:  molecPDBxyz, etc.,          
c
c OUT: atomPWatTag()                     
c
c        implicit none
        include 'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBcrd.h'
        include 'xyzPDBinfo.h'
	include 'enForce_PW.h'
c local
        character*4 watResName,watAtName(4)
        integer i,kanalp
	logical CONTROL
c
        kanalp=kanalRunOut 
c
        CONTROL  = .false. !.true. 
c
c init/correct
c
        watResName = 'HOH '
        watAtName(1)='O   '
        watAtName(2)='H1  '
        watAtName(3)='H2  '
c
        do i=1,natomMAX
	atomPWatTag(i) = 0
	end do!
c
        if(CONTROL)then
	write(kanalp,*)'initPWatAtomTag: '
        write(kanalp,*)'ia,atomName,resName,atomPWatTag(ia)'
        end if!C
c
        natomSolutMol=0
	natomWSolv=0
        do i=1,natom
        if(resName(i) .eq. watResName)then
	atomPWatTag(i) = 2
	natomWSolv=natomWSolv+1
	else 
	atomPWatTag(i) = 1
	natomSolutMol=natomSolutMol+1
	end if
c
        if(CONTROL)then
	write(kanalp,*)i,atomName(i),resName(i),atomPWatTag(i)
	end if!C
	end do!i
c
        if(CONTROL)then
	write(kanalp,*)'initPWatAtomTag:natomSolutMol,natomWSolv:',
     &   natomSolutMol,natomWSolv
        end if
c
        if(CONTROL)STOP
c
        return 
c
        end
