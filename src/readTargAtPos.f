c read targetAtomPositions from file
c
	subroutine readTargAtPos(pdbfile,natomMX,nresMX,
     &             natom,atomNumb,atomName,resName,resNumb,
     &             atomXYZ)
c
c natom - MX number of atoms, ATOM number as it is in the PDB file
c nres  - MX number of RESiduesNUMB line
c data for TargetAtPosition: atomNumb,atomName,resName,resNumb,atomXYZ
c                            as it is in the molec.pdb file
c
        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        character*(*) pdbfile
        integer natomMX,nresMX,natom
        character*4 atomName(*)
        character*4 resName(*)
        integer resNumb(*),atomNumb(*)
        real atomXYZ(*) 
c
        character*80 line
        character*6 head
        integer i,j,k,k3
        integer kanalXYZ,kanalp
        logical readAtom
        logical CONTROL,CONTROL0
c
        kanalXYZ=kanalInPdb
        kanalp=kanalRunOut
c
        CONTROL0 = .true.
        CONTROL  = .false.
c
       write(kanalp,*)'In readTargAtPos:ReadAtTargFile:',pdbfile
c
       open(unit=kanalXYZ,file=pdbfile,form='formatted',
     &      status='old' )
 
c READ PDB file and assign to variables
           k=0
c read XYZ
       rewind kanalXYZ
 
400     read(kanalXYZ,'(a80)',end=401)line
c
         readAtom=.false.
         if(line(1:3).eq.'END' .or. line(1:3).eq.'end')goto 401
         if(line(1:4).eq.'ATOM' )then   !ATOM
         readAtom=.true.
         if(CONTROL)then
           write(kanalp,'(a68,a12,a2,L1)')
     &     line(1:68),' OrigPDBline','R=',readAtom
           end if
         end if
c
         if(readAtom)then
         k=k+1
         k3=3*k-3 
c control
           if(k.gt.natomMX)then
           write(kanalp,*)'ERROR!: Too much atoms in PDB file'
           write(kanalp,*)'       Allowed atomMX = ',natomMX
           write(kanalp,*)'       increase parameter (atomMAX)...'
c
           write(kanalPStat,*)mError,
     &    ' (-c file) has to much atoms, increase parameter (atomMX)',
     &    ' in xyzPDBsize.h file '
c
          stop
          end if
c
           read(line,7071)head,
     &     atomNumb(k),atomName(k),resName(k),resNumb(k),
     &     (atomXYZ(k3+i),i=1,3) 
c
cx          goto 400
          end if !atomRead 
c
          goto 400
c end of file is reached:
401        natom=k
c
        close(kanalXYZ)
c
        if(CONTROL0)then
        write(kanalp,*)'readTargAtPos: targetAtPositions:'
        do k=1,natom
        k3=3*k-3
        write(kanalp,7071)'ATOM  ',
     &     atomNumb(k),atomName(k),resName(k),resNumb(k),
     &     (atomXYZ(k3+i),i=1,3)
        end do !k
        end if!C
c
	return
7071    format(a6,i5,2x,a4,a4,1x,i4,4x,3f8.3) !orig PDB 
	end
c
