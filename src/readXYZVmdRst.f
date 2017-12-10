c READ in PDBXYZV  to do mdRestart 
c atomXYZ() atomVel0 ()  : real*4
c      read BAD PDB, some residues (atoms) are missing, no XYZ
c
        subroutine readPDBxyzV01(pdbfile,natom,
     &       atomXYZ,atomVel0,nRecPdb,ntime0)
c
        implicit none
        character*(*) pdbfile
        integer natom
        real*4 atomXYZ(*),atomVel0(*)
        integer nRecPdb,ntime0
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local variables
        character*6 head
        character*30 head30
c
        integer i,j,k,k3  
        character*80 line
        integer kanalXYZV,kanalp
        logical readAtom
        logical CONTROL
c
        kanalXYZV=31
        kanalp=kanalRunOut
        CONTROL  = .true. 
c
        write(kanalp,*)'In readXYZVmdRst in file:',pdbfile 
c
       open(unit=kanalXYZV,file=pdbfile,form='formatted',
     &      status='old')
c
c READ PDB file and assign to variables
           k=0
c read XYZV
          rewind kanalXYZV

400       read(kanalXYZV,'(a80)',end=401)line
c     
           if(CONTROL)then
           write(kanalp,'(a79,a8)')line(1:79),' Origline'
           end if
c
           if(line(1:9) .eq. '$nRecPDB:')then
           read(line,'(9x,i7)')nRecPdb
           end if
c
           if(line(1:7) .eq. '$nstep:')then
           read(line,'(9x,i7)')ntime0
           end if
c
           readAtom=.false.
           if(line(1:4).eq.'ATOM' )then   !ATOM
           readAtom=.true. 
           end if !ATOM
c
           if(readAtom)then
           k=k+1
           k3=3*k-3
c
           if(k.gt.natom)then
           write(kanalp,*)
     &     'ERROR!: Too much atoms in PDBV mdRestart file'
           write(kanalp,*)'       Allowed atomMAX = ',natom
           write(kanalp,*)' Wrong mdXYZin.pdb file ! '
c
           write(kanalPStat,*)mError,' Wrong mdXYZin.pdb file ! '
c
           stop
           end if
  
           read(line,7072)head30,
     &     (atomXYZ(k3+j),j=1,3),(atomVel0(k3+j),j=1,3)
ccontrol
           if(CONTROL)then
           write(kanalp,7072)head30,
     &     (atomXYZ(k3+j),j=1,3),(atomVel0(k3+j),j=1,3) 
           end if !CONTROL

           end if !readatom
c
           goto 400

c end of file is reached:
401        continue
c
           if(natom .ne. k)then
           write(kanalp,*)
     &    'ERROR!! mdRestXYZVfile is notEqual molec.pdb records'
           stop
           end if
c
          close(kanalXYZV)
c
7071    format(a6,i5,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB
7072    format(a30,3f8.3,1x,3f8.3) !(i5,1x) !PDBVel
c
            return
            end
