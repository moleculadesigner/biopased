c
c Yurii Vorobjev  2002
c
c initialization of the Gauss Shell Implicit solvation model
	subroutine initSolvatGSmod

        include 'xyzPDBsize.h'
        include 'filedat.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'pair1234array.h'
        include 'solvGSarray.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer i
        integer nld,kanalp
        real atRadH_OPT
        logical CONTROL
c
        kanalp = kanalRunOut
        CONTROL = .true.
        atRadH_OPT=1.1
c
cx        nld = ndirLett
cx         solvAtmTypeFile = dir_loopMod(1:nld)//
cx     &                     '/dat/solvGSPar_aa_amb.dat'
cx         solvAtmDatFile = dir_loopMod(1:nld)//
cx     &                     '/dat/solvGSPar.dat'
c
c define solvation atomic parameters
        call solvAtmParDef
     &             (solvAtmTypeFile,solvAtmDatFile,
     &              natom,atomNameEx,resName,nres,
     &              resNameRes,startAtInRes,
     &              chNameRes,atomQ,atmChargeNeutr,atomSolPar)
c
        do i=1,natom
        atomRad(i) = atomSolPar(6*i-1)
        if(atomName(i)(1:1) .eq. 'H' )atomRad(i) = atRadH_OPT
        end do!i
c
        write(kanalp,*)'initSolvatGSmod: Mod:ihydrPhSolv:',ihydrPhSolv
        write(kanalp,*)'initSolvatGSmod:atomHydPhSolv:',atomHydPhSolv
c
       call solvAtmParModif(ihydrPhSolv,
     &               natom,atomNameEx,resName,
     &               atomHydPhSolv,atomSolPar)
c
         return
         end
