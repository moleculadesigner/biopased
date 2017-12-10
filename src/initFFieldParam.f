c* initialize: loop Modeller                       
c
	subroutine initFFieldParam
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include 'xyzPDBcrd.h'
        include 'pair1234array.h'
        include 'nbondPairVCS.h'
        include 'filedat.h'
        include 'vdw12Par.h'
        include 'hbond128.h'
c
        integer nld,kanalp
        integer i,n,np12s
        logical OPT_defFFatomName
        logical CONTROL
        
c replica of array sizes for nonBondPairLists
        nnbpLVMX = nnbpLVMAX
        nnbpLCMX = nnbpLCMAX
        nnbpLSMX = nnbpLSMAX
c
        kanalp =  kanalRunOut
        CONTROL = .true.
c
c define file PATH/names 
c prepare amber ff parameters: allGeomDef, VDW, COUL (atomQ)
c
cx ffAtomTypeFile = dir_loopMod(1:nld)//'/dat/atmAAambff.dat'!defined in inputMDSApar.f
c
c get ff-atom code from atomNames
          OPT_defFFatomName = .false.  ! the ffAtomName() and atomQ() are defined
c                                      ! in the initMolecTopSeq01
         if(OPT_defFFatomName)then
	 call defFFatomName (ffAtomTypeFile,
     &              natom,atomNameEx,ResName,chName,
     &              ffAtomName,atomQ)
c
         end if !OPT_defFFatomName
c       
c define bondDef parameters for pair12List()
c
cx        ffParFile = dir_loopMod(1:nld)//'/dat/bsparBATV.dat'
c
        call getBondDefPar(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              bond12List,nbond12,bond12ParL)
c
c define valence angles def parameters
        call getVangDefPar(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              trip123List,nTrip123,ang123ParL)

c define Improper angle def parameters
        call getImpDefPar(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              quarImp1234L,nImp1234,impAng1234ParL)  
c
c define torsion parameters
        call getTorsPar(ffParFile,
     &        natom,atomNameEx,ResName,chName,ffAtomName,
     &        quar1234List,nQuar1234,quar1234ParL,quar1234nPar)
c
c assign atomMass and vdwParameters
        nVDWtypeMX = nVDWtypeMAX
        call  getVDWatMass(ffParFile,
     &              natom,atomNameEx,ResName,chName,ffAtomName,
     &              nVDWtypeMX,nVDWtype,atomVDWtype,
     &                                     atomVDW12ab,atomMass) 
c        
c hBond128
        nHXhb128MX = nHXhb128MAX
        nYYhb128MX = nYYhb128MAX
	call getHbond128Par(ffParFile,nHXhb128MX,nYYhb128MX,
     &             nHXhb128,nYYhb128,hb128HxList,hb128YYList,
     &             nHbond128Types,hB128pariHjY,hB128parRE,
     &             hb128parAB)
c all FField Parameters are defined
	call getHb128AtType
c
        write(kanalp,*)'initFFieldParam : DONE!'
c
	return
        end 
