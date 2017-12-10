c expLSASPotential data
c surface atoms of molecule P01=protein
c calculated by call  surf_SAS04 ( )
c 
cx        include 'ligInfoSiz.h'
c
        integer nsurfSASP01            ! number of SAS atoms of ProtMolec
c the next data are taken from dataSASdr.h data
        real atSurfArP01(natomMAX)     !SAS ProteinAtoms ias=1,nsurfSASP01
        real atSurfNrmP01(3*natomMAX)
        real atSurfXYZP01(3*natomMAX)
	integer nP01AtVsSASAtList(natomMAX)  !sas atom List points to GLobIndex ia of Prot atom List
c                                     ia=nP01AtVsSASAtList(isas)
c
        common/exWatSASP01/nsurfSASP01,atSurfArP01,atSurfNrmP01,
     &                    atSurfXYZP01,nP01AtVsSASAtList
c
c 
        integer ncall_exWSASPot     ! numbers of call of subroutine initExWatShRestrSASPot 
	real rcut_exWSASPot
c nbPairSAS03.h data structure 
c calculated by expLpotDataCalculate( )  
        integer natomExWatSh      !number of water atoms in exWatSell
	real atomXYZExWatSh(3*nAtomWatMAX)  !coord of atoms in exWatSell
        integer indxSASatNeighExWatSh(nAtomWatMAX) ! glob index=isas of SAS ProtAt nearest to exWatAtoms
c                                            from atSurfXYZP01( )
        real d2SASatNeighExWatSh(nAtomWatMAX) ! dist D2 from exWatAtom to nearest SASatom
c
        common/exWatSASP02/ ncall_exWSASPot,natomExWatSh,
     &    rcut_exWSASPot,atomXYZExWatSh
	common/exWatSASP03/indxSASatNeighExWatSh,d2SASatNeighExWatSh
c energy/forces of expL from protein SAS
        real KhSASexWatSh,ddSASexWatShMX
	common/exWatSASP04/KhSASexWatSh,ddSASexWatShMX
        real ddexWatSASat(nAtomWatMAX)
        real engSASexWatSh
	real fatSASexWatSh(3*nAtomWatMAX)
	common/exWatSASP05/ddexWatSASat,engSASexWatSh,fatSASexWatSh
	real dSOLVSheLLvar,dSOLVSHELLVolFix,areaProtSAS04
	common/exWatSASP06/dSOLVSheLLvar,dSOLVSHELLVolFix,areaProtSAS04
c
