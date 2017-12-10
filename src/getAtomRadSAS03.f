c
	subroutine getAtomRadSAS03(natom,atomName,atomRad)
c
c assign approximate atomRad(*) for SAS03 calculation
c
        integer natom
        character*4 atomName(*)
        real atomRad(*)
c
        integer i
        real ratS03(8)
        data ratS03/1.90,0.75,2.05,1.65,1.65,2.0,2.0,3.50/
c identical to solvGSPar.dat
c
        do i=1,natom
        atomRad(i) = ratS03(1)    ! unknownAtom
        if(atomName(i)(1:1) .eq. 'H')atomRad(i) = ratS03(2)
        if(atomName(i)(1:1) .eq. 'C')atomRad(i) = ratS03(3) 
        if(atomName(i)(1:1) .eq. 'N')atomRad(i) = ratS03(4) 
        if(atomName(i)(1:1) .eq. 'O')atomRad(i) = ratS03(5) 
        if(atomName(i)(1:1) .eq. 'P')atomRad(i) = ratS03(6) 
        if(atomName(i)(1:1) .eq. 'S')atomRad(i) = ratS03(7) 
        if(atomName(i)(1:3) .eq. 'NA+')atomRad(i) = ratS03(8)
        end do !i
c
        return
        end
c
