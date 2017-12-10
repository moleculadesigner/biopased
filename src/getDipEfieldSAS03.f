c calculate  eField on SAS dots correction from the WBridgeDipol
c            dot_eField_occFlag(*)-sterically occupied dots
c Y Vorobjev 2005
c
	subroutine getDipEfieldSAS03(wDipXYZ,wDipV,
     &               ndot,dotXYZ,dot_eField,dot_eField_occFlag,
     &               dotExclRad)
c
      	real wDipXYZ(*)
      	real wDipV(*)
        integer ndot
        real dotXYZ(*)
        real dot_eField(*)
        integer dot_eField_occFlag(*)
        real dotExclRad
c
c local
        real rdw(3),corEfield(3)
        real rd
        real wBrgMRad_OPT
	real rMX,rMX2
        real rMN,rMN2
        real eLectScaleC
c
        eLectScaleC=332.0/80.0 !regular
        eLectScaleC=0.0      !EfieldCorrfromWater = 0.0 
c
        wBrgMRad_OPT = dotExclRad + 0.2  ! 2.8 + 0.2
        rMN = wBrgMRad_OPT
        rMN2 = rMN**2
        rMX = rMN + 1.25*wBrgMRad_OPT
        rMX2=rMX**2
c
        do i=1,ndot
        i3=3*i-3
        d2 = 0.0
        rd = 0.0
        do k=1,3
        rdw(k)=dotXYZ(i3+k) - wDipXYZ(k)
        rd = rd + rdw(k)*wDipV(k)
        d2= d2 + rdw(k)**2  
        end do !k
        if( d2.le. rMN2)dot_eField_occFlag(i)=1  !occupated dots
c        
        if( d2.gt. rMN2 .and. d2.le. rMX2)then
        d = sqrt(d2)
        do k=1,3
        corEfield(k)=(3*rd*rdw(k)/d2 - wDipXYZ(k))/(d*d2)
        dot_eField(i3+k) = dot_eField(i3+k) + corEfield(k)*eLectScaleC
        end do !k 
        end if !do correction
        end do !i
c
	return
	end 

