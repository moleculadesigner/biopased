c calculate  eLectrPot on IonAS dots correction from the  new ion          
c            dot_ePot_occFlag(*)-sterically occupied dots
c Y Vorobjev 2005
c
	subroutine getIonEpotSAS03(ionXYZ,ionZZ,
     &               ndot,dotXYZ,dot_ePot,dot_ePot_occFlag)
c
c dot_ePot(*) - corrected ElectPotential on dot IonAS
c
      	real ionXYZ(*)
      	real ionZZ   
        integer ndot
        real dotXYZ(*)
        real dot_ePot(*)
        integer dot_ePot_occFlag(*)
c local
        real rdw(3),corEPot
        real rd
        real ionMRad_OPT
	real rMX,rMX2
        real rMN,rMN2
        real eLectScaleC
c
        eLectScaleC=332.0*ionZZ/80.0 !regular
c
        ionMRad_OPT = 4.2  ! ion-w-ion = ionExclusionRad
c
        rMN = ionMRad_OPT
        rMN2 = rMN**2
        rMX = 14.0         ! max electrostatic CUTOFF
        rMX2=rMX**2
c
        do i=1,ndot
        i3=3*i-3
        d2 = 0.0
        rd = 0.0
        do k=1,3
        rdw(k)=dotXYZ(i3+k) - ionXYZ(k)
        rd = rd + rdw(k)*ionXYZ(k)
        d2= d2 + rdw(k)**2  
        end do !k
        if(d2 .le. rMN2)dot_ePot_occFlag(i)=1  !occupated dots
c        
        if(d2 .gt. rMN2 .and. d2 .le. rMX2)then
        d = sqrt(d2)
        corEPot=eLectScaleC/d 
        dot_ePot(i) = dot_ePot(i) + corEPot
        end if !do correction
        end do !i
c
	return
	end 

