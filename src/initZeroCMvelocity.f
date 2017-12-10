c make zeroVelocity for CM                   
c
c Yu. Vorobjev 2002
c
	subroutine initZeroCMVelocity(natom,atomMass,atomVel)
c
	implicit none
	integer natom
        real atomMass(*)
	real atomVel(*)
c
        include "charStringSiz.h"
        include "output.h"
c
        integer i,i3,k
        integer kanalp
        logical CONTROL
        real V(3),mtot
	real ta,ta0,sc
c
        kanalp = kanalRunOut
        CONTROL = .false. 
c
c check C.M. move and adjust
        V(1)=0.0
        V(2)=0.0
        V(3)=0.0
        mTot = 0.0
        ta0=0.0
	do i = 1,natom
        i3 = 3*i-3
        do k=1,3
        V(k) = V(k) + atomMass(i)*atomVel(i3+k)
        ta0 = ta0 + atomMass(i)*atomVel(i3+k)**2 
        end do !k
        mTot = mTot +  atomMass(i)
        end do!i
c
         do k=1,3
         V(k) = V(k)/mTot
         end do !
c
        do i = 1,natom
        i3 = i*3-3 
        do k=1,3
        atomVel(i3+k)=atomVel(i3+k) - V(k)                        
        end do !k       
        end do !i
c check velocities
        if(CONTROL)then
        write(kanalp,*)'initZeroVelocity:  VcM:',V 
        write(kanalp,*)'initZeroVelocity: corrected on cMVel: '
        write(kanalp,*)'   ia    vx     vy      vz    [A/ps]'
        end if!C
c scale CM corrected V to initial Temp
        V(1)=0.0
        V(2)=0.0
        V(3)=0.0
        ta = 0.0
        do i = 1,natom
        i3 = i*3-3
        do k=1,3
        V(k) = V(k) + atomMass(i)*atomVel(i3+k)
        ta = ta + atomMass(i)*atomVel(i3+k)**2
        end do !k
        end do !i
c
        do k=1,3
        V(k) = V(k)/mTot
        end do!k
c
        sc = sqrt(ta0/ta)
c
        if(CONTROL)then
        write(kanalp,*)'initZeroVelocity: final: VcM:',V
        write(kanalp,*)'initZeroVelocity: eKin0, eKin:', ta0,ta,
     &   ' ScailCf:',sc
        write(kanalp,*)'initZeroVelocity: VCmass: ', V    
        end if !C
c correct V by scaling
        do i = 1,3*natom
        atomVel(i) = atomVel(i)*sc
        end do !i
c
	return
	end
