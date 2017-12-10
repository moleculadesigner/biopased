c MD atomicCoord and Velocities
        real atomXYZ0(3*natomMAX)        ! ti: current time
        real atomXYZp(3*natomMAX)        !  ti + dt
c        real atomXYZm(3*natomMAX)        !  ti - dt
        real atomXYZp05(3*natomMAX)       !  ti + dt/2
c atomVelocities
	real atomVel0(3*natomMAX)         !  ti        ! V0 = (V1+V2)/2
	real atomVelm(3*natomMAX)        !  ti - dt/2 
	real atomVelp(3*natomMAX)        !  ti + dt/2
c
	common/mdAtXYZ/atomXYZ0,atomXYZp
c
c MD  atomVelocities
	common/mdAtVel/atomVel0,atomVelm,atomVelp
c
cend
