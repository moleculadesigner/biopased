c SAS: ePotField,ePot
c
	subroutine ePotEfieldSAS(axyz,qe,dxyz,rcut,ep,ef)
c
c ePot, eField : on neutralbackground in sphere=rcut, eps=1
c
        implicit none
        real axyz(*)
        real dxyz(*)
        real qe
        real eP,eF(3)
        real rcut(3)
        real rda(3)
        real ds2,ds
        integer k
c
           ds2=0.0
           do k=1,3
           rda(k) = dxyz(k) - aXYZ(k)
           ds2=ds2+rda(k)**2
           end do!k
           if(ds2 .lt. rcut(2))then
           ds = sqrt(ds2)
           ep = qe*(1.0/ds+
     &        0.5*(ds2/rcut(2)-3.0)/rcut(1))
           do k=1,3
             ef(k)= qe*rda(k)*(1.0/(ds2*ds) - 1.0/rcut(3))
           end do !k
           end if !ds2 .lt. rcut(2)
c
         return
         end
c
        subroutine ePotEfieldSAS03(axyz,qe,dxyz,rcut,cReps,ep,ef)
c
c ePot, eField : on neutralbackground in sphere=rcut, eps=r*cReps
c
        implicit none
        real axyz(*)
        real dxyz(*)
        real qe
        real eP,eF(3)
        real rcut(3),cReps
c
        real rda(3)
        real ds2,ds,qep
        integer k
c
           ds2=0.0
           do k=1,3
           rda(k) = dxyz(k) - aXYZ(k)
           ds2=ds2+rda(k)**2
           end do!k
           if(ds2 .lt. rcut(2))then
           ds = sqrt(ds2)
           qep=qe/(ds*cReps)
           ep = qep*(1.0/ds+
     &        0.5*(ds2/rcut(2)-3.0)/rcut(1))
           do k=1,3
             ef(k)= rda(k)*qep*(1.0/(ds2*ds) - 1.0/rcut(3)) 
     &              + rda(k)*ep/ds2
           end do !k
           end if !ds2 .lt. rcut(2)
c
         return
         end
c
