c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c*  correct atomic force to be zero summ over all atoms                      *
c*                                                                           *
c*     Yury Vorobjev 2004                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
        subroutine zeroAllAtomForce(nmoveatom,moveAtomList,
     &                              atomForce)
c
c correct all atomic Forces to give ZERO total summ
c calcuted for moveAtomList(im)=ia, im=1,nmoveatom  
c
	implicit none
        integer nmoveatom  
        integer moveAtomList(*)
        real atomForce(*)
c
        include "charStringSiz.h"
        include "output.h"
c
        real fsum(3),fcsum(3)
        integer i,i3,ia,ia3,k
        integer kanalp
        logical CONTROL
c
        kanalp =  kanalRunOut
        CONTROL = .false. 

c initialize
        do i=1,3
        fsum(i)=0.0
        end do!i
c
        if(nmoveatom .ge. 1) then
c
	do i = 1,nmoveatom
        i3 = 3*i - 3
        ia = moveAtomList(i)
        ia3 = 3*ia - 3
c
        fsum(1) = fsum(1) + atomForce(ia3+1)
        fsum(2) = fsum(2) + atomForce(ia3+2)
        fsum(3) = fsum(3) + atomForce(ia3+3)
c
        end do !i
c
        if(CONTROL)then 
        write(kanalp,'(a28,3f8.3)')'getZeroAllAtomForce01: fsum: ',fsum
        write(kanalp,*)'Corrected forces:'
        end if !C
c correct forces
        fsum(1)=fsum(1)/nmoveatom
        fsum(2)=fsum(2)/nmoveatom
        fsum(3)=fsum(3)/nmoveatom
c
        fcsum(1)=0.0
        fcsum(2)=0.0
        fcsum(3)=0.0
c
        do i = 1,nmoveatom
        i3 = 3*i - 3
        ia = moveAtomList(i)
        ia3 = 3*ia - 3
        atomForce(ia3+1) = atomForce(ia3+1)-fsum(1) 
        atomForce(ia3+2) = atomForce(ia3+2)-fsum(2) 
        atomForce(ia3+3) = atomForce(ia3+3)-fsum(3) 
c
        fcsum(1) = fcsum(1) + atomForce(ia3+1)
        fcsum(2) = fcsum(2) + atomForce(ia3+2)
        fcsum(3) = fcsum(3) + atomForce(ia3+3)
c
        if(CONTROL)then
        write(kanalp,'(2i5,1x,3f8.3,2x,3f8.3)')
     &  i,ia,(atomForce(ia3+k),k=1,3),fcsum
        end if !CONTROL
c
        end do ! i
c
        end if ! nmoveatom .ge. 1
c
	return
	end
c
