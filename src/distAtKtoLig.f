c 
c nearest distances from atom kP to Ligand atoms
c
	subroutine distAtKtoLigAt(kP,rKLig,iLm)
c
c rKlig - min dist; iLm - ligAt number
c
        include 'xyzPDBsize.h'
        include'xyzPDBcrd.h'
        include 'ligInfo.h'
c
        integer kP
        real rKLig
        integer iLm
c
        real d2,dm
	integer i,j,i3,k3      
        integer iL  
        integer kanalp
c
        kanalp = 6
c
        iLm=1
        dm=1.0e10
        k3=3*kP-3
c
c        write(kanalp,*)'distAtKtoLig: kP:',kP
c        
        do i=1,nAtomInLig
        iL=atomInLigList(i)
        i3=3*iL-3
c       
        d2=0.0
        do j=1,3
        d2=d2 + (atomXYZ(k3+j)-atomXYZ(i3+j))**2 
        end do !j
        if(d2 .lt. dm) then
        dm=d2
        iLm=i
        end if
c
        end do!i
c
        rKlig = sqrt(dm)
c
cx        write(kanalp,*)'distAtKtoLig: kP,iLm,rAtKLig:',kP,iLm,rKLig
c
	return
	end
