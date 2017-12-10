c random generator for fi,psi pair of aaResidues
c
c Yurii Vorobjev, 2002
c
	subroutine getFiPsi01(nresLoop,FiPsiLoop,iseed)
c
	integer nresLoop
        real FiPsiLoop(*)
        integer iseed
c local
        integer nRamMAX
        parameter (nRamMAX=10)
        real ramFiPsi(2*nRamMAX)
c
        real rv
        integer i,irv,nr
	data ramFiPsi /-120.0,120.0, -90.0,150.0, -60.0,120.0,
     &                 -60.0,-60.0,  -90.0,-30.0, -90.0,-15.0,
     &                 -60.0, 120.0, -90.0, 0.0,  -90.0,-15.0,
     &                 -90.0, -15.0 /
c
         nr = nresLoop + 2
	 do i = 1,nr
         CALL RANDOM(rv,iseed)
         irv = rv*nRamMAX + 1
         FiPsiLoop(2*i-1) = ramFiPsi(2*irv-1)
         FiPsiLoop(2*i) = ramFiPsi(2*irv)
         end do !i
c
	FiPsiLoop(1) = 0.0
c
	return
	end  

