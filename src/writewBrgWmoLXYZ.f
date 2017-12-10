c
	subroutine writewBrgWmoLXYZ
c
c write waterBridge molecules on molecule Surface
c
	include'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"

        include 'solvWBrg01.h'
c
        integer k,i,iw,iw3,i3
        integer kanalp,kanalXYZwbr
c
        kanalp = kanalRunOut
c
        write(kanalp,*)'WaterBridge molecules: strong bound water'
        write(kanalp,*)'            X Y Z        bindingEng(kcal/mol)'
c
        kanalXYZwbr = kanalwWatBrg
c
        open(unit=kanalXYZwbr,file=fileEFSwbrXYZ, form='formatted',
     &       status='unknown')
c
         write(kanalXYZwbr,*)'writeEFSolModWBrgXYZ: '
         do iw=1,3*nWBrgNow
          iw3=3*iw-3
          i3=(iw+2)/3
c
          write(kanalXYZwbr,7003)
     &   'ATOM  ',iw,wBrgAtomName(iw),wBrgResName(i3),i3,
     &   (wBrgWmoLXYZ(iw3+k),k=1,3),wBrgAtomQ(iw),wBrgEpotQ(i3)
cx     &   wBrgEpotQ(i3)*scaleWBrgEng(i3)
          end do !iw
c
        close(kanalXYZwbr)
c
	return
7003   format(a4,1x,i6,2x,a4,1x,a4,i4,4x,3f8.3,1x,f8.5,1x,f7.2,1x,f7.2) ! PDB
	end
