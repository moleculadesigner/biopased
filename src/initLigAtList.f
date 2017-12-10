c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2003                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c    InPut:  nLigd, resStEndLig()
c    OUTput:
c           ligInfo.h data : nLigd,nAtomInLig,atomInLigList(ia) 
c                            atomLigXYZ()
c
	subroutine initLigAtList
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'movingAtom.h'
        include 'ligInfo.h' 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
clocal
        integer i,i3,i2,i1,j
        integer ia,ia3,ilg
        integer kanalp
        logical CONTROL,eStop   
c        
c initialize:
        kanalp = kanalRunOut
        CONTROL = .true.
        eStop = .false.
c init
        do i=1,natomMAX
        atInLigFlag(i) = 0
        end do !i
c
c define ligand atom list
c
         nAtomInLig = 0
         do ilg = 1,nLigd        
         i1 = resStEndLig(2*ilg-1)
         i2 = resStEndLig(2*ilg)
         if( i1 .le. i2 )then
         do ia = startAtInRes(i1),stopAtInRes(i2)
         nAtomInLig = nAtomInLig + 1
c
         if(nAtomInLig .gt. nAtomInLigMAX)then
	 write(kanalp,*)
     &  'ERROR! initLigAtLis: nAtomInLigMAX too LOW:',nAtomInLigMAX,
     &  ' increase nAtomInLigMAX in ligInfo.h file!'
         stop
         end if
c
         atomInLigList(nAtomInLig)=ia
         atInLigFlag(ia) = ilg
c 
         end do !ia
         end if
         end do !ilg
c
c define ligGeoCentr
        atomLigXYZcGeo0(1)=0.0
        atomLigXYZcGeo0(2)=0.0
        atomLigXYZcGeo0(3)=0.0
        if(nAtomInLig .ge. 1) then
        do i = 1,nAtomInLig
        i3 = i*3-3
        ia3 =  3*atomInLigList(i)-3
        atomLigXYZ(i3+1) = atomXYZ(ia3+1)
        atomLigXYZ(i3+2) = atomXYZ(ia3+2)
        atomLigXYZ(i3+3) = atomXYZ(ia3+3)
c
	atomLigXYZcGeo0(1)=atomLigXYZcGeo0(1)+atomLigXYZ(i3+1)
        atomLigXYZcGeo0(2)=atomLigXYZcGeo0(2)+atomLigXYZ(i3+2)
        atomLigXYZcGeo0(3)=atomLigXYZcGeo0(3)+atomLigXYZ(i3+3)
c
        end do!i
c
	atomLigXYZcGeo0(1)=atomLigXYZcGeo0(1)/nAtomInLig
        atomLigXYZcGeo0(2)=atomLigXYZcGeo0(2)/nAtomInLig
        atomLigXYZcGeo0(3)=atomLigXYZcGeo0(3)/nAtomInLig
        atomLigXYZcGeo(1)=atomLigXYZcGeo0(1)
        atomLigXYZcGeo(2)=atomLigXYZcGeo0(2)
        atomLigXYZcGeo(3)=atomLigXYZcGeo0(3)
c
        do i = 1,nAtomInLig
        i3 = i*3-3
        atomLigXYZ0c(i3+1) = atomLigXYZ(i3+1) - atomLigXYZcGeo0(1)
        atomLigXYZ0c(i3+2) = atomLigXYZ(i3+2) - atomLigXYZcGeo0(2)
        atomLigXYZ0c(i3+3) = atomLigXYZ(i3+3) - atomLigXYZcGeo0(3)
        end do !i
        end if !nAtomInLig .ge.1
c
        if(CONTROL)then
        write(kanalp,*)
     &  'initLigAtList: nAtomInLig:',nAtomInLig
        write(kanalp,'(a30)')'ATOM  Lig  XYZ         mvFlag:'
c
        do i = 1,nAtomInLig
        ia = atomInLigList(i)
        i3=3*i-3
c
        write(kanalp,7073)
     &     "ATOM  ",
     &     atomNumb(ia),atomName(ia),resName(ia),chName(ia),resNumb(ia),
     &     (atomLigXYZ(i3+j),j=1,3),atomQ(ia),moveFlag(ia)
c
        if( moveFlag(ia) .lt. 1)then  
        write(kanalp,*)'ERROR!! LigAtom is not in moveAtList!! ia:',ia
        eStop = .true.
        end if !
        end do ! i
        write(kanalp,'(a6,1x,i4,2x,a4,13x,3f8.3)')
     &  'ATOM  ',(ia+1),'LGct',atomLigXYZcGeo0
c
        write(kanalp,*)'Lig: atomLigXYZ0c: relative LigCnt:'
        do i = 1,nAtomInLig
        ia = atomInLigList(i)
        i3=3*i-3
c
        write(kanalp,7073)
     &     "ATOM  ",
     &     atomNumb(ia),atomName(ia),resName(ia),chName(ia),resNumb(ia),
     &     (atomLigXYZ0c(i3+j),j=1,3),atomQ(ia),moveFlag(ia)
c 
         end do !i
       
c
        end if !CONTROL
c
        write(kanalp,*)'initLigAtList: Finish:'
        if(eStop)then
        write(kanalp,*)'initLigAtList:eStop: due to ERROR!!'
        stop
        end if
c
	return
7073    format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,7x,f8.5,i5) ! PDB rqmv 
        end
