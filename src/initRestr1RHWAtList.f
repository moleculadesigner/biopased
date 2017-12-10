c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2005                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c defines:  restrain1RHW.h  data
c harmonic wall restraints on REsidue CMass
c restr1MHWatomList(ia) for a given list of nRestr1MHWSeg,resEndRestr1MHW(2*nRestr1Seg)
c
	subroutine initRestr1RHWAtList(atomXYZs,
     &        nRestr1Seg,resEndRestr1Resid,
     &        nRestr1Resid,restr1ResidList,
     &        nRestr1RHWat,restr1RHWatList,
     &        restr1RrefResidPos,restr1ResidMass)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'filedat.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        real atomXYZs(*)
        integer nRestr1Seg,nRestr1Resid
        integer resEndRestr1Resid(*)
        integer restr1ResidList(*)
        real restr1ResidMass(*)
        real restr1RrefResidPos(*)
        integer nRestr1RHWat
        integer restr1RHWatList(*)
clocal
        integer kanalp
        integer i,ia,k
        integer ip,ip2,ips1,ips2
        integer i3,ia3,ir,n3
        real ma
        logical CONTROL
c        
c initialize:
        kanalp = kanalRunOut
        CONTROL = .true.
c
c init:                                             
         nRestr1Resid=0 
         nRestr1RHWat=0
c
         if(nRestr1Seg .ge. 1)then
         do ip = 1,nRestr1Seg
         ip2 = ip*2-1
         ips1 = resEndRestr1Resid(ip2)
         ips2 = resEndRestr1Resid(ip2+1)
         if(CONTROL)then
         write(kanalp,*)' resNst, resNfin:', ips1,ips2 
         if(ips1 .lt. 0 .or. ips1 .gt. nres 
     &      .or. ips2 .lt. 0 .or. ips2 .gt. nres )then
         write(kanalp,*)
     &  'ERROR!initRestr1MHWAtList: WRONG! resNumb of restSegm:',ip,
     &  ' resNst, resNfin:', ips1,ips2
         write(kanalp,*)' correct file: restr1.inp ',
     &  ' are out of 1-nREs range !!! nres = ',nres
         stop
          end if !
          end if !Control
c
         do ir = ips1,ips2
         nRestr1Resid = nRestr1Resid + 1   
         n3=3*nRestr1Resid-3      
         restr1ResidList(nRestr1Resid)=ir
         restr1ResidMass(nRestr1Resid)=0.0 
         restr1RrefResidPos(n3+1)=0.0 
         restr1RrefResidPos(n3+2)=0.0
         restr1RrefResidPos(n3+3)=0.0
         do ia = startAtInRes(ir),stopAtInRes(ir)
         nRestr1RHWat = nRestr1RHWat + 1
         restr1RHWatList(nRestr1RHWat) = ia
         ia3=3*ia-3
         ma = atomMass(ia)
         restr1ResidMass(nRestr1Resid)= restr1ResidMass(nRestr1Resid)+
     &        + ma            
         do k=1,3
         restr1RrefResidPos(n3+k)=restr1RrefResidPos(n3+k)+
     &      atomXYZs(ia3+k)*ma 
         end do !k 
         end do !ia
         do k=1,3
         restr1RrefResidPos(n3+k)=restr1RrefResidPos(n3+k)
     &                      /restr1ResidMass(nRestr1Resid)
         end do !k
c
         end do !ir
         end do !ip
c
        if(CONTROL)then
        write(kanalp,*)'initRestr1HWAtList: nRestr1Seg:',nRestr1Seg
        write(kanalp,*)'nRestr1HWResid :',nRestr1Resid
        write(kanalp,*)'restr1HWRefResidPos:'
        write(kanalp,*)'iRes   iGlob  Mass   X Y Z C.M. '
        do i = 1,nRestr1Resid
        i3=3*i-3
        write(kanalp,'(2i6,f7.3,1x,3f7.2)')i,restr1ResidList(i),
     &  restr1ResidMass(i),
     &  restr1RrefResidPos(i3+1),restr1RrefResidPos(i3+2),
     &  restr1RrefResidPos(i3+3)
        end do !i
c
        end if !CONTROL
c
        write(kanalp,*)'initRestr1RHWAtList: Finish:'
c
        end if !nRestr1Seg .ge. 1
c
	return
        end
