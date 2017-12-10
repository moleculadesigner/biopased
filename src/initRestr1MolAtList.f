c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2004                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c defines:  restrain1MHW.h  data
c restr1MHWatomList(ia) for a given list of nRestr1MHWSeg,resEndRestr1MHW(2*nRestr1Seg)
c
	subroutine initRestr1MHWAtList(nRestr1Seg,resEndRestr1,
     &                 nRestr1atom,restr1atomList,refAtomPos,
     &                 atMs1MHWatList)   
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
        integer nRestr1Seg,nRestr1atom
        integer resEndRestr1(*),restr1atomList(*)
        real refAtomPos(*)
        real atMs1MHWatList(*)
clocal
        integer kanalp
        integer i,ia
        integer ip,ip2,ips1,ips2
        integer i3,ia3
        real msf,ma
        logical CONTROL
c        
c initialize:
        kanalp = kanalRunOut
        CONTROL = .true.
c
c all realAtom(ia) [protCORE + LOOP] are moving 
         nRestr1atom=0 
         msf = 0.0   
         refAtomPos(1) = 0.0
         refAtomPos(2) = 0.0
         refAtomPos(3) = 0.0
c
         if(nRestr1Seg .ge. 1)then
         do ip = 1,nRestr1Seg
         ip2 = ip*2-1
         ips1 = resEndRestr1(ip2)
         ips2 = resEndRestr1(ip2+1)
         if(CONTROL)then
         if(ips1 .lt. 0 .or. ips1 .gt. nres 
     &      .or. ips2 .lt. 0 .or. ips2 .gt. nres )then
         write(kanalp,*)
     &  'ERROR!initRestr1MHWAtList: WRONG! resNumb of restSegm:',ip,
     &  ' resNst, resNfin:', ips1,ips2
         write(kanalp,*)' correct file: restr1Mol.inp ',
     &  ' are out of 1-nREs range !!! nres = ',nres
         stop
          end if !
          end if !Control
c         
         do ia = startAtInRes(ips1),stopAtInRes(ips2)
         nRestr1atom = nRestr1atom + 1
         restr1atomList(nRestr1atom) = ia            
         end do !ia
         end do !ip
c
c define refAtomPos()
        do i = 1,nRestr1atom
        i3 = i*3-3
        ia3 = 3*restr1atomList(i) - 3
        ma = atomMass(restr1atomList(i))
        refAtomPos(1) = refAtomPos(1)+atomXYZ(ia3+1)*ma
        refAtomPos(2) = refAtomPos(2)+atomXYZ(ia3+2)*ma
        refAtomPos(3) = refAtomPos(3)+atomXYZ(ia3+3)*ma
c
        msf = msf + ma
        atMs1MHWatList(i) = ma
c
        end do !i
c
        msf = 1.0/msf
        refAtomPos(1) = refAtomPos(1)*msf             
        refAtomPos(2) = refAtomPos(2)*msf           
        refAtomPos(3) = refAtomPos(3)*msf               
c
        do i = 1,nRestr1atom
        atMs1MHWatList(i) = atMs1MHWatList(i)*msf
        end do !i
c
        if(CONTROL)then
        write(kanalp,*)'initRestr1MolAtList: nRestr1Seg:',nRestr1Seg
        write(kanalp,*)'nRestr1atom :',nRestr1atom
        write(kanalp,*)'ref1MolPos :',
     &  refAtomPos(1),refAtomPos(2),refAtomPos(3)
        write(kanalp,*)'i    restr1MHWatList atMs1MHWatList: '
        do i = 1,nRestr1atom
        write(kanalp,'(2i6,f7.3)')i,restr1atomList(i),atMs1MHWatList(i)
        end do !i
c
        end if !CONTROL
c
        write(kanalp,*)'initRestr1MolAtList: Finish:'
c
        end if !nRestr1Seg .ge. 1
c
	return
        end
