c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2004                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c defines:
c startAtInRGB(*) stopAtInRGB(*) list of list of atoms in rigBodySegments
c
	subroutine initRigBodyAtList(nRigBody,rigBodyStEndRes,
     &      nAtRigBody,nAtRigBodySeg,startAtInRGB,stopAtInRGB,
     &      rigBodyFlag)
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
c
        integer nRigBody
        integer rigBodyStEndRes(*)
        integer nAtRigBody
        integer nAtRigBodySeg(*)
        integer startAtInRGB(*),stopAtInRGB(*)
        integer rigBodyFlag(*)
clocal
        integer kanalp
        integer i,ia
        integer ip,ip2,ips1,ips2
        integer i3,ia3
        integer nAtRigBodyAv
        logical CONTROL
c        
c initialize:
        kanalp = 6
        CONTROL = .true.
c
c all realAtom(ia) [protCORE + LOOP] are moving 
         do i=1,natomMAX
         rigBodyFlag(i)=0
         end do !i
c
         nAtRigBody=0    
c
         if(nRigBody .ge. 1)then
         nAtRigBodyAv = 0
         do ip = 1,nRigBody
         ip2 = ip*2-1
         ips1 = rigBodyStEndRes(ip2)
         ips2 = rigBodyStEndRes(ip2+1)
         if(CONTROL)then
         if(ips1 .lt. 0 .or. ips1 .gt. nres 
     &      .or. ips2 .lt. 0 .or. ips2 .gt. nres )then
         write(kanalp,*)
     &  'ERROR!initRigBodyAtList WRONG! resNumb of rigBody:',ip,
     &  ' resNst, resNfin:', ips1,ips2
         write(kanalp,*)' correct file: rigBodyRes.inp',
     &  ' are out of 1-nREs range !!! nres = ',nres
         stop
          end if !
          end if !Control
c        
         nAtRigBodySeg(ip)=stopAtInRes(ips2)-startAtInRes(ips1)+1 
         nAtRigBodyAv = nAtRigBodyAv + nAtRigBodySeg(ip)
         startAtInRGB(ip) = startAtInRes(ips1)
         stopAtInRGB(ip) = stopAtInRes(ips2)
c
         do ia = startAtInRGB(ip),stopAtInRGB(ip)
         nAtRigBody = nAtRigBody + 1
c rigBodyFlag
         rigBodyFlag(ia) = ip
         end do !ia
         end do !ip
c
         nAtRigBodyAv = nAtRigBodyAv/nRigBody
c correct atomMass(ia) out of RigBody to get comparable thermal move!
         do ia = 1,natom
         if(rigBodyFlag(ia) .eq. 0) 
     &   atomMass(ia) = atomMass(ia)*nAtRigBodyAv
        end do !ia
c
        if(CONTROL)then
        write(kanalp,*)'initRigBodyAtList: nRigBody:',nRigBody
        write(kanalp,*)'nAtRigBody Total :',nAtRigBody
        write(kanalp,*)'nAtRigBodyAv per RigBody:',nAtRigBodyAv
        write(kanalp,*)'i nAtRigBodySeg,startAtInRGB,stopAtInRGB:'
        do i = 1,nRigBody  
        write(kanalp,'(2i6,6x,2i6)') i,nAtRigBodySeg(i),
     &  startAtInRGB(i),stopAtInRGB(i)
        end do !i
c
        do ia = 1,natom
        write(kanalp,*)'initRigBodyAtL: rgbFlag = ',ia,rigBodyFlag(ia)
        end do !ia
c
        end if !CONTROL
c
        write(kanalp,*)'initRigBodyAtList: Finish:'
c
        end if !nRigBody .ge. 1
c
	return
        end
