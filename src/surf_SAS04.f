c SAS:
c Solvent Acessible Surface         
c  + atomMinDistfromSAS + SASelectrostaticPot/Efield
c Yury N Vorobjev 2005
c
 
       subroutine surf_SAS04(atomXYZs,atomRads,
     &                 natoms,dprobe,dotden,
     &                 dotXYZ,dotnrm,dotarea,
     &                 dot_IATNUM,ndot,ndotMX,
     &                 nsurfAt,nsurfAtList,
     &                 atSurfAr,atSurfNrm,
     &                 atSurfXYZ,atSurfProj,bindSiteAt,
     &                 dot_jATNUM,
     &                 head_dotNum,linkListDotAt,  
     &                 atDistToSAS,atSASexp,
     &                 dot_eField,dot_ePot)
c
c nsurfAt : number of atoms onSAS
c nsurfAtList(is=1,nsurfAt) : list of atoms on SAS
c atSurfAr(ia) : SASarea for all atoms
c atSurfNrm(3*ia): SASarea normal for all atoms
c 
c INPUT:
c
cx	implicit none
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'pair1234array.h'
        include 'solvWBrg01.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        real atomXYZs(*)
        real atomRads(*)
        integer natoms
        real dprobe,dotden
        real dotXYZ(*),dotnrm(*),dotarea(*)
        integer dot_IATNUM(*),ndot,ndotMX
        integer nsurfAt,nsurfAtList(*)
        integer dot_jATNUM(*)         ! dot_jATNUM(id)=ja atomSphera intersects with ia sphera
        integer head_dotNum(*),linkListDotAt(*)
        real atSurfAr(*),atSurfNrm(*)
        real atSurfProj(*),atSurfXYZ(*)
        real bindSiteAt(*)
        real atDistToSAS(*)
        real atSASexp(*)
        real dot_eField(*),dot_ePot(*)
c Local:
        logical find12,find13
        integer ncallPL,make
        real rcut(3),rbuff
        integer nnbpLT,nnbpLMX 
c
        include 'sasDataSiz.h'
c
        include 'nbPairSAS03.h'
c local pairList dataStructure
c        integer nnbpLSAS03MAX
c        parameter (nnbpLSAS03MAX = natomsMAX*400)    ! 
c        integer nnbPairLSAS03(natomsMAX)
c        integer startnbPairLSAS03(natomsMAX)
c        integer nbpairListSAS03(nnbpLSAS03MAX)
c        real nbpairListSAS03D2(nnbpLSAS03MAX) 
c
c dots on the atom sphere
        real dotLocXYZ(3*ndotAtMAX)
        real dotLocAr(ndotAtMAX)
        real dotLocTet(ndotAtMAX)
        integer ndotLoc
        integer ndotAtMX
        integer i,ia,ia3
        integer id,id1,id3
        integer nd3,k,is
        integer ja,jn,ja3
        integer jAtNumbOL
        real diffdij2rj2,diffdij2rj2new,big
        real dixyz(3),aixyz(3),rda(3)
        real pi,rj2,qe
        real dij2,aij2
        real Ratav,pr,ri,rj
        real rjp,rjp2,dotSize
        real snrm,rimrj
        real ds,ds2,dSASmax 
        real cExp1_OPT
        real dotSizeScale_OPT
        real dePot,deField(3)
        real cReps
        integer epotVar_OPT
        logical atOnSurf,dotExist,atomijInside
        logical CONTROL,CONTROL0,CONTROL00
        integer kanalp
c
        CONTROL00 = .true. 
        CONTROL0 = .false.   ! write SAS to kanalpRunOut
        CONTROL = .false. 
        kanalp = kanalRunOut
c
        big = 1.0e10
        eLectScaleEFSMod=332.0
        epsEFSMod = 1.0  ! 80.0  ! dielConst WATer
        eLectScaleEFSMod=eLectScaleEFSMod/epsEFSMod
        cExp1_OPT = 4.0   ! 3.5     ! 2.0 huzhe
        cReps = 3.3   ! eps(r)=r*cReps ! SAS04 pot/efield
cx        epotVar_OPT = 1
        epotVar_OPT = 2
        pi = 3.1415927
        Ratav = 2.0        ! average atomic radius
        pr = dprobe*0.5
c        
        dotSizeScale_OPT = 1.25  ! 1.1=malo !
        dotSize = dotSizeScale_OPT/sqrt(dotden)
c        
        RcutEFSMAX=9.0   !  12.0       !maxRcut to calculate ePot/eField
c
        if(CONTROL0)then
        write(kanalp,*)'surf_SAS04: atXYZ,atRad'
        do ia = 1,natoms
        ia3=3*ia-3
        write(kanalp,7071)
     &  'ATOM  ',ia,'ATOM','RES ',' ',ia,(atomXYZs(ia3+k),k=1,3),
     &   atomRads(ia)
         end do!ia
        end if
c init
        do i=1,ndotMX
        dot_IATNUM(i) = 0
        dot_jATNUM(i) = 0
        end do !i
c
c neighbours list
c
        ncallPL = 0
        nnbpLMX = nnbpLSAS03MAX 
cx        rcut = 2.0*Ratav + dprobe
        rcut(1) = RcutEFSMAX
        rcut(2) = rcut(1)**2
        rcut(3) = rcut(2)*rcut(1) !
        dSASmax = rcut(1) 
c         
        call nonbondListFull(ncallPL,
     &          natoms,atomXYZs,rcut(1),
     &          nbpairListSAS03,startnbPairLSAS03,nnbPairLSAS03,
     &          nbpairListSAS03D2,nnbpLT,nnbpLMX)
c
c calculate SAS in long loop
c 
        ndotAtMX = ndotAtMAX
        ndot = 0
        nsurfAt = 0
c
        do ia = 1,natoms
c
        ia3 = 3*ia
        atSurfAr(ia) = 0.0
        head_dotNum(ia) = 0
c
        atSurfNrm(ia3-2) = 0.0
        atSurfNrm(ia3-1) = 0.0
        atSurfNrm(ia3)   = 0.0
c
        atSurfXYZ(ia3-2) = 0.0
        atSurfXYZ(ia3-1) = 0.0
        atSurfXYZ(ia3) = 0.0
c
        bindSiteAt(ia3-2) = 0.0
        bindSiteAt(ia3-1) = 0.0
        bindSiteAt(ia3)   = 0.0
c
        atOnSurf = .false.
c
        ri = atomRads(ia) + pr
        aixyz(1) = atomXYZs(ia3-2)
        aixyz(2) = atomXYZs(ia3-1)
        aixyz(3) = atomXYZs(ia3)
c
        ndotLoc = 4.0*pi*ri**2*dotden 
c
cx        write(kanalp,*)'surf_SAS:ia, pi,ri,dotden:'
cx     *  ,ia,atomRads(ia),pr,ri,dotden,ndotLoc
cx        write(kanalp,*)'surf_SAS: ia : ', ia,' xyz:',aixyz
c
	call genun01dot(ri,dotLocXYZ,dotLocAr,dotLocTet,
     &                  ndotLoc,ndotAtMX)
c
        do id1 = 1,ndotLoc 
            id3 = id1*3
            dixyz(1)=dotLocXYZ(id3-2) + aixyz(1)
            dixyz(2)=dotLocXYZ(id3-1) + aixyz(2)
            dixyz(3)=dotLocXYZ(id3) + aixyz(3)
c
            if(CONTROL)then
            write(kanalp,*)' id1:',id1,' dixyz:',dixyz
            end if
c loop over neihbour atoms
            dotExist = .true.
            jAtNumbOL = 0
            diffdij2rj2 = big
c
          if(CONTROL)then
        write(kanalp,*)
     &    'surf_SAS: ia,startnbPairLSAS03(ia),nnbPairLSAS03(ia):',
     &    ia,startnbPairLSAS03(ia),nnbPairLSAS03(ia)
         end if !C
c
        if(nnbPairLSAS03(ia) .ge. 1) then
       do jn=startnbPairLSAS03(ia),
     &              startnbPairLSAS03(ia)+nnbPairLSAS03(ia)-1        
           ja = nbpairListSAS03(jn) 
           ja3 = ja*3
           rj = atomRads(ja) + pr
           rj2 = rj**2
           rjp = rj + dotSize
           rjp2 = rjp**2
c
           if(CONTROL)then
           write(kanalp,*)'surf_SAS: ia,id1:',ia,id1,' jn,ja:',jn,ja
           end if
c
           dij2 = (dixyz(1)-atomXYZs(ja3-2))**2 +
     &            (dixyz(2)-atomXYZs(ja3-1))**2 +
     &            (dixyz(3)-atomXYZs(ja3))**2
c
           aij2 = (aixyz(1)-atomXYZs(ja3-2))**2 +
     &            (aixyz(2)-atomXYZs(ja3-1))**2 +
     &            (aixyz(3)-atomXYZs(ja3))**2
c
           rimrj=ri-rj
           if(aij2 .lt. rimrj**2)then
           atomijInside = .true.
           else
           atomijInside = .false.
           end if !inSide
c
           if(atomijInside )then
           if(ri .lt. rj)then
           goto 201  ! atom ia Inside ja
c                                     goto new atom
           else 
           goto 301   ! atom ja inside ia  goto new ja
           end if !ri .lt. rj
           end if !atomijInside
c
           if(CONTROL)then
           write(kanalp,*)'surf_SAS: ia,id1,jn,ja:',ia,id1,jn,ja
           write(kanalp,*)' aij2,dij2,ri-rj:',aij2,dij2,rimrj,
     &     ' atomijInside:',atomijInside
           end if
c
            if(dij2 .lt. rj2) then
            dotExist = .false.

            if(CONTROL)write(kanalp,*)'dot deleted'
c
            goto 101            ! take next dot
            end if ! dij2 .lt. rj2
c
c mark dot on intersection ia/ja spheres
c atoms ia,ja are not 12, 13 bonded
            if(dij2 .gt. rj2 .and. dij2 .lt. rjp2)then
            diffdij2rj2new = dij2 - rj2
            if(diffdij2rj2new .lt. diffdij2rj2)then !
c take nearest dot to intersection ia/ja spheres
            find12 = .false.
            find13 = .false.
            call compAtomPairList(ja,ia,
     &             pair12List,startPairL12,npairL12,find12)
            if (.not. find12)
     &      call compAtomPairList(ja,ia,
     &             pair13List,startPairL13,npairL13,find13)
c
             if(.not. find12 .and. .not. find13)then
c
            jAtNumbOL = ja
            diffdij2rj2 = diffdij2rj2new
c
            end if ! not in the 12 13 list
            end if !diffdij2rj2new .lt. diffdij2rj2
            end if !dij2 .gt. rj2 .and. dij2 .lt. rjp2
c
301     continue
        end do !jn
        end if ! jn
c collect surviving dots
        if(dotExist)then
        ndot = ndot + 1
c
        if(ndot .gt. ndotMX)then
        write(kanalp,*)
     &  'ERROR!: surf_SAS: ndotMAX is low',
     &  '  increase parameter sasDataSiz.h ',ndot
c
        write(kanalPStat,*)mError,
     & ' surf_SAS: ndotMAX is low increase parameter ',
     & ' in sasDataSiz.h '
c
        stop
        end if
c
        nd3=ndot*3
        dotXYZ(nd3-2) = dixyz(1)
        dotXYZ(nd3-1) = dixyz(2)                        
        dotXYZ(nd3)   = dixyz(3)                        
c
        dotnrm(nd3-2) = dotLocXYZ(id3-2)/ri
        dotnrm(nd3-1) = dotLocXYZ(id3-1)/ri
        dotnrm(nd3)   = dotLocXYZ(id3)/ri 
c
        dotarea(ndot) = dotLocAr(id1)
        dot_IATNUM(ndot) = ia
        dot_jATNUM(ndot) = jAtNumbOL
c
        atSurfAr(ia) = atSurfAr(ia) + dotLocAr(id1)
        atSurfNrm(ia3-2)=atSurfNrm(ia3-2) + 
     &                   dotnrm(nd3-2)*dotarea(ndot)
        atSurfNrm(ia3-1)=atSurfNrm(ia3-1) + 
     &                   dotnrm(nd3-1)*dotarea(ndot)
        atSurfNrm(ia3)=atSurfNrm(ia3) + 
     &                   dotnrm(nd3)*dotarea(ndot)
c
        atSurfXYZ(ia3-2) = atSurfXYZ(ia3-2) + 
     &                     dotXYZ(nd3-2)*dotarea(ndot)
        atSurfXYZ(ia3-1) = atSurfXYZ(ia3-1) +
     &                     dotXYZ(nd3-1)*dotarea(ndot)
        atSurfXYZ(ia3) = atSurfXYZ(ia3) +
     &                       dotXYZ(nd3)*dotarea(ndot)
c
        
        atOnSurf = .true.
c
        end if !dotExist
c
101     continue
c
        if(CONTROL)then
        write(kanalp,*)'dot id1:',id1,' EXIST:',dotExist, ' ndot;',ndot
        end  if !C
c
        end do !id1
c
        if(atOnSurf)then
        nsurfAt = nsurfAt+1
        nsurfAtList(nsurfAt) = ia
        atSurfNrm(ia3-2)=atSurfNrm(ia3-2)/atSurfAr(ia)
        atSurfNrm(ia3-1)=atSurfNrm(ia3-1)/atSurfAr(ia)
        atSurfNrm(ia3)=atSurfNrm(ia3)/atSurfAr(ia)
c
        atSurfXYZ(ia3-2) = atSurfXYZ(ia3-2)/atSurfAr(ia)
        atSurfXYZ(ia3-1) = atSurfXYZ(ia3-1)/atSurfAr(ia)
        atSurfXYZ(ia3)   = atSurfXYZ(ia3)/atSurfAr(ia)
c
        snrm = sqrt(atSurfNrm(ia3-2)**2+atSurfNrm(ia3-1)**2
     &              +atSurfNrm(ia3)**2)
c
        atSurfProj(ia) = snrm*atSurfAr(ia)  ! 
c     atomSurf area progected on plane perpendicular to average Normal vect
c
        bindSiteAt(ia3-2)=aixyz(1) + atSurfNrm(ia3-2)*rj/snrm
        bindSiteAt(ia3-1)=aixyz(2) + atSurfNrm(ia3-1)*rj/snrm
        bindSiteAt(ia3)  =aixyz(3) + atSurfNrm(ia3)*rj/snrm
c
        end if !atOnSurf
c
201     continue
 	end do !ia
c
c all SAS dots are collected
c
c linked List of dots on atom
      do id=1,ndot
      linkListDotAt(id) =  head_dotNum(dot_IATNUM(id))
      head_dotNum(dot_IATNUM(id))=id
      end do!id 
c
c atDistToSAS
c
       do ia = 1,natoms 
       ia3=3*ia-2
       atSASexp(ia) = 0.0
       if(atSurfAr(ia) .gt. 0.0 ) then
       atDistToSAS(ia) = 0.0 
c atSASexposure: 
       ri = atomRads(ia) + pr 
       atSASexp(ia)=sqrt(cExp1_OPT*atSurfAr(ia)/(4.0*pi*ri**2))
       if(atSASexp(ia) .gt. 1.0)atSASexp(ia) = 1.0
c
       else
c       
       atDistToSAS(ia) = dSASmax  ! MAXdist to SAS(xyz)
       if(nnbPairLSAS03(ia) .ge. 1) then
       do jn = startnbPairLSAS03(ia),
     &         startnbPairLSAS03(ia) + nnbPairLSAS03(ia)-1
          ja = nbpairListSAS03(jn)
          ja3=3*ja-2 
       if(atSurfAr(ja) .gt. 0.0 )then
       ds = sqrt((atomXYZs(ia3)-bindSiteAt(ja3))**2+
     &       (atomXYZs(ia3+1)-bindSiteAt(ja3+1))**2+
     &       (atomXYZs(ia3+2)-bindSiteAt(ja3+2))**2) - ri
       if(ds .lt. 0.0)ds = 0.0
c
       if(atDistToSAS(ia) .gt. ds)atDistToSAS(ia) = ds
c
       end if !atSurfAr(ja) .gt. 0.0       
       end do !jn
       end if !nnbPairLSAS03(ia) .ge. 1 
       end if !atSurfAr(ia) .gt. 0.0      
       end do ! ia
c
c calculate ePot/eField on dots
c
       do id1=1,ndot
        id3=3*id1-3
        ia = dot_IATNUM(id1)
        ia3=3*ia-3
c
        dixyz(1) = dotXYZ(id3+1)
        dixyz(2) = dotXYZ(id3+2)
        dixyz(3) = dotXYZ(id3+3)
c field from atom ia:
         qe=eLectScaleEFSMod*atomQ(ia)
         if(epotVar_OPT .eq. 1)
     &   call ePotEfieldSAS(atomXYZs(ia3+1),qe,dixyz,rcut,
     &                   dot_ePot(id1),dot_eField(id3+1))
c
         if(epotVar_OPT .eq. 2)
     &   call ePotEfieldSAS03(atomXYZs(ia3+1),qe,dixyz,rcut,cReps,
     &                   dot_ePot(id1),dot_eField(id3+1))
c
c field from neighbours of ia:
c
        if(nnbPairLSAS03(ia) .ge. 1) then
        do jn = startnbPairLSAS03(ia),
     &          startnbPairLSAS03(ia) + nnbPairLSAS03(ia)-1
           ja = nbpairListSAS03(jn)   ! atom j
           ja3 = ja*3-3
c
c
           if(atomQ(ja) .ne. 0.0 )then
           qe = eLectScaleEFSMod*atomQ(ja)
c
           if(epotVar_OPT .eq. 1)
     &     call ePotEfieldSAS(atomXYZs(ja3+1),qe,dixyz,rcut,
     &                   dePot,deField)
c
           if(epotVar_OPT .eq. 2)
     &      call ePotEfieldSAS03(atomXYZs(ja3+1),qe,dixyz,rcut,cReps,
     &                   dePot,deField)
c
           dot_ePot(id1)=dot_ePot(id1)+dePot          
c
           do k=1,3
           dot_eField(id3+k)=dot_eField(id3+k)+deField(k)
           end do !k 
c       
            if(CONTROL)then
            write(kanalp,*)'ja,qa:',ja,atomQ(ja),
     &      ' xyz:',(atomXYZs(ja3+k),k=1,3)
            write(kanalp,*)'id1:',id1,
     &      ' eField:',(dot_eField(id3+k),k=1,3)
            end if
           end if !atomQ(ja) .ne. 0.0 
c
           end do !jn
           end if !nnbPairLSAS03(ia) .ge. 1
c
           end do !id1 
c
c WRITE SAS:
c
        if(CONTROL0)then
        write(kanalp,*)'surf_SAS04:  '
        write(kanalp,*)'iDOT  iAt   XYZ  aRea nrmXYZ dot_jATNUM'
c
        do id1=1,ndot
        id3=3*id1
        write(kanalp,7004)
     &  'DOT   ',id1,'ATOM ',dot_IATNUM(id1),
     &  dotXYZ(id3-2),dotXYZ(id3-1),dotXYZ(id3),dotarea(id1),
     &  dotnrm(id3-2),dotnrm(id3-1),dotnrm(id3),dot_jATNUM(id1)
        end do!id1
c
        write(kanalp,*)'surf_SAS04: '
        write(kanalp,*)'iDOT  iAt   XYZ  ePot  eField'
c
        do id1=1,ndot
        id3=3*id1
        write(kanalp,7002)
     &  'ATOM  ',id1,'DOT ',dot_IATNUM(id1),
     &  dotXYZ(id3-2),dotXYZ(id3-1),dotXYZ(id3),dot_ePot(id1),
     &  dot_eField(id3-2),dot_eField(id3-1),dot_eField(id3)
        end do!id1
c
        write(kanalp,*)'SAS04: Atom on surface:'
        write(kanalp,*)
     &  '      isa       iat         atXYZ   atSurf(A^2)  atNrm atRad' 
        is = 0
	do ia = 1,natoms
c
        if(atSurfAr(ia) .gt. 0.0)then
        ia3=3*ia-3
        is = is + 1
c
        write(kanalp,7002)
     &  'ATOM  ',is,'ATOM',nsurfAtList(is),(atomXYZs(ia3+k),k=1,3),
     &   atSurfAr(ia),(atSurfNrm(ia3+k),k=1,3),atomRads(ia)
c
        end if
c
        end do !ia
c
        write(kanalp,*)'BindSiteAtom on surface:'
        write(kanalp,*)
     &  '      isa       iat         bindStXYZ   ' 
        is = 0
	do ia = 1,natoms
        ia3=3*ia-3
c
        if(atSurfAr(ia) .gt. 0.0) then
        is=is+1
c
        write(kanalp,7003)
     &  'ATOM  ',is,'BsAt',ia,(bindSiteAt(ia3+k),k=1,3),
     &   atSurfAr(ia),(atSurfNrm(ia3+k),k=1,3),
     &   'd,e:', atDistToSAS(ia),atSASexp(ia)
        else
c
         write(kanalp,7003)
     &  'ATOM  ',is,'SrAt',ia,(atomXYZs(ia3+k),k=1,3),
     &   atSurfAr(ia),(atSurfNrm(ia3+k),k=1,3),
     &   'd,e:', atDistToSAS(ia),atSASexp(ia)
c
        end if
        end do !ia
        end if !CONTROL0
c
       if(CONTROL00) write(kanalp,*)'surf_SAS04: finish ! ' 
c
7001   format(a6,i6,a18,3f8.3,f8.3,1x,f8.3)
7002   format(a4,1x,i6,2x,a4,5x,i4,4x,3f8.3,1x,f6.2,3f6.2,f5.2) ! 
7004   format(a4,1x,i6,2x,a4,5x,i4,4x,3f8.3,1x,f6.2,3f6.2,1x,i5)  
7003   format(a4,1x,i6,2x,a4,5x,i4,4x,3f8.3,1x,f6.2,3f6.3,1x,a4,1x,
     &        f5.1,1x,f4.2) ! PDB+
7071   format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3,2x,f6.3,f8.3) ! PDB
c
	return
	end
