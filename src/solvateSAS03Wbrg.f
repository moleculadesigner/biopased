c put ImplicitWater molec around SoluteMolecule 
c to the SAS points of Efield MAX
c Yury Vorobjev 2005
c
	subroutine solvateMoL03(atomXYZ)
c
c READ in:  molecPDBxyz, etc.,          
c
c dSOLVSheLL: solvationShell thickness
c
        include 'xyzPDBsize.h'
        real atomXYZ(*)
        include 'xyzPDBinfo.h'
        include 'coulEnPar.h'
        include 'dataSASdr.h'
        include 'solvWBrg01.h'
        include 'optionPar.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        real  eFieldMIN
        real efM2
        integer ip,iw,iw3,i3
        integer id3mx,k
        integer id1,id3
        integer ndotMX
        logical efMXfound
        integer kanalp
        logical CONTROL,CONTROL0
c
        kanalp=kanalRunOut 
        CONTROL  = .false.
        CONTROL0 = .false.
c OPT_ variables:
        dipMWEFS_OPT = 0.477  ! TIP3/SPC watDipMom
        eBindMIN_OPT = 0.60 !kcal,porogWBR(eps=1)! 4.0 OK'!
        eFieldMIN = eBindMIN_OPT/dipMWEFS_OPT 
cx        wBrgEpotQMAX_OPT = -0.60 !kcal/mol Wbrg ! MAX eng for WatBrg 
        wBrgEpotQMAX_OPT = -2.50 ! strong WatBrg in proteins
c
c calculate SAS and SAS dielectric model parameters
        ndotMX = ndotMAX
        do id1=1,ndotMX
        dot_eField_occFlag(id1)=0
        end do !id
c
	call getAtomRadSAS03(natom,atomName,atomRad) 
c
        ip = 1
c
cx        dotden = 4.0  !inputMDpar.f
        dprobe(ip) = solMolRad_OPT  ! water
c
c dot_SAS and dot_ePOT/eField
        if(.not.OPT_SolvGBorn .or. (nSAScall .eq. 0))then
c the SAS is alredy calculated in the bornRad subroutine
c 
               call  surf_SAS04(atomXYZ,atomRad,
     &                 natom,dprobe(ip),dotden,
     &                 dotXYZ(1,ip),dotnrm(1,ip),dotarea(1,ip),
     &                 dot_IATNUM(1,ip),ndot(ip),ndotMX,
     &                 nsurfAt(ip),nsurfAtList(1,ip),
     &                 atSurfAr(1,ip),atSurfNrm(1,ip),
     &                 atSurfXYZ(1,ip),atSurfProj(1,ip),
     &                 bindSiteAt01(1,ip),
     &                 dot_jATNUM(1,ip),
     &                 head_dotNum,linkListDotAt,
     &                 atDistToSAS,atSASexp,
     &                 dot_eField,dot_ePot)
c
          nSAScall = nSAScall + 1
c
          write(kanalp,*)
     &    'solvateSAS03Wbrg: surf_SAS04 call: nSAScall=',nSAScall
c
          end if !
c
cx         write(kanalp,*)
cx     & 'solvateMoL03: doneInitializ SAS ePotential for Solute MOLecule:'
c
c 2) calculate MAX of efield on the SAS from MOLECatoms
c 3) put water molecule to the MAX of e-field if efMAX > efPOROG
c 3b) calculate ePot in 4-point to orient water molecule
c 4) correct potential/efield on SAS taking the new WaterMolecule
c 4) goto 2
c
               nWBrgNow = 0
               iw = 1
1001           eFieldMX(iw) = 0.0
               ideFMX(iw)= 0
c
               if(CONTROL)then
               write(kanalp,*)
     &         'solvateMoL03: ndot(ip):',ip,ndot(ip),' tryNext iw:',iw
               end if!C
c
               efMXfound=.false.
c find MAX eField
               do id1=1,ndot(ip)
cx               if(dot_eField_occFlag(id1) .eq. 0)then
               if(dot_eField_occFlag(id1) .eq. 0 .and. 
     &            dot_jATNUM(id1,ip) .gt. 0) then
               id3=3*id1-3
               efM2=dot_eField(id3+1)**2+
     &              dot_eField(id3+2)**2+dot_eField(id3+3)**2
c         
               if(efM2 .gt. eFieldMX(iw))then
               eFieldMX(iw) = efM2
               ideFMX(iw) = id1
               efMXfound=.true.
               end if
c
               end if !occFlag
               end do!id1
c
	       if( .not. efMXfound)goto 1002
c               
               eFieldMX(iw) = sqrt(eFieldMX(iw))
c
               if(CONTROL)then
               write(kanalp,*)'solvateMoL03:nWBrgNow=',nWBrgNow,
     &         ' ideFMX(iw)=',ideFMX(iw),' eFieldMX(iw)=',eFieldMX(iw),
     &         ' eFieldMIN :',eFieldMIN
               end if !C
c
c put the implicit water molecule
               if(eFieldMX(iw) .gt. eFieldMIN)then   ! newWater
cx               nWBrgNow = iw
               iw3=iw*3-3
               id3mx=3*ideFMX(iw)-3
c
               do k=1,3
             wBrgXYZ(iw3+k)=dotXYZ(id3mx+k,ip)
c unitDipMomVect dipMWEFS_OPT
             wBrgDip(iw3+k)=dot_eField(id3mx+k)
     &                      /eFieldMX(iw)
               end do !k
c
               if(CONTROL)then
               write(kanalp,*)'solvateMoL03:nWBrgNow=',nWBrgNow,
     &         ' ideFMX(iw)=',ideFMX(iw),' eFieldMX(iw)=',eFieldMX(iw),
     &         ' eFieldMIN :',eFieldMIN, 
     &         ' ePotwbDip:',-eFieldMX(iw)*dipMWEFS_OPT 
               write(kanalp,7002)
     &   'ATOMwb',ideFMX(iw),'DOTw',dot_IATNUM(ideFMX(iw),ip),
     &   (wBrgXYZ(iw3+k),k=1,3),(wBrgDip(iw3+k),k=1,3) 
               write(kanalp,7002)
     &   'ATOMwb',ideFMX(iw),'H   ',dot_IATNUM(ideFMX(iw),ip),
     &   ((wBrgXYZ(iw3+k)+wBrgDip(iw3+k)*0.5),k=1,3),
     &    (wBrgDip(iw3+k),k=1,3)
              write(kanalp,7002)
     &   'ATOMwb',ideFMX(iw),'O   ',dot_IATNUM(ideFMX(iw),ip),
     &   ((wBrgXYZ(iw3+k)-wBrgDip(iw3+k)*0.5),k=1,3),
     &    (wBrgDip(iw3+k),k=1,3)
               end if !C
c
c do H1,H2 orientation
         
           call getWatHOrientSAS03(atomXYZ,ideFMX(iw),
     &           dotXYZ(id3mx+1,ip),wBrgDip(iw3+1),
     &           wBrgWmoLXYZ(9*iw-8),wBrgAtomQ(iw3+1),
     &           wBrgAtomName(iw3+1),wBrgResName(iw),wBrgEpotQ(iw))
c
c correct eField: 
c
               if(iw .lt.nWBrgMAX)then
c correct eField and calculate dot_eField_occFlag(*)
               call getDipEfieldSAS03(wBrgXYZ(iw3+1),wBrgDip(iw3+1),
     &                              ndot(ip),dotXYZ(1,ip),dot_eField,
     &         dot_eField_occFlag,dprobe(ip))
c
               if(wBrgEpotQ(iw) .lt. wBrgEpotQMAX_OPT)then
               nWBrgNow = iw 
               iw = iw+1
               if(CONTROL0)then
               write(kanalp,*)'solvateSAS03: accepted nWBrgNow=',
     &         nWBrgNow,' wBrgEpotQ(iwBr):',wBrgEpotQ(nWBrgNow)
               end if !C
               goto 1001
               else

               if(CONTROL0)then
               write(kanalp,*)'solvateSAS03: doNotAccepted nWBrgNow=',
     &         nWBrgNow, ' wBrgEpotQ(iw):',wBrgEpotQ(iw), 
     &         'wBrgEpotQMAX_OPT:',wBrgEpotQMAX_OPT
               end if !C
c
               goto 1001  ! gotoNextCandidate
cx             goto 1002  ! stop WB search
               end if  !wBrgEpotQ(iw) .lt. wBrgEpotQMAX_OPT
c
               else 
               write(kanalp,*)
     &         'solvateMoL03: nWBrgNow > nWBrgMAX:',nWBrgMAX
               write(kanalp,*)'ERROR:nWBrgMAX (solvWBrg01.h) is low !'
c
               write(kanalPStat,*)mError, 
     &           ' nWBrgMAX (solvWBrg01.h) is low'
               stop
               end if ! iw .lt.nWBrgMAX: goto newWBrgmolecule
c
               end if !eFieldMX(iw) .gt. eFieldMIN
c
1002           continue
c
               if(CONTROL0)then
          write(kanalp,*)'solvateMoL03: final watBridges XYZ:'
          do iw=1,3*nWBrgNow
          iw3=3*iw-3
          i3=(iw+2)/3
          write(kanalp,7003)
     &   'ATOMwb',iw,wBrgAtomName(iw),wBrgResName(i3),iw,
     &   (wBrgWmoLXYZ(iw3+k),k=1,3),wBrgAtomQ(iw)
               end do !iw
               end if!C
c
cx               write(kanalp,*)'solvateMoL03: Finish!'
c
        return
7002   format(a4,1x,i6,2x,a4,5x,i4,4x,3f8.3,1x,f6.2,3f6.2,f5.2) ! PDB
7003   format(a4,1x,i6,2x,a4,1x,a4,i4,4x,3f8.3,1x,f8.5) ! PDB
        end
