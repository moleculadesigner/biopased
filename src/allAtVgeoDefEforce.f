c calculate Val Bond defEn forces for all atoms
c
	subroutine allAtVBondEForce(atomXYZ,
     &           natom,bond12List,nbond12,nbond12noH,iShake,
     &           bond12ParL,eVbondDef,vbdefForce )

c defines total Evbonddef and forces on the set of moving atoms
c natom - total atoms in system
c bond12List(2*ib) - bondPairList12; ib=1,nbond12
c                    bonds are collected if atoms 1or2 are in the movingAtList:
c                    nmoveatom,moveAtomList()
c bond12ParL(2*ib) - vbond Parameters :Kb,b0

c	implicit none
        include 'xyzPDBsize.h'
        real atomXYZ(*)
        integer natom
        integer bond12List(*),nbond12 
        integer nbond12noH,iShake     !nonH bonds; iShake flag
        real bond12ParL(*)
        real  eVbondDef, vbdefForce(*)
c
c        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        include "enForce_PW.h"
	include "optionPar.h"
c local
        integer iPWtag1,iPWtag2
        integer i,i2,i3,k,ia13,ia23
        integer nbsum
        real ed,f1(3),f2(3),fs(3)
        logical CONTROL,CONTROLx
        integer kanalp

c initialize 
         CONTROL = .false.
         CONTROLx = .false.
c
         kanalp = kanalRunOut
c
         eVbondDef = 0.0
         do i = 1,3*natom
         vbdefForce(i) = 0.0
         end do
c PWat:
         do k=1,3
	 eGeoDef_PWat(k) = 0.0
	 eVbondDef_PWat(k) = 0.0
	 eVangDef_PWat(k) = 0.0
	 end do!k
c
         if(iShake .ge. 2) return
c
         nbsum = nbond12
         if(iShake .eq. 1)nbsum = nbond12noH
c
         if(CONTROL)then
         write(kanalp,*)'allAtVgeoDefEforce: atomXYZ:'
         do i=1,natom
         write(kanalp,'(a8,i6,3f12.3)')
     &   'atom:',i,(atomXYZ(3*i-3+k),k=1,3)
         end do!i
         end if !Cx
c loop over all bonds
         do i = 1,nbsum  
         i2 = 2*i
         ia13=3*bond12List(i2-1)-3
         ia23=3*bond12List(i2)-3
c
         call vbonddefenf(atomXYZ(ia13+1),atomXYZ(ia23+1),
     &                         bond12ParL(i2-1),ed,f1,f2)
c
         eVbondDef = eVbondDef + ed
c
c PWat:define eVbondDef_PWat():
         if(OPT_SolvateExWat)then 
	 iPWtag1=atomPWatTag(bond12List(i2-1))
	 iPWtag2=atomPWatTag(bond12List(i2))
	 call engPWatDistribute(iPWtag1,iPWtag2,ed,eVbondDef_PWat)
	 end if!OPT_SolvateExWat
c	 
	 do k = 1,3
         vbdefForce(ia13+k)=vbdefForce(ia13+k) + f1(k) 
         vbdefForce(ia23+k)=vbdefForce(ia23+k) + f2(k) 
         end do !k
c
          if(CONTROL)then
         write(kanalp,
     &   '(a9,3i4,1x,a5,f6.1,1x,f5.3,2(a3,3f6.1),a6,f8.3,f8.1)')
     &   'nb,i1,i2:',i,bond12List(i2-1),bond12List(i2),
     &   'bPar:',bond12ParL(i2-1),bond12ParL(i2),
     &   ' x1',(atomXYZ(ia13+k),k=1,3),' x2',(atomXYZ(ia23+k),k=1,3),
     &    ' e,eT:',ed,eVbondDef
          end if
c
         end do !i
c
	if(CONTROL)then
        write(kanalp,*)'vBondDefF :ed,eds, f1,f2,f2 fs1,fs2,fs2 '
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do i = 1,natom
        i3=3*i-3
        fs(1) = fs(1) + vbdefForce(i3+1)
        fs(2) = fs(2) + vbdefForce(i3+2)
        fs(3) = fs(3) + vbdefForce(i3+3)
        write(kanalp,'(i6,2x,3f8.3,2x,3f8.3)')i,
     &  vbdefForce(i3+1),
     &  vbdefForce(i3+2), vbdefForce(i3+3), fs(1),fs(2),fs(3)
        end do!i
	end if
c
	return
	end 
c calculate Val Angl defEn forces for all atoms
c
	subroutine allAtVangEForce(atomXYZ,
     &           natom,trip123List,nTrip123,ang123ParL,
     &           eVangDef,vAngdefForce )

c defines total Evbonddef and forces on the set of moving atoms
c natom - total atoms in system
c trip123List(3*ib) - tripletList12; it=1,nTrip123
c ang123ParL(2*ib) - vAng Parameters :Ka,a0
c
cx	implicit none
        include "xyzPDBsize.h"
c
        real atomXYZ(*)
        integer natom
        integer trip123List(*),nTrip123 
        real ang123ParL(*)
        real  eVangDef, vAngdefForce(*)
c
cx        include "charStringSiz.h" 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
        include "enForce_PW.h"
	include 'optionPar.h'
c local
        integer iPWtag1,iPWtag2
        integer i,i2,k,ia13,ia23,ia33,i3
        real ed,f1(3),f2(3),f3(3),fs(3)
        logical CONTROL,CONTROLx
        integer kanalp

c initialize 
         CONTROL = .false.
         CONTROLx = .false.
         kanalp = kanalRunOut
c initialize 
         eVangDef = 0.0
         do i = 1,3*natom
         vAngdefForce(i) = 0.0
         end do
c loop over all bonds
         do i = 1,nTrip123
         i2 = 2*i
         i3 = 3*i
         ia13=3*trip123List(i3-2)-3
         ia23=3*trip123List(i3-1)-3
         ia33=3*trip123List(i3)-3
c
         if(CONTROLx)then
         write(kanalp,*)'vAnglDefF: ia1,ia2,ia3 :',
     &    trip123List(i3-2),trip123List(i3-1),trip123List(i3)
         write(kanalp,'(a8,3f12.3)')'atom1:',(atomXYZ(ia13+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom2:',(atomXYZ(ia23+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom3:',(atomXYZ(ia33+k),k=1,3)
         end if
c
         call vangldefenf(atomXYZ(ia13+1),atomXYZ(ia23+1),
     &        atomXYZ(ia33+1),ang123ParL(i2-1),ed,f1,f2,f3)
c
         eVangDef = eVangDef + ed
         do k = 1,3
         vAngdefForce(ia13+k)=vAngdefForce(ia13+k) + f1(k) 
         vAngdefForce(ia23+k)=vAngdefForce(ia23+k) + f2(k) 
         vAngdefForce(ia33+k)=vAngdefForce(ia33+k) + f3(k) 
         end do !k
c PWat: 
c PWat:define eVangDef_PWat():
         if(OPT_SolvateExWat)then
         iPWtag1=atomPWatTag(trip123List(i3))
	 iPWtag2=atomPWatTag(trip123List(i3-1))
	 call engPWatDistribute(iPWtag1,iPWtag2,ed,eVangDef_PWat)
	 end if!OPT_SolvateExWat
c
         end do !i
c
	if(CONTROL)then
        write(kanalp,*)'vAngDefF : f1,f2,f3  fs1,fs2,fs3 :'
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do i = 1,natom
        i3=3*i-3
        fs(1) = fs(1) + vAngdefForce(i3+1)
        fs(2) = fs(2) + vAngdefForce(i3+2)
        fs(3) = fs(3) + vAngdefForce(i3+3)
        write(kanalp,'(i4,2x,3f8.3,2x,3f8.3)')i,vAngdefForce(i3+1),
     &  vAngdefForce(i3+2),vAngdefForce(i3+3),fs(1),fs(2),fs(3)
        end do!i
	end if
c
	return
	end 
c
c calculate ImpAng defEn and forces for all atoms
c
	subroutine allAtImpTEForce(atomXYZ,
     &           natom,quarImp1234L,nImp1234,impAng1234ParL,
     &           eImpDef,impDefForce )

c defines total impAngdef and forces on the set of moving atoms
c natom - total atoms in system
c quarImp1234L(4*ib) - quartetsList12; ib=1,nImp1234
c impAng1234ParL(2*ib) - vbond Parameters :Kb,b0

	implicit none
        real atomXYZ(*)
        integer natom
        integer quarImp1234L(*),nImp1234 
        real impAng1234ParL(*)
        real  eImpDef, impDefForce(*)
c
        include "charStringSiz.h" 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer i,i2,k,ia13,ia23,i4,ia33,ia43,i3
        integer errorImp
        real ed,f1(3),f2(3),f3(3),f4(3),fs(3)
        logical CONTROL,CONTROLx
        integer kanalp

c initialize 
         CONTROL = .false.
         CONTROLx = .false.
         kanalp = kanalRunOut

c initialize 
         errorImp = 0
         eImpDef = 0.0
         do i = 1,3*natom
         impDefForce(i) = 0.0
         end do
c loop over all bonds
         do i = 1,nImp1234
         i2 = 2*i
         i4 = 4*i
         ia13=3*quarImp1234L(i4-3)-3
         ia23=3*quarImp1234L(i4-2)-3
         ia33=3*quarImp1234L(i4-1)-3
         ia43=3*quarImp1234L(i4)-3
         call imprtorsanglenf(atomXYZ(ia13+1),atomXYZ(ia23+1),
     &                    atomXYZ(ia33+1),atomXYZ(ia43+1),
     &                 impAng1234ParL(i2-1),ed,f1,f2,f3,f4,errorImp)
c
         if(errorImp .eq. 1)then
         write(kanalp,*)'WARNING imprtorsanglenf ERROR! '
         write(kanalp,*)'atom 1,2,3,4:',
     &    quarImp1234L(i4-3),quarImp1234L(i4-2),
     &   quarImp1234L(i4-1),quarImp1234L(i4)
cX         stop
         end if
c
         eImpDef = eImpDef + ed
         do k = 1,3
         impDefForce(ia13+k)=impDefForce(ia13+k) + f1(k) 
         impDefForce(ia23+k)=impDefForce(ia23+k) + f2(k) 
         impDefForce(ia33+k)=impDefForce(ia33+k) + f3(k) 
         impDefForce(ia43+k)=impDefForce(ia43+k) + f4(k) 
         end do !k
c
         end do !i
c
	if(CONTROL)then
        write(kanalp,*)'impAngDefF : f1,f2,f3  fs1,fs2,fs3 :'
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do i = 1,natom
        i3=3*i-3
        fs(1) = fs(1) + impDefForce(i3+1)
        fs(2) = fs(2) + impDefForce(i3+2)
        fs(3) = fs(3) + impDefForce(i3+3)
        write(kanalp,'(i4,2x,3f8.3,2x,3f8.3)')i, impDefForce(i3+1),
     &  impDefForce(i3+2),impDefForce(i3+3),fs(1),fs(2),fs(3)
        end do!i
	end if
c
	return
	end 
c
	subroutine allAtTorsEForce(atomXYZ,
     &           natom,quar1234L,nQuar1234,
     &           quar1234ParL,quar1234nPar,
     &           eTorsDef,torsAngForce )

c defines total impAngdef and forces on the set of moving atoms
c natom - total atoms in system
c quar1234L(4*iq) - quartetsList12; iq=1,nQuar1234
c impAng1234ParL(16*iq) - torsParam, 4 param for 4 torsHarm

	implicit none
        real atomXYZ(*)
        integer natom
        integer quar1234L(*),nQuar1234 
        real quar1234ParL(*)
        integer quar1234nPar(*)
        real  eTorsDef,torsAngForce(*)
c
        include "charStringSiz.h" 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        integer i,i2,i3,i4,i16
        integer alarm
        integer ia13,ia23,ia33,ia43,k 
        real ed,f1(3),f2(3),f3(3),f4(3),fs(3)
        logical CONTROL,CONTROLx
        integer kanalp

c initialize 
         CONTROL = .false. 
         CONTROLx = .false.
         kanalp = kanalRunOut
c initialize 
         eTorsDef = 0.0
         do i = 1,3*natom
         torsAngForce(i) = 0.0
         end do
c loop over all atom quartets : torsAngles 
         do i = 1,nQuar1234
         i2 = 2*i
         i4 = 4*i
         i16 = 16*i-16
         ia13=3*quar1234L(i4-3)-3
         ia23=3*quar1234L(i4-2)-3
         ia33=3*quar1234L(i4-1)-3
         ia43=3*quar1234L(i4)-3
c
         if(CONTROL)then
         write(kanalp,*)'allAtVgeoDefEforce:torsEng:'
         write(kanalp,'(a10,i6,1x,a8,4i6)')
     &   'quartetN: ',i,'i-j-k-l:',(quar1234L(i4-4+k),k=1,4)
         end if
c
         if(CONTROLx)then
         write(kanalp,*)'vTorsDefF: ia1,ia2,ia3,ia4 :',
     &    quar1234L(i4-3),quar1234L(i4-2),quar1234L(i4-1),quar1234L(i4)
         write(kanalp,'(a8,3f12.3)')'atom1:',(atomXYZ(ia13+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom2:',(atomXYZ(ia23+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom3:',(atomXYZ(ia33+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom4:',(atomXYZ(ia43+k),k=1,3)
         end if
c
         call torsanglenf(atomXYZ(ia13+1),atomXYZ(ia23+1),
     &                    atomXYZ(ia33+1),atomXYZ(ia43+1),
     &         quar1234nPar(i),quar1234ParL(i16+1),ed,f1,f2,f3,f4,alarm)
c
         if(alarm .ge. 1)then
c
cx        iStatus = iStatus + iMolTopoMX
cx        write(kanalPStat,*)mStatus,iStatus,messageStatus
cx     &  ,' WARNING: torsForce: i-j-k-l in line', '...'
c
         write(kanalp,*)'allAtVgeoDefEforce:torsEng:Warning!:'
         write(kanalp,'(a10,i6,1x,a8,4i6)')
     &   'quartetN: ',i,'i-j-k-l:',(quar1234L(i4-4+k),k=1,4)
         write(kanalp,'(a8,3f12.3)')'atom1:',(atomXYZ(ia13+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom2:',(atomXYZ(ia23+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom3:',(atomXYZ(ia33+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom4:',(atomXYZ(ia43+k),k=1,3)
         end if !alarm
         
         eTorsDef = eTorsDef + ed
         do k = 1,3
         torsAngForce(ia13+k)=torsAngForce(ia13+k) + f1(k) 
         torsAngForce(ia23+k)=torsAngForce(ia23+k) + f2(k) 
         torsAngForce(ia33+k)=torsAngForce(ia33+k) + f3(k) 
         torsAngForce(ia43+k)=torsAngForce(ia43+k) + f4(k) 
         end do !k
c
         end do !i
c
	if(CONTROL)then
        write(kanalp,*)'torsAngDefF : f1,f2,f3  fs1,fs2,fs3 :'
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do i = 1,natom
        i3=3*i-3
        fs(1) = fs(1) + torsAngForce(i3+1)
        fs(2) = fs(2) + torsAngForce(i3+2)
        fs(3) = fs(3) + torsAngForce(i3+3)
        write(kanalp,'(i4,2x,3f8.3,2x,3f8.3)')i, torsAngForce(i3+1),
     &  torsAngForce(i3+2),torsAngForce(i3+3),fs(1),fs(2),fs(3)
        end do!i
	end if
c
	return
	end 
c
