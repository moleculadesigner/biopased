c calculate Val Bond defEn forces for LigList atoms
c
	subroutine allAtVBondEForceLig(atomXYZ,
     &           vFlag,bond12LigFlag,
     &           natom,bond12List,nbond12,nbond12noH,iShake,
     &           bond12ParL,eVbondDef,vbdefForce )
c
c defines total Evbonddef and forces on the set of moving atoms
c 
c natom - total atoms in system
c bond12List(2*ib) - bondPairList12; ib=1,nbond12
c                    bonds are collected if atoms 1or2 are in the movingAtList:
c                    nmoveatom,moveAtomList()
c bond12ParL(2*ib) - vbond Parameters :Kb,b0
c vFlag - flag defines details of calculations
c bond12LigFlag(*) - flag defines this bond12 will be processed if=1
c

	implicit none
        real atomXYZ(*)
        integer natom
        integer vFlag
        integer bond12LigFlag(*)
        integer bond12List(*),nbond12 
        integer nbond12noH,iShake     !nonH bonds; iShake flag
        real bond12ParL(*)
        real  eVbondDef, vbdefForce(*)
c local
        integer i,i2,i3,k,ia13,ia23
        integer nbsum
        real ed,f1(3),f2(3),fs(3)
        logical CONTROL,CONTROLx
        logical doProces
        integer kanalp
c
c initialize 
         CONTROL = .false.
         CONTROLx = .false.
         kanalp = 6
c
         eVbondDef = 0.0
         do i = 1,3*natom
         vbdefForce(i) = 0.0
         end do
c
         if(iShake .ge. 2) return
c
         nbsum = nbond12
         if(iShake .eq. 1)nbsum = nbond12noH
c
         if(CONTROLx)then
         write(kanalp,*)'allAtVgeoDefEforceLig: atomXYZ:'
         do i=1,natom
         write(kanalp,'(a8,i6,3f12.3)')
     &   'atom:',i,(atomXYZ(3*i-3+k),k=1,3)
         end do!i
         end if !Cx
c
c process bonds
c loop over all bonds
c
         if(vFlag .eq. 0)then
         doProces = .true.
         else
         doProces = .false.
         end if
c
         do i = 1,nbsum  
c
         doProces = .false.
         if(vFlag .eq. 0)doProces = .true.
c
         if(vFlag .eq. 1 .and.
     &     bond12LigFlag(i) .eq. 1) then
         doProces = .true.
         end if
c
         if(doProces)then
         i2 = 2*i
         ia13=3*bond12List(i2-1)-3
         ia23=3*bond12List(i2)-3
c
         if(CONTROLx)then
         write(kanalp,*)'vBondDefF: ia1,ia2 :',
     &    bond12List(i2-1),bond12List(i2)
         write(kanalp,'(a8,3f12.3)')'atom1:',(atomXYZ(ia13+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom2:',(atomXYZ(ia23+k),k=1,3)
         end if
c
         call vbonddefenf(atomXYZ(ia13+1),atomXYZ(ia23+1),
     &                         bond12ParL(i2-1),ed,f1,f2)
c
         eVbondDef = eVbondDef + ed
         do k = 1,3
         vbdefForce(ia13+k)=vbdefForce(ia13+k) + f1(k) 
         vbdefForce(ia23+k)=vbdefForce(ia23+k) + f2(k) 
         end do !k
c
c          if(CONTROL)then
c        write(kanalp,'(a25,3i6,2f8.3)')'vBondDefFLig :ed,eds:',
c     & i,ia13,ia23, ed,eVbondDef
c          end if
c
         end if ! doProces
c
         end do !i
c
	if(CONTROL)then
        write(kanalp,*)'vBondDefFLig :ed,eds, f1,f2,f2 fs1,fs2,fs2 '
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
c
c calculate Val Angl defEn forces for all atoms
c
	subroutine allAtVangEForceLig(atomXYZ,
     &           vFlag,trip123LigFlag,      
     &           natom,trip123List,nTrip123,ang123ParL,
     &           eVangDef,vAngdefForce )

c defines total Evbonddef and forces on the set of moving atoms
c natom - total atoms in system
c trip123List(3*ib) - tripletList12; it=1,nTrip123
c ang123ParL(2*ib) - vAng Parameters :Ka,a0
c
	implicit none
        real atomXYZ(*)
        integer natom
        integer vFlag
        integer trip123LigFlag(*)
        integer trip123List(*),nTrip123 
        real ang123ParL(*)
        real  eVangDef, vAngdefForce(*)
c local
        integer i,i2,k,ia13,ia23,ia33,i3
        real ed,f1(3),f2(3),f3(3),fs(3)
        logical CONTROL,CONTROLx
        logical doProces
        integer kanalp

c initialize 
         CONTROL = .false.
         CONTROLx = .false.
         kanalp = 6

c initialize 
         eVangDef = 0.0
         do i = 1,3*natom
         vAngdefForce(i) = 0.0
         end do
c
c loop over all triplets
c
         do i = 1,nTrip123
c
         doProces = .false.
         if(vFlag .eq. 0)doProces = .true.
c
         if(vFlag .eq. 1 .and.
     &      trip123LigFlag(i) .eq. 1) then
         doProces = .true.
         end if
c
         if(doProces)then
c
         i2 = 2*i
         i3 = 3*i
         ia13=3*trip123List(i3-2)-3
         ia23=3*trip123List(i3-1)-3
         ia33=3*trip123List(i3)-3
c
         if(CONTROLx)then
         write(kanalp,*)'vAnglDefFLig: ia1,ia2,ia3 :',
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
c
         end if ! doProces 
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
	subroutine allAtImpTEForceLig(atomXYZ,
     &           vFlag,qImp1234LigFlag, 
     &           natom,quarImp1234L,nImp1234,impAng1234ParL,
     &           eImpDef,impDefForce )

c defines total impAngdef and forces on the set of moving atoms
c natom - total atoms in system
c quarImp1234L(4*ib) - quartetsList12; ib=1,nImp1234
c impAng1234ParL(2*ib) - vbond Parameters :Kb,b0

	implicit none
        real atomXYZ(*)
        integer natom
        integer vFlag
        integer qImp1234LigFlag(*)
        integer quarImp1234L(*),nImp1234 
        real impAng1234ParL(*)
        real  eImpDef, impDefForce(*)
c local
        integer i,i2,k,ia13,ia23,i4,ia33,ia43,i3
        real ed,f1(3),f2(3),f3(3),f4(3),fs(3)
        logical CONTROL,CONTROLx
        logical doProces
        integer errorImp
        integer kanalp

c initialize 
         CONTROL = .false.
         CONTROLx = .false.
         kanalp = 6

c initialize 
         errorImp=0
         eImpDef = 0.0
         do i = 1,3*natom
         impDefForce(i) = 0.0
         end do
c
c loop over all bonds
         do i = 1,nImp1234
c
         doProces = .false.
         if(vFlag .eq. 0)doProces = .true.
c
         if(vFlag .eq. 1 .and.
     &     qImp1234LigFlag(i) .eq. 1) then
         doProces = .true.
         end if
c
         if(doProces)then
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
         end if !doProces
c
         end do !i
c
	if(CONTROL)then
        write(kanalp,*)'impAngDefFLig : f1,f2,f3  fs1,fs2,fs3 :'
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
	subroutine allAtTorsEForceLig(atomXYZ,
     &           vFlag,quar1234LigFlag,
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
        integer vFlag
        integer quar1234LigFlag(*)
        integer quar1234L(*),nQuar1234 
        real quar1234ParL(*)
        integer quar1234nPar(*)
        real  eTorsDef,torsAngForce(*)
c local
        integer i,i2,i3,i4,i16
        integer alarm
        integer ia13,ia23,ia33,ia43,k 
        real ed,f1(3),f2(3),f3(3),f4(3),fs(3)
        logical CONTROL,CONTROLx
        logical doProces 
        integer kanalp

c initialize 
         CONTROL = .false. 
         CONTROLx = .false.
         kanalp = 6

c initialize 
         eTorsDef = 0.0
         do i = 1,3*natom
         torsAngForce(i) = 0.0
         end do
c loop over all atom quartets : torsAngles 
         do i = 1,nQuar1234
c        
         doProces = .false.
         if(vFlag .eq. 0)doProces = .true.
c
         if(vFlag .eq. 1 .and.
     &     quar1234LigFlag(i) .eq. 1) then
         doProces = .true.
         end if
c
         if(doProces)then 
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
         write(kanalp,*)'allAtVgeoDefEforce:torsEng:'
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
         end if ! doProces
c
         end do !i
c
	if(CONTROL)then
        write(kanalp,*)'torsAngDefFLig : f1,f2,f3  fs1,fs2,fs3 :'
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
