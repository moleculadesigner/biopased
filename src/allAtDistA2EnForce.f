c calculate A2 dist RestrEng  forces for all atoms
c
	subroutine allAtDistA2EnForce(atomXYZ,
     &           natom,nRestrDistA2,restrDistA2List,
     &           restrDistA2DistHK,restrDistA2Eng,restrDistA2Force)

c defines DistRestr Eng  and forces on the set of moving atoms
c natom - total atoms in system
c nRestrDistA2: idist=1,nRestrDistA2
c restrDistA2List(2*nRestrDistA2) : at1,at2 
c restrDistA2DistHK(2*nRestrDistA2): dist^2  hKonst

	implicit none
        real atomXYZ(*)
        integer natom
        integer nRestrDistA2,restrDistA2List(*)
        real restrDistA2DistHK(*)
        real restrDistA2Eng,restrDistA2Force(*)
c local
        integer i,i2,i3,k
        integer ia1,ia2,ia13,ia23
        integer nbsum
        real ed,f1(3),f2(3),fs(3)
        real restrBondPar(2)
        real restrDistMin
        logical CONTROL,CONTROLx
        integer kanalp

c initialize 
         CONTROL = .false.
         CONTROLx = .false.
         kanalp = 6
c
         restrDistA2Eng = 0.0
         do i = 1,3*natom
         restrDistA2Force(i) = 0.0
         end do
c
         if(CONTROL)then
         write(kanalp,*)'allAtDistRestrEforce:nRestrDistA2:',
     &   nRestrDistA2
         end if !C
         
         if(nRestrDistA2 .ge. 1)then
c
c make calculation of force 
c
         restrDistMin = 0.1    !minimal a1-a2 restr distance
c
         if(CONTROLx)then
         write(kanalp,*)'allAtDistRestrEforce:atomXYZ:'
         do i=1,natom
         write(kanalp,'(a8,i6,3f12.3)')
     &   'atom:',i,(atomXYZ(3*i-3+k),k=1,3)
         end do!i
         end if !Cx

c loop over all restrDist
         do i = 1,nRestrDistA2
         i2 = 2*i
         ia1 = restrDistA2List(i2-1)
         ia2 = restrDistA2List(i2)
         ia13=3*ia1-3
         ia23=3*ia2-3
c
         if(CONTROL)then
         write(kanalp,*)'allAtDistA2EnForce: ia1,ia2:',
     &   ia1,ia2
         write(kanalp,'(a8,3f12.3)')'atom1:',(atomXYZ(ia13+k),k=1,3)
         write(kanalp,'(a8,3f12.3)')'atom2:',(atomXYZ(ia23+k),k=1,3)
         write(kanalp,*)'restrDistA2DistHK:',
     &   restrDistA2DistHK(i2-1),restrDistA2DistHK(i2)
         end if
c
         restrBondPar(2) = restrDistA2DistHK(i2-1)         ! dist0
         restrBondPar(1) = restrDistA2DistHK(i2)           ! harmConst
         if(restrBondPar(2) .lt. restrDistMin)
     &    restrBondPar(2) = restrDistMin
         restrBondPar(1) = restrBondPar(1)/(2.0*restrBondPar(2)**2) ! Kb
c
         call vbonddefenf(atomXYZ(ia13+1),atomXYZ(ia23+1),
     &        restrBondPar,ed,f1,f2)
c
         restrDistA2Eng = restrDistA2Eng + ed
         do k = 1,3
         restrDistA2Force(ia13+k)=restrDistA2Force(ia13+k) + f1(k) 
         restrDistA2Force(ia23+k)=restrDistA2Force(ia23+k) + f2(k) 
         end do !k
c
         end do !i
c
	if(CONTROLx)then
        write(kanalp,*)'allAtDistA2EnForce:restrDistA2Eng:',
     &   restrDistA2Eng
        write(kanalp,*)
     &  'allAtDistRestrEforce:f11,f2,f2 fs1,fs2,fs2 '
        fs(1)= 0.0
        fs(2) = 0.0
        fs(3) = 0.0
        do i = 1,natom
        i3=3*i-3
        fs(1) = fs(1) + restrDistA2Force(i3+1)
        fs(2) = fs(2) + restrDistA2Force(i3+2)
        fs(3) = fs(3) + restrDistA2Force(i3+3)
c
        write(kanalp,'(i6,2x,3f8.3,2x,3f8.3)')i,
     &  restrDistA2Force(i3+1),
     &  restrDistA2Force(i3+2),
     &  restrDistA2Force(i3+3),
     &  fs(1),fs(2),fs(3)
        end do!i
	end if ! Cx
c
        end if ! nRestrDistA2 .ge. 1
c
        if(CONTROL)then
        write(kanalp,*)'allAtDistA2EnForce: finis OK: '
        end if
c
	return
	end 
