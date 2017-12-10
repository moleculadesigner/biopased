c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2005                                               *
c                2005                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine initRigBodyVel(atXYZ,dtt)
c
c calculation of initial                  
c      Velocities due to rigBody rotation 
c IN:  atXYZ(*)[t=t],atVelM(*)[t=t-dtt], dtt=dt/2=MdtimeStep/2
c      dOmegaDt[t=t]  from rigBody01.h
c
c Out: atomVelP(*)[t=t+dt/2], absolute coords atXYZrot(*)[t=t+dt] 
c      atomVel0(*)[t=t]                       due to rigBrotation
c      
c        real rgbCmXYZ(3*nRigBodyMAX)        !rigBodycenter of Mass XYZ
c        real rgbMass(nRigBodyMAX)           !rigBodyMass
c        real rgbTensorInrt(9*nRigBodyMAX)   !tensorInercii
c        real rgbTensorInrt01(9*nRigBodyMAX) !invertedTensorInert
c        real rgbCmVel(3*nRigBodyMAX)        ! rigBodyCmassVelocity
c        real rgbCmForce(3*nRigBodyMAX)      ! rigBodyCmassForce=totForce
c        real rgbPmomVel(3*nRigBodyMAX)       ! momImpulse
c        real rgbOmegaRot(3*nRigBodyMAX)     ! rgbVectorOmega
c             rgbAtForceRot0(natomMAX)       ! centrostremitel'naja sila
cx:
cIn get02real rgbKmomForce(3*nRigBodyMAX)     ! momForce
cIn get02  real rgbOmegaRotDt(3*nRigBodyMAX)   ! rgb dOmega/dt05
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'rigBody01.h'
        include 'mdAtomXYZvel.h'
c
        real atXYZ(*)
        real dtt             !  md dtSmallstep/2
c local
        real Omega(3)
        real rrc(3),rrot(3),drrot(3)
        real vrot(3),vrotM(3)
        real x1(3),y1(3),z1(3)
        real v1(3),v2(3)
        real t01(9),t10(9)
        real t11(9),t00(9)
        real trLoc(9),trotF(9)
        real alfaRot,sc
        logical zeroRot,defLocXYZ1
c
        integer ir,ir3
        integer ia,ia3,k
        integer k1,k2 
c
        logical CONTROL
        integer kanalp
c
        CONTROL = .true.
        kanalp = 6
c
        if(nRigBody .eq. 0) return
c
         do ir = 1,nRigBody
         ir3 = 3*ir-3
c  omegaM= rgbOmega(ir)-dtt*dOmegaDt(ir)
         sc=0.0
         do k=1,3
         rgbOmegarotP(ir3+k)=rgbOmegaRot(ir3+k)+ 
     &                       dtt*rgbOmegaRotDt(ir3+k)
         omega(k)=rgbOmegaRot(ir3+k)-dtt*rgbOmegaRotDt(ir3+k) 
         rgbOmegarotM(ir3+k)=omega(k)
         sc=sc + omega(k)**2
         end do!k
c
         if(CONTROL)then
         write(kanalp,'(a25,f7.4)')'initRigBodyRotV:dtt=',dtt
         write(kanalp,'(a40,3f7.3)')
     &   'initRigBodyRotV: rgbOmegaRot:', (rgbOmegaRot(ir3+k),k=1,3) 
         write(kanalp,'(a40,3f7.3)')
     &   'initRigBodyRotV: rgbOmegaRotM:', (rgbOmegaRotM(ir3+k),k=1,3)
          write(kanalp,'(a40,3f7.3)')
     &   'initRigBodyRotV: rgbOmegaRotDt:',(rgbOmegaRotDt(ir3+k),k=1,3)
         end if
c
         sc=sqrt(sc)
         alfaRot = - sc*dtt       !angle of rotation= - dtt*Omega
         if(sc .gt. 0.0)then
         zeroRot = .false.
         do k=1,3
         z1(k)=omega(k)/sc        !normalized Omega
         end do!k
         else 
         zeroRot = .true.
         end if
c
         if(CONTROL)then
         write(kanalp,'(a40,i5,f7.3,2x,3f7.3)')
     &   'initRigBodyRotV: ir,OmegaMod,omgVect:',
     &   ir, sc, z1 
         end if!C
c
         defLocXYZ1 = .true.
         do ia = startAtInRGB(ir),stopAtInRGB(ir)
         ia3 = ia*3-3
c define rigBody Loc xyz1 system
c
         if(.not. zeroRot)then
          do k=1,3
          rrc(k) = atXYZ(ia3+k) - rgbCmXYZ(ir3+k)
          rrot(k) = rrc(k)
          drrot(k)=0.0
          vrot(k) = 0.0
          end do !k
c
         if(defLocXYZ1)then
c define local system
         call  vectorp(rrc,z1,x1)
c
         call vectorNrm1(x1)
c
         call vectorp(z1,x1,y1)
c rotationMatrix T01, T10
         if(CONTROL)then
         write(kanalp,*)'initRigBodyRotVel: x1:',x1
         write(kanalp,*)'initRigBodyRotVel: y1:',y1
         write(kanalp,*)'initRigBodyRotVel: z1:',z1
         end if
c
         do k=1,3
         t01(k)=x1(k)
         t01(3+k)=y1(k)
         t01(6+k)=z1(k)
c
         t10(3*k-2)=x1(k)
         t10(3*k-1)=y1(k)
         t10(3*k)  =z1(k)
         end do!k
c 
         if(CONTROL)then
         call matrix2p(t01,t10,t11)
         write(kanalp,*)'initRigBodyRotV: locXYZ: ir,ia:',ir,ia 
         write(kanalp,*)'initRigBodyRotV: RotMatrTest:'
         write(kanalp,'(a6,9f6.2)')'t01*t10:',t11
         end if
c
c  TrotLocM
         trLoc(9)=1.0
         trLoc(6)=0.0
         trLoc(3)=0.0
         trLoc(8)=0.0
         trLoc(7)=0.0
         trLoc(1)= cos(alfaRot)
         trLoc(2)=-sin(alfaRot)
         trLoc(4)=-trLoc(2)
         trLoc(5)= trLoc(1)
c
         call matrix2p(trLoc,t01,t00)
         call matrix2p(t10,t00,trotF)          !trotF(9) = fullRotMatrix
c
         if(CONTROL)then
         do k1=1,3
         v1(k1)=0.0
         do k2=1,3
         v1(k1)=v1(k1) + trotF(3*k1-3+k2)**2
         end do!k2
         end do !k1
         write(kanalp,*)'test:trotF(line2):',v1
         end if! C
c 
         defLocXYZ1 = .false.
         end if ! defLocXYZ1
c rotationMatrix is done!
c rotate coords
         call vectMatrx3Prod(rrc,trotF,rrot)
c
         end if !.not. zeroRot
c
         call vectorp(rgbOmegaRot(ir3+1),rrc,vrot)
         call vectorp(omega,rrot,vrotM)       ! rigBrotational Velocity
c
         do k=1,3
         atVelRigBodyRotM(ia3+k) = vrotM(k)         
         atVelRigBodyRot(ia3+k) = vrot(k)
         atVelRigBodyRotP(ia3+k) = vrot(k)
         atVelRigBodyTraM(ia3+k) = atVelRigBodyTra(ia3+k) 
     &     - dtt*rgbTransForce(ir3+k)/rgbMass(ir)
c
         atomVelm(ia3+k) = atVelRigBodyTraM(ia3+k)
     &                     + atVelRigBodyRotM(ia3+k)
         atomVel0(ia3+k) = atVelRigBodyTra(ia3+k) + 
     &                     atVelRigBodyRot(ia3+k)
         end do !k
c
         end do !ia
c end RigBodyLoop
         end do ! ir
c
         return
         end
