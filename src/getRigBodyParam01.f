c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2004                                               *
c                2005                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine getRigBodyParam01(doOmega,atXYZ,atomVel)
c
c INput: atomXYZ(*), atomVel(*)
c     
c 1)calculation ot the RigBody parameters: i.e. 
c M, RCM, VCM. inertiaTensor, Iinvert, OmegaRotVect(from vi)
c 2)convertion of InPut atomVel(*) --> atomVel(*)[rigidBody translation ]
c elimination of internal movement of atoms
c
c Out: atomVelRigTra(*)--> rigBody01.h                   
c      rigBody01.h  data structure
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
cIn get02  real rgbOmegaRotDt(3*nRigBodyMAX)   ! rgb dOmega/dt
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'rigBody01.h'
        include 'enForce.h'
        real atXYZ(*)
        real atomVel(*)
        integer doOmega
c
        integer ir,ir3,ir9
        integer ia,ia3,k
c
        real mInrt(9),mInrt01(9)
        real mInrtd(3)
        real detmInrt
        real v1(3),v2(3),v3(3)
        real mmp(9)
        real sc
        logical CONTROL
        integer kanalp
c
        CONTROL = .false.
        kanalp = 6
c
c rgbCenter, force and forceMom
        if(nRigBody .eq. 0) return
c
         do ir = 1,nRigBody
cinitialization
         if(CONTROL)then
         write(kanalp,*)'getRigBodyParam01: ir=',ir
         end if
c
         rgbMass(ir)=0.0
c
         ir3 = 3*ir-3
         do k=1,3
         rgbCmXYZ(ir3+k)=0.0
         rgbCmVel(ir3+k)=0.0
         rgbPmomVel(ir3+k)=0.0
         rgbOmegaRot(ir3+k)=0.0
cx         rgbOmegaRotDt(ir3+k)=0.0
         mInrtd(k)=0.0
         rgbTransForce(ir3+k)=0.0
         end do !k
c
         ir9 = 9*ir-9
         do k=1,9
         rgbTensorInrt(ir9+k)=0.0
         rgbTensorInrt01(ir9+k)=0.0  
         mInrt(k) = 0.0
         mInrt01(k) = 0.0
         end do!k9
c
c mass and cMassXYZ   
         do ia = startAtInRGB(ir),stopAtInRGB(ir)
         ia3 = ia*3-3
         rgbMass(ir) = rgbMass(ir) + atomMass(ia)
c
         do k=1,3
         rgbCmXYZ(ir3+k) = rgbCmXYZ(ir3+k)+
     &                     atomMass(ia)*atXYZ(ia3+k)
         rgbCmVel(ir3+k) = rgbCmVel(ir3+k) +
     &                     atomMass(ia)*atomVel(ia3+k)
c
         rgbTransForce(ir3+k) = rgbTransForce(ir3+k) +
     &                     atForceTot(ia3+k)
         end do !k
c
         end do !ia
c
         do k=1,3 
         rgbCmXYZ(ir3+k) = rgbCmXYZ(ir3+k)/rgbMass(ir)
         rgbCmVel(ir3+k) = rgbCmVel(ir3+k)/rgbMass(ir)       
         end do !k
c
c rgbMass and rgbCmXYZ, rgbCmVel(*) are done 
c     
c moment inertiaTensor
c
         do ia = startAtInRGB(ir),stopAtInRGB(ir)
         ia3 = ia*3-3
c
         mInrtd(1) = mInrtd(1) + atomMass(ia)*atXYZ(ia3+2)**2
c
         mInrtd(2) = mInrtd(2) + atomMass(ia)*atXYZ(ia3+3)**2
c
         mInrtd(3) = mInrtd(3) + atomMass(ia)*atXYZ(ia3+1)**2
c
         mInrt(2) = mInrt(2) - atomMass(ia)*
     &              atXYZ(ia3+1)*atXYZ(ia3+2)
         mInrt(3) = mInrt(3) - atomMass(ia)*
     &              atXYZ(ia3+1)*atXYZ(ia3+3)
         mInrt(6) = mInrt(6) - atomMass(ia)*
     &              atXYZ(ia3+2)*atXYZ(ia3+3)
c
c momVelocity
         do k=1,3
         v1(k) = atXYZ(ia3+k)-rgbCmXYZ(ir3+k)
         v2(k) = atomVel(ia3+k)-rgbCmVel(ir3+k)
         atomXYZrigBodyFrame(ia3+k) = v1(k)
         atVelRigBodyFrame(ia3+k) = v2(k)
         atVelRigBodyTra(ia3+k) = rgbCmVel(ir3+k)
         end do!k
c
         atomXYZrigBMod(ia)=sqrt(v1(1)**2+v1(2)**2+v1(3)**2)
c
         call vectorp(v1,v2,v3)
         do k=1,3
         rgbPmomVel(ir3+k) = rgbPmomVel(ir3+k) + 
     &   atomMass(ia)*v3(k)
         end do !k
c         
         end do !ia
c
         mInrt(1) = mInrtd(1) + mInrtd(2)
         mInrt(5) = mInrtd(2) + mInrtd(3)
         mInrt(9) = mInrtd(1) + mInrtd(3)
         mInrt(4) = mInrt(2)
         mInrt(7) = mInrt(3)
         mInrt(8) = mInrt(6)
c
c invert tensorInertia for RigBody      
c
         call invertM33(mInrt,mInrt01,detmInrt)
c
c keep tensorInertia
c
         do k=1,9
         rgbTensorInrt(ir9+k)=mInrt(k)
         rgbTensorInrt01(ir9+k)=mInrt01(k)
         end do !k
c
         if(CONTROL)then
      mmp(1)=mInrt(1)*mInrt01(1)+mInrt(2)*mInrt01(4)+mInrt(3)*mInrt01(7)
      mmp(5)=mInrt(4)*mInrt01(2)+mInrt(5)*mInrt01(5)+mInrt(6)*mInrt01(8)
      mmp(9)=mInrt(7)*mInrt01(3)+mInrt(8)*mInrt01(6)+mInrt(9)*mInrt01(9)
      mmp(2)=mInrt(1)*mInrt01(2)+mInrt(2)*mInrt01(5)+mInrt(3)*mInrt01(8)
      mmp(3)=mInrt(1)*mInrt01(3)+mInrt(2)*mInrt01(6)+mInrt(3)*mInrt01(9)
      mmp(8)=mInrt(7)*mInrt01(2)+mInrt(8)*mInrt01(5)+mInrt(9)*mInrt01(8)
      mmp(4)=mInrt(4)*mInrt01(1)+mInrt(5)*mInrt01(4)+mInrt(6)*mInrt01(7)
         mmp(6)=mmp(8)
         mmp(7)=mmp(3)
c
         write(kanalp,*)'getRigBodyParam01: mInrt*mInrt01:',mmp
         end if
c
       if(doOmega .eq. 0) goto 3001
c
c calculate OmegaVector for the rigBody number = ir
c 
         call vectMatrx3Prod(rgbPmomVel(ir3+1),mInrt01,
     &                              rgbOmegaRot(ir3+1))
c
         rgbOmegaRotMod(ir) = sqrt(rgbOmegaRot(ir3+1)**2+
     &       rgbOmegaRot(ir3+2)**2+rgbOmegaRot(ir3+3)**2)
c
         if(rgbOmegaRotMod(ir) .gt. 0.0)then
         sc=1.0/rgbOmegaRotMod(ir)
         else
         sc = 0.0
         end if
c
         do k=1,3
         rgbOmegaRot1(ir3+k)=rgbOmegaRot(ir3+k)*sc
         end do!k
c
         if(CONTROL)then
         write(kanalp,'(a32,f7.3,1x,a14,3f7.3)')
     &   'getRigBodyParam01: rgbOmegaRotMod:',rgbOmegaRotMod(ir), 
     &   ' rgbOmegaRot1:',(rgbOmegaRot1(ir3+k),k=1,3) 
         end if !C
c
c atomXYZrigBpar[allel]Omg(*) and atomXYZrigBppdOmg1(*) vektors
cx         do ia = startAtInRGB(ir),stopAtInRGB(ir)
cx         ia3 = ia*3-3
cx         sc = 0.0
cx         do k=1,3
cx         sc = sc + rgbOmegaRot1(ir3+k)*atomXYZrigBodyFrame(ia3+k)
cx         end do !k
c       
cx         do k = 1,3
cx         atomXYZrigBppdOmg1(ia3+k)=atomXYZrigBodyFrame(ia3+k) 
cx     &    - sc*rgbOmegaRot1(ir3+k)
cx         end do !k
c
cx         sc=sqrt(atomXYZrigBppdOmg1(ia3+1)**2+atomXYZrigBppdOmg1(ia3+2)**2+
cx     &           atomXYZrigBppdOmg1(ia3+3)**2)
c
cx         if(sc .gt. 0.0)then
cx         sc = 1.0/sc
cx         else
cx         sc = 0.0
cx         end if
cx         do k=1,3
cx         atomXYZrigBppdOmg1(ia3+k)=atomXYZrigBppdOmg1(ia3+k)*sc
cx         end do !k
c
cx         end do !ia Omega
c
3001     continue
c END of rigidBodyParameters: inerTiaTensor(), omega(*) etc.
c end RigBodyLoop
         end do ! ir
c
         return
         end
