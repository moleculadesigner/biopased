c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2004                                               *
c                2005                                               * 
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine getRigBdOmegaDt(atXYZ,atForce)
c
c calculate from detailed atomic forces atForce(*) dOmegaDt(ir) vectors
c atForce(rigBodyAtomList) = atForceRigBrot()+atForceRigBtra()
c
c USE data:
c        rigBody01.h  data structure
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
        real atForce(*)
        real atXYZ(*)
c local
cx        real rgbMomForce(3*nRigBodyMAX)
cx        real atForceRigBrot(3*natomMAX)   
cx        real atForceRigBtra(3*natomMAX)
c
        real scf
        integer ir,ir3,ir9
        integer ia,ia3,k
c
        real v1(3),v2(3),v3(3),vx(3)
        real sc
c
        scf = 418.40    ! scaling coef if force[Kcal/mol*A]
c                                   to get dOmg/dt in [rad/ps**2]
c rgbCenter, force and forceMom
        if(nRigBody .eq. 0) return
c
         do ir = 1,nRigBody
cinitialization
         ir3 = 3*ir-3
         ir9 = 3*ir3
c
         do k=1,3
         rgbTransForce(ir3+k)=0.0
         rgbMomForce(ir3+k)=0.0  
         end do !k
c    
         do ia = startAtInRGB(ir),stopAtInRGB(ir)
         ia3 = ia*3-3
         do k=1,3
         rgbTransForce(ir3+k) = rgbTransForce(ir3+k) +
     &                     atForce(ia3+k)
         end do !k
         end do !ia
c translationalForce= SUM of atForces = totalForce  
         sc = 1.0/rgbMass(ir)
         do k=1,3 
         v1(k) = rgbTransForce(ir3+k)*sc                     
         end do !k
c
c moment of force         
c
         do ia = startAtInRGB(ir),stopAtInRGB(ir)
         ia3 = ia*3-3
c 
         do k=1,3
         v2(k) = atForce(ia3+k) - v1(k)*atomMass(ia)             
         vx(k) = atXYZ(ia3+k) - rgbCmXYZ(ir3+k)    
         end do!k
c
 
         call vectorp(vx,v2,v3)
c
         do k=1,3
         rgbMomForce(ir3+k)=rgbMomForce(ir3+k)+v3(k)
         end do !k
         end do !ia
c momForce is done
c
c calculate dOmegaVector/dt for RigBody=ir 
c 
         call vectMatrx3Prod(rgbMomForce(ir3+1),
     &                       rgbTensorInrt01(ir9+1),
     &                       rgbOmegaRotDt(ir3+1))
c end RigBodyLoop
c        scf = 418.4       !to get dOmg/dt in rad/ps**2
         do k=1,3
         rgbOmegaRotDt(ir3+k) = rgbOmegaRotDt(ir3+k)*scf
         end do !k
c
         end do ! ir
c
         return
         end
