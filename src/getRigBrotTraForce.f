c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2004                                               *
c                2005                                               * 
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine getRigBRotTraForce(atForce)
c 
c calculation: rgbOmegaRotDt(*)
c convertion of detailed atomic forces atForce(*) to atForceRigBody(*)
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
c global
cx        real rgbTransForce(3*nRigBodyMAX)
cx        real rgbMomForce(3*nRigBodyMAX)
cx        real atForceRigBrot(3*natomMAX)   
cx        real atForceRigBtra(3*natomMAX)
c local
        integer ir,ir3,ir9
        integer ia,ia3,k
c
        real v0(3),v1(3),v2(3),v3(3)
        real sc
        logical CONTROL
        integer kanalp
c
        CONTROL = .false. 
        kanalp = 6
c rgbCenter, traforce and forceMom
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
c translationalForce=totalForce  
         sc = 1.0/rgbMass(ir)
         do k=1,3 
         v0(k) = rgbTransForce(ir3+k)*sc                     
         end do !k
c
c moment of force         
c
         do ia = startAtInRGB(ir),stopAtInRGB(ir)
         ia3 = ia*3-3
c 
         do k=1,3
         v1(k) = v0(k)*atomMass(ia)
         v2(k) = atForce(ia3+k) - v1(k)                 
         atForceRigBtra(ia3+k) = v1(k) 
         atForceRigBrot(ia3+k) = v2(k)
         end do!k
c
         call vectorp(atomXYZrigBodyFrame(ia3+1),v2,v3)
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
c
         if(CONTROL)then
         write(kanalp,*)'getRigBrotTraForce: ir:',ir
         write(kanalp,*)'getRigBrotTraForce: rgbTransForce:',
     &   (rgbTransForce(ir3+k),k=1,3) 
         write(kanalp,*)'getRigBrotTraForce:rgbMomForce:',
     &   (rgbMomForce(ir3+k),k=1,3)
         write(kanalp,*)'getRigBrotTraForce:rgbOmegaRotDt:',
     &   (rgbOmegaRotDt(ir3+k),k=1,3)
         write(kanalp,*)'getRigBrotTraForce:rgbTensorInrt:',
     &   (rgbTensorInrt(ir9+k),k=1,9)
          write(kanalp,*)'getRigBrotTraForce:rgbTensorInrt01:',
     &   (rgbTensorInrt01(ir9+k),k=1,9)
         end if !C
c end RigBodyLoop
         end do ! ir
c
         return
         end
