c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c                                                                   *
c  Yury Vorobjev 2005                                               *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c modify velocities to get zero Translation and Rotation as a RigidBody
c
	subroutine initZeroRotationRigB(natoms,atXYZ,atVel)
c
c INput: atXYZ(*), atomVel(*)
c      : natoms,atomMass(*)  : from xyzPDBinfo.h
c     
c 1)calculation ot the RigBody parameters for the whole MOLECULE: i.e. 
c M, RCM, VCM. inertiaTensor, Iinvert, OmegaRotVect(from vi)
c 2)convertion of InPut atomVel(*) --> atomVel(*) - [rigidBody Rotation ]
c elimination of overal rotation of atoms
c
c Out: atomVel(*)--> atomVel(*) modified           
c
cIn get02real rgbKmomForce(3*nRigBodyMAX)     ! momForce
cIn get02  real rgbOmegaRotDt(3*nRigBodyMAX)   ! rgb dOmega/dt
c
        include 'xyzPDBsize.h'
        include 'xyzPDBinfo.h'
        include 'enForce.h'
        integer natoms
        real atXYZ(*)
        real atVel(*)
c
        include "output.h"
c local
        real atXYZcm(3*natomMAX)
c
cx        integer ir,ir3,ir9
        integer ia,ia3,k
c
        real rgbMass,rgbCmXYZ(3),rgbCmVel(3)
        real rgbPmomVel(3),rgbOmegaRot(3)
        real rgbTransForce(3)
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
        kanalp = kanalRunOut
c
c rgbCenter, force and forceMom
c
         rgbMass=0.0
c
         do k=1,3
         rgbCmXYZ(k)=0.0
         rgbCmVel(k)=0.0
         rgbPmomVel(k)=0.0
         rgbOmegaRot(k)=0.0
         mInrtd(k)=0.0
         end do !k
c
         do k=1,9
         mInrt(k) = 0.0
         mInrt01(k) = 0.0
         end do!k9
c
c mass and cMassXYZ   
         do ia = 1,natoms
         ia3 = ia*3-3
         rgbMass = rgbMass + atomMass(ia)
c
         do k=1,3
         rgbCmXYZ(k) = rgbCmXYZ(k)+
     &                     atomMass(ia)*atXYZ(ia3+k)
         rgbCmVel(k) = rgbCmVel(k) +
     &                     atomMass(ia)*atVel(ia3+k)
c
         rgbTransForce(k) = rgbTransForce(k) +
     &                     atForceTot(ia3+k)
         end do !k
c
         end do !ia
c
         do k=1,3 
         rgbCmXYZ(k) = rgbCmXYZ(k)/rgbMass
         rgbCmVel(k) = rgbCmVel(k)/rgbMass       
         end do !k
c
c rgbMass and rgbCmXYZ, rgbCmVel(*) are done 
c     
c moment inertiaTensor
c
         do ia = 1,natoms
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
         v1(k) = atXYZ(ia3+k)-rgbCmXYZ(k)
         atXYZcm(ia3+k) = v1(k)
         v2(k) = atVel(ia3+k)-rgbCmVel(k)     ! in CMass system
         atVel(ia3+k) = v2(k)
         end do!k
c
         call vectorp(v1,v2,v3)
c
         do k=1,3
         rgbPmomVel(k) = rgbPmomVel(k) + 
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
         write(kanalp,*)'initZeroRotationRigB:Inrt*mInrt01:',mmp
         end if
c
c calculate OmegaVector for the rigBody = ALL MOlecule
c 
         call vectMatrx3Prod(rgbPmomVel,mInrt01,
     &                              rgbOmegaRot)
c
         if(CONTROL)then
         write(kanalp,'(a32,f7.3,1x,a14,3f7.3)')
     &   'initZeroRotationRigB: rgbOmegaRot:',(rgbOmegaRot(k),k=1,3) 
         end if !C
c
c modify velocities to eliminate Wholl MOlecule rotation
c
          do ia = 1,natoms
          ia3=3*ia-3
          call vectorp(rgbOmegaRot,atXYZcm(ia3+1),v3)          
          do k=1,3
          atVel(ia3+k) = atVel(ia3+k) - v3(k)
          end do !k
          end do!ia
c
         return
         end
