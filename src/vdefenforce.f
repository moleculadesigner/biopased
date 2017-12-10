c energy and force of deformation of valence structure
c Yury Vorobjev 2002
c
c gromos expression: E = (Kb/4)[b^2-b0^2]^2
c 
c amber Kharm can be transformed to Kb=Kh/(2*b0^2)
c
c bondPar(2)=b0,Kb
c
	subroutine vbonddefenf(xyz1,xyz2,bondPar,edef,f1,f2)      

	implicit none
        include "charStringSiz.h"
        include "output.h"
        real xyz1(3),xyz2(3)
        real f1(3),f2(3)
        real bondPar(2)
c        real b0,Kb
        real edef,b2,Kbb2
        integer kanalp
        logical CONTROL
c
c        kanalp = kanalRunOut
c        CONTROL = .false.
c
c        b0 = bondPar(2)
c        Kb = bondPar(1)
      
        f1(1)=xyz1(1)-xyz2(1) 
        f1(2)=xyz1(2)-xyz2(2) 
        f1(3)=xyz1(3)-xyz2(3) 
c
c        b2 = f1(1)**2+f1(2)**2+f1(3)**2 - b0**2
c
        b2 = f1(1)**2+f1(2)**2+f1(3)**2 - bondPar(2)**2
        Kbb2 = bondPar(1)*b2
        edef = Kbb2*b2*0.25
        f1(1)=-Kbb2*f1(1) 
        f1(2)=-Kbb2*f1(2) 
        f1(3)=-Kbb2*f1(3) 
        f2(1)=-f1(1)
        f2(2)=-f1(2)
        f2(3)=-f1(3)
c
c	if(CONTROL)then
c        write(kanalp,*)'vdefenforce: xyz1,xyz2, Kb,b0,ed,f:'
c        write(kanalp,'(3f8.3,2x,3f8.3,2x,2f8.3,2x,f8.2,1x,3f8.3)')
c     &  xyz1(1),xyz1(2),xyz1(3),xyz2(1),xyz2(2),xyz2(3),
c     &  bondPar(1),bondPar(2),
c     &  edef,f1(1),f1(2),f1(3)
c        end if

	return
        end
c
c angle deformation
c gromos function: edef = 0.5*Kang*(cosA-cos0)^2
c
c                  Kharm=Kang*(sinAng0)^2
c
	subroutine vangldefenf(xyz1,xyz2,xyz3,angPar,
     &                         edef,f1,f2,f3)

        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        real xyz1(3),xyz2(3),xyz3(3)
        real f1(3),f2(3),f3(3)
        real angPar(2),fs(3)
        real cos0,Kang,edef
c
        real rij(3),rkj(3)
        real rij1(3),rkj1(3)
        real dij,dkj,cost,cost0,kcos 
        real aij,akj
        integer k
        integer kanalp
        logical CONTROL

        kanalp = kanalRunOut
        CONTROL = .false.
c
c init
        edef = 0.0
        do k=1,3
        f1(k) = 0.0
        f2(k) = 0.0
        f3(k) = 0.0
        end do !k
c
        cos0 = angPar(2)
        Kang = angPar(1)
c
	rij(1)=xyz1(1)-xyz2(1)
	rij(2)=xyz1(2)-xyz2(2)
	rij(3)=xyz1(3)-xyz2(3)
	rkj(1)=xyz3(1)-xyz2(1)
	rkj(2)=xyz3(2)-xyz2(2)
	rkj(3)=xyz3(3)-xyz2(3)
        dij = sqrt(rij(1)**2+rij(2)**2+rij(3)**2)
        dkj = sqrt(rkj(1)**2+rkj(2)**2+rkj(3)**2)
c
        if(dij .eq. 0.0 .or. dkj .eq. 0.0 ) then
        write(kanalp,*)'vangdefenforce: WARNING! dij,dkj = zero:'
     &  ,dij,dkj
        write(kanalp,*)'atom1:',xyz1
        write(kanalp,*)'atom2:',xyz2
        write(kanalp,*)'atom3:',xyz3
        return
        end if
c
        do k=1,3
        rij1(k)=rij(k)/dij
        rkj1(k)=rkj(k)/dkj
        end do !k
c
        cost=(rij(1)*rkj(1)+rij(2)*rkj(2)+rij(3)*rkj(3))/(dij*dkj)
c
        cost0=cost-cos0
        kcos=Kang*cost0
	edef = 0.5*Kang*cost0**2
c
        aij = -kcos/dij
        akj = -kcos/dkj
c
        do k=1,3
        f1(k) = aij*(rkj1(k)-rij1(k)*cost)
        f3(k) = akj*(rij1(k)-rkj1(k)*cost)
        f2(k) = - f1(k) - f3(k)
        end do !k
c
	if(CONTROL)then
        do k=1,3
        fs(k) = f1(k)+f2(k)+f3(k)
        end do !k
c
        write(kanalp,*)'vAngenforce: x1,x2,x3:'
        write(kanalp,'(3f8.3,2x,3f8.3,2x,3f8.3)')
     &  xyz1,xyz2,xyz3
        write(kanalp,'(a21,2f7.3,2x,f7.3)')
     &  'Par2,ed,f1,f2,f3:', angPar, edef
        write(kanalp,'(a4,4(3f7.3,1x))')
     &  'f:',f1(1),f1(2),f1(3),f2(1),f2(2),f2(3),f3(1),f3(2),f3(3),
     &   fs(1),fs(2),fs(3)
        end if !C
c
	return
        end
c
c torsion energy forces
c
	subroutine torsanglenf(xyz1,xyz2,xyz3,xyz4,nTorsH,
     &                          torsPar,eTors,f1,f2,f3,f4,alarm)


c torsPar(6) = K1,delt1,K2,delt2,K3,delt3
c  eTors = sum{ Ki*[1+cos(delti)cos(i*Ftors)] }
c
   	implicit none
        include "charStringSiz.h"
        include "output.h" 
        real xyz1(3),xyz2(3),xyz3(3),xyz4(3)
        real f1(3),f2(3),f3(3),f4(3)
        integer nTorsH
        real torsPar(*)
        real eTors
        integer alarm
c
        real rij(3),rkj(3),rkl(3)
        real rim(3),rln(3)
        real dij,dkj,dim,dln 
        real aij,akj
        real cosfi,delt
	real cosfi2,kb
        real scijkj,scklkj
        real det,epsi
        real epsR
        integer nFi
        integer k,m,m4
        integer kanalp
        logical CONTROL
c
        CONTROL = .false.
        kanalp =  kanalRunOut
        alarm = 0
c
        epsi = 0.10
        epsR = 0.10
c
        eTors = 0.0
        do k=1,3
        f1(k) = 0.0
        f2(k) = 0.0
        f3(k) = 0.0
        f4(k) = 0.0
        end do !k
c
        dij=0.0
        dkj=0.0
        scijkj=0.0
        scklkj=0.0
	do k =1,3
        rij(k) = xyz1(k) - xyz2(k) 
        rkj(k) = xyz3(k) - xyz2(k) 
        rkl(k) = xyz3(k) - xyz4(k)
        scijkj=scijkj +  rij(k)*rkj(k)
	scklkj=scklkj +  rkl(k)*rkj(k)
        dkj=dkj +  rkj(k)**2
        end do !k
c
        if(dkj .le. epsR ) then
        write(kanalp,*)
     &  'torsanglenf: WARNING!! dkj=zero, return ZERO force!'
        alarm = 1
        return
        end if
c
        scijkj=scijkj/dkj
        scklkj=scklkj/dkj 
c
        cosfi = 0.0
        dim = 0.0
        dln = 0.0
	do k =1,3
        rim(k) = rij(k) - rkj(k)*scijkj
        rln(k) = rkj(k)*scklkj - rkl(k)
        cosfi = cosfi + rim(k)*rln(k)
        dim = dim + rim(k)**2
	dln = dln + rln(k)**2
        end do!k
c
        dim = sqrt(dim)
        dln = sqrt(dln)
c
        if(dim .le. epsR .or. dln .le. epsR ) then
        write(kanalp,*)
     &  'torsanglenf: WARNING!! dim,dln=zero, return ZERO force!',
     &  dim,dln
        alarm = 1
        return
        end if
c
        do k=1,3
        rim(k) = rim(k)/dim 
        rln(k) = rln(k)/dln
        end do!
c
	cosfi = cosfi/(dim*dln)
        cosfi2=cosfi**2
c
        if(CONTROL)then
        write(kanalp,*)
     & 'in torsanglenf: x1,x2,x3,x4:'
        write(kanalp,'(4(3f8.3,1x))')
     &  (xyz1(k),k=1,3),(xyz2(k),k=1,3),(xyz3(k),k=1,3),(xyz4(k),k=1,3)
        write(kanalp,*)'torsAng: cosfi=',cosfi
        end if
c
        do m = 1,nTorsH
        m4 = 4*m-4
        kb = torsPar(m4+2) 
        delt = torsPar(m4+3)
        nFi = int(torsPar(m4+4) + epsi)
c
        if(kb .ne. 0.0)then 
        if(nFi .eq. 1)then
        det = kb*(1.0 + delt*cosfi) 
        eTors = eTors + det
c
        if(CONTROL)then
        write(kanalp,'(a23,f7.3,f5.1,i3,a18,2f8.3)')
     &  'torsPar:Vt/2,delt,nFi:',kb,delt,nFi,' cosFi, eTors:', cosfi,det
        end if
c
        do k=1,3
        f1(k) = f1(k) - kb*delt*(rln(k)-rim(k)*cosfi)/dim 
        f4(k) = f4(k) - kb*delt*(rim(k)-rln(k)*cosfi)/dln
        end do!k
        end if !m=1
        if( nFi .eq. 2) then
        det = kb*(1.0 + delt*(2.0*cosfi2-1.0))
        eTors = eTors + det
c
        if(CONTROL)then
        write(kanalp,'(a23,f7.3,f5.1,i3,a18,2f8.3)')
     &  'torsPar:Vt/2,delt,nFi:',kb,delt,nFi,' cosFi, eTors:', cosfi,det
        end if
c
        do k=1,3
        f1(k) = f1(k) - kb*delt*4.0*cosfi*(rln(k)-rim(k)*cosfi)/dim
        f4(k) = f4(k) - kb*delt*4.0*cosfi*(rim(k)-rln(k)*cosfi)/dln
        end do!k
        end if !nfi=2
c
        if( nFi .eq. 3) then
        det = kb*(1.0 + delt*(4.0*cosfi2-3.0)*cosfi)
        eTors = eTors + det
c
        if(CONTROL)then
        write(kanalp,'(a23,f7.3,f5.1,i3,a18,2f8.3)')
     &  'torsPar:Vt/2,delt,nFi:',kb,delt,nFi,' cosFi, eTors:', cosfi,det
        end if
c
        do k=1,3
        f1(k) = f1(k) - kb*delt*(12.0*cosfi2-3.0)
     &                         *(rln(k)-rim(k)*cosfi)/dim
        f4(k) = f4(k) - kb*delt*(12.0*cosfi2-3.0)
     &                         *(rim(k)-rln(k)*cosfi)/dln
        end do!k
        end if !nfi=3
c      
        if( nFi .eq. 4) then
        det = kb*(1.0 + delt*(8.0*cosfi2*(cosfi2-1.0) + 1.0))
        eTors = eTors + det
c
        if(CONTROL)then
        write(kanalp,'(a23,f7.3,f5.1,i3,a18,2f8.3)')
     &  'torsPar:Vt/2,delt,nFi:',kb,delt,nFi,' cosFi, eTors:', cosfi,det
        end if
c
        do k=1,3
        f1(k) = f1(k) - kb*delt*16.0*cosfi2*(cosfi - 1.0)
     &                         *(rln(k)-rim(k)*cosfi)/dim
        f4(k) = f4(k) - kb*delt*16.0*cosfi2*(cosfi - 1.0)
     &                         *(rim(k)-rln(k)*cosfi)/dln
        end do!k
        end if !nfi=4
c      
        end if!kb.ne.0
        end do!m
c
        do k =1,3
        f2(k) =  (scijkj - 1.0)*f1(k) - scklkj*f4(k)
        f3(k) =  - f1(k) - f2(k) - f4(k)
        end do!k
        
	return
	end
c
c improper torsion energy force
c
        subroutine imprtorsanglenf(xyz1,xyz2,xyz3,xyz4,impPar,
     &                         eImpt,f1,f2,f3,f4,error)

c
c ImptPak(2) = K1, ksi0
c  eImpt =  K1*[(ksi - ksi0)**2] }
c
        implicit none
        include "charStringSiz.h"
        include "output.h" 
        real xyz1(3),xyz2(3),xyz3(3),xyz4(3)
        real f1(3),f2(3),f3(3),f4(3)
        real impPar(2)
        real eImpt
        integer error
c
        real rij(3),rkj(3),rkl(3)
        real rmj(3),rnk(3)
        real dij,dkj,dmj,dnk
        real dkj1,dmj1,dnk1
        real aij,akj,a1,a4
        real cosksi,sksi,ksi,dksi
        real scijkj,scklkj,scmjnk
        real dmjnk1
        real epsR  
        integer k
        integer kanalp
        logical CONTROL

        kanalp = kanalRunOut
        CONTROL = .false.
c
        epsR = 0.10
        error = 0
        eImpt = 0.0
        do k=1,3
        f1(k) = 0.0
        f2(k) = 0.0
        f3(k) = 0.0
        f4(k) = 0.0
        end do !k
c
        dkj=0.0
        scijkj=0.0
        scklkj=0.0
        do k =1,3
        rij(k) = xyz1(k) - xyz2(k)
        rkj(k) = xyz3(k) - xyz2(k)
        rkl(k) = xyz3(k) - xyz4(k)
        scijkj=scijkj +  rij(k)*rkj(k)
        scklkj=scklkj +  rkl(k)*rkj(k)
        dkj=dkj +  rkj(k)**2
        end do !k
c
        if(dkj .le. epsR) then
        write(kanalp,*)
     &  'impAnglenf: WARNING!! dkj=zero, return ZERO force!'
        error = 1
        return
        end if
c
        scijkj=scijkj/dkj
        scklkj=scklkj/dkj
c
        call vectorp(rij,rkj,rmj)
        call vectorp(rkj,rkl,rnk)
        call scalarp(rmj,rnk,scmjnk)
        call scalarp(rmj,rmj,dmj)
        call scalarp(rnk,rnk,dnk)
        dkj1=sqrt(dkj)
        dmjnk1=sqrt(dmj*dnk)
c
	call scalarp(rij,rnk,sksi)

        if(sksi .ge. 0.0 ) then
         sksi = 1.0
         else
         sksi = -1.0
        end if
c
        if(dmjnk1 .le. epsR .or. dmj .le. epsR
     &                     .or. dnk .le. epsR ) then
        write(kanalp,*)
     &  'impAnglenf: WARNING!! dmjnk1,dmj,dnk=zero, return ZERO force!'
     &  ,dmjnk1,dmj,dnk
        write(kanalp,*)'atom1:',xyz1
        write(kanalp,*)'atom2:',xyz2
        write(kanalp,*)'atom3:',xyz3
        write(kanalp,*)'atom4:',xyz4
        error = 1
        return
        end if
c
        cosksi = scmjnk/(dmjnk1)
        if(cosksi .lt. -1.0)cosksi = -1.0
        if(cosksi .gt. 1.0)cosksi = 1.0
        ksi = sksi*acos(-cosksi)
c
        dksi = ksi - ImpPar(2)
        eImpT = 0.5*ImpPar(1)*dksi**2
c
        a1 = dkj1/dmj
        a4 = dkj1/dnk
        do k=1,3
        f1(k) =  ImpPar(1)*dksi*a1*rmj(k)
        f4(k) = -ImpPar(1)*dksi*a4*rnk(k)
        end do!k
c
        do k =1,3
        f2(k) =  (scijkj - 1.0)*f1(k) - scklkj*f4(k)
        f3(k) =  - f1(k) - f2(k) - f4(k)
        end do!k
c
	if(CONTROL)then
        write(kanalp,*)'impEforce: x1,x2,x3,x4:'
        write(kanalp,'(3f8.3,2x,3f8.3,2x,3f8.3,2x,3f8.3)')
     &  xyz1,xyz2,xyz3,xyz4
        write(kanalp,*)'impEforce: dmj,dmk,cosksi,ksi,dksi:'
        write(kanalp,'(5f12.7)') dmj,dnk,cosksi,ksi,dksi
        write(kanalp,'(a21,2f7.3,2x,f7.3)')
     &  'Par2,ed:', impPar,eImpT
        write(kanalp,'(a4,4(3f7.3,1x))')
     &  'f:',f1(1),f1(2),f1(3),f2(1),f2(2),f2(3),f3(1),f3(2),f3(3),
     &   f4(1),f4(2),f4(3)
        end if
	return
        end
c
