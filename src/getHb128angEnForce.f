c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* Hx...Y Hb12-8 angle  H-bond potential                                     *
c* soft VDW potential 
c*                                                                           *
c*     Yury Vorobjev 2004                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
	subroutine getHb128angEforce(xyz1,xyz2,xyz3,hbAngPar,
     &                         eHbAng,f1,f2,f3)
c
c xyz3,2,1 = X,H,Y  in the X-H...Y h-bond
c
        implicit none
        real xyz1(3),xyz2(3),xyz3(3)
        real f1(3),f2(3),f3(3)
        real hbAngPar(*)
        real eHbAng
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
        real rm12,em12,A12,B12,ehb0
        real ctm0,sg2,fAng,dFang  
        logical OPT_AngleDep
        real f0i(3),fs(3)
        real cos0,Kang
        real rij(3),rkj(3)
        real rij1(3),rkj1(3)
        real dij,dkj,cost
        real dij2
        real aij(3),akj(3)
        integer k
        integer kanalp
        logical CONTROL

        kanalp = kanalRunOut
        CONTROL = .false.
c
        OPT_AngleDep = .true.
c        OPT_AngleDep = .false.
c
c init
        eHbAng = 0.0
        do k=1,3
        f1(k) = 0.0
        f2(k) = 0.0
        f3(k) = 0.0
        end do !k
c
        rm12=hbAngPar(3)
        em12=hbAngPar(4)
        A12=hbAngPar(1)
        B12=hbAngPar(2)
        cos0 = hbAngPar(5)
        sg2 = hbAngPar(6)
c
        rij(1)=xyz1(1)-xyz2(1)
        rij(2)=xyz1(2)-xyz2(2)
        rij(3)=xyz1(3)-xyz2(3)
        rkj(1)=xyz3(1)-xyz2(1)
        rkj(2)=xyz3(2)-xyz2(2)
        rkj(3)=xyz3(3)-xyz2(3)
        dij2 = rij(1)**2+rij(2)**2+rij(3)**2
        dij = sqrt(dij2)
        dkj = sqrt(rkj(1)**2+rkj(2)**2+rkj(3)**2)
c
        if(dij .eq. 0.0 .or. dkj .eq. 0.0 ) then
        write(kanalp,*)'Hb128efenforce: WARNING! dij,dkj = zero:'
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
        cost=rij1(1)*rkj1(1)+rij1(2)*rkj1(2)+rij1(3)*rkj1(3)
c
        ctm0=cost-cos0
c
        call getHb128EFij00
     &               (dij2,dij,rij,rm12,em12,A12,B12,ehb0,f0i)
c
c force10 = -fi0 
c initilize anglDependence
        fAng = 1.0
        dFang = 0.0
c
        if(OPT_AngleDep)then
        call eHb128angleDep(ctm0,sg2,fAng,dFang)
        end if 

c
        do k=1,3
        aij(k) = (rkj1(k)-rij1(k)*cost)/dij
        akj(k) = (rij1(k)-rkj1(k)*cost)/dkj
c
        f1(k) = -f0i(k)*fAng - ehb0*dFang*aij(k)
        f3(k) = -ehb0*dFang*akj(k)
        f2(k) = - f1(k) - f3(k)
        end do !k
c toal energy
        eHbAng = ehb0*fAng
c
        if(CONTROL)then
        do k=1,3
        fs(k) = f1(k)+f2(k)+f3(k)
        end do !k
c
        write(kanalp,*)'vAngenforce: x1,x2,x3:'
        write(kanalp,'(3f8.3,2x,3f8.3,2x,3f8.3)')
     &  xyz1,xyz2,xyz3
        write(kanalp,'(a6,6f7.3,a8,f8.3)')
     &  'Par6:', (hbAngPar(k),k=1,6), ' eHbAng:',eHbAng
        write(kanalp,'(a4,4(3f7.3,1x))')
     &  'f:',f1(1),f1(2),f1(3),f2(1),f2(2),f2(3),f3(1),f3(2),f3(3),
     &   fs(1),fs(2),fs(3)
        end if !C
c
        return
	end

        subroutine getHb128EFij00
     &               (dij2,dij1,rij,rm12,em12,A12,B12,ehb,fi)
c
c  rij = rj - ri
c - 12-8 potential - !? is not good, too shortRange
c - 12-6 - ?
c  rm12 - Rmin of potential
c  em12 - Emin
c  OUT:
c       fi = force on atom i
c
        implicit none
        real dij2,dij1
        real rm12,em12,A12,B12
        real rij(3),ehb,fi(3)
c local
        real d4,d6,d8,fv
c        integer m,n
        integer kanalp
        logical CONTROL
c
c potential A/r^n - B/r^m
c
        kanalp = 6
        CONTROL = .false.
c
c        n=12
c        m=6
c
        d4=dij2**2
c        d8=d4**2
        d6=d4*dij2
 
        if( dij1 .gt. rm12 ) then
        ehb = (A12/d6 - B12)/d6
c
        fv = (6.0*B12 - 12.0*A12/d6)/d6/dij2
c
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
c
        else
c
        ehb = -em12 
        fi(1)=0.0          
        fi(2)=0.0           
        fi(3)=0.0          
        end if
c
        if(CONTROL)then
        write(kanalp,'(a35,3f8.4,a5,f8.3,a5,3f8.3)')
     &  'getPairHb128EFij: dij1,rm12,em12:',dij1,rm12,em12,' ehb:',ehb,
     &  ' fi:',fi
        end if
c
        return
        end
c
c hB128angle dependence
	subroutine eHb128angleDep(ctm0,sg2,fAng,dFang)
c
c ctmo = cosTeta  cost0, teta=angle(xHy)
c sg2 = sigma**2
c fAng - angle dependence function
c dFang = d fAng/d costeta
c
	implicit none
        real ctm0,b
        real fAng,dFang,sg2
c
         b = ctm0/sg2
	fAng = exp(-ctm0*b)
ce
	dFang = -2.0*b*fAng
c
	return
	end
c

