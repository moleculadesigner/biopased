c get rotated molecularCoordinates XYZ via X,Y,Z rotation+translation
c rotation are done around GeoCenter
c
c Y.Vorobjev 2003
c
c
c
	subroutine rotateMolZXZ(n,XYZ,tet,vtra,XYZout)
c
c xNEW = tz2*tx*tz1*(xO - xGeoC) + xGeoC + vtra 
c n - number of atoms to be rotated in atXYZ()
c XYZ() - initial and final(rotated translated) XYZ
c tet(1,2,3) - tz1,tX,tz2 : Eiler Angles: tz1,tX=Fi: defines axis orientation
c                           tz2 = rotation around ritationAxis
c vtra(3) - vX,vY,vZ
c
	implicit none
	integer n
	real XYZ(*)
        real XYZout(*)
        real tet(3),vtra(3)
c
	integer i,i3  
        real tt1(3,3),tt2(3,3),tt3(3,3)
        real c(3),c1(3),c2(3)
        real geoC(3)
c
        geoC(1)=0.0
        geoC(2)=0.0
	geoC(3)=0.0
        do i=1,n
        i3=3*i-3
        geoC(1)=geoC(1)+XYZ(i3+1)
        geoC(2)=geoC(2)+XYZ(i3+2)
        geoC(3)=geoC(3)+XYZ(i3+3)
	end do !i
        geoC(1)=geoC(1)/n        
        geoC(2)=geoC(2)/n          
        geoC(3)=geoC(3)/n          
c
        call gen_tz(tet(3),tt3)
        call gen_tx(tet(2),tt2)
        call gen_tz(tet(1),tt1)
c
        do i = 1,n
        i3=3*i-3
        c(1)=XYZ(i3+1) - geoC(1)
        c(2)=XYZ(i3+2) - geoC(2)
        c(3)=XYZ(i3+3) - geoC(3)
c
        call matZvb(tt1,c,c1)
        call matXvb(tt2,c1,c2)
        call matZvb(tt3,c2,c)
c
        XYZout(i3+1) = c(1)+geoC(1)+vtra(1)
        XYZout(i3+2) = c(2)+geoC(2)+vtra(2)
        XYZout(i3+3) = c(3)+geoC(3)+vtra(3)
c
        end do!i
c
	return
	end
c
	subroutine rotateMolXYZ(n,XYZ,tet,vtra)
c
c xNEW = tx*ty*tz*(xO - xGeoC) + vtra 
c n - number of atoms to be rotated in atXYZ()
c XYZ() - initial and final(rotated translated) XYZ
c tet(3) - tX,tY,tZ
c vtra(3) - vX,vY,vZ
c
	implicit none
	integer n
	real XYZ(*)
        real tet(3),vtra(3)
c
	integer i,i3  
        real ttx(3,3),tty(3,3),ttz(3,3)
        real c(3),c1(3),c2(3)
        real geoC(3)
c
        geoC(1)=0.0
        geoC(2)=0.0
	geoC(3)=0.0
        do i=1,n
        i3=3*i-3
        geoC(1)=geoC(1)+XYZ(i3+1)
        geoC(2)=geoC(2)+XYZ(i3+2)
        geoC(3)=geoC(3)+XYZ(i3+3)
	end do !i
        geoC(1)=geoC(1)/n        
        geoC(2)=geoC(2)/n          
        geoC(3)=geoC(3)/n          
c
        call gen_tz(tet(3),ttz)
        call gen_ty(tet(2),tty)
        call gen_tx(tet(1),ttx)
c
        do i = 1,n
        i3=3*i-3
        c(1)=XYZ(i3+1) - geoC(1)
        c(2)=XYZ(i3+2) - geoC(2)
        c(3)=XYZ(i3+3) - geoC(3)
c
        call matZvb(ttz,c,c1)
        call matYvb(tty,c1,c2)
        call matXvb(ttx,c2,c)
c
        XYZ(i3+1) = c(1)+vtra(1)
        XYZ(i3+2) = c(2)+vtra(2)
        XYZ(i3+3) = c(3)+vtra(3)
c
        end do!i
c
	return
	end
c

c XYZ rotation set of subroutines
c
c tx(3,3) rotation matrix
        subroutine gen_tx(tet,ttm)
        implicit none
        real tet
        real ttm(3,3)
C generate tx
        ttm(1,1)=1.0
        ttm(1,2)=0.0
        ttm(1,3)=0.0
        ttm(2,1)=0.0
        ttm(2,2)=cos(tet)
        ttm(2,3)=-sin(tet)
        ttm(3,1)=0.0
        ttm(3,2)=-ttm(2,3)
        ttm(3,3)=ttm(2,2)
        return
        end
c
        subroutine gen_ty(tet,ttm)
        implicit none
        real tet
        real ttm(3,3)
C generate ty
        ttm(2,1)=0.0
        ttm(2,2)=1.0
        ttm(2,3)=0.0
        ttm(1,2)=0.0
        ttm(1,1)=cos(tet)
        ttm(1,3)=-sin(tet)
        ttm(3,2)=0.0
        ttm(3,1)=-ttm(1,3)
        ttm(3,3)=ttm(1,1)
        return
        end
c
        subroutine gen_tz(tet,ttm)
        implicit none
        real tet
        real ttm(3,3)
C generate tz
        ttm(3,3)=1.0
        ttm(3,2)=0.0
        ttm(3,1)=0.0
        ttm(1,3)=0.0
        ttm(1,1)=cos(tet)
        ttm(1,2)=-sin(tet)
        ttm(2,3)=0.0
        ttm(2,1)=-ttm(1,2)
        ttm(2,2)=ttm(1,1)
        return
        end
c
        subroutine gen_grtx(tet,ttm)
        implicit none
        real tet
        real ttm(3,3)
C generate grad tx
        ttm(1,1)=0.0
        ttm(1,2)=0.0
        ttm(1,3)=0.0
        ttm(2,1)=0.0
        ttm(2,2)=-sin(tet)
        ttm(2,3)=-cos(tet)
        ttm(3,1)=0.0
        ttm(3,2)=-ttm(2,3)
        ttm(3,3)=ttm(2,2)
        return
        end
c
        subroutine gen_grty(tet,ttm)
        implicit none
        real tet
        real ttm(3,3)
C generate ty
        ttm(2,1)=0.0
        ttm(2,2)=0.0
        ttm(2,3)=0.0
        ttm(1,2)=0.0
        ttm(1,1)=-sin(tet)
        ttm(1,3)=-cos(tet)
        ttm(3,2)=0.0
        ttm(3,1)=-ttm(1,3)
        ttm(3,3)=ttm(1,1)
        return
        end
c
        subroutine gen_grtz(tet,ttm)
        implicit none
        real tet
        real ttm(3,3)
C generate tz
        ttm(3,3)=0.0
        ttm(3,2)=0.0
        ttm(3,1)=0.0
        ttm(1,3)=0.0
        ttm(1,1)=-sin(tet)
        ttm(1,2)=-cos(tet)
        ttm(2,3)=0.0
        ttm(2,1)=-ttm(1,2)
        ttm(2,2)=ttm(1,1)
        return
        end
c vc = ma*vb
        subroutine mavb(ma,vb,vc)
        implicit none
        real ma(3,3)
        real vb(3),vc(3)
        vc(1)=ma(1,1)*vb(1)+ma(1,2)*vb(2)+ma(1,3)*vb(3)
        vc(2)=ma(2,1)*vb(1)+ma(2,2)*vb(2)+ma(2,3)*vb(3)
        vc(3)=ma(3,1)*vb(1)+ma(3,2)*vb(2)+ma(3,3)*vb(3)
        return
        end
c
c special for rotation mtrix
        subroutine matXvb(ma,vb,vc)
c ma - mTx() matrix TX
        implicit none
        real ma(3,3)
        real vb(3),vc(3)
        vc(1)=vb(1)
        vc(2)=ma(2,2)*vb(2)+ma(2,3)*vb(3)
        vc(3)=ma(3,2)*vb(2)+ma(3,3)*vb(3)
        return
        end
c
        subroutine matYvb(ma,vb,vc)
        implicit none
        real ma(3,3)
        real vb(3),vc(3)
        vc(1)=ma(1,1)*vb(1)+ma(1,3)*vb(3)
        vc(2)=vb(2)
        vc(3)=ma(3,1)*vb(1)+ma(3,3)*vb(3)
        return
        end
c
        subroutine matZvb(ma,vb,vc)
        implicit none
        real ma(3,3)
        real vb(3),vc(3)
        vc(1)=ma(1,1)*vb(1)+ma(1,2)*vb(2)
        vc(2)=ma(2,1)*vb(1)+ma(2,2)*vb(2)
        vc(3)=vb(3)
        return
        end
c
