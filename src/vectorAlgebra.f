c vector Algebra
        subroutine scalarp(v1,v2,sp)
        implicit none
        real v1(3),v2(3)
        real sp
        sp = v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
        return
        end
c
        subroutine vectorNrm1(z)
c return Normalized to 1.0 vector z
        implicit none
        real z(3)
        integer i
        real a
c
        call scalarp(z,z,a)
        a = sqrt(a)
        z(1) = z(1)/a
        z(2) = z(2)/a
        z(3) = z(3)/a
c
        return
        end
c
        subroutine vectorNrm2(v,vn)
c return Normalized to 1.0 vector vn
        implicit none
        real v(3),vn(3)
        integer i
        real a
c
        call scalarp(v,v,a)
        a = sqrt(a)
        vn(1) = v(1)/a
        vn(2) = v(2)/a
        vn(3) = v(3)/a
c
        return
        end
c vector product
c v3=v1*v2
        subroutine vectorp(v1,v2,v3)
        implicit none
        real v1(3),v2(3),v3(3)
c
        v3(1)=v1(2)*v2(3)-v1(3)*v2(2)
        v3(2)=v1(3)*v2(1)-v1(1)*v2(3)
        v3(3)=v1(1)*v2(2)-v1(2)*v2(1)
         return
         end
c
c local coord system on three points
	subroutine getLocRefSyst(r1,r2,r3,r0,x,y,z,p0)
c
c points r1,r2,r3 ;  r0-origin, define
c local Decart coordinate system: x,y,z, p0 - origin in global
c x0y0z0 system
c
	implicit none
        real r1(3),r2(3),r3(3),r0(3)
        real x(3),y(3),z(3),p0(3)
c
        real r12(3),r32(3),r31(3)
        integer i
c
        do i=1,3
        r12(i) = r1(i) - r2(i)
        r32(i) = r3(i) - r2(i)
        r31(i) = r3(i) - r1(i)
        p0(i) = r0(i)
        end do!i
c
        call vectorp(r32,r12,z)
        call vectorNrm1(z)
c
        call vectorNrm2(r12,x)
c
	call vectorp(z,x,y)
        call vectorNrm1(y)
c
        return
        end         
c
c rotation Matrix from x0,y0,z0 TO x1,y1,z1
c x0,y0,z0 - ortogonal unit vectors
c x1,y1,z1 -ortogonal unit vectors of the new Decart system
c
	subroutine getRotMatrix01(x0,y0,z0,p0,
     &                    x1,y1,z1,p1,t01,p01)
c
	implicit none
        real x0(3),y0(3),z0(3),p0(3)
        real x1(3),y1(3),z1(3),p1(3)
        real t01(3,3),p01(3)
        integer i,j
c
        p01(1) = p1(1) - p0(1)
        p01(2) = p1(2) - p0(2)
        p01(3) = p1(3) - p0(3)
c
        call scalarp(x1,x0,t01(1,1))
        call scalarp(x1,y0,t01(1,2))
        call scalarp(x1,z0,t01(1,3))
c
        call scalarp(y1,x0,t01(2,1))
        call scalarp(y1,y0,t01(2,2))
        call scalarp(y1,z0,t01(2,3))
c
        call scalarp(z1,x0,t01(3,1))
        call scalarp(z1,y0,t01(3,2))
        call scalarp(z1,z0,t01(3,3))
c
        return
        end
c
c rotation Matrix from x0,y0,z0 TO x1,y1,z1
c x0,y0,z0 - ortogonal unit vectors
c x1,y1,z1 -ortogonal unit vectors of the new Decart system
c
        subroutine getRotMatrix02(x0,y0,z0,
     &                    x1,y1,z1,t01)
c
        implicit none
        real x0(3),y0(3),z0(3)
        real x1(3),y1(3),z1(3)
        real t01(9)
        integer i,j
c
        call scalarp(x1,x0,t01(1))
        call scalarp(x1,y0,t01(4))
        call scalarp(x1,z0,t01(7))
c
        call scalarp(y1,x0,t01(2))
        call scalarp(y1,y0,t01(5))
        call scalarp(y1,z0,t01(8))
c
        call scalarp(z1,x0,t01(3))
        call scalarp(z1,y0,t01(6))
        call scalarp(z1,z0,t01(9))
c
        return
        end
c
c	rotate MolecFragment 
	subroutine rotMolBy3P(r01,r02,r03,rc0,nat,xyz0,
     &                        r11,r12,r13,rc1,xyz1)
c rotate+translate r01,r02,r03 and rc0(center) points to the NEW position
c                  r11,r12,r13 and rc1
c transform molecule coords xyz0 to a new position xyz1(*)
c
	implicit none
        real r01(3),r02(3),r03(3)
        real rc0(3),rc1(3)
        real xyz0(*)
        integer nat
        real r11(3),r12(3),r13(3)
        real xyz1(*)
c
        integer i,j,ia,ia3
        real t01(3,3),t20(3,3)
        real x0(3),y0(3),z0(3),p0(3) 
        real x1(3),y1(3),z1(3),p1(3) 
        real x2(3),y2(3),z2(3),p2(3) 
        real x(3),y(3)
        real p01(3),p02(3)
c x0,y0,z0 - unitvectors of LAB system
       do i=1,3
       x0(i) = 0.0
       y0(i) = 0.0
       z0(i) = 0.0
       p0(i) = 0.0
       end do !i
c
       x0(1) = 1.0
       y0(2) = 1.0
       z0(3) = 1.0
c get INITial local system
	call getLocRefSyst(r01,r02,r03,rc0,x1,y1,z1,p1)  
c
c get rotation matrix t01
        call getRotMatrix01(x0,y0,z0,p0,
     &                      x1,y1,z1,p1,t01,p01)
c get NEW local syatem
        call getLocRefSyst(r11,r12,r13,rc1,x2,y2,z2,p2)
c
c get  rotation matrix t20
       call getRotMatrix01(x2,y2,z2,p2,
     &                     x0,y0,z0,p0,t20,p02)  
c
c transform 
       do ia = 1,nat
       x(1) = 0.0
       x(2) = 0.0
       x(3) = 0.0
       ia3 = 3*ia-3
       do i=1,3
       do j=1,3
       x(i) = x(i) + t01(i,j)*(xyz0(ia3+j)-p01(j))
       end do!j
       end do!i
c
       y(1) = 0.0
       y(2) = 0.0
       y(3) = 0.0
       do i=1,3
       do j=1,3
       y(i) = y(i) + t20(i,j)*x(j)
       end do!j
       xyz1(ia3+i) = y(i) - p02(i)
       end do!i
       end do !ia
c
	return
	end
c
c Matrix*vector  Product
        subroutine vectMatrx3Prod(v1,Mx,v2)
c
c v2 = Mx*v1
c     1 2 3
c M = 4 5 6
c     7 8 9
c
         implicit none
         real v1(3)
         real Mx(9),v2(3)
c
        v2(1) = v1(1)*Mx(1) + v1(2)*Mx(2) + v1(3)*Mx(3)
        v2(2) = v1(1)*Mx(4) + v1(2)*Mx(5) + v1(3)*Mx(6)
        v2(3) = v1(1)*Mx(7) + v1(2)*Mx(8) + v1(3)*Mx(9)
c
        return
        end
c
c invert Matrix 3*3
c
        subroutine invertM33(mx,mx1,dmx)
c
c     1 2 3
c M = 4 5 6
c     7 8 9
c mx - inMatrix(3*3)
c mx1 - inverted
c dmx - det mx
        implicit none
        real mx(9)
        real mx1(9)
        real dmx
c
        integer i
c
        mx1(1) = mx(5)*mx(9) - mx(6)*mx(8)
        mx1(2) = -mx(2)*mx(9) + mx(3)*mx(8)
        mx1(3) = mx(4)*mx(8) - mx(5)*mx(7)
        mx1(4) = -mx(4)*mx(9) + mx(6)*mx(7)
        mx1(5) = mx(1)*mx(9) - mx(3)*mx(7)
        mx1(6) = -mx(1)*mx(6) + mx(3)*mx(4)
        mx1(7) = mx(4)*mx(8) - mx(5)*mx(7)
        mx1(8) = -mx(1)*mx(8) + mx(2)*mx(7)
        mx1(9) = mx(1)*mx(5) - mx(2)*mx(4)
c
        dmx = mx(1)*mx1(1)+mx(2)*mx1(4)+mx(3)*mx1(7)
c
        do i=1,9
        mx1(i) = mx1(i)/dmx
        end do !i
c
        return
	end
c
        subroutine vector2p(v1,v2,v3,v4)
c double vector product v4 = [v1*[v2*v3]] = v2*(v1*v3) - v3*(v1*v2)
        real v1(3),V2(3),v3(3),v4(3)
        real a,b
        b=v1(1)*v2(1)+v1(2)*v2(2)+v1(3)*v2(3)
        a=v1(1)*v3(1)+v1(2)*v3(2)+v1(3)*v3(3)
        v4(1) = v2(1)*a - v3(1)*b
        v4(2) = v2(2)*a - v3(2)*b
        v4(3) = v2(3)*a - v3(3)*b
        return
        end
c
        subroutine matrix2p(m1,m2,m3)
c
c matrix m3 = m1*m2
c      1  2  3
c  m = 4  5  6
c      7  8  9
        real m1(9),m2(9),m3(9)
c
        m3(1)=m1(1)*m2(1)+m1(2)*m2(4)+m1(3)*m2(7)
        m3(2)=m1(1)*m2(2)+m1(2)*m2(5)+m1(3)*m2(8)
        m3(3)=m1(1)*m2(3)+m1(2)*m2(6)+m1(3)*m2(9)
c
        m3(4)=m1(4)*m2(1)+m1(5)*m2(4)+m1(6)*m2(7)
        m3(5)=m1(4)*m2(2)+m1(5)*m2(5)+m1(6)*m2(8)
        m3(6)=m1(4)*m2(3)+m1(5)*m2(6)+m1(6)*m2(9)
c
        m3(7)=m1(7)*m2(1)+m1(8)*m2(4)+m1(9)*m2(7)
        m3(8)=m1(7)*m2(2)+m1(8)*m2(5)+m1(9)*m2(8)
        m3(9)=m1(7)*m2(3)+m1(8)*m2(6)+m1(9)*m2(9)
c
        return
        end 
