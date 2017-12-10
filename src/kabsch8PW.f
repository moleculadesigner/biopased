c *******************************************************************
c *  Determination of the best rotation of a given vector set into  *
c * a second vector set by minimazing the sum of squared deviations *
c *             E = 1/N * SUM(w(k) * |U * x(k) - y(k)|**2)          *   
c *                        k                                        *
c ***                                                             ***
c *** Required functions : jacobi() for diagonalization,          ***
c ***                      eigsrt() for sorting eigenvalues       ***  
c ***                                                             ***
c *** x(3, n) -- rotated vector set, on output x = U*x            ***
c *** y(3, n) -- fixed vector set                                 ***
c *** n       -- number of vectors                                ***
c *** w(n)    -- weight of each vector                            ***
c *** u(3, 3) -- returned matrix of proper rotation               ***
c ***                                                             ***
c *** WARNING: x & y must have physical length = 3 in 1st dim     ***
c *******************************************************************

      SUBROUTINE kabsch8(x, y, n, w,u,tx,ty)
c	implicit none
        include 'xyzPDBsize.h'
	integer n
	real x(*), y(*), w(*)
	real*8 u(3, 3)
	real*8 tx(*),ty(*)
c
        real*8 buf(3*natomMAX)
        real*8 x8(3*natomMAX)
        real*8 y8(3*natomMAX)
        real*8 w8(natomMAX) 
c ............. local variables ..................
	real*8 r(3, 3), rt(3, 3), rr(3, 3)
	real*8 a(3, 3), b(3, 3)
	real*8 mu(3)
	real*8 s, b3(3)
cx        real*8 tx(3),ty(3)
        real*8 sw
        real*8 zero,one
        logical jdone
	integer i,j,k,nrot
        integer i3,k3
c initiate
        zero = 0.0d0
        one = 1.0d0
        do i=1,3*n
         x8(i) = dble(x(i))
         y8(i) = dble(y(i))
        end do!i
        do i=1,n
         w8(i)=dble(w(i))
        end do !i
c eliminate overal translation
         do k=1,3
          tx(k)=zero
          ty(k)=zero 
         end do !k
c
         do i = 1, 3
            do j = 1, 3
                r(i, j) = zero
            end do
          end do !i
c
          do i=1,n
          i3 = 3*i-3
          do k=1,3
          tx(k) = tx(k) + x8(i3+k)*w8(i)
          ty(k) = ty(k) + y8(i3+k)*w8(i)
          end do!k
          end do!i
c
        sw = zero
        do k=1,n   
         sw = sw + w8(k)    
        end do!k
c 
        if(sw .ne. zero)sw = one/sw
        do k=1,3
        tx(k)=tx(k)*sw
        ty(k)=ty(k)*sw
        end do!k
c
        do i=1,n
         i3=3*i-3
          do k=1,3
          x8(i3+k) = x8(i3+k) - tx(k)  
          y8(i3+k) = y8(i3+k) - ty(k) 
          end do!k
        end do!i
c x8,y8 translated to the origin CM
c .............. determine matrix R ............
     	do k = 1, n
        k3=3*k-3
         do i = 1, 3
            do j = 1, 3
	      r(i, j) = r(i, j) + w8(k)*y8(k3+i)*x8(k3+j)
	    end do
	   end do
	end do
c
c ............. transpon R .................

	do i = 1, 3
            do j = 1, 3
		rt(i, j) = r(j, i)
	    end do
        end do

c ............ compute matrix Rt * R ...........

	do i = 1, 3
	    do j = 1, 3
		rr(i, j) = zero
		do k = 1, 3
		    rr(i, j) = rr(i, j) + rt(i, k) * r(k, j)
		end do
            end do
        end do

c ........ determine eigenvalues & vectors of Rt * R .....
         i=3
	call jacobi8(rr, i, i, mu, a, nrot,jdone)
        if( .not. jdone)return   !! yv
c
	call eigsrt8(mu, a, i, i)               !!! sort mu

c ........... determine a(3) = [a(1) * a(2)] ..............

	a(1, 3) = a(2, 1) * a(3, 2) - a(3, 1) * a(2, 2)
	a(2, 3) = - a(1, 1) * a(3, 2) + a(3, 1) * a(1, 2)
	a(3, 3) = a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2)

c ......... determine b(k) = R * a(k)  .........

	do i = 1, 3
            do k = 1, 3
		b(i, k) = zero
		do j = 1, 3
		    b(i, k) = b(i, k) + r(i, j) * a(j, k) 
		end do
	        if(k .ne. 3) b(i, k) = b(i, k) / sqrt(mu(k))
	    end do
        end do

c ......... determine b3 = [b(1) * b(2)] .............

	b3(1) = b(2, 1) * b(3, 2) - b(3, 1) * b(2, 2)
        b3(2) = - b(1, 1) * b(3, 2) + b(3, 1) * b(1, 2)
        b3(3) = b(1, 1) * b(2, 2) - b(2, 1) * b(1, 2)

c ......... determine proper rotation ..............

	s = zero
	do i = 1, 3
	    s = s + b3(i) * b(i, 3)
	end do
	if(s .lt. 0) then
	    s = -one
	else
	    s = one
	end if
	do i = 1, 3
	    b(i, 3) = s * b3(i)
	end do

c ......compute matrix of rotation U ............

	do i = 1, 3
            do j = 1, 3
                u(i, j) = zero
                do k = 1, 3
                    u(i, j) = u(i, j) + b(i, k) * a(j, k)
                end do
            end do
           end do
c
        do i=1,3*n
        buf(i) = x8(i)  !!! buffer
        end do !i
c rotate given vector set x = Ux
	do i = 1, n
         i3=3*i-3
            do j = 1, 3
		x8(i3+j) = zero
	      do k = 1, 3
		x8(i3+j) = x8(i3+j) + buf(i3+k) * u(j, k)
	      end do	
	    end do
        end do
c return translation to reference y
        do i=1,n
         i3=i*3-3
        do k=1,3
         x(i3+k) = x8(i3+k)+ty(k)
        end do !k 
        end do !i
c
        return
        end
c
        subroutine xyz_rmsd(x,y,n,sd)
c
        implicit none
        real x(*)
        real y(*)
        integer n
c
        real sd
        integer i,i3,k
        real*8 s,ss,zero
c
        zero = 0.0d0
        ss=zero
        do i = 1, n      
          i3=3*i-3
            s = zero
            do k = 1, 3
                s = s + (x(i3+k) - y(i3+k))**2
            end do
            ss = ss + s
        end do
        sd = dsqrt(ss/n)
c
        return
        end
c rotation+ translation:
c rotate x by rotMatrx u : x = x*U
c tranclate x CMass to y CMass
c
      SUBROUTINE rotateMOLEC8(x, n, u,tx,ty)
c       implicit none
        include 'xyzPDBsize.h'
        integer n
        real x(*)
        real*8 u(3, 3),tx(3),ty(3)
c
        real*8 buf(3*natomMAX)
        real*8 x8(3*natomMAX)
c ............. local variables ..................
        real*8 zero,one
        logical jdone
        integer i,j,k,nrot
        integer i3,k3
c initiate
        zero = 0.0d0
        one = 1.0d0
        do i=1,3*n
        x8(i) = dble(x(i))
        end do!i
c
        do i=1,n
        i3=3*i-3
        do k=1,3
        x8(i3+k) = x8(i3+k) - tx(k)
        end do!k
        end do!i
c x8,y8 translated to the origin CM
c
        do i=1,3*n
        buf(i) = x8(i)  !!! buffer
        end do !i
c rotate given vector set x = Ux
          do i = 1, n
        i3=3*i-3
           do j = 1, 3
        x8(i3+j) = zero
            do k = 1, 3
        x8(i3+j) = x8(i3+j) + buf(i3+k) * u(j, k)
           end do
         end do
       end do
c return translation to reference y
        do i=1,n
        i3=i*3-3
        do k=1,3
        x(i3+k) = x8(i3+k)+ty(k)
        end do !k
        end do !i
c
        return
        end
