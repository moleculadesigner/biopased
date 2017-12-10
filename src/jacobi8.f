c jacobi: computes all egenvalues and egenvectors of real symmetric matrix A.
c of size n*n, with maxdimens(parameter) np*np: when call in main progr.
c On output elements of A above diaginal are destroed.
c D returns egenvalues
c V(k,i) is matrix with colomns equal to normalized egenvectors egv-i(k)
c nrot  returns the number of Jacobi rotations
c reguire N**3 poperations
c
      SUBROUTINE jacobi8(a,n,np,d,v,nrot,done)
      implicit none
      INTEGER n,np,nrot,NMAX
      REAL*8 a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=3)     !local MAX matrix A size
      INTEGER i,ip,iq,j
      REAL*8 c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      logical done
c
      include "charStringSiz.h"
      include "output.h"  
c
      real*8 zero,one
      integer jiterMX
      integer kanalp
c
      done = .false.
      kanalp = kanalRunOut
c
      jiterMX = 500
      zero = 0.0d0
      one = 1.0d0
      if(n .gt. NMAX)then
      write(kanalp,*)'jacobi:ERROR:! local matrix size NMAX is low!'
      stop
      end if
c
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=zero
11      continue
        v(ip,ip)=one
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=zero
13    continue
      nrot=0
      do 24 i=1,jiterMX
        sm=zero
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.zero)then
        done = .true.
        return
        end if !sm.eq.zero
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=zero
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.0d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=zero
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5*h/a(ip,iq)
                t=one/(abs(theta)+dsqrt(one+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=one/dsqrt(one+t**2)
              s=t*c
              tau=s/(one+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=zero
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
cx      pause 'too many iterations in jacobi'
      write(kanalp,*)'too many iterations in jacobi'
      done = .false.
      return
      end
c
c sorting of egenvalues in ascending order
c D - given eigenvalue, V(,) - eigenvectors
c returns sorted D and V
c require N**2 operations
c
      SUBROUTINE eigsrt8(d,v,n,np)
      INTEGER n,np
      REAL*8 d(np),v(np,np)
      INTEGER i,j,k
      REAL*8 p
      do 13 i=1,n-1
        k=i
        p=d(i)
        do 11 j=i+1,n
          if(d(j).ge.p)then
            k=j
            p=d(j)
          endif
11      continue
        if(k.ne.i)then
          d(k)=d(i)
          d(i)=p
          do 12 j=1,n
            p=v(j,i)
            v(j,i)=v(j,k)
            v(j,k)=p
12        continue
        endif
13    continue
      return
      END
