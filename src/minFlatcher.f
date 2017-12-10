c energy optimization by local optimization
c
c Yurii Vorobjev  2003
c 
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
C func_tomin(p) - minimized function      - must be added
C dfunc_tomin(p,xi) - grad=xi, at point p - must be added 
c
c INCLUDES:
C linmin.f
C mnbrak_min.f
C brent.f
C
C minimization,optimization,minima
C at starting point P, vectordim n,
C Fletcher-Reeves-Polak-Ribiere minimization is performed on function FUNC
C convergence tolerance for FUnc = ftol,
C return - p = minimum, fret - value in the min
C
      SUBROUTINE flatcher_Optim(n,p,ftol,iterMX,fret)
c
c      implicit none
      include 'xyzPDBsize.h' 
c
      include "output.h"
      include "kanalUse.h"
      include "statusMessg_mDyn.h"
c
      INTEGER n
      REAL fret,ftol,p(n)
c local
      INTEGER iter,ITMAX,iterMX
      real TOL,EPS
      integer NMAX
      PARAMETER (NMAX=natomMAX*3)    ! from xyzPDBsize.h
      PARAMETER(ITMAX=100,EPS=1.e-10)
      parameter (TOL = 1.e-2)   !for XYZshift in energyMinimization
C NMAX - max numb of variables
C
CU    USES dfunc,func_tomin,linmin
      INTEGER its,j
      REAL dgg,fp,gam,gg,g(NMAX),h(NMAX),xi(NMAX)
C loc
      logical CONTROL,CONTROL1
      integer k,kanalp
      integer ie,ix,i0,i1
      integer iterDO
      real sumgrad,shift
      real sumgradmin,shiftmin
      real p1(NMAX)

      CONTROL = .true. 
      CONTROL1 = .false.
      kanalp =   kanalRunOut
c
      sumgradmin = TOL  
      shiftmin = TOL*0.0001
      iterDO = iterMX
      if(iterDO .gt. ITMAX)iterDO=ITMAX
c
      i0 = 0  ! write flag
      i1 = 1
      ie = 0
c     
      write(kanalp,*)'fletcherOptimizer : starts! nOptXYZ:',n
c  
      if(n.gt.nmax)then
      write(kanalp,*)'fletcher : too low NMAX,n:',NMAX,n,
     & '  stop'
      stop
      end if
C!start:
      ix = i0
      call engOptFunct(p,fp,xi,ie,ix)
      ie = ie + 1
c returns function F(p) and gradF=xi(p)
c
      sumgrad=0.0
      do k=1,n
      sumgrad=sumgrad+abs(xi(k))
      end do
      sumgrad=sumgrad/n
c
      if(CONTROL)then
      write(kanalp,*)
     & 'Flatcher: initial:func  =',fp,' sumgrad=',sumgrad
      end if
      if(CONTROL1)then 
      write(kanalp,*)'Xinit: '
      call print9r(kanalp,n,p)
      write(kanalp,*)'GradInit:'
      call print9r(kanalp,n,xi)
      end if
C
      do 11 j=1,n
        g(j)=-xi(j)
        h(j)=g(j)
        xi(j)=h(j)
11    continue
        do 14 its=1,iterDO   ! opimization steps
        iter=its
C! variable shift
        do k=1,n
        p1(k)=p(k)
        end do !k

        call linmin(p,xi,n,fret)
      
        shift=0.0
        
        do k=1,n
        shift=shift+abs(p1(k)-p(k))
        end do !k
        shift=shift/n
ctt
	if(CONTROL)then
        write(kanalp,*)'Fletcher: iter=',its
        write(kanalp,*)'Fletcher: finit,fncLmin=',fp,fret
        write(kanalp,*)'Fletcher: dX = ',shift 
        write(kanalp,*)'sumgrad = ',sumgrad
        end if ! Cntrol

        if((shift.lt.shiftmin))then       
c   
      call engOptFunct(p,fp,xi,ie,i1)
c
        if(CONTROL)then
        write(kanalp,'(a14,i4,a18,f8.6)')
     &  'Fletchet:iter=',iter,' reach XshiftMin',shift  
        write(kanalp,*)
     &  'func  =',fp,' Xshift=',shift,' sumgrad:',sumgrad   
        end if!C
c
        return
        end if
c
Cyv    if((2.*abs(fret-fp).le.ftol*(abs(fret)+abs(fp)+EPS)))return
c    
        call engOptFunct(p,fp,xi,ie,i0)
        ie = ie + 1
c        
        sumgrad=0.0
        do k=1,n
        sumgrad=sumgrad+abs(xi(k))
        end do !k
        sumgrad=sumgrad/n
C! 
	 if(CONTROL1)then
         write(kanalp,*)'Fletcher: p()=',(p1(k),k=1,n)
         write(kanalp,*)'Fletcher: pm()=',(p(k),k=1,n)
         write(kanalp,*)'Fletcher: fret,func=',fret,fp
         write(kanalp,*)'Fletcher: grad=',(xi(k),k=1,n)
         write(kanalp,*)'Fletcher: sumgrad = ', sumgrad
         end if ! CONTROL

        if(sumgrad.lt.sumgradmin)then
c
      call engOptFunct(p,fp,xi,ie,i1)
c
        if(CONTROL)then
        write(kanalp,*)
     & 'Fletchet:iter=',iter,' reach gradmin=', sumgrad
        write(kanalp,*)'          func  =',fp,' Xshift=',shift  
        end if!C
        if(CONTROL1)then 
        write(kanalp,*)'Xfin: '
        call print9r(kanalp,n,p)
        write(kanalp,*)'GradFin:'
        call print9r(kanalp,n,xi)
        end if !Control1
c
        return      
        end if 
C!
        gg=0.
        dgg=0.
        do 12 j=1,n
          gg=gg+g(j)**2
C!        dgg=dgg+xi(j)**2            !this for Fletcher-Reeves
          dgg=dgg+(xi(j)+g(j))*xi(j)  !this for Polak-Ribiere, recomended
12      continue
        if(gg.eq.0.)return
        gam=dgg/gg
        do 13 j=1,n
          g(j)=-xi(j)
          h(j)=g(j)+gam*h(j)
          xi(j)=h(j)
13      continue
c
        iStatus = iStatus + 1            
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' engOptimization step = ',its, ' is done ...'
c
14    continue
c
      write(kanalp,*)'fletcher-rprmn maximum OptSteps exceeded'
      write(kanalp,*)'fletcher OptSt= ', its
c
      call engOptFunct(p,fp,xi,ie,i1)
c
      write(kanalp,*)'fletcher-optimization Finish!:'
c
      return
      END
C
C
C powel.f
C linmin.f
C mnbrak_min.f
C brent.f
C
c      SUBROUTINE powell(p,xi,n,np,ftol,iter,fret)
cC minimization, minimum,optimization.
cC minimization of func_tomin by Powell method
cC p(n) - initial variable, np-max dimesion in parameter(  )
cC xi(,) - initial matrix
cC ftol  - tolerance for function to stop minimiz
cC fret  - value of func in min point
cC func_tomin - minimized function
cC
c      implicit none
c      INTEGER iter,n,np,NMAX,ITMAX
c      REAL fret,ftol,p(np),xi(np,np),func_tomin
c      EXTERNAL func_tomin
c      PARAMETER (NMAX=50,ITMAX=200)
cCU    USES func_tomin,linmin
c      INTEGER i,ibig,j
c      REAL del,fp,fptt,t,pt(NMAX),ptt(NMAX),xit(NMAX)
c      fret=func_tomin(p)
c      do 11 j=1,n
c        pt(j)=p(j)
c11    continue
c      iter=0
c1     iter=iter+1
c      fp=fret
c      ibig=0
c      del=0.
c      do 13 i=1,n
c        do 12 j=1,n
c          xit(j)=xi(j,i)
c12      continue
c        fptt=fret
c        call linmin(p,xit,n,fret)
c        if(abs(fptt-fret).gt.del)then
c          del=abs(fptt-fret)
c          ibig=i
c        endif
c13    continue
c      if(2.*abs(fp-fret).le.ftol*(abs(fp)+abs(fret)))return
c      if(iter.eq.ITMAX) pause 'powell exceeding maximum iterations'
c      do 14 j=1,n
c        ptt(j)=2.*p(j)-pt(j)
c        xit(j)=p(j)-pt(j)
c        pt(j)=p(j)
c14    continue
c      fptt=func_tomin(ptt)
c      if(fptt.ge.fp)goto 1
c      t=2.*(fp-2.*fret+fptt)*(fp-fret-del)**2-del*(fp-fptt)**2
c      if(t.ge.0.)goto 1
c      call linmin(p,xit,n,fret)
c      do 15 j=1,n
c        xi(j,ibig)=xi(j,n)
c        xi(j,n)=xit(j)
c15    continue
c       goto 1
c      END
cC
      SUBROUTINE linmin(p,xi,n,fret)
C minimization, optimization, minimum.
c p() - point, xi(p) - grad(p) 
c      implicit none 
      include 'xyzPDBsize.h'
      INTEGER n
      REAL fret,p(n),xi(n)
c local
C isolates MIN with acc=TOL
      real TOL
      INTEGER NMAX
      PARAMETER (NMAX=3*natomMAX) ! from xyzPDBsize.h
      PARAMETER (TOL=1.e-5)   !endMin in XYZ
CU    USES brent,f1dim_min,mnbrak_min
      INTEGER j,ncom
      REAL ax,bx,fa,fb,fx,xmin,xx,pcom(NMAX),xicom(NMAX),brent
      COMMON /f1com_min/ pcom,xicom,ncom
      real lambinit, xim2
c USE function f1dim_min()
      real f1dim_min
      EXTERNAL f1dim_min
c
c      write(*,*)'linmin: starts!:'
c
c define initial step in 1-dim lambda minimization
      lambinit = 1.0e-1  ! dXYZ angstren
      ncom=n
      xim2 = 0.0
      do 11 j=1,n
        pcom(j)=p(j)
        xicom(j)=xi(j)
	xim2 = xim2 + xi(j)**2
11    continue
      ax=0.
C      xx=1.   !yv should be adjasted for given problem
      if(xim2 .gt. 1.0 )lambinit = lambinit/sqrt(xim2)
      xx = lambinit

c!c
cx      call mnbrak_min(ax,xx,bx,fa,fx,fb,f1dim_min)
      call mnbrak_min(ax,xx,bx,fa,fx,fb)
c!OK
      fret=brent(ax,xx,bx,f1dim_min,TOL,xmin)
c
      do 12 j=1,n
        xi(j)=xmin*xi(j)
        p(j)=p(j)+xi(j)
12    continue
c
c	write(*,*)'linmin finish:!'
c
      return
      end
C
cx      SUBROUTINE mnbrak_min(ax,bx,cx,fa,fb,fc,func_tomin)
       SUBROUTINE mnbrak_min(ax,bx,cx,fa,fb,fc)
C minimize givenNAME function f1dim_min, on interval: ax,bx;
Creturns : ax,bx,cx points which bracket the minimum: ax<bx<cx
        implicit none
cx        include 'xyzPDBsize.h'
C mininization, optimization, minimum.
        REAL ax,bx,cx,fa,fb,fc
        real f1dim_min
c local
c        integer NMAX
c        parameter (NMAX=natomMAX*3)    ! from xyzPDBsize.h
c        real gfn(NMAX)                 ! gradFunc 
        real GOLD,GLIMIT,TINY 
      PARAMETER (GOLD=1.618034, GLIMIT=100., TINY=1.e-20)
      REAL dum,fu,q,r,u,ulim
c
      fa=f1dim_min(ax)
      fb=f1dim_min(bx)
c
      if(fb.gt.fa)then
        dum=ax
        ax=bx
        bx=dum
        dum=fb
        fb=fa
        fa=dum
      endif
      cx=bx+GOLD*(bx-ax)
c
      fc=f1dim_min(cx) 
c
1     if(fb.ge.fc)then
        r=(bx-ax)*(fb-fc)
        q=(bx-cx)*(fb-fa)
        u=bx-((bx-cx)*q-(bx-ax)*r)/(2.*sign(max(abs(q-r),TINY),q-r))
        ulim=bx+GLIMIT*(cx-bx)
        if((bx-u)*(u-cx).gt.0.)then
c
       fu=f1dim_min(u) 
c
          if(fu.lt.fc)then
            ax=bx
            fa=fb
            bx=u
            fb=fu
            return
          else if(fu.gt.fb)then
            cx=u
            fc=fu
            return
          endif
          u=cx+GOLD*(cx-bx)
c
          fu=f1dim_min(u)
c
        else if((cx-u)*(u-ulim).gt.0.)then
c
          fu=f1dim_min(u)
cx          fu=func_tomin(u)
c
          if(fu.lt.fc)then
            bx=cx
            cx=u
            u=cx+GOLD*(cx-bx)
            fb=fc
            fc=fu
c
           fu=f1dim_min(u)
cx            fu=func_tomin(u)
c
          endif
        else if((u-ulim)*(ulim-cx).ge.0.)then
          u=ulim
c
         fu=f1dim_min(u)
cx          fu=func_tomin(u)
c
        else
          u=cx+GOLD*(cx-bx)
c
          fu=f1dim_min(u)
cx        fu=func_tomin(u)
c
        endif
        ax=bx
        bx=cx
        cx=u
        fa=fb
        fb=fc
        fc=fu
        goto 1
      endif
      return
      END
C
      FUNCTION brent(ax,bx,cx,f,tol,xmin)
C minimization,optimization,minimum
c  f(x) oneDim function
C given: points ax<bx<cx
C isolates the minium with accuracy = tol
      implicit none
c
      INTEGER ITMAX
      REAL brent,ax,bx,cx,tol,xmin
      REAL CGOLD,ZEPS
c
      REAL f
      EXTERNAL f
      PARAMETER (ITMAX=100)
      PARAMETER (CGOLD=.3819660,ZEPS=1.0e-10)
      INTEGER iter
      REAL a,b,d,e,etemp,fu,fv,fw,fx
      REAL p,q,r,tol1,tol2,u,v,w,x,xm
cyv
      real tolerba
c
      tolerba = tol
cyv
      a=min(ax,cx)
      b=max(ax,cx)
      v=bx
      w=v
      x=v
      e=0.
c
      fx=f(x)
c
      fv=fx
      fw=fx
c
      do 11 iter=1,ITMAX
        xm=0.5*(a+b)
        tol1=tol*abs(x)+ZEPS
        tol2=2.*tol1
        if(abs(x-xm).le.(tol2-.5*(b-a))) goto 3
cyv 
        if(abs(b-a).le.tolerba) goto 3
cyv
        if(abs(e).gt.tol1) then
          r=(x-w)*(fx-fv)
          q=(x-v)*(fx-fw)
          p=(x-v)*q-(x-w)*r
          q=2.*(q-r)
          if(q.gt.0.) p=-p
          q=abs(q)
          etemp=e
          e=d
          if(abs(p).ge.abs(.5*q*etemp).or.p.le.q*(a-x).or.p.ge.q*(b-x))
     *goto 1
          d=p/q
          u=x+d
          if(u-a.lt.tol2 .or. b-u.lt.tol2) d=sign(tol1,xm-x)
          goto 2
        endif
1       if(x.ge.xm) then
          e=a-x
        else
          e=b-x
        endif
        d=CGOLD*e
2       if(abs(d).ge.tol1) then
          u=x+d
        else
          u=x+sign(tol1,d)
        endif
c
        fu=f(u)
c
        if(fu.le.fx) then
          if(u.ge.x) then
            a=x
          else
            b=x
          endif
          v=w
          fv=fw
          w=x
          fw=fx
          x=u
          fx=fu
        else
          if(u.lt.x) then
            a=u
          else
            b=u
          endif
          if(fu.le.fw .or. w.eq.x) then
            v=w
            fv=fw
            w=u
            fw=fu
          else if(fu.le.fv .or. v.eq.x .or. v.eq.w) then
            v=u
            fv=fu
          endif
        endif
11    continue
C!yv      pause 'brent exceed maximum iterations'
3     xmin=x
      brent=fx
      return
      END
C
      FUNCTION f1dim_min(x)
c      implicit none
      include 'xyzPDBsize.h'
c
      INTEGER NMAX
      REAL f1dim_min,x
cx    REAL func_tomin
      PARAMETER (NMAX=natomMAX*3)    ! from xyzPDBsize.h
CU    USES func
      integer i0
      INTEGER j,ncom
      REAL pcom(NMAX),xicom(NMAX),xt(NMAX)
      real gfn(nMAX)
      real f1d
      COMMON /f1com_min/ pcom,xicom,ncom
c
      i0 = 0
c
      do 11 j=1,ncom
        xt(j)=pcom(j)+x*xicom(j)
11    continue
c
      call engOptFunct(xt,f1d,gfn,i0,i0)  
      f1dim_min = f1d
c
      return
      END
C
