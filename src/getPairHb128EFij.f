c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* Hx...Y Hb12-8 H-bond potential                                            *
c* soft VDW potential 
c*                                                                           *
c*     Yury Vorobjev 2004                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine getPairHb128EFij
     &               (dij2,dij1,rij,rm12,em12,A12,B12,ehb,fi)
c
c  rij = rj - ri
c - 12-8 potential - !? is not good, too shortRange
c - 12-6 - ?
c  rm12 - Rmin of potential
c  em12 - Emin
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
