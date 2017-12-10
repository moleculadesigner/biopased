c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* modified 2003 aug.11                                                      *
c* soft VDW potential by option
c* subroutines for pair VWD   energy and Forces                              *
c*                                                                           *
c*     Yury Vorobjev 2002                                                    *
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine vdwenforceij(dij2,dij1,rij,A12,B12,evdw,fi)
c
c  rij = rj - ri
c  force on atom i
c
        implicit none
        real dij2,dij1,A12,B12
        real rij(3),evdw,fi(3)
c local
        real d6,fv
        real sc,sc2,sc6,evsc,fvsc
        integer kanalp
c
        sc = 1.2  ! softCore
        if(B12 .gt. 0.0 )then
        sc = 0.50*(2.0*A12/B12)**0.16667
        end if
c
        d6=dij2**3
 
        if( dij1 .gt. sc ) then
        evdw = (A12/d6 - B12)/d6
c
        fv = 6.0*(B12 - 2.0*A12/d6)/d6/dij2
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
c
        else
c
        sc2 = sc**2
        sc6 = sc2**3
        evsc = (A12/sc6 -  B12)/sc6               !vinesti v predrasch Table
        fvsc =  6.0*(B12 - 2.0*A12/sc6)/sc6/sc2   !vinesti v predrasch Table
c linear interpolation:
        evdw = evsc + fvsc*(dij1 - sc)
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
        return
        end
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
        subroutine vdwenforceijSC
     &               (dij2,dij1,rij,A12,B12,rsc12,evdw,fi)
c
c  rij = rj - ri
c  rs12 - soft core radius
c
        implicit none
        real dij2,dij1,A12,B12
        real rij(3),evdw,fi(3)
        real rsc12
c local
        real d6,fv
        real sc2,sc6,evsc,fvsc,A12dsc6
        integer kanalp
c
        d6=dij2**3
 
        if( dij1 .gt. rsc12 ) then
        evdw = (A12/d6 - B12)/d6
c
        fv = 6.0*(B12 - 2.0*A12/d6)/d6/dij2
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
c
        else
c
        sc2 = rsc12**2
        sc6 = sc2**3
        A12dsc6 = A12/sc6
        evsc = (A12dsc6 -  B12)/sc6               !vinesti v predrasch Table
c*      fvsc =  6.0*(B12 - 2.0*A12dsc6)/sc6/sc2   !vinesti v predrasch Table
        fvsc = -6.0*(evsc + A12dsc6/sc6)/sc2
c linear interpolation:
        evdw = evsc + fvsc*(dij1 - rsc12)
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
        return
        end
c
