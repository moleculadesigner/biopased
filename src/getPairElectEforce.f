c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c* modified 2003 aug.11                                                      *
c* soft VDW potential by option                                              * 
c* subroutines for pair  ELECTrostatic energy and Forces                     *
c*                                                                           *
c*     Yury Vorobjev 2002                                                    *
c*                   2005             GenBorn                                * 
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c standart Coulon law
c
        subroutine
     &  coulenforceij0(dij2,dij1,rij,qi,qj,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real dij2,dij1,qi,qj
        real rij(3),ecoul,fi(3)
c local
        real qij,fv
        real sc,fvsc,ecsc,sc2
        sc = 1.20    ! min dji1(A) to approximate Coul by linear
c
        qij = qi*qj
c simple COUL
        if(dij1 .gt. sc )then
        ecoul = qij/dij1
        fv = -ecoul/dij2
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small dist
        ecsc = qij/sc
        fvsc = -ecsc/sc**2
        ecoul = ecsc + (dij1 - sc)*fvsc*sc
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
        return
        end
c
        subroutine
     &  coulenforceij0SC(dij2,dij1,rij,qi,qj,rsc,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real dij2,dij1,qi,qj,rsc
        real rij(3),ecoul,fi(3)
c local
        real qij,fv
        real fvsc,ecsc
c        rsc = rscVDW ! min dji1(A) to approximate Coul by linear
c
        qij = qi*qj
c simple COUL
        if(dij1 .gt. rsc )then
        ecoul = qij/dij1
        fv = -ecoul/dij2
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small dist
        ecsc = qij/rsc
        fvsc = -ecsc/rsc**2
        ecoul = ecsc + (dij1 - rsc)*fvsc*rsc
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
        return
        end
c
c charge on the uniform background + eps=1
        subroutine
     &  coulenforceij01(rc1,rc2,dij2,dij1,rij,qi,qj,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real rc1,rc2
        real dij2,dij1,qi,qj
        real rij(3),ecoul,fi(3)
c local
        real qij,fv,e1
        real sc,fvsc,ecsc,e1sc,sc2
c
        sc = 1.20    ! min dji1(A) to approximate Coul by linear
c
        qij = qi*qj
c
c charge on the uniform background
        if(dij1 .gt. sc)then
        e1 = 1.0/dij1
        ecoul = e1 + (dij2*rc2 - 3.0)*rc1*0.5
        ecoul =  ecoul*qij
        fv = -qij*(e1/dij2 - rc1*rc2)
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small dist
        e1sc = 1.0/sc
        sc2 = sc**2
        ecsc = e1sc + (sc2*rc2 - 3.0)*rc1*0.5
        ecsc = ecsc*qij
        fvsc = -qij*(e1sc/sc2 - rc1*rc2)
        ecoul = ecsc + (dij1 - sc)*fvsc*sc
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
        return
        end
c
c charge on the uniform background + eps=1
        subroutine
     &  coulenforceij01SC(rc1,rc2,dij2,dij1,rij,qi,qj,rsc,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real rc1,rc2
        real dij2,dij1,qi,qj,rsc
        real rij(3),ecoul,fi(3)
c local
        real qij,fv,e1
        real fvsc,ecsc,e1sc,sc2
c
c       rsc = 1.20    ! min dji1(A) to approximate Coul by linear
c
        qij = qi*qj
c
c charge on the uniform background
        if(dij1 .gt. rsc)then
        e1 = 1.0/dij1
        ecoul = e1 + (dij2*rc2 - 3.0)*rc1*0.5
        ecoul =  ecoul*qij
        fv = -qij*(e1/dij2 - rc1*rc2)
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small dist
        e1sc = 1.0/rsc
        sc2 = rsc**2
        ecsc = e1sc + (sc2*rc2 - 3.0)*rc1*0.5
        ecsc = ecsc*qij
        fvsc = -qij*(e1sc/sc2 - rc1*rc2)
        ecoul = ecsc + (dij1 - rsc)*fvsc*rsc
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
        return
        end
c
c charge on the uniform background + eps=(r*cReps)
c
        subroutine
     &  coulenforceij02(cReps,rc1,rc2,dij2,dij1,rij,qi,qj,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real cReps
        real rc1,rc2
        real dij2,dij1,qi,qj
        real rij(3),ecoul,fi(3)
c local
        real qij,fv,e1
        real sc,fvsc,ecsc,e1sc,sc2
c
        sc = 1.2  ! linear appriximation
        qij = qi*qj
c
c charge on the uniform background + eps=(r*cReps)
        if(dij1 .gt. sc)then
        e1 = 1.0/dij1
        ecoul = e1 + (dij2*rc2 - 3.0)*rc1*0.5
        ecoul =  ecoul*qij/dij1/cReps
c
        fv = -qij*(e1/dij2 - rc1*rc2)*e1/cReps
        fv = fv - ecoul/dij2
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small distance
        e1sc = 1.0/sc
        sc2 = sc**2
        ecsc = e1sc + (sc2*rc2 - 3.0)*rc1*0.5
        ecsc = ecsc*qij*e1sc/cReps
c
        fvsc = -qij*(e1sc/sc2 - rc1*rc2)*e1sc/cReps
        fvsc = fvsc - ecsc/sc2
        ecoul = ecsc + (dij1 - sc)*fvsc*sc
c
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
c        write(*,*)'cReps,rc1,rc2,dij2,dij1,rij,qi,qj:',
c     &            cReps,rc1,rc2,dij2,dij1,rij,qi,qj 
c        write(*,*)'ecoul:',ecoul,' fi',fi
c        write(*,*)'coulenforceij02: end'
c	stop
c
        return
        end
c
c charge on the uniform background + eps=(r*cReps)
c
        subroutine  coulenforceij02SC
     &  (cReps,rc1,rc2,dij2,dij1,rij,qi,qj,rsc,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real cReps
        real rc1,rc2
        real dij2,dij1,qi,qj,rsc
        real rij(3),ecoul,fi(3)
c local
        real qij,fv,e1
        real fvsc,ecsc,e1sc,sc2
c
c       rsc = 1.2  ! linear appriximation
        qij = qi*qj
c
c charge on the uniform background + eps=(r*cReps)
        if(dij1 .gt. rsc)then
        e1 = 1.0/dij1
        ecoul = e1 + (dij2*rc2 - 3.0)*rc1*0.5
        ecoul =  ecoul*qij/dij1/cReps
c
        fv = -qij*(e1/dij2 - rc1*rc2)*e1/cReps
        fv = fv - ecoul/dij2
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small distance
        e1sc = 1.0/rsc
        sc2 = rsc**2
        ecsc = e1sc + (sc2*rc2 - 3.0)*rc1*0.5
        ecsc = ecsc*qij*e1sc/cReps
c
        fvsc = -qij*(e1sc/sc2 - rc1*rc2)*e1sc/cReps
        fvsc = fvsc - ecsc/sc2
        ecoul = ecsc + (dij1 - rsc)*fvsc*rsc
c
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
c        write(*,*)'cReps,rc1,rc2,dij2,dij1,rij,qi,qj:',
c     &            cReps,rc1,rc2,dij2,dij1,rij,qi,qj 
c        write(*,*)'ecoul:',ecoul,' fi',fi
c        write(*,*)'coulenforceij02: end'
c	stop
c
        return
        end
c
        subroutine
     &  coulenforceij03(eij,deij,rc1,rc2,dij2,dij1,rij,qi,qj,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real rc1,rc2
        real eij,deij                  !epsIJ, de(r)/dr
        real dij2,dij1,qi,qj
        real rij(3),ecoul,fi(3)
c local
        real qij,fv,e1
        real sc,fvsc,ecsc,e1sc,sc2
c        logical CONTROL
c        CONTROL = .true.
c
        sc = 1.2  ! linear appriximation
        qij = qi*qj
c
c charge on the uniform background + atomDconst(r), r- dependent
        if(dij1 .gt. sc)then
        e1 = 1.0/dij1
        ecoul = e1 + (dij2*rc2 - 3.0)*rc1*0.5
        ecoul =  ecoul*qij/eij
c
        fv = -qij*(e1/dij2 - rc1*rc2)
        fv = (fv - ecoul*deij*e1)/eij
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small distance
        e1sc = 1.0/sc
        sc2 = sc**2
        ecsc = e1sc + (sc2*rc2 - 3.0)*rc1*0.5
        ecsc = ecsc*qij/eij
c
        fvsc = -qij*(e1sc/sc2 - rc1*rc2)
        fvsc = (fvsc - ecsc*deij*e1sc)/eij
        ecoul = ecsc + (dij1 - sc)*fvsc*sc
c
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
c       if(CONTROL)then
c       write(*,*)'coulenforceij03:'
c       write(*,*)'eij,deij : ',eij,deij
c       write(*,*)'rc1,rc2,dij2,dij1:',rc1,rc2,dij2,dij1
c       write(*,*)'rij,qi,qj:',rij,qi,qj
c       write(*,*)'ecoul,fi:',ecoul,fi
cx       stop
c       end if !C
c
        return
        end
c
        subroutine  coulenforceij03SC
     &  (eij,deij,rc1,rc2,dij2,dij1,rij,qi,qj,rsc,ecoul,fi)
c
c  rij = rj - ri
c
        implicit none
        real rc1,rc2
        real eij,deij                  !epsIJ, de(r)/dr
        real dij2,dij1,qi,qj,rsc
        real rij(3),ecoul,fi(3)
c local
        real qij,fv,e1
        real fvsc,ecsc,e1sc,sc2
c        logical CONTROL
c        CONTROL = .true.
c
c       rsc = 1.2  ! linear appriximation
        qij = qi*qj
c
c charge on the uniform background + atomDconst(r), r- dependent
        if(dij1 .gt. rsc)then
        e1 = 1.0/dij1
        ecoul = e1 + (dij2*rc2 - 3.0)*rc1*0.5
        ecoul =  ecoul*qij/eij
c
        fv = -qij*(e1/dij2 - rc1*rc2)
        fv = (fv - ecoul*deij*e1)/eij
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
        else
c small distance
        e1sc = 1.0/rsc
        sc2 = rsc**2
        ecsc = e1sc + (sc2*rc2 - 3.0)*rc1*0.5
        ecsc = ecsc*qij/eij
c
        fvsc = -qij*(e1sc/sc2 - rc1*rc2)
        fvsc = (fvsc - ecsc*deij*e1sc)/eij
        ecoul = ecsc + (dij1 - rsc)*fvsc*rsc
c
        fi(1)=fvsc*rij(1)
        fi(2)=fvsc*rij(2)
        fi(3)=fvsc*rij(3)
        end if
c
c       if(CONTROL)then
c       write(*,*)'coulenforceij03:'
c       write(*,*)'eij,deij : ',eij,deij
c       write(*,*)'rc1,rc2,dij2,dij1:',rc1,rc2,dij2,dij1
c       write(*,*)'rij,qi,qj:',rij,qi,qj
c       write(*,*)'ecoul,fi:',ecoul,fi
cx       stop
c       end if !C
c
        return
        end
c pair term for GBorn approximation
c + charge on the uniform background + epsIN=1,espW=espW
c rc1=Rcut, rc2=rc1**2
        subroutine
     &  electPforceijGBSC(rc1,rc2,dij2,dij1,
     &              rij,qi,qj,rsc1,RBi,RBj,eng,fi)
c
c  rij = rj - ri
c
        implicit none
        include 'coulEnPar.h'
        real rc1,rc2
        real dij2,dij1,qi,qj,rsc1
        real RBi,RBj
        real rij(3),eng,fi(3)
c local
        real epsSOL,epsMOL
        real qij
        real fv,fv1,fv2,rsc2
        real gbex,e1,engsc,rbp
        real fGB1,fGB2,nbg,dnbg,dfGB,ksi,smal
        logical CONTROL
        integer kanalp
c
        CONTROL = .false.
        kanalp = 6
c
c       rsc = 1.20    ! min dji1(A) to approximate Coul by linear
c
        smal = 0.000001
        epsSOL = epsSol_OPT
        epsMOL = epsMol_OPT
c
        qij = qi*qj
        ksi = 1.0/epsSOL - 1.0/epsMOL
        rbp = RBi*RBj
c
c charge on the uniform background
c
        if(dij1 .gt. rsc1)then
c GB function
        gbex = exp(-dij2/(4.0*rbp))
        fGB2 = dij2 + rbp*gbex
        fGB1 = sqrt(fGB2)
        dfGB = (1.0 - 0.25*gbex)*dij1/fGB1
        nbg = 1.0/dij1 + 0.5*(dij2*rc2 - 3.0)*rc1
        dnbg = -1.0/dij2 +dij1*rc2*rc1
c
        e1 = ksi*dij1/fGB1 + 1.0/epsMOL
        eng = e1*nbg*qij
c
        fv1 = e1*dnbg
        fv2 = ksi*(1.0/fGB1 - dij1*dfGB/fGB2)*nbg
c
        fv = qij*(fv1 + fv2)/dij1
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
c
        else
c small dist
        rsc2 = rsc1**2
        gbex = exp(-rsc2/(4.0*rbp))
        fGB2 = rsc2 + rbp*gbex
        fGB1 = sqrt(fGB2)
        dfGB = (1.0 - 0.25*gbex)/fGB1
        nbg = 1.0/rsc1 + 0.5*(rsc2*rc2 - 3.0)*rc1
        dnbg = -1.0/rsc2 +rsc1*rc2*rc1
c
        e1 = ksi*dij1/fGB1 + 1.0/epsMOL
        engsc = qij*e1*nbg
c
        fv1 = engsc*dnbg
cx        fv2 = ksi*(1.0/fGB1 - dij1*dfGB/fGB2)*nbg
cx        fv = qij*(fv1 + fv2)/(dij1 + smal)
        fv2 = ksi*(1.0/fGB1 - rsc1*dfGB/fGB2)*nbg
        fv = qij*(fv1 + fv2)/rsc1             
c
        eng = engsc + (dij1 - rsc1)*fv*rsc1
        fi(1)=fv*rij(1)
        fi(2)=fv*rij(2)
        fi(3)=fv*rij(3)
c
        end if
c
        if(CONTROL)then
        write(kanalp,*)
     &  'getPelectEforceGB:rc1,rc2,dij2,dij1,qi,qj,rsc1,epsSOL,RBi,RBj:'
        write(kanalp,*)
     &  rc1,rc2,dij2,dij1,qi,qj,rsc1,epsSOL,RBi,RBj
         write(kanalp,*)'rij:', rij 
        write(kanalp,*)'eng:',eng, ' fi:', fi
        end if !C
c
        return
        end
c
c selfatom ai GB energy  and force
c
        subroutine selfEPolEngGBorn(qi,Rbi,dRbi,eng,fi)
c
        implicit none        
        include 'coulEnPar.h'
        real qi,Rbi
        real dRbi(3)
        real eng,fi(3)
c
        real epsSOL,epsMOL
        real bornSelfEngScale_OPT
        real ksi,a,sckal
        integer k
c
        epsSOL = epsSol_OPT
        epsMOL = epsMol_OPT
c
        bornSelfEngScale_OPT = 1.00   !OK'
c
        sckal = 166.0*bornSelfEngScale_OPT  !332.0 
        ksi = 1.0/epsSOL - 1.0/epsMOL
c
        a = sckal*ksi*qi**2
        eng = a/Rbi
c
        do k =1,3
        fi(k) = -a*dRbi(k)
        end do !k
c
        return
        end 
