c inititValecity in [A/ps] for a given T [K]
c
c Yu. Vorobjev 2002
c
	subroutine initVelocity(temp,natom,
     &                    nmoveatom,moveAtomList,atomMass,atomVel)
c
	implicit none
	integer natom
        integer nmoveatom,moveAtomList(*)
        real temp
        real atomMass(*)
c
	real atomVel(*)
c
        include "charStringSiz.h"
        include "output.h" 
c
        integer IG
        integer n,i,ia,iam
        integer i3,k,ir,ntot
        integer kanalp
        logical CONTROLcm,CONTROL
        real V(3),va,mtot
	real sc,sdi,sd
        real ta,sa,ra,low
        real sct,scs,sccal
        integer Ndf
c
        kanalp = kanalRunOut 
        CONTROLcm = .true.
        CONTROL = .false.
        n = 3
        IG = 1007
c init
        low = 0.00005
        sct = 1.2028     ! scaling coeff T[k]=sct*MwV^2/Nsfredom
        scs = 1.0/sct
        sccal = 2.390    ! sccal*Mw*V^2/2 = KinEn[cal]
        Ndf = 3*nmoveatom
c
        va = 0.0
        do i=1,3*natom
        atomVel(i) = 0.0
        end do
c
        if (temp .le. low )return
c
        iam = 0
100     iam = iam +1 
c choos atom randomly
        call random(ra,IG)
        ir = ra*nmoveatom + 1
        if (ir .gt. nmoveatom ) ir = nmoveatom
        ia = moveAtomList(ir)
        i3 = ia*3
c
        sd = sqrt(temp*scs/atomMass(ia))

        call gauss (va,sd,n,V,IG) 
        atomVel(i3-2)=V(1)
        atomVel(i3-1)=V(2)
        atomVel(i3)  =V(3)

        if(iam .lt. nmoveatom) goto 100
c
c some atoms my have not been assigned
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3 = ia*3
        if(atomVel(i3-2) .eq. 0.0 .and. 
     &     atomVel(i3-1) .eq. 0.0 .and. 
     &       atomVel(i3) .eq. 0.0 ) then

        call gauss (va,sd,n,V,IG)
        atomVel(i3-2)=V(1)
        atomVel(i3-1)=V(2)
        atomVel(i3)  =V(3)
        end if
        end do !iam
c
c check C.M. move and adjust
        V(1)=0.0
        V(2)=0.0
        V(3)=0.0
        mtot = 0.0
	do iam = 1,nmoveatom
        ia = moveAtomList(iam) 
        i3 = ia*3-3
        do k=1,3
        V(k) = V(k) + atomMass(ia)*atomVel(i3+k)
        end do !k
        mTot = mTot +  atomMass(ia)
        end do!iam
c
         do k=1,3
         V(k) = V(k)/mTot
         end do !
c
         write(kanalp,*)'initVelocity: initial: VcM:',V
c
        do iam = 1,nmoveatom
        ia = moveAtomList(iam)
        i3 = ia*3-3 
        do k=1,3
        atomVel(i3+k)=atomVel(i3+k) - V(k)                        
        end do !k       
        end do !iam
c check velocities
        if(CONTROLcm)then
        write(kanalp,*)'initVelocity: corrected on cMVel: T[K]:',temp
        write(kanalp,*)'   ia    vx     vy      vz    [A/ps]'
        V(1)=0.0
        V(2)=0.0
        V(3)=0.0
        ta = 0.0
	do i = 1,natom
        i3 = i*3-3
        do k=1,3
        V(k) = V(k) + atomMass(i)*atomVel(i3+k)
        ta = ta + atomMass(i)*atomVel(i3+k)**2
        end do !k
        if(CONTROL)then
        write(kanalp,'(i6,3f8.3)')
     &  i, (atomVel(i3+k), k=1,3)
        end if !C
        end do!i
c
        do k=1,3
        V(k) = V(k)/mtot
        end do!k
c
        write(kanalp,*)'initVelocity: final: VcM:',V
c
	ta = sct*ta/Ndf
        sc = sqrt(temp/ta)
        write(kanalp,*)'initVelocity: Tinit[K], targT:', ta,temp,
     &   ' ScailCf:',sc
        write(kanalp,*)'initVelocity: VCmass: ', V    
c correct V by scaling
        do i = 1,natom*3
        atomVel(i) = atomVel(i)*sc
        end do 
        end if!CONTROLcm
c
	return
	end
