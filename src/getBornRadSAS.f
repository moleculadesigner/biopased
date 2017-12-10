c calculation of the BornRad from SAS03
c
	subroutine getBornRadSAS(atomXYZ,atomRad,
     &              natom,
     &              dotXYZ,dotnrm,dotarea, 
     &              ndot,
     &              nsurfAt,nsurfAtList,
     &              atSurfAr,atSurfNrm,atSurfXYZ,
     &              head_dotNum,linkListDotAt,
     &              atBornRad,atBornRadDr) 
c
c
         implicit none
         real atomXYZ(*)
         real atomRad(*)
         integer natom,ndot,nsurfAt
         real dotXYZ(*),dotnrm(*),dotarea(*) 
         real atSurfAr(*),atSurfNrm(*)
         real atSurfXYZ(*)
         real atBornRad(*),atBornRadDr(*)
         integer nsurfAtList(*)
         integer head_dotNum(*),linkListDotAt(*)
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c local
         real probeRad_OPT
         real bornRadCorr_OPT
         real bornRadCorr
         logical atomNB,atomFA
         logical CONTROL,CONTROL1
         integer kanalp
c
         real rBnb,rbOut,rbTot
         real pi4,rcutRB2
         real rcutRB_OPT
         real aixyz(3),ajixyz(3)
         real rBnbDr(3),rBOutDr(3)
         real djai(3),dji2,p1,djai2,sj
         real djai4,a0,a1
cx         real scaleBornRad_OPT
         integer i,j,k,ia,ia3
         integer jd,ja,ja3,jas,jd3
c
         probeRad_OPT = 1.4  !dataSASdr.h
         bornRadCorr_OPT = 0.5
         bornRadCorr = bornRadCorr_OPT - probeRad_OPT
c         
cx         scaleBornRad_OPT = 1.0          
         pi4 = 12.5663708
c
cx         kanalp = 6
         kanalp = kanalRunOut
c
         CONTROL = .false.
         CONTROL1 = .false.
c
cx         rcutRB_OPT=6.0
         rcutRB_OPT=1000.0
         rcutRB2 = rcutRB_OPT**2
c
         do ia = 1,natom
         rBnb =0.0
         rbOut=0.0
         do k=1,3
         rBnbDr(k)=0.0
         rbOutDr(k)=0.0
         end do!k
         atomNB = .false.
         atomFA = .false.
         ia3=3*ia-3
         aixyz(1) = atomXYZ(ia3+1)
         aixyz(2) = atomXYZ(ia3+2)
         aixyz(3) = atomXYZ(ia3+3)
c
         atBornRadDr(ia3+1) = 0.0
         atBornRadDr(ia3+2) = 0.0
         atBornRadDr(ia3+3) = 0.0
c
         do jas=1,nsurfAt
         ja=nsurfAtList(jas)
         ja3=3*ja-3
c
         ajixyz(1)= atomXYZ(ja3+1) - aixyz(1) 
         ajixyz(2)= atomXYZ(ja3+2) - aixyz(2)
         ajixyz(3)= atomXYZ(ja3+3) - aixyz(3)
c
         dji2 = ajixyz(1)**2+ajixyz(2)**2+ajixyz(3)**2
c
         if(dji2 .le. rcutRB2)then
         atomNB = .true.
         else 
         atomNB = .false.
         end if !dji2
         atomFA = (.not. atomNB)
c
         sj = 0.0
c
         if(atomNB)then
c sum over neighbour ja atomSurf
c start from head id0t
         jd = head_dotNum (ja)
c
100      if(jd .gt. 0) then
         jd3 = 3*jd-3
         djai(1)=dotXYZ(jd3+1)-aixyz(1)       
         djai(2)=dotXYZ(jd3+2)-aixyz(2)
         djai(3)=dotXYZ(jd3+3)-aixyz(3)
c
         djai2 = djai(1)**2+djai(2)**2
     &           + djai(3)**2
c
         djai4 = djai2**2
c
         p1 = djai(1)*dotnrm(jd3+1)+
     &        djai(2)*dotnrm(jd3+2)+
     &        djai(3)*dotnrm(jd3+3)  
c
         a0 = dotarea(jd)/djai4
         a1 = p1*a0
         rBnb = rBnb + a1                     
c
c rBnbDr: 
c
         do k=1,3
         rBnbDr(k) = rBnbDr(k) - dotnrm(jd3+k)*a0
     &                 + 4.0*a1*djai(k)/djai2
         end do !k
c
         if(CONTROL1)then
         sj = sj + dotarea(jd)
         write(kanalp,*)'jd,dotarea(jd),djai2,p1,sj,rBnb:',
     &   jd,dotarea(jd),djai2, p1,sj,rBnb
         end if !C
c
         else 
         goto 101
         end if !jd .gt. 0
c next jd
         jd = linkListDotAt(jd)
         goto 100
c
101      continue  
c goto next jas Atom
         if(CONTROL1)then
         write(kanalp,*)
     &   'getBornRadSAS: Sja(extracted): ',sj,atSurfAr(ja)
         end if !C
c
         end if!atomNB
c         
         if(atomFA)then
c use atomPatch as surface element: 
         djai(1)=atSurfXYZ(ja3+1)-aixyz(1)
         djai(2)=atSurfXYZ(ja3+2)-aixyz(2)
         djai(3)=atSurfXYZ(ja3+3)-aixyz(3)
c
         djai2 = djai(1)**2+djai(2)**2
     &           + djai(3)**2
c
         djai4 = djai2**2
c
         p1 = djai(1)*atSurfNrm(ja3+1)+
     &        djai(2)*atSurfNrm(ja3+2)+
     &        djai(3)*atSurfNrm(ja3+3)
c
         a0 = 1.0/djai4
         a1 = p1*a0
c
         rbOut = rbOut + a1   
c
c rBOutDr:
c
         do k=1,3
         rBOutDr(k) = rBOutDr(k) - atSurfNrm(ja3+k)*a0
     &                          +4.0*a1*djai(k)/djai2
         end do !k
c
         end if! atomFA
c
         end do!jas
c
         rbTot = rBnb + rbOut 
c
         atBornRad(ia) = pi4/rbTot + bornRadCorr
c
         do k=1,3
         atBornRadDr(ia3+k) = (rBnbDr(k) + rBOutDr(k))/pi4
         end do!k    
c
         if(CONTROL)then
         write(kanalp,*)'getBornRadSAS: ia,atomRad(ia),atBornRad(ia):',
     &   ia,atomRad(ia),atBornRad(ia),
     &   ' rbTot,rBnb,rbOut:',rbTot,rBnb,rbOut
         write(kanalp,*)'getBornRadSAS: BrDR: nb=',rBnbDr
         write(kanalp,*)'getBornRadSAS: BrDR: Out=',rBOutDr
         write(kanalp,*)
     &   'getBornRadSAS: BrDR: BornRadDrTot=',(atBornRadDr(ia3+k),k=1,3)
         end if!C
c
         end do !ia

         return
         end
