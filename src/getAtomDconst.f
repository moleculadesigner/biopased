C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
C*     getAtomDconst()                                                        * 
C* Computes the local atomic dielectric constant                              *
C* solvent excluded approximation:
C* eps(ai) = [esp(mol)*Vmol + eps(solv)*Vsol]/(Vmol+Vsol)                     *
C* Vmol(at) - volume of molecular atoms                                       *
C* Vsol     - solvent volume in the first hydration sphera                    *
C* solute-solute Coul energy and solute-solute coulombic forces               *
C*                                                                            *
C*     Yury Vorobjev 2003                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

	subroutine getAtomDconst(epsMol,epsSol,natom,
     &          atomXYZ,atomSolPar,solMolRad,atomBlockName,
     &          startPairL,nPairL,pairList,atomDconst)

C* epsMol - molecular Dielectric constant (inside molecular volume)
C* epsSol - dielectric constant in solvent volume
C* natom - number of atoms                                                           
c* pairList()  -   consequtive List of neighbours for all atoms ia=1,...,natom  
C* startPairL(ia) -  startPos in pairList() which are neighbours to ia */
C* npairL(ia) - numbers of neighbours in the pairList() file for atom ia 
C* atomSolPar[] - atomSolvationPar: delV,delGref,lambd,alfa,R,qsmod
C* solMolRad - radius of solvent molecules
C*RESULT:
C*   atomDconst() - local Dielectric constant for atom ia
C*
        implicit none
c
        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        real epsMol,epsSol
        integer  natom 
        real atomXYZ(*)
        real atomSolPar(*)
        real solMolRad
        character*4 atomBlockName(*)
        integer startPairL(*)
        integer nPairL(*)
        integer pairList(*)
        real atomDconst(*)    
c local
       integer naPhGr,naPhnew
       real dCsum
       integer ia,ja,j
       integer nia,i3,j3,j6
       real pi,scpack
       real ri,rj,ra,rav
       real di,dj,dvi,dvj
       real Vih,Vim,Vis
       real vectorij(3)
       real ijdist,ijdist2
       logical CONTROL,CONTROL1
       integer kanalp
c
        CONTROL = .false. ! .true.  ! control print*/
        CONTROL1 = .false.
        kanalp = kanalRunOut                    
c
        pi = 3.14159
c       scpack = 0.8 ! packing koeff for hard spheres
        scpack = 0.75 
        rav = 1.5  ! average at rad
c
        if(CONTROL)then
        write(kanalp,*)'start:getAtomDconst '
        end if
C/*initialize */
     	do ia=1,natom 
        atomDconst(ia) = epsMol
        end do! ia
c
          do  ia=1,natom   
c
          i3 = 3*ia-3
          ra = (3.0*atomSolPar(6*ia-5)/(4.0*pi))**(1./3.)
c
c          if(ra .eq. 0.0 )ra = rav   !atom Hydr
           if(ra .lt. rav )ra = rav 
c
          ri = ra  + solMolRad      ! radius of sphere arround ia
          Vih = 4.0*pi*ri**3/3.0  - atomSolPar(6*ia-5) !totV
          Vim = 0.0
          nia = 0
c
          if(CONTROL1)then
          write(kanalp,*)'getAtomDconst:ia,dVa,ra,ri,Vih: ',ia,
     &    atomSolPar(6*ia - 5),ra,ri,Vih
          write(kanalp,*)'startPairL(ia),npairL(ia):',
     &                    startPairL(ia),npairL(ia)
          end if!C 
          if(nPairL(ia) .ne. 0)then  !nonZero List
c 
           do j=startPairL(ia),startPairL(ia)+npairL(ia)-1
           ja = pairList(j)
c
           if(CONTROL1)then
           write(kanalp,*)'getAtomDconst:j, pairListS(j):',
     &     j,pairList(j)
           end if
c
           j6 = 6*(ja-1)
           rj = (3.0*atomSolPar(6*ja - 5)/(4.0*pi))**(1./3.)
c
                 j3 = 3*(ja-1)   
                 vectorij(1)=(atomXYZ(i3+1) - atomXYZ(j3+1)) 
                 vectorij(2)=(atomXYZ(i3+2) - atomXYZ(j3+2)) 
                 vectorij(3)=(atomXYZ(i3+3) - atomXYZ(j3+3)) 

                ijdist2  = vectorij(1)**2 +  vectorij(2)**2  
     &                     +  vectorij(3)**2  

                ijdist = sqrt(ijdist2)  
c         
                dvi = 0.0
                dvj = 0.0 
c
                if(ijdist .le. (ri+rj) )then 
                nia = nia + 1
                if(ijdist .gt. (ri-rj))then
c calculate intersection between  sphera ri and rj
                di = (ri**2-rj**2 + ijdist2)/(2.0*ijdist)
                dj = (rj**2-ri**2 + ijdist2)/(2.0*ijdist)
c
                dvi = (2*ri**3-3.0*ri**2*di+di**3)*pi/3.0
                dvj = (2*rj**3-3.0*rj**2*dj+dj**3)*pi/3.0 
                else
                dvi = 0.0
                dvj = 4.0*rj**3/3.0
                end if
                end if            
c
           if(CONTROL1)then
           write(kanalp,*)'getAtomDconst:j, pairListS(j):',
     &     j,pairList(j), ' rj:',rj,' rij:',ijdist,' dvi,dvj:',dvi,dvj
           end if!C1
c
               Vim = Vim + dvi + dvj
c
            end do !  ja
            if(CONTROL1)then
            write(kanalp,*)'getAtomDconst: Vim,nia:',Vim,nia
            end if  !C
c correct excluded volume
            Vim = Vim/scpack
            if (Vim .gt. Vih)Vim = Vih
            Vis = Vih - Vim
c
            atomDconst(ia) = (epsMol*Vim + epsSol*Vis)/Vih
c    
          end if ! 
          if(CONTROL1)then
          write(kanalp,*)'getAtomDconst:ia,atomDconst:',
     &                    ia,atomDconst(ia)
          end if !C
c
          end do!ia
c
c calculate average D group for NA Ph group
          naPhGr = 4
          naPhnew = 0
          dCsum = 0.0
c
          do ia=1,natom
           if(atomBlockName(ia)(1:1) .eq. 'P')then
c add new atom
           naPhnew = naPhnew + 1
           dCsum = dCsum + atomDconst(ia) 
c distribute
           if(naPhnew .eq. naPhGr)then ! PhGroup is complete
           dCsum = dCsum/naPhGr
           atomDconst(ia) = dCsum
           atomDconst(ia-1) = dCsum
           atomDconst(ia-2) = dCsum
           atomDconst(ia-3) = dCsum
           dCsum = 0.0
           naPhnew = 0
           end if! .eq. naPhGr
           end if! P
           end do!ia
c
          if(CONTROL)then
          write(kanalp,*)'getAtomDconst:ia, Dcost:'
          do ia=1,natom
          write(kanalp,*)ia,atomDconst(ia)
          end do !ia
          write(kanalp,*)'final:getAtomDconst '
          end if !CONTROL
c
	  return 
          end
