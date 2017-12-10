C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
C*     SolventEnForce()                                                       * 
C* Computes the solvent-solute Free energy and                                *
C* solvent-solute atomic forces                                               *
C* via Lazaridis & Karplus Gaussian shell                                     *
C* solvation model PROTEINS 1999, 35, 133-152                                 *
C* solute-solute coul energy and solute-solute coulombic forces               *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 

	subroutine SolventEnForces(natom, atomXYZ, 
     &         atomName,startPairL12,nPairL12,pair12List,
     &         nbpairList,startnbPairL,nnbPairL,
     &         atomSolPar, molSolEn, atomSolEn, atomSolFr) 

C*atomXYZ[3*natom] -  x,y,z = atomXYZ[3(i-1)+1],atomXYZ[3(i-1)+2],atomXYZ[3(i-1)+3]
c* nbpairList(i) -  consequtive List of neighbours for ia=1,...,natom  
C* startnbPairL(ia) -  startPos in nbpairList() which are 
C* neighb to ia */
C* nnbPairL(ia) - numbers of neighbours in the pairList() file nbpairList()  
C* atomSolPar[] - atomSolvationPar: delV,delGref,lambd,alfa,R,qsmod
C*RESULT:
C*   molSolEn - total solvation energy
C*   atomSolEn[ia]-atomic Solvation Energies
C*   atomSolFr[3*ia-2],[3*ia-1],[3*ia] - atomicSolvation Forces for atom ia

c        implicit none
        include "xyzPDBsize.h"
c
        integer  natom 
        real    atomXYZ(*) 
        character*4 atomName(*)
        integer startPairL12(*)
        integer nPairL12(*)
        integer pair12List(*)
        integer nbpairList(*) 
        integer startnbPairL(*) 
        integer nnbPairL(*) 
        real  atomSolPar(*) 
        real molSolEn, atomSolEn(*), atomSolFr(*) 
c
c        include "charStringSiz.h"
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c PWat:
        include 'enForce_PW.h'
        include 'optionPar.h'
c local
        integer i,i3,ia,i6,k
        integer j,ja,j3,j6
        real xi,fisol, xj,fjsol, forcei 
        real position_i(3), position_j(3) 
        real ijdist2,ijdist1,ijdist,ijdist3
        real frsol_i(3)  
        real vectorij(3)  
        real fsum(3)
        real esum
        real de,scale,scale1
        real dijm,dijm2,dijm3
        real cijd2
        real dr2inv,r2inv
c
        real mHyd,mHeavy
        real epsdist,sEbyG
        real scaleV
        integer nHMAX
        parameter (nHMAX = 6)
        integer atomH(nHMAX),nH
        logical CONTROL,CONTROL0
        logical OPT_Solv
        integer kanalp

        CONTROL = .false.  ! .true.  ! control print*/
        CONTROL0 = .false.  ! printTot En/force and STOP
cx	CONTROL0 = .true.
c
        OPT_Solv = .true.  ! .true. = include GSmodel 
c
        kanalp = kanalRunOut ! kanalp for outPut = std
        epsdist = 0.010    ! A prevent 1/dist error
c
        dijm = 1.00      ! A min ijdist to correct 1/r^2
        dijm2 = dijm**2
        dijm3 = dijm*dijm2
c
cx        scaleV = 0.75
cx         scaleV = 0.50    ! YV 
cx        scaleV = 0.60
cx          scaleV = 0.65
          scaleV = 1.0
c
C* atomSolPar( ) should be defined. */
C* atmChargeSol() - are neutralized (modified) set of charges 
C* to calculate Coul En */ 
C/*initialize */
        molSolEn = 0.0  
     	do i=1,natom 
        i6 = 6*(i-1)
     	atomSolEn(i) = atomSolPar(i6+2) 
     	i3=3*(i-1)  
        atomSolFr(i3+1) = 0.0 
        atomSolFr(i3+2) = 0.0 
        atomSolFr(i3+3) = 0.0 
        if(CONTROL)then
        write(kanalp,'(a6,2x,i5,6f8.3)')
     &  'SolPar:',i,(atomSolPar(i6+k),k=1,6)
        end if !control
        end do !i=1,natom
c
        molSolEn_PWat(1)=0.0
	molSolEn_PWat(2)=0.0
	molSolEn_PWat(3)=0.0
c
        if(.not.OPT_Solv)return

C* loop over all atoms ia  */   
          if(CONTROL)then
          fsum(1)=0.0
          fsum(2)=0.0
          fsum(3)=0.0
          write(kanalp,*)'in solvEnGShellSB:'
          end if
c
          do  ia=1,natom   
c
           i3 = 3*(ia-1) 
           i6 = 6*(ia-1)
c
           frsol_i(1)=0.0
           frsol_i(2)=0.0
           frsol_i(3)=0.0
c loop over neighbour ja atoms
          if(atomName(ia)(1:1) .ne. 'H' .and. 
     &        atomSolPar(i6+3) .ne. 0.0  )then  !nonHydrogen
c
          if(nnbPairL(ia) .ne. 0)then  !nonZero List
c 
           do j=startnbPairL(ia),startnbPairL(ia)+nnbpairL(ia)-1
           ja = nbpairList(j)
           j6 = 6*(ja-1)
c
          if(atomName(ja)(1:1) .ne. 'H' .and.
     &       atomSolPar(j6+3) .ne. 0.0 )then  !ja nonHydrogen
c
C/*ja - from special neigbour list which includes 12 and 13 vbonded atoms
C/* calculate contribution of atom ja to the solvation of iatom */
                 j3 = 3*(ja-1)   
                 vectorij(1)=(atomXYZ(i3+1) - atomXYZ(j3+1)) 
                 vectorij(2)=(atomXYZ(i3+2) - atomXYZ(j3+2)) 
                 vectorij(3)=(atomXYZ(i3+3) - atomXYZ(j3+3)) 

                ijdist2  = vectorij(1)**2 +  vectorij(2)**2  
     &                     +  vectorij(3)**2  + epsdist

                ijdist = sqrt(ijdist2)  
                ijdist3 = ijdist2*ijdist 

C* solvation energy */
                xi = (ijdist- atomSolPar(i6+5))/atomSolPar(i6+3) 
                if ( xi .lt.  0.0 ) xi = 0.0  
c function r2inv(r,d,d2) 
                cijd2 = r2inv(ijdist,ijdist2,dijm,dijm2)
c
                fisol = atomSolPar(i6+4)*exp(-xi*xi)*cijd2*scaleV 
c
                atomSolEn(ia) = atomSolEn(ia) - atomSolPar(j6+1)*fisol
c - - -
c check for overHydration
cx                if(atomSolEn(ia)/atomSolPar(i6+2) .lt. 0.0) goto 1001
cx ploho goto 1001
c stop the loop because all volume is excluded                   
c - - - 
c solvation forces */
                if(ijdist .gt. epsdist)then
                ijdist1 = 1.0/ijdist 
                else 
                ijdist1 = 2.0/(ijdist+epsdist)
                end if
c
                xj = (ijdist- atomSolPar(j6+5))/atomSolPar(j6+3) 
                if ( xj .lt.  0.0 ) xj = 0.0  
                fjsol = atomSolPar(j6+4)*exp(-xj*xj)*cijd2*scaleV
c
                forcei = ( fisol*(xi/atomSolPar(i6+3) 
     &             + dr2inv(ijdist,dijm) )*atomSolPar(j6+1)
     &             + fjsol*(xj/atomSolPar(j6+3) 
     &             + dr2inv(ijdist,dijm) )*atomSolPar(i6+1))*ijdist1

c collect forces on atom ia
            do   k=1, 3 
            frsol_i(k) = frsol_i(k) - forcei*vectorij(k)
            end do !k 
                
            if(CONTROL)then
         de = - atomSolPar(j6+1)*fisol
         write(kanalp,
     & '(a7,i4,1x,f8.3,a4,f8.3,a4,i4,a4,f8.3,2(a4,f8.3),4f8.3)')
     & 'SE:ia',ia,atomSolEn(ia),' de:',de,' ja:',ja,
     & ' ijd:',ijdist,' fis:',fisol,' fjs:',fjsol,forcei,
     & xi,xj
            end if!control
c
            end if ! ja nonHydr
            end do !} /*ja for*/
c
1001        continue
c     
c* collect total solvEn */
c overhydration: all volume excluded
             if(atomSolEn(ia)/atomSolPar(i6+2) .lt. 0.0)
     &          atomSolEn(ia) = 0.0
             molSolEn = molSolEn + atomSolEn(ia)  
c PWat: distribute
          if(OPT_SolvateExWat)
     &  call engPWatDistribute(atomPWatTag(ia),atomPWatTag(ia),
     &          atomSolEn(ia),molSolEn_PWat)
c
          end if !nonZero List
          end if !ia nonHydrogen
c
          if(CONTROL)then
          fsum(1) = fsum(1)+frsol_i(1)
          fsum(2) = fsum(2)+frsol_i(2)
          fsum(3) = fsum(3)+frsol_i(3)

          write(kanalp,'(i6,a3,3f8.3,a5,f8.3,a5,3f8.3)')
     &    ia,' f:',(frsol_i(k),k=1,3),' sEn:',atomSolEn(ia),
     &    'fsum:',fsum 
          end if !control 
c
          atomSolFr(i3+1) = frsol_i(1)          
          atomSolFr(i3+2) = frsol_i(2)          
          atomSolFr(i3+3) = frsol_i(3)          
          end do !} /* ia for*/
c
c for MD distribute forces over H connected to C,N,O
          mHeavy = 14 ! average C,N,O
          mHyd = 1
          do ia = 1,natom
          if(atomName(ia)(1:1) .ne. 'H')then
c find H atom attached to ia
          nH = 0
          do i=startPairL12(ia),startPairL12(ia)+nPairL12(ia)-1
          ja = pair12List(i)
          if( atomName(ja)(1:1) .eq. 'H')then
          nH = nH+1
          if(CONTROL)then
          if(nH .gt. nHMAX)then
          write(kanalp,*)'ERROR!:solvEnGShell: nHMAX is low'
          end if
          end if!CONTR
          atomH(nH) = ja
          end if! H 
          end do !i
          if (nH .ne. 0)then
c distribute forces
          i3 = 3*(ia-1)
          frsol_i(1) = atomSolFr(i3+1)
          frsol_i(2) = atomSolFr(i3+2)
          frsol_i(3) = atomSolFr(i3+3)
          scale = mHyd/(nH*mHyd + mHeavy)
          scale1 = mHeavy/(nH*mHyd + mHeavy)
c
          atomSolFr(i3+1) = frsol_i(1)*scale1   
          atomSolFr(i3+2) = frsol_i(2)*scale1   
          atomSolFr(i3+3) = frsol_i(3)*scale1
c
          do i = 1,nH
          ja=atomH(i)
          j3=3*(ja-1)  
          atomSolFr(j3+1) = frsol_i(1)*scale    
          atomSolFr(j3+2) = frsol_i(2)*scale    
          atomSolFr(j3+3) = frsol_i(3)*scale    
          end do! i
          end if ! nH >0

          end if !non H
          end do!ia
c
          if(CONTROL0)then
          write(kanalp,*)'solvGSmodEnForce:'
          esum = 0.0
          fsum(1) = 0.0 
          fsum(2) = 0.0 
          fsum(3) = 0.0 
          do ia=1,natom
c
          esum = esum + atomSolEn(ia)
c
          i3 = 3*(ia-1)
          frsol_i(1) = atomSolFr(i3+1)    
          frsol_i(2) = atomSolFr(i3+2)    
          frsol_i(3) = atomSolFr(i3+3)    
c
          fsum(1) = fsum(1)+frsol_i(1)
          fsum(2) = fsum(2)+frsol_i(2)
          fsum(3) = fsum(3)+frsol_i(3)
c
          if(atomSolPar(ia*6-4) .ne. 0.0)then
          sEbyG = atomSolEn(ia)/atomSolPar(ia*6-4)
          else 
          sEbyG = 0.0
          end if!
c
          write(kanalp,
     &    '(i6,a4,1x,a3,3f7.2,a5,f8.3,a7,f6.3,a6,3f7.3,a6,f8.2)')
     &    ia,atomName(ia),' f:',(frsol_i(k),k=1,3),
     &    ' sEn:',atomSolEn(ia), ' sE/G:', 
     &    sEbyG, ' fsum:',fsum,' eSt:',esum 

          end do !ia
c
          write(kanalp,'(33x,a10,f8.3)')'eSolvTot=',esum
          stop
c
          end if !CONTROL0
c
	  return 
          end
c
          real function r2inv(r,r2,d,d2)
c r2inv = f(r) = 1/r**2
          implicit  none
          real r,r2,d,d2
c
          if (r .ge. d) then
          r2inv = 1.0/r2
          return
          else
          r2inv = (3.0 - 2.0*r/d)/d2
          end if
c
          return
          end
c
          real function dr2inv(r,d)
c dr2inv = -(1/2)(df/dr)/f(r)
          implicit  none
          real r,d,d3
c
          if (r .ge. d) then
          dr2inv = 1.0/r
          return
          else
          dr2inv = 1.0/(3.0*d-2.0*r)
          end if
c 
          return
          end

