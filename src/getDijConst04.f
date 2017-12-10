C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     getDijConst()                                                          *
C* Computes the dielectric constant  for atomPair i,j                         *
C*                                                                            *
C*     Yury Vorobjev 2003                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
        subroutine getDijConst04(atomDconst,D0,Ds,rMol,aPar,nPar,
     &                           ia,ja,dij1,DCij,dDCij)
c
c IN: ia,ja - atom N
c     atomDconst() - atomDconst(ia)
c     D0,Ds - molCore, bulkSolvent D constants
c     rMol -  radii of molecular sphera ~ rGyration
c     aPar,nPar - parameters of dielectric model Lazaridis-JCC2002 
c     dij1 - at-at distant
c OUT:
c      DCij - dielectric constant for pair i,j
c      dDCij - dD(r)/dr
c
	implicit none
        real atomDconst(*)
        real D0,Ds,rMol
        real aPar
        integer nPar 
        integer ia,ja
        real dij1	
        real DCij,dDCij
c local
        real DC1,dDC1
        real DC2,dDC2
        real Dcmax
        real a,b
        integer dType
        integer kanalp
        logical CONTROL
c
        kanalp = 6
        CONTROL = .false.           
c
        DCij = D0
        dDCij = 0.0
c
c estimate Dcmax
        Dcmax = max(atomDconst(ia),atomDconst(ja)) 
c
c surfaceAtoms 
        dType = 1
        a = (dij1/aPar)**nPar
        b = a/(1.0+a)
        DC1 = D0 + (Dcmax-D0)*b               ! dij const
        dDC1 = (Dcmax-D0)*b*nPar/(1.0+a)/dij1 ! dDij/dr
c 
c coreAtoms
        dType = 2
        a = 1.0/D0
        b = (a-1.0/Ds)/rMol
        if(dij1 .lt. rMol)then
	DC2 = 1.0/(a - b*dij1)              ! dij
        dDC2 = b*DCij**2                   !dDij/dr
        else 
        DC2 = Ds
        dDC2 = 0.0
        end if
c take a maximal dij for the pair
        if(DC2 .lt. DC1)then
         DCij = DC1
         dDCij = dDC1
         dType = 1
         else
         DCij = DC2
         dDCij = dDC2
         dType = 2
        end if
c
c
c         if(CONTROL)then
c         write(kanalp,'(a29,2i5,2f7.2)')
c     &  'getDijConst:i,j,atomDconsti:',
c     &   ia,ja,atomDconst(ia),atomDconst(ja)
c        write(*,*)' D0,R1,R2:',D0,R1,R2
c         write(kanalp,'(a22,2i5,7f7.1,1x,i4)')
c     &  'DijConst dij,Dij,dDij:',
c     &   ia,ja,dij1,DC1,dDC1,DC2,dDC2,DCij,dDCij,dType
c         end if !C
c
c        stop
c
        return
        end
