C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*     SolventAtmParModif()                                                   *
C* MoDify   the solvent atomic parameters                                     *
C* for Lazaridis & Karplus Gaussian shell                                     *
C* solvation model PROTEINS 1999, 35, 133-152                                 *
C*                                                                            *
C*     Yury Vorobjev 2002                                                     *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
	subroutine solvAtmParModif(imodif,
     &                natom,atmNameEx,atmResName,
     &                atmHydphSolv,atmSolPar)
c*
c* InPut:
c  imodif - flag to do modification 0/1
C* natom
C* atmNameEx(ia) - atomNames Extended
C* atmResName(ia) - resName for atom ia
c  atmHydphSolv - value of hydroPhobic solvation to be added
c  atmSolPar() - initial atomic Solvation Parameters
c
c* OUTPut: modified solvaton parameters
C* atmSolPar(6*natom): ia: delV,delGref,lamb,alfa(*),R,qmod
C*                            qmod-neutralized set of atmic Q
        implicit none
        integer imodif
        integer natom
        real atmSolPar(*)
        real atmHydphSolv
        character*8 atmNameEx(*)     !extendedName atm+termTag
        character*4 atmResName(*)
c
        include "charStringSiz.h"
        include "output.h"
c*      
        real v23,r42
        integer ia,i6,i
        integer kanalp
        logical CONTROL,CONTROL0

c
        kanalp =  kanalRunOut
c
        if(imodif .eq. 0 ) then
        write(kanalp,*)'solvAtmParModif : NO '
        return
        end if
c
c CONTROL printOut
        CONTROL0 = .true.  ! .true.
        CONTROL = .false.
c
        if(CONTROL)then
       write(kanalp,*)
       write(kanalp,*)'In solvAtmParModif: initial solvPar:'
       write(kanalp,*)
c 
       write(kanalp,'(a6,1x,a8,a4,6(a5,3x))')
     & 'atomNb','atomName','res','   dV  ','dGref ','lamb ','alfa',
     & 'Ra  ','qat '
       do ia = 1,natom
        i6 = 6*ia-6
        write(kanalp,'(i6,1x,a8,a4,6f8.3)')
     &  ia,atmNameEx(ia),atmResName(ia),
     &  (atmSolPar(i6+i), i=1,6)
        end do
       end if !C
c do modification
       do ia = 1,natom
       i6=(ia-1)*6
       r42 = 1.0
       if(atmSolPar(i6+2) .ne. 0.0 )
     & r42 = atmSolPar(i6+4)*atmSolPar(i6+3)/(0.0898*atmSolPar(i6+2))
       v23 = 0.0
       if(atmSolPar(i6+1) .ne. 0.0)
     & v23 = atmSolPar(i6+1)**0.66667
c
       atmSolPar(i6+2) = atmSolPar(i6+2) + 
     &           atmHydphSolv*v23   
c       
       if(atmSolPar(i6+3) .ne. 0.0)
     & atmSolPar(i6+4) = atmSolPar(i6+4) + 
     &           atmHydphSolv*v23*r42*0.0898/atmSolPar(i6+3)
c
       end do !ia       
c
c control print
       if(CONTROL0)then
       write(kanalp,*)
       write(kanalp,*)'solvAtmParModif: modified SolvPar:'
       write(kanalp,*)
     & 'SolvationGSmodel Param + HydrPhSolv:',atmHydphSolv*6.5
c 
       write(kanalp,'(a6,1x,a8,a4,6(a5,3x))')
     & 'atomNb','atomName','res','   dV  ',
     & '    dGref ','lamb ','alfa',
     & 'Ra  ','qat '
       do ia = 1,natom
        i6 = 6*ia-6
        write(kanalp,'(i6,1x,a8,a4,6f8.3)')
     &  ia,atmNameEx(ia),atmResName(ia),
     &  (atmSolPar(i6+i), i=1,6)
        end do
c       
       write(kanalp,*)'End of solvAtmParModif: '
       end if !C
c
	return
        end
