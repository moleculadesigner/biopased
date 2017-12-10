c calculate SAS restrainWatShell Potential/Forces
c
c Yury Vorobjev 2009
c
c restrain water molecules in SAS solvation shell
c
	subroutine get_expLSASPot02(natomLg,atomLgXYZ,
     &                                  engExpL,fatExpL)
c
c INP: natomLg, atomLgXYZ( )
c      indxSASatNeigh(iaLg) = indxNearestProtAt in atSurfXYZP01(*),atSurfNrmP01(*)
c      atSurfXYZP01(*)  : protSAS atoms XYZ
c      atSurfNrmP01(*)  : protSAS atoms NormalVectors
c OUT: ddexWatSASat(*) : dist of Lig atom penetration into SAS
c      engExpL       : energy of expulsion Pot
c      fatExpL(*)    : atomExpuLForces
c      ddSASexWatShMX : OPTion = max ddPenetration after that expLPot=linear function
c      KexLoc    : OPTion = harmonic constant
c
cx       implicit none
       include 'xyzPDBsize.h'
cx       include "charStringSiz.h"
       include "output.h"
       include "kanalUse.h"
       include "statusMessg_mDyn.h"
       include "solvate01.h"
       include "exWatShRestrSASPot.h"
       include "enForce_PW.h"
c
       integer natomLg
       real atomLgXYZ(*)
       real engExpL,fatExpL(*)
c local
       integer iSLT
       integer kanalp
       real KexLoc
       real dWatShLoc
       real ddSShL_OPT
       integer i,i3,k,iL,iL3,ip,ip3
       real es,fex(3)
       real ex,ds
       logical OPT_WrotatSymm
       integer nWat
       logical CONTROL
c
       kanalp=kanalRunOut
       CONTROL = .false. ! .true.
c
       if(CONTROL)then
       write(*,*)'get_expLSASpot02: STARt!'
       end if
       if(CONTROL)then
       write(kanalp,*)'natomLg,atomLgXYZ:'
       do i=1,natomLg
       i3=i*3-3
       write(kanalp,*)'ATOMlg',i,i,(atomLgXYZ(i3+k),k=1,3),
     & indxSASatNeighExWatSh(i)
       end do!i
       write(kanalp,*)'nsurfSASP01:',nsurfSASP01,
     &   ' atSurfXYZP01() atSurfNrmP01():'
       do i=1,nsurfSASP01
       i3=i*3-3
       write(kanalp,*)'ATOMPr',i,i,(atSurfXYZP01(i3+k),k=1,3),
     & (atSurfNrmP01(i3+k),k=1,3)
       end do!i
       do i=1,nsurfSASP01
       i3=i*3-3
       write(*,*)'ATOMPr',i,i,(atSurfXYZP01(i3+k),k=1,3),
     & (atSurfNrmP01(i3+k),k=1,3)
       end do!i
       end if!C
c
       ddSShL_OPT = 5.0 ! 4.5 !  3.5   !
       dWatShLoc = dSoLVSHELLvar - ddSShL_OPT
       KexLoc=0.5*KhSASexWatSh
       ddSASexWatShMX = 1.5 ! A
       es=0.0
c
       iSLT=3*natomSolutMol
c init ZERO ff_exWatSASP()
       do i=1,iSLT+3*natomLg
       fatExpL(i) = 0.0
       end do !i
c
       OPT_WrotatSymm=.true.
       nWat=3
       do iL=1,natomLg
c water molec: O, H1,H2 - fixed atom order
c forceces on O are not eq ZERO: to have rotational symmetry
c
       if(OPT_WrotatSymm.and.(nWat*int((iL-1)/nWat).eq.
     &      (iL-1)))then
c
       iL3=iL*3-3
       ip=indxSASatNeighExWatSh(iL)
       ip3=3*ip-3
       ds=0.0
       do k=1,3
       ds=ds + (atomLgXYZ(iL3+k)-atSurfXYZP01(ip3+k))
     &       *atSurfNrmP01(ip3+k)
       end do!k
       ds = ds - dWatShLoc
c
       if(ds .le. 0.0)then   ! watMol is inside solvationSHELL
       ex = 0.0
       fex(1)=0.0
       fex(2)=0.0
       fex(3)=0.0
       else
       if(ds .le. ddSASexWatShMX)then
       ex = ds*ds*KexLoc
       fex(1)=-KexLoc*ds*atSurfNrmP01(ip3+1)
       fex(2)=-KexLoc*ds*atSurfNrmP01(ip3+2)
       fex(3)=-KexLoc*ds*atSurfNrmP01(ip3+3)
       else
       ex = ddSASexWatShMX*KexLoc*(2.0*ds - ddSASexWatShMX) 
       fex(1) = -KexLoc*ddSASexWatShMX*atSurfNrmP01(ip3+1)
       fex(2) = -KexLoc*ddSASexWatShMX*atSurfNrmP01(ip3+2)
       fex(3) = -KexLoc*ddSASexWatShMX*atSurfNrmP01(ip3+3)
       end if!ds .le. ddSASexWatShMX
       end if !ds .le. 0.0
       es = es +ex
       fatExpL(iSLT+iL3+1) = fex(1)
       fatExpL(iSLT+iL3+2) = fex(2)
       fatExpL(iSLT+iL3+3) = fex(3)
c
       ddexWatSASat(iL) = ds
c
       if(CONTROL)then
       write(kanalp,'(a15,2i5,2x,f6.2,1x,3f8.3)')
     &       'iL,iwp,ds,ff:',iL,ip,ds,fex
       write(kanalp,'(a6,3f7.2,a6,3f7.2,a4,f7.2,a4,3f7.2)')
     &  'LgXYZ:',(atomLgXYZ(iL3+k),k=1,3),
     & ' pSAS:',(atSurfXYZP01(ip3+k),k=1,3),' ar:',atSurfArP01(ip),
     & 'Sn:',(atSurfNrmP01(ip3+k),k=1,3)
       end if
       end if! OPT_WrotatSymm
       end do !iL
c
       engExpL = es
c
       if(CONTROL)then
       write(*,*)'get_expLSASpot02: FINish!'
       end if
c
       return
       end
