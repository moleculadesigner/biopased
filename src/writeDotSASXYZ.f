c
	subroutine writeDotSASXYZ(iv)
c
c iv ; 0/1 : no/yes indepenedent SAS calculation
c
c write SolventAccessibleSurface SAS of solute molecule
c
	include'xyzPDBsize.h'
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
        include 'xyzPDBinfo.h'
        include 'xyzPDBcrd.h'
        include 'dataSASdr.h'
c
        integer iv
c
        real ss
        integer k,i,ia,ia3,id1,id3
        integer ip,ndotMX
        integer kanalp,kanalXYZsas
c
        kanalp = kanalRunOut
        kanalXYZsas = kanalwWatBrg
c
        ip=1
        ndotMX = ndotMAX
        dotden = 4.0
        dprobe(ip) = solMolRad_OPT  ! water
c
c the SAS is alredy calculated in the bornRad subroutine
c
        if(iv .ge. 1)then
               call  surf_SAS04(atomXYZ,atomRad,
     &                 natom,dprobe(ip),dotden,
     &                 dotXYZ(1,ip),dotnrm(1,ip),dotarea(1,ip),
     &                 dot_IATNUM(1,ip),ndot(ip),ndotMX,
     &                 nsurfAt(ip),nsurfAtList(1,ip),
     &                 atSurfAr(1,ip),atSurfNrm(1,ip),
     &                 atSurfXYZ(1,ip),atSurfProj(1,ip),
     &                 bindSiteAt01(1,ip),
     &                 dot_jATNUM(1,ip),
     &                 head_dotNum,linkListDotAt,
     &                 atDistToSAS,atSASexp,
     &                 dot_eField,dot_ePot)
c
       end if
c
        open(unit=kanalXYZsas,file=fileSASdotXYZ,form='formatted',
     &       status='unknown')
c
        write(kanalXYZsas,'(a23,a19,f5.3,a14,f5.3)')
     &  '#SolvAccessibleSurface:',
     &  '    dotden per A^2=',dotden,'   dprobe(A)=',dprobe(ip)
c
        write(kanalXYZsas,'(a30,a55)')'         id            iat     ',
     &  '   dotXYZ                dotArea     ELcPot(kcal/mol) '
c
	ss=0.0
        do id1=1,ndot(ip)
        id3=3*id1
        write(kanalXYZsas,7003)
     &  'ATOM ',id1,'DOT ',
     &  atomName(dot_IATNUM(id1,ip)),dot_IATNUM(id1,ip),
     &  dotXYZ(id3-2,ip),dotXYZ(id3-1,ip),dotXYZ(id3,ip),dotarea(id1,ip)
     &  ,dot_ePot(id1)
     	ss=ss + dotarea(id1,ip)
        end do!id1
c
         write(kanalXYZsas,'(a4)')'TER '
         write(kanalXYZsas,'(a4)')'END '
	 write(kanalXYZsas,'(a18,f8.2)')'#SASarea(A**2)MOL:',ss
c
c
        write(kanalXYZsas,'(a38)')
     &  '#SolvAccessibleSurface Area per atoms:'
        write(kanalXYZsas,'(a65)')
     & '#PDB:                                             atSurfAr(A^2)' 
        do ia = 1,natom
        ia3=3*ia-3
c
        write(kanalXYZsas,7003)
     &  'ATOM ',ia,atomName(ia),resName(ia),
     &   resNumb(ia),(atomXYZ(ia3+k),k=1,3),atSurfAr(ia,ip)
c
        end do!ia
c
        write(kanalXYZsas,'(a4)')'TER '
        write(kanalXYZsas,'(a4)')'END '
c
        close(kanalXYZsas)
c
	return
7003   format(a5,1x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2x,1x,2f8.3) ! PDB+s
	end
