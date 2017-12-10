c add H-atoms to PDB atoms to PDB file,  required: data atmhtag.dat 
c                                              file to be processed
c                            in to the assign_htg.h data structure 
c	
        subroutine add_hatom
     &   (natm0,head0,atmnam0,resnam0,resnumb0,chnam0,atmxyz0,
     &    nres0,startA0,stopA0,resnam_r0,
     &    natm1,head1,atmnam1,resnam1,resnumb1,chnam1,atmxyz1,
     &    nres1,startA1,stopA1,resnam_r1)
c
c..................................................................
        implicit none
        include "charStringSiz.h"
        include "assign_htg.h"
c
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c assign_htg.h : parameters
       integer natm0,nres0
       integer resnumb0(*)
       character*6 head0(*)
       character*4 atmnam0(*)
       character*4 resnam0(*)
       character*4 resnam_r0(*)
       character*1 chnam0(*)
       integer startA0(*),stopA0(*)
       real*8 atmxyz0(3,*)
       integer natm1,nres1
       integer resnumb1(*)
       character*6 head1(*)
       character*4 atmnam1(*)
       character*4 resnam1(*)
       character*4 resnam_r1(*)
       character*1 chnam1(*)
       integer startA1(*),stopA1(*)
       real*8 atmxyz1(3,*)
c...................................................................
        integer naddh,naddhtot
	integer hvar
	integer nJKL
        integer nline,natmtot
	character*4 atmJKLai(3)
	character*4 hnameloc(3)
	real*8 atmJKLxyz(3,3)
	real*8  hatxyz(3,3)
        logical CONTROL
        logical CONTROL0
        integer i,j,k
	integer iiprint0,iiprint1
	integer kanalpl
c initialize
	natm1 = 0
	naddhtot = 0
	natmtot = 0
c
cx        kanalpl = 6	
        kanalpl = kanalRunOut
c
c       CONTROL = .false.    ! regular
        CONTROL = .true.     !local print
        CONTROL0 = .true.
c
	if(CONTROL)then
	write(kanalpl,*)'add_hatom : start: natm0:', natm0
        write(kanalpl,*)'nres0:', nres0
	write(kanalpl,*)'          i,resnam_r0(i), startA0(i),stopA0(i):'
	do i = 1,nres0
	write(kanalpl,*)i,resnam_r0(i),startA0(i),stopA0(i)
	end do!i

	do i=1,natm0
	if(atmrecn(i).gt.0)then
	write(kanalpl,'(a25,3i5,1x,a9,3a5,a7,3a5)')
     &	'in addhat:atN,resN,atrec:',i,resnumb0(i),atmrecn(i),
     &  'hname_tb:',(hname_tb(k,atmrecn(i)),k=1,3),
     &  'atmJKL:',(atmJKL(k,atmrecn(i)), k=1,3)
	end if
	end do !i
	end if !CONTROL
c
	nres1 = nres0
	do i=1,nres0
	resnam_r1(i)=resnam_r0(i)
        startA1(i) = 0
        stopA1(i) = 0
	end do !i
c add Hatoms
	do i = 1,natm0
c
c write atom to result file
        natm1 = natm1 + 1
	head1(natm1) = head0(i)
	atmnam1(natm1) = atmnam0(i)
	resnam1(natm1) = resnam0(i)
	resnumb1(natm1) = resnumb0(i)
	chnam1(natm1) = chnam0(i)
	atmxyz1(1,natm1) = atmxyz0(1,i)
	atmxyz1(2,natm1) = atmxyz0(2,i)
	atmxyz1(3,natm1) = atmxyz0(3,i)

	nline = atmrecn(i)

c h-add section
        hvar = 0
	if(atmrecn(i).ne.0)then
	hvar = htag_tb(nline)
	hnameloc(1)=hname_tb(1,nline)
	hnameloc(2)=hname_tb(2,nline)
	hnameloc(3)=hname_tb(3,nline)
	end if

	if(CONTROL)then
	write(kanalpl,*)
	write(kanalpl,*)
     &	'in addhatom:AtN,resN,nline:',i,resnumb0(i),atmrecn(i)
	if(nline.gt.0)then
	write(kanalpl,*)' hnameloc:',hnameloc, 'atmJKL:',
     &  (atmJKL(k,nline), k=1,3)
	end if
	end if
c
	if(hvar.ne.0)then
c
c define coordinates atmJKLxyz(): atom I and J,K,L defines hatxyz attached to atom I
c processing atmJKL() from data assign_htg.h
        naddh = 0  
        if(hvar .eq. 1 .or. hvar .eq. 2 .or. hvar .eq. 3 .or.
     &     hvar .eq. 4 .or. hvar .eq. 5 .or. hvar .eq. 6 .or.
     &     hvar .eq. 7) then
	do k = 1,3
	atmJKLai(k) = atmJKL(k,nline)
	end do
c
        call def_atmJKLxyz(i,nline,natm0,atmnam0,resnam0,resnumb0,
     &         	 chnam0,atmxyz0,nres0,startA0,stopA0,
     &           atmJKLai,nJKL,atmJKLxyz)

c test
	if(CONTROL0)then
	write(kanalpl,*)'add_hatom:'
	write(kanalpl,*)'hvar:', hvar,' hnameloc:', hnameloc
	write(kanalpl,*)'PDBfile OrigiN atom:'
        write(kanalpl,7071)'ATOM  ',natm1,
     &  	atmnam1(natm1),resnam1(natm1),chnam1(natm1),
     &          resnumb1(natm1),
     &          atmxyz1(1,natm1),atmxyz1(2,natm1),atmxyz1(3,natm1)
        end if !C
c calculate h-atoms xyz       
        call hatm_xyz(hvar,atmnam0(i),resnam0(i),resnumb0(i),chnam0(i),
     &          atmxyz0(1,i),atmJKLxyz,hnameloc,naddh,hatxyz)
        end if ! hvar.ne.8
c
        if(hvar .eq. 8)then
        call hatmWat_xyz(atmnam0(i),resnam0(i),resnumb0(i),chnam0(i),
     &         atmxyz0(1,i),hnameloc,naddh,hatxyz)
        end if !hvar.eq.8
	
	if(naddh.gt.0)then
	naddhtot = naddhtot + naddh
	do j=1,naddh
	natm1 = natm1 + 1
	head1(natm1) = head0(i)
	atmnam1(natm1) = hnameloc(j)
	resnam1(natm1) = resnam0(i)
	resnumb1(natm1) = resnumb0(i)
	chnam1(natm1) = chnam0(i)
	atmxyz1(1,natm1) = hatxyz(1,j)
	atmxyz1(2,natm1) = hatxyz(2,j)
	atmxyz1(3,natm1) = hatxyz(3,j)
c test
	if(CONTROL0)then
	write(kanalpl,*)'Atom i:',i,' Added Hatom:',naddh,naddhtot
        write(kanalpl,7071)'ATOMh:',natm1,
     &	    atmnam1(natm1),resnam1(natm1),chnam1(natm1),
     &      resnumb1(natm1),
     &      atmxyz1(1,natm1),atmxyz1(2,natm1),atmxyz1(3,natm1)
        end if !C
c
	end do !j
	end if !naddh.gt.0
c
	end if !hvar ne 0
c
c startA1(), stopA1()
        stopA1(resnumb1(natm1)) = stopA0(resnumb0(i)) + naddhtot	
c
        write(kanalpl,*)'addhatom: iat:',i,
     & ' stopA1(resnumb1(natm1)):',stopA1(resnumb1(natm1))
c
	end do ! i

c correct 
	startA1(1) = startA0(1) 
	do i = 2,nres1+1
	startA1(i)=stopA1(i-1) + 1
	end do !i

c
	if(CONTROL0)then
	write(kanalpl,*)'add_hatom : finish: TotAtom natm_h:', natm1
        end if
c
        if(CONTROL)then
        do i = 1,natm1
        write(kanalpl,7071)'ATOM  ',i,
     &          atmnam1(i),resnam1(i),chnam1(i),
     &          resnumb1(i),
     &          atmxyz1(1,i),atmxyz1(2,i),atmxyz1(3,i)
        end do!i
c
        write(kanalpl,*)'add_hatom: nres1:', nres1
	write(kanalpl,*)
     &	'  i, resnam_r1(i),startA1(i),stopA1(i) startA0(i),stopA0(i):'
	do i = 1,nres1
	write(kanalpl,*)i,resnam_r1(i),startA1(i),stopA1(i)
     &	                 ,startA0(i),stopA0(i)
	end do!i
	end if !control

c
707    format(a6,1x,i4,2x,a4,a8,a1,i4,4x,3f8.3) !orig PDB
7071   format(a6,1x,i4,2x,a4,a4,a1,i4,4x,3f8.3) !orig PDB
c
        return
	end
c----------------------------------------------------------------------71
c
       subroutine def_atmJKLxyz(ia,recn,natm0,atmnam0,resnam0,resnumb0,
     &         chnam0,atmxyz0,nres0,startA0,stopA0,
     &         atmJKLai,nJKL,atmJKLxyz)

c defines JKLxyz for JKL atoms , neighbours of atomI = ia, to define local geometry
c recn - linenumb in tables of the data structure assign_htg.h
c natm0,atmnam0,resnam0,resnumb0,chnam0,atmxyz0 : original PDB(noHatom) info 
c atmJKLai(3): names of atomIJK related with line recn for atomI
c ..................................................................
       implicit none
c
        include "charStringSiz.h" 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
       integer ia,recn,natm0
       integer resnumb0(natm0)
       character*4 atmnam0(natm0)
       character*4 resnam0(natm0)
       character*1 chnam0(natm0)
       integer nres0
       integer startA0(nres0),stopA0(nres0)
       real*8 atmxyz0(3,natm0)
       character*4 atmJKLai(3)
c result:
       integer nJKL
       real*8 atmJKLxyz(3,3)
c......................................................................
c local
       character*4 compaund,prevres,pusto,prot
       character*4 atnamloc
       integer i,j,k,nr,nrm1,n
       logical CONTROL
       integer kanalpl

       compaund = "PROT"
       prot = "PROT"
c
       CONTROL = .true.
cx     kanalpl = 6
       kanalpl = kanalRunOut   

       if(CONTROL)then
       write(kanalpl,*)'start def_atmJKL :'
       end if
       
       if(compaund.eq.prot)then
       pusto   = "    "
       prevres = "C-1 "
       nr = resnumb0(ia)
       nrm1 = nr-1
c
c nJKL - number of JKL atoms: 2 or 3
       nJKL=0
       do k = 1,3
       if(atmJKLai(k).ne.pusto)then
       n = nr
       if(atmJKLai(k).eq.prevres .and. nr .gt. 1) n = nrm1
c
       if(CONTROL)then
       write(kanalpl,*)
     & 'In def_atmJKLxyz:addHtabAt(lookFor):k,atmJKLai(k):',
     &  k,atmJKLai(k) 
       end if
c define atNumbers in this residue or res-1
        do i = startA0(n),stopA0(n)

c special for resn=1,to add to CA,  reduce NT->N ; C-1 -> C 
       atnamloc = atmnam0(i)
       
       if(atmJKLai(k) .eq. prevres .and. 
     &       atnamloc .eq. 'C    ' ) atnamloc = prevres

c compare Ia atom with the tableList
        if(CONTROL)then
        write(kanalpl,*)'In def_atmJKLxyz: PDBati;',i, atnamloc
        end if
c
        if(atnamloc .eq. atmJKLai(k) )then 
        nJKL=nJKL+1
        atmJKLxyz(1,k) = atmxyz0(1,i)
        atmJKLxyz(2,k) = atmxyz0(2,i)
        atmJKLxyz(3,k) = atmxyz0(3,i)
ctest
	if(CONTROL)then
	write(kanalpl,*)'In def_atmJKLxyz:'
        write(kanalpl,'(a22,i6,2a6)')'i,atnamloc,atmnam0(i):',
     &  i,atnamloc,atmnam0(i)
	write(kanalpl,'(a15,i6, a6)')'k, atmJKLai(k):',
     &	k,atmJKLai(k)
	end if

	end if
	end do !i

	end if ! ne.pusto
	end do !k
     
ctest
	if(CONTROL)then
	write(kanalpl,*)'In def_atmJKLxyz: END:'
	write(kanalpl,*)'nJKL: ', nJKL
        end if 
      
        end if !prot

       if(CONTROL)then
       write(kanalpl,*)'END def_atmJKL :'
       write(kanalpl,*)
       end if
       
       return
	end
c ----------------------------------------------------------------------------
c defines H atoms xyz attached to the atom: atmnam,resnam,resnumb,chnam,atmxyz(3)
        subroutine hatm_xyz(hvar,atmnam,resnam,resnumb,chnam,
     &          atmxyz,atmJKLxyz,hnameloc,nhadd,hatxyz)

c hvar - flag defines method of attachment
c atmnam,resnam,resnumb,chnam - atom I name, residuename, resnumb,chain name
c atmxyz(3) - coord of atom I
c atmJKLxyz(3,3) - coords of atom J,K,L
c naddh - number of H atoms added
c hname_tb(3) - names of H atoms
c hatxyz(3,3) - coords of H atoms added to structure
	implicit none
c
        include "charStringSiz.h" 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
        integer hvar 
	character*4 atmnam
cx	character*8 resnam
        character*4 resnam
	character*1 chnam
	real*8 atmxyz(3),atmJKLxyz(3,3)
	integer nhadd,resnumb
	character*4 hnameloc(3)
	real*8 hatxyz(3,3)
c local
	real*8 bondh,sint,cost,cos30
	real*8 hxyzloc(3,3)
	real*8 nulV(3),vH1(3),vH2(3),vH3(3)
	real*8 hxyzlc(3,3),s,PI
	real*8 vIJ(3), vIK(3), VIL(3), vJK(3)
	real*8 vX(3),vY(3),vZ(3)
	real*8 trot(3,3)
	integer i,j,k,kanalpl
	logical CONTROL

cx	kanalpl = 6
        kanalpl = kanalRunOut
c
	CONTROL = .false.
c        CONTROL = .true.
	bondh = 1.04d0
	PI = 3.141592654d0
	sint = 0.8166415d0*bondh   !sin(109.5/2)
	cost = 0.5771452d0*bondh
	cos30 = dcos(PI/6.0d0)
	nulV(1)=0.0d0
	nulV(2)=0.0d0
	nulV(3)=0.0d0
c
        if(CONTROL)then
        write(kanalpl,*)'in  hatm_xyz: hvar: ',hvar,' atomI:',
     &   atmnam,resnam,resnumb,chnam
        write(kanalpl,*)'atmxyzI:',atmxyz
        write(kanalpl,*)'atomxyzJ:',(atmJKLxyz(k,1),k=1,3)
        write(kanalpl,*)'atomxyzK:',(atmJKLxyz(k,2),k=1,3)
        write(kanalpl,*)'atomxyzL:',(atmJKLxyz(k,3),k=1,3)
        end if !C
c
c H-atoms positions for different cases=hvar
	if(hvar.eq.1)then
	nhadd = 1
        call vectNrmh(atmxyz,atmJKLxyz(1,1),vIJ)
        call vectNrmh(atmxyz,atmJKLxyz(1,2),vIK)
        call vectNrmh(atmxyz,atmJKLxyz(1,3),vIL)
        vIJ(1)=(vIJ(1)+vIK(1)+vIL(1))/3.0
        vIJ(2)=(vIJ(2)+vIK(2)+vIL(2))/3.0
        vIJ(3)=(vIJ(3)+vIK(3)+vIL(3))/3.0
        call vectNrmh(nulV,vIJ,vH1)
        hatxyz(1,1)=atmxyz(1) - vH1(1)*bondh
        hatxyz(2,1)=atmxyz(2) - vH1(2)*bondh
        hatxyz(3,1)=atmxyz(3) - vH1(3)*bondh
c
	return
	end if ! hvar=1

	if(hvar.eq.2)then
	nhadd = 2
        call vectNrmh(atmxyz,atmJKLxyz(1,1),vIJ)
        call vectNrmh(atmxyz,atmJKLxyz(1,2),vIK)
        vH1(1)=(vIJ(1)+vIK(1))/2.0
        vH1(2)=(vIJ(2)+vIK(2))/2.0
        vH1(3)=(vIJ(3)+vIK(3))/2.0
        call vectNrmh(vH1,nulV,vH1)

	call vectPrh(VIJ,VIK,vH2)
	call vectNrmh(nulV,vH2,vH2)

        hatxyz(1,1)=atmxyz(1) + (vH1(1)*cost + vH2(1)*sint)
        hatxyz(2,1)=atmxyz(2) + (vH1(2)*cost + vH2(2)*sint)
        hatxyz(3,1)=atmxyz(3) + (vH1(3)*cost + vH2(3)*sint)

        hatxyz(1,2)=atmxyz(1) + (vH1(1)*cost - vH2(1)*sint)
        hatxyz(2,2)=atmxyz(2) + (vH1(2)*cost - vH2(2)*sint)
        hatxyz(3,2)=atmxyz(3) + (vH1(3)*cost - vH2(3)*sint)
	
	return
	end if ! hvar=2

	if(hvar.eq.3)then
	nhadd = 3
c hxyz in local system
	sint = 0.8166415d0*bondh   !sin(109.5/2)
	cost = 0.5771452d0*bondh
	hxyzloc(1,1)= sint/cos30
	hxyzloc(2,1)= 0.0d0
	hxyzloc(3,1)= dsqrt(bondh**2-hxyzloc(1,1)**2)
	hxyzloc(3,2)= hxyzloc(3,1)
	hxyzloc(3,3)= hxyzloc(3,1)
	hxyzloc(1,2)= hxyzloc(1,1)-2.0*sint*cos30
	hxyzloc(1,3)= hxyzloc(1,2)
	hxyzloc(2,2)= sint
	hxyzloc(2,3)= -sint

c test 
	if(CONTROL)then
	write(kanalpl,'(a10,3f8.3)')
     &	'hxyzloc1: ',hxyzloc(1,1),hxyzloc(2,1),hxyzloc(3,1) 
	write(kanalpl,'(a10,3f8.3)')
     &	'hxyzloc1: ',hxyzloc(1,2),hxyzloc(2,2),hxyzloc(3,2) 
	write(kanalpl,'(a10,3f8.3)')
     &	'hxyzloc1: ',hxyzloc(1,3),hxyzloc(2,3),hxyzloc(3,3) 
	end if
c define loc XYZ
        call vectNrmh(atmJKLxyz(1,1),atmxyz,vZ)     

        call vectNrmh(atmJKLxyz(1,1),atmJKLxyz(1,2),vX)
c ortogonalize to vZ
	s = vZ(1)*vX(1)+vZ(2)*vX(2)+vZ(3)*vX(3)
	vX(1) = vX(1)-s*Vz(1)
	vX(2) = vX(2)-s*Vz(2)
	vX(3) = vX(3)-s*Vz(3)
        call vectNrmh(vX,nulV,vX)
	call vectPrh(vZ,vX,vY)

c rotation matrix
	trot(1,1) = vX(1)
	trot(2,1) = vY(1)
	trot(3,1) = vZ(1)
	trot(1,2) = vX(2)
	trot(2,2) = vY(2)
	trot(3,2) = vZ(2)
	trot(1,3) = vX(3)
	trot(2,3) = vY(3)
	trot(3,3) = vZ(3)
c 
	do j = 1,nhadd

	do i =1,3
	hatxyz(i,j) = 0.0d0
	end do
c rotate        
	do i = 1,3
	do k = 1,3
	hatxyz(i,j)=hatxyz(i,j) + trot(k,i)*hxyzloc(k,j)
	end do
	hatxyz(i,j) = hatxyz(i,j) + atmxyz(i)
	end do !i
	end do !j

	return
	end if ! hvar=3

	if(hvar.eq.4 .or. hvar .eq. 40)then
        nhadd=1
	cost = bondh*0.33333333d0    !cos(70.5= 180 -109.5)
	sint = bondh*0.94280904d0 
	if(hvar.eq.40)sint = -sint

c local vectors
        call vectNrmh(atmxyz,atmJKLxyz(1,1),vIJ)
	call vectNrmh(atmJKLxyz(1,1),atmJKLxyz(1,2),vJK)

	s = vIJ(1)*vJK(1) + vIJ(2)*vJK(2) + vIJ(3)*vJK(3)
	vY(1) = vJK(1) - s*vIJ(1)
	vY(2) = vJK(2) - s*vIJ(2)
	vY(3) = vJK(3) - s*vIJ(3)
        call vectNrmh(nulV,vY,vY)

        hatxyz(1,1) = atmxyz(1) - cost*vIJ(1) - sint*vY(1)
        hatxyz(2,1) = atmxyz(2) - cost*vIJ(2) - sint*vY(2)
        hatxyz(3,1) = atmxyz(3) - cost*vIJ(3) - sint*vY(3)

        return
	end if !hvar=4

	if(hvar.eq.5)then
	nhadd = 1
        call vectNrmh(atmxyz,atmJKLxyz(1,1),vIJ)
        call vectNrmh(atmxyz,atmJKLxyz(1,2),vIK)
        vIJ(1)=(vIJ(1)+vIK(1))/2.0
        vIJ(2)=(vIJ(2)+vIK(2))/2.0
        vIJ(3)=(vIJ(3)+vIK(3))/2.0
        call vectNrmh(nulV,vIJ,vH1)
        hatxyz(1,1)=atmxyz(1) - vH1(1)*bondh
        hatxyz(2,1)=atmxyz(2) - vH1(2)*bondh
        hatxyz(3,1)=atmxyz(3) - vH1(3)*bondh

	return
	end if ! hvar=5

	if(hvar.eq.6)then
        nhadd = 2
	cost = bondh*0.5d0    !cos(60.0= 180 -120.0)
	sint = bondh*0.86602540d0 
c local vectors
        call vectNrmh(atmxyz,atmJKLxyz(1,1),vIJ)
	call vectNrmh(atmJKLxyz(1,1),atmJKLxyz(1,2),vJK)

	s = vIJ(1)*vJK(1) + vIJ(2)*vJK(2) + vIJ(3)*vJK(3)
	vY(1) = vJK(1) - s*vIJ(1)
	vY(2) = vJK(2) - s*vIJ(2)
	vY(3) = vJK(3) - s*vIJ(3)
        call vectNrmh(nulV,vY,vY)

        hatxyz(1,1) = atmxyz(1) - cost*vIJ(1) - sint*vY(1)
        hatxyz(2,1) = atmxyz(2) - cost*vIJ(2) - sint*vY(2)
        hatxyz(3,1) = atmxyz(3) - cost*vIJ(3) - sint*vY(3)
        hatxyz(1,2) = atmxyz(1) - cost*vIJ(1) + sint*vY(1)
        hatxyz(2,2) = atmxyz(2) - cost*vIJ(2) + sint*vY(2)
        hatxyz(3,2) = atmxyz(3) - cost*vIJ(3) + sint*vY(3)
 
	return
	end if !hvar=6

	if(hvar.eq.7)then
        nhadd = 2
	cost = bondh*0.5d0    !cos(60)
	sint = bondh*0.86602540d0 
c define loc XYZ
        call vectNrmh(atmJKLxyz(1,1),atmxyz,vZ)     
        call vectNrmh(atmJKLxyz(1,1),atmJKLxyz(1,2),vX)
c ortogonalize to vZ
	s = vZ(1)*vX(1)+vZ(2)*vX(2)+vZ(3)*vX(3)
	vX(1) = vX(1)-s*Vz(1)
	vX(2) = vX(2)-s*Vz(2)
	vX(3) = vX(3)-s*Vz(3)
        call vectNrmh(vX,nulV,vX)
	call vectPrh(vZ,vX,vY)

        hatxyz(1,1) = atmxyz(1) + cost*vZ(1) - sint*vY(1)
        hatxyz(2,1) = atmxyz(2) + cost*vZ(2) - sint*vY(2)
        hatxyz(3,1) = atmxyz(3) + cost*vZ(3) - sint*vY(3)
        hatxyz(1,2) = atmxyz(1) + cost*vZ(1) + sint*vY(1)
        hatxyz(2,2) = atmxyz(2) + cost*vZ(2) + sint*vY(2)
        hatxyz(3,2) = atmxyz(3) + cost*vZ(3) + sint*vY(3)

        return
	end if !hvar=7
c

	return
	end

c calculate unit vector between points
	subroutine vectNrmh(xyzI,xyzJ,vIJ)
        implicit none	
c
        include "charStringSiz.h" 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c
	real*8 xyzI(3),xyzJ(3),vIJ(3)
	real*8 a
        integer kanalpl
c
        kanalpl = kanalRunOut
c
	vIJ(1)=xyzJ(1)-xyzI(1)
	vIJ(2)=xyzJ(2)-xyzI(2)
	vIJ(3)=xyzJ(3)-xyzI(3)
	a = dsqrt(vIJ(1)**2 + vIJ(2)**2 + vIJ(3)**2)
	
	if(a.gt.0.0d0)then
	vIJ(1) = vIJ(1)/a
	vIJ(2) = vIJ(2)/a
	vIJ(3) = vIJ(3)/a
	else
	write(kanalpl,*)'ERROR: vectorNrm: norm=0 !'
        iStatus = iStatus + iMolTopoMX
        write(kanalPStat,*)mStatus,iStatus,messageStatus
     &  ,' ERROR: vectorNrm: norm=0 ! '
	stop
	endif

	return
	end
c calculate vector product
	subroutine vectPrh(v1,v2,v3)
	implicit none
	real*8 v1(3),v2(3),v3(3)

	v3(1)= v1(2)*v2(3) - v1(3)*v2(2)
	v3(2)= -v1(1)*v2(3) + v1(3)*v2(1)
        v3(3)= v1(1)*v2(2) - v1(2)*v2(1)

	return
	end
c
c defines H atoms xyz attached to the WatO: atmnam,resnam,resnumb,chnam,atmxyz(3)
        subroutine hatmWat_xyz(atmnam,resnam,resnumb,chnam,
     &          atmxyz,hnameloc,nhadd,hatxyz)
 
c atmnam,resnam,resnumb,chnam - atom I name, residuename, resnumb,chain name
c atmxyz(3) - coord of atom I
c naddh - number of H atoms added
c hname_tb(3) - names of H atoms
c hatxyz(3,3) - coords of H atoms added to structure
        implicit none
        include "charStringSiz.h" 
        include "output.h"
        include "kanalUse.h"
        include "statusMessg_mDyn.h"
c        
        integer hvar
        character*4 atmnam
        character*4 resnam
        character*1 chnam
        real*8 atmxyz(3)
        integer nhadd,resnumb
        character*4 hnameloc(3)
        real*8 hatxyz(3,3)
c local
        real*8 bondh,sint,cost
        real*8 vH1(3),vH2(3)
        integer i,j,k,kanalpl
        logical CONTROL
        data vH1/1.,0.,0./
        data vH2/0.,1.,0./     ! ex,ey
c
cx        kanalpl = 6
          kanalpl = kanalRunOut 
c
        CONTROL = .false.

        bondh = 0.975d0
        sint = 0.8166415d0*bondh   !sin(109.5/2)
        cost = 0.5771452d0*bondh
c WatH-atoms positions
        nhadd = 2
c
        hatxyz(1,1)=atmxyz(1) + (vH1(1)*cost + vH2(1)*sint)
        hatxyz(2,1)=atmxyz(2) + (vH1(2)*cost + vH2(2)*sint)
        hatxyz(3,1)=atmxyz(3) + (vH1(3)*cost + vH2(3)*sint)
 
        hatxyz(1,2)=atmxyz(1) + (vH1(1)*cost - vH2(1)*sint)
        hatxyz(2,2)=atmxyz(2) + (vH1(2)*cost - vH2(2)*sint)
        hatxyz(3,2)=atmxyz(3) + (vH1(3)*cost - vH2(3)*sint)

        if(CONTROL)then
        write(kanalpl,*)'hatmWat_xyz:O:',(atmxyz(k),k=1,3)
        write(kanalpl,*)'hatmWat_xyz:H1:',(hatxyz(k,1),k=1,3)
        write(kanalpl,*)'hatmWat_xyz:H2:',(hatxyz(k,2),k=1,3) 
        end if 
        return
        end
