C* ** * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C* getDijConstLaZmodel()                                                      *
C* Computes the dielectric constant  for atomPair i,j                         *
C* T.Lazaridis model                                                          *
C* JCC.23: 1090-1099,2002
C* Yury Vorobjev 2005                                                         *
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c 
        subroutine getDijConst05(parDrLZmodel,
     &                atDistToSAS,atSASexp,ia,ja,dij1,DCij,dDCij)
c
c IN: ia,ja - atom N
c     parDrLZmodel(*) = D0,Ds - molCore, bulkSolvent D constants,+
c                       A1,A2,A3,A4,A5 constants
c     atDistToSAS(ia),atSASexp(ia) - 
c     dij1 - at-at distant
c OUT:
c      DCij - dielectric constant for pair i,j
c      dDCij - dD(r)/dr
c
        include 'xyzPDBsize.h'
cx	implicit none
        real parDrLZmodel(*)
        real atDistToSAS(*),atSASexp(*)
        real D0,Ds
        integer ia,ja
        real dij1	
        real DCij,dDCij
        include 'controlCall.h'
        include 'vdw12Par.h'
        include 'xyzPDBinfo.h'
c local
        real dsij,esij,fsij
        real r,rn,rn1,a,aA,dMin
        integer nP,iP,ip7 
        integer t1,t2,t12
        real rm12,rm12MIN_OPT
        logical sasDielModVDWrad_OPT
        integer kanalp
        logical CONTROL
c
        np = 1
        ncallSubTot(np)= ncallSubTot(np) + 1
c
        kanalp = 6
        CONTROL = .false.          
c
        sasDielModVDWrad_OPT=.true.   
        rm12MIN_OPT=2.2  !min Rvdwij
c defineAtomcVDWrad
        if(sasDielModVDWrad_OPT)then
        t1 = atomVDWtype(ia)
        t2 = atomVDWtype(ja)
c get pointer to the AB table
        call vdw12TablePos(nVDWtype,t1,t2,t12)
        rm12 = atomVDW12ab(4*t12-1)
        if(rm12MIN_OPT .gt. rm12)rm12 = rm12MIN_OPT
c redefine parDrLZmodel(4)
        parDrLZmodel(4)=exp(rm12/parDrLZmodel(3))
        end if ! sasDielModVDWrad_OPT
c
        dMin = 0.5   ! A
        DCij = parDrLZmodel(1)
        dDCij = 0.0
        if(dij1 .lt. dMin)return
c
c estimate dsij  
        dsij = min(atDistToSAS(ia),atDistToSAS(ja)) 
        esij = max(atSASexp(ia),atSASexp(ja))
        fsij = 1.0 - esij
c
        a = parDrLZmodel(4)+parDrLZmodel(5)*dsij+parDrLZmodel(6)*fsij
        aA = parDrLZmodel(3)*alog(a)
        r = dij1/aA
        rn = r**parDrLZmodel(7)
        rn1 = rn/(1.0+rn)
        DCij = parDrLZmodel(1)+parDrLZmodel(2)*rn1
        dDCij = parDrLZmodel(2)*parDrLZmodel(7)/dij1
        dDCij = dDCij*(rn1**2)/rn
c
         if(CONTROL)then
         ip = ncallSubTot(np)/nPrintNcall(np)
         if(ip*nPrintNcall(np) .eq. ncallSubTot(np))then
c
         write(kanalp,'(a22,2i5,2f7.1,f7.3,f7.1,1x,2f5.1)')
     &  'DijConst05: dij,Dij,dDij:',
     &   ia,ja,dij1,DCij,(1.0/DCij),dDCij,dsij,esij 
         end if ! print
         end if !C
c
        return
        end
