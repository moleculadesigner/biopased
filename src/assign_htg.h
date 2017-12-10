c---------------------------------------------------
c parameters for assign_htg subroutine
c---------------------------------------------------
        integer ntbmax,natmx
c
	parameter (ntbmax = 1007)
	parameter (natmx = 50000)
c---------------------------------------------------
c
c     atom names for htg table
       character*4 atnam_tb(ntbmax)		

c residue Terminal tag
      character*4 term_tb(ntbmax)       

c     residue names for htg table
       character*4 rnam_tb(ntbmax)		
c
c     residue numbers for htg table
       character*4 rnumb_tb(ntbmax)
c
c     chain name
       character*1 chnam_tb(ntbmax)
c
c     htag id
       integer htag_tb(ntbmax)
c
c     H atoms names to be added
       character*4 hname_tb(3,ntbmax)
c
c     names of iatoms J K L to meke geometry for H position
       character*4 atmJKL(3,ntbmax)
c 
c     atom ---> record link table 
       integer atmrecn(natmx)

       common/htg_tbl1/atnam_tb,term_tb,rnam_tb,rnumb_tb,chnam_tb,
     &               hname_tb, atmJKL
       common/htg_tbl2/htag_tb
       common/htg_atmlink/atmrecn
c end
