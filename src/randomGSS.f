C
C     gauss (AM,SD,NUM,V,IG)
C
C     gauss will supply NUM normally distributed random numbers with
C     a given mean and standard deviation. It uses a polar Box-Muller
C     method to compute normal random numbers.
C     (Black, S.C., Computers in Physics, 1989)
C
C     AM    = the desired mean of the normal distribution
C     SD    = the desired standard deviation of the normal distribution
C     NUM   = the desired number of elements in V to produce
C     V(NUM)= delivered with NUM normal random variables
C     IG    = random number generator seed. This variable is modified
C             by subroutine RANDOM.
C
C     gauss uses a random number generator L<RANDOM>.
C
C     This can be replaced by subroutine RANDU(IG,IY,Y) on IBM machines,
C     function RANDU(IG) on VAX machines, or intrinsic function RAN(IG)
C     on CDC machines.

      subroutine gauss(AM,SD,NUM,V,IG)
C args
      INTEGER NUM,IG
      REAL AM,SD,V(NUM)
C local 
      INTEGER I
      LOGICAL LNEW
      REAL R,STORE,SQLN,W1,W2
      SAVE LNEW,STORE
      DATA LNEW/.FALSE./

C begin
      DO 20 I=1,NUM

         IF (LNEW) THEN
            V(I) = AM + SD * STORE
            LNEW = .FALSE.
         ELSE
C repeat until (r<1.0 and r<>0.0)
 10         CONTINUE

               CALL RANDOM(W1,IG)
               CALL RANDOM(W2,IG)

               W1 = 2.0*W1 - 1.0
               W2 = 2.0*W2 - 1.0

               R = W1**2 + W2**2

            IF ((R.GE.1.0) .OR. (R.EQ .0.0)) GO TO 10
C end of repeat until loop

            SQLN = SQRT(-2.0E0*LOG(R)/R)
            V(I) = AM + SD* (W1*SQLN)
            STORE = W2*SQLN
            LNEW = .TRUE.
         ENDIF
 20   CONTINUE
C end gauss
      return
      END
C
C     subroutine RANDOM (RAND,IG)
C
C     RANDOM GENERATES A RANDOM NUMBER RAND, USING A LINEAR
C     CONGRUENTIAL METHOD. THE RECURSION FORMULA
C
C         IRAND = MOD(IRAND * B + 1, A)
C
C     IS USED WITH  B = 31415821  AND  A = 100000000. THE LAST DIGIT
C     FROM THE RANDOM INTEGER IRAND IS CHOPPED OF, AND THE NUMBER
C     IS SCALED TO A REAL VALUE RAND BETWEEN 0 AND 1, INCLUDING 0 BUT
C     EXCLUDING 1.
C
C     RAND = DELIVERED WITH RANDOM NUMBER BETWEEN 0 AND 1
C     IG = RANDOM NUMBER GENERATOR SEED, IS DELIVERED WITH A RANDOM
C          INTEGER.
C
C     THIS VERSION A  NUMBER OF RANDOM NUMBERS LOCALLY

      subroutine RANDOM(RAND,IG)
C args
      REAL RAND
      INTEGER IG
C local params
      INTEGER M,M1,MULT
      PARAMETER (M=100000000,M1=10000,MULT=31415821)
C     how many random numbers to store locally
      INTEGER MAXLOC
      PARAMETER (MAXLOC = 20)
C local vars
      LOGICAL LNEW
      INTEGER IRAND
      SAVE LNEW,IRAND
C
      INTEGER IRANDH,IRANDL,MULTH,MULTL,IUSED,I
      INTEGER ISAVE(MAXLOC)
      REAL R,RSAVE(MAXLOC)
      SAVE IUSED,ISAVE,RSAVE
C data
      DATA IRAND/0/
      DATA LNEW/.TRUE./
C begin
      IF (LNEW) THEN
         LNEW = .FALSE.
         IRAND = MOD(IABS(IG),M)
         IUSED = MAXLOC
      ENDIF

      IF (IUSED .EQ. MAXLOC) THEN
C     have to generate new numbers which we store
C     into ISAVE and RSAVE
C
         DO 10 I=1,MAXLOC
C*****MULTIPLY IRAND BY MULT, BUT TAKE INTO ACCOUNT THAT OVERFLOW MUST
C*****BE DISCARDED, AND DO NOT GENERATE AN ERROR.
C
            IRANDH = IRAND/M1
            IRANDL = MOD(IRAND,M1)
            MULTH = MULT/M1
            MULTL = MOD(MULT,M1)
C
            IRAND = MOD(IRANDH*MULTL+IRANDL*MULTH,M1)*M1 + IRANDL*MULTL
            IRAND = MOD(IRAND+1,M)
C
C*****CONVERT IRAND TO A REAL RANDOM NUMBER BETWEEN 0 AND 1.
C
            R = REAL(IRAND/10)*10/REAL(M)
            IF ((R .LE. 0.0) .OR. (R .GT. 1.0)) R = 0.0
            RSAVE(I) = R
            ISAVE(I) = IRAND
 10      CONTINUE
         IUSED = 0
      ENDIF
C     now give out a number
      IUSED = IUSED + 1
      RAND = RSAVE(IUSED)
      IG   = ISAVE(IUSED)
C end random
      END
