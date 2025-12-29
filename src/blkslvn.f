************************************************************************
************************************************************************
*
*  Version:   0.1
*
*  Authors:   Esmond G. Ng, Lawrence Berkeley National Laboratory
*             Barry W. Peyton, Oak Ridge National Laboratory
*
*  Created:   11-15-2000
*  Modified:  01-30-2001
*
************************************************************************
************************************************************************
*****   BLKSLVN ... Block multiple triangular solutions            *****
************************************************************************
************************************************************************
*
*  --------
*  Purpose:
*  --------
*
*     Given the Cholesky factorization of a sparse symmetric matrix,
*     this subroutine performs the triangular solutions of multiple
*     linear system.  It uses output from BLKLDL.
*
*  ----------
*  Arguments:
*  ----------
*
*     NPANEL      (input) integer
*                 Number of panels.
*                 [N = XPANEL(NPANEL+1)-1 is the number of rows and
*                 columns.]
*
*     XPANEL      (input) integer array, dimension NPANEL+1
*                 The panel partition.
*
*     XLINDX      (input) integer array, dimension NPANEL+1
*                 Pointers to column structure of the Cholesky factor.
*
*     LINDX       (input) integer array, dimension XSUPER(NPANEL+1)-1
*                 Row indices of nonzero entries in the Cholesky
*                 factor, stored by columns, in a compressed
*                 representation.
*
*     XLNZ        (input) integer array, dimension N+1
*                 Pointers to nonzero entries in the Cholesky factor.
*
*     LNZ         (input) double precision array, dimension XLNZ(N+1)-1
*                 Numerical values of nonzero entries in the Cholesky
*                 factor, stored by columns.
*
*     DEFBLK      (input) integer
*                 The size of the last supernode, which contains the
*                 known deficiency.
*
*     NDEF        (input) integer
*                 Rank deficiency of the matrix.
*
*     LBNDEF      (input) integer
*                 If DEFBLK is nonzero, LNBDEF is the deficiency
*                 of the last block.
*
*     DEF         (input) integer array, dimension NDEF
*                 It identifies the columns of the matrix that are
*                 linearly dependent.
*
*     IPROW       (input) integer array, dimension DEFBLK
*                 If DEFBLK is nonzero, IPROW contains the row
*                 pivoting sequence for the last block.
*
*     IPCOL       (input) integer array, dimension DEFBLK
*                 If DEFBLK is nonzero, IPCOL contains the column
*                 pivoting sequence for the last block.
*
*     PERM        (input) integer array, dimension N
*                 The permutation vector.
*
*     INVP        (input) integer array, dimension N
*                 The inverse of the permutation vector.
*
*     NRHS        (input) integer
*                 Number of right hand side vectors.
*
*     LRHS        (input) integer
*                 Leading dimension of RHS.
*
*     RHS         (input) double precision array, dimension LRHS by NRHS
*                 The right hand side array.
*
*     LSOLN       (input) integer
*                 Leading dimension of SOLN.
*
*     SOLN        (output) double precision array, dimension LRHS by NRHS
*                 The solution array; can be the same as RHS.
*
*     LTEMP       (input) integer
*                 Leading dimension of TEMP.
*
*     TEMP        (temporary) double precision array, dimension
*                 LTEMP by NRHS
*                 Temporary work array.
*
************************************************************************
*
      SUBROUTINE  BLKSLVN (
     &                        NPANEL, XPANEL, XLINDX, LINDX , XLNZ  ,
     &                        LNZ   , DEFBLK, NDEF  , LBNDEF, DEF   ,
     &                        IPROW , IPCOL , PERM  , INVP  , LRHS  ,
     &                        NRHS  , RHS   , LSOLN , SOLN  , LTEMP ,
     &                        TEMP
     &                    )
*
************************************************************************
*
*     --------------------
*     Scalar Arguments ...
*     --------------------
      INTEGER                 NPANEL
      INTEGER                 LRHS  , NRHS
      INTEGER                 LSOLN , LTEMP
      INTEGER                 DEFBLK, LBNDEF, NDEF
*
*     -------------------
*     Array Arguments ...
*     -------------------
      INTEGER                 DEF(*)
      INTEGER                 INVP(*)       , PERM(*)
      INTEGER                 IPCOL(*)      , IPROW(*)
      INTEGER                 LINDX(*)      , XPANEL(*)
      INTEGER                 XLINDX(*)     , XLNZ(*)
      DOUBLE PRECISION        LNZ(*)
      DOUBLE PRECISION        RHS(LRHS,*)   , SOLN(LSOLN,*)
      DOUBLE PRECISION        TEMP(LTEMP,*)
*
************************************************************************
*
*     --------------
*     Parameters ...
*     --------------
      DOUBLE PRECISION        ONE   , ZERO
      PARAMETER       (       ONE  = 1.0D0,
     &                        ZERO = 0.0D0 )
*
*     --------------------------
*     Local Scalar Variables ...
*     --------------------------
      INTEGER                 FJ    , INFO  , ISUB  , J     , JDPTR,
     &                        JJ    , JLEN  , JLPNT , JPANEL, JXPNT ,
     &                        K     , LJ    , LSTPNL, N     , NJ
      DOUBLE PRECISION        T
*
*     ------------------------
*     External Subroutines ...
*     ------------------------
      EXTERNAL                DGEMM , DGERS , DTRSM
*
************************************************************************
*
      IF  ( NPANEL .LE. 0 )  RETURN
*
*     -----------------------------
*     N is the number of equations.
*     -----------------------------   
      N = XPANEL(NPANEL+1) - 1
*
*     ------------------
*     Initialization ...
*     ------------------
      DO  K = 1, NRHS
          DO  J = 1, N
              T         = RHS(PERM(J),K)
              TEMP(J,K) = T
              SOLN(J,K) = T
          END DO
      END DO
*
*     ---------------------------------------------------------
*     Some key local variables ...
*         FJ    is the first column in the current panel.
*         LJ    is the last column in the current panel.
*         NJ    is the number of columns in the current panel.
*         JLEN  is the length of the column FJ.
*         JXPNT points to the first row index for column FJ.
*         JLPNT points to the first nonzero entry in column FJ.
*     ---------------------------------------------------------
*
      LSTPNL = NPANEL
      IF  ( DEFBLK .NE. 0 )  LSTPNL = NPANEL - 1
*
*     ------------------------
*     Forward substitution ...
*     ------------------------
      DO  JPANEL = 1, LSTPNL
*
          FJ    = XPANEL(JPANEL)
          NJ    = XPANEL(JPANEL+1) - FJ
          JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
          JXPNT = XLINDX(JPANEL)
          JLPNT = XLNZ(FJ)
*
*         ---------------------------------------------------
*         Lower triangular solution using the diagonal block.
*         ---------------------------------------------------
          CALL  DTRSM (
     &            'LEFT', 'LOWER', 'NO TRANSPOSE', 'UNIT',
     &            NJ, NRHS,
     &            ONE,
     &            LNZ(JLPNT), JLEN,
     &            TEMP(FJ,1), LTEMP
     &        )
*
*         ----------------------------------
*         Compute update to right-hand side.
*         ----------------------------------
          CALL  DGEMM (
     &            'NO TRANSPOSE', 'NO TRANSPOSE',
     &            JLEN-NJ, NRHS, NJ,
     &            ONE,
     &            LNZ(JLPNT+NJ), JLEN,
     &            TEMP(FJ,1), LTEMP,
     &            ZERO,
     &            SOLN, LSOLN
     &        )
*
*         ------------------------------------
*         Apply update to the right-hand side.
*         ------------------------------------
          DO  K = 1, NRHS
              JJ = JXPNT + NJ - 1
              DO  J = 1, JLEN-NJ
                  JJ = JJ + 1
                  ISUB = LINDX(JJ)
                  TEMP(ISUB,K) = TEMP(ISUB,K) - SOLN(J,K)
              END DO
          END DO
*
      END DO
*
      IF  ( DEFBLK .NE. 0 )  THEN
*         ------------------------------------------------
*         If the last block is potentially rank deficient,
*         then the corresponding segment of the right-hand
*         side should be handled differently.
*         ------------------------------------------------
          FJ = XPANEL(NPANEL)
          LJ = XPANEL(NPANEL+1) - 1
          NJ = LJ - FJ + 1
          JLEN = XLNZ(FJ+1) - XLNZ(FJ)
          JLPNT = XLNZ(FJ)
          CALL  DGERS (
     &            NJ, NRHS, LNZ(JLPNT), JLEN, LBNDEF, IPROW, IPCOL,
     &            TEMP(FJ,1), LTEMP, INFO
     &        )
      END IF
*
*     -------------------------------------
*     Handle zero diagonal entries, if any.
*     -------------------------------------
      DO  K = 1, NRHS
          DO  JJ = 1, NDEF-LBNDEF
              J = DEF(JJ)
              TEMP(J,K) = ZERO
          END DO
      END DO
*
*     ----------------------------------------------------
*     Scale the right-hand side using the diagonal matrix.
*     ----------------------------------------------------
      DO  JPANEL = 1, LSTPNL
          FJ = XPANEL(JPANEL)
          LJ = XPANEL(JPANEL+1) - 1
          JLEN = XLNZ(FJ+1) - XLNZ(FJ)
          JDPTR = XLNZ(FJ)
          DO  J = FJ, LJ
              DO  K = 1, NRHS
                  TEMP(J,K) = TEMP(J,K)/LNZ(JDPTR)
              END DO
              JDPTR = JDPTR + JLEN + 1
          END DO
      END DO
*
*     -------------------------
*     Backward substitution ...
*     -------------------------
      DO  JPANEL = LSTPNL, 1, -1
*
          FJ    = XPANEL(JPANEL)
          NJ    = XPANEL(JPANEL+1) - FJ
          JLEN  = XLNZ(FJ+1) - XLNZ(FJ)
          JXPNT = XLINDX(JPANEL)
          JLPNT = XLNZ(FJ)
*
*         ---------------------------
*         Gather the right-hand size.
*         ---------------------------
          DO  K = 1, NRHS
              JJ = JXPNT + NJ - 1
              DO  J = 1, JLEN-NJ
                  JJ = JJ + 1
                  ISUB = LINDX(JJ)
                  SOLN(J,K) = TEMP(ISUB,K)
              END DO
          END DO
*
*         ---------------------------
*         Update the right-hand size.
*         ---------------------------
          IF  ( JLEN .GT. NJ )  THEN
              CALL  DGEMM (
     &                'TRANSPOSE', 'NO TRANSPOSE',
     &                NJ, NRHS, JLEN-NJ,
     &                -ONE,
     &                LNZ(JLPNT+NJ), JLEN,
     &                SOLN, LSOLN,
     &                ONE,
     &                TEMP(FJ,1), LTEMP
     &            )
          ENDIF
*
*         ---------------------------------------------------
*         Upper triangular solution using the diagonal block.
*         ---------------------------------------------------
          CALL  DTRSM (
     &            'LEFT', 'LOWER', 'TRANSPOSE', 'UNIT',
     &            NJ, NRHS,
     &            ONE,
     &            LNZ(JLPNT), JLEN,
     &            TEMP(FJ,1), LTEMP
     &        )
*
      END DO
*
*     ------------------------------------
*     Apply the permutation appropriately.
*     ------------------------------------
      DO  K = 1, NRHS
          DO  J = 1, N
              SOLN(J,K) = TEMP(INVP(J),K)
          END DO
      END DO
*
      RETURN
*
*     ---------------
*     End of BLKSLVN.
*     ---------------
      END
