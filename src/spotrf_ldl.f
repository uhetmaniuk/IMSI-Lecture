      SUBROUTINE  SPOTRF_LDL (
     &              UPLO, N, A, LDA, NDEF, IDEF, TOL, WORK, INFO
     &          )
*
*  -- LAPACK-like routine --
*     Esmond G. Ng, Oak Ridge National Laboratory
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N, NDEF
      REAL TOL
*     ..
*     .. Array Arguments ..
      INTEGER            IDEF(*)
      REAL A(LDA,*)
      REAL WORK(*)
*     ..
*
*  Purpose
*  =======
*
*  SPOTRF_LDL computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**T * U,   if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the block version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) SINGLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  NDEF    (output) INTEGER
*          The rank deficiency of the matrix A.  NDEF <= N.
*
*  IDEF    (output) INTEGER array, dimension NDEF
*          Indices of columns for which zero pivots are encountered.
*
*  TOL     (input) REAL 
*          Tolerance for checking if a pivot is zero.
*
*  WORK    (input) REAL array, dimension N.
*          Temporary work array.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*
*  =====================================================================
*
*     .. Parameters ..
      REAL ONE, ZERO
      PARAMETER         ( ONE = 1.0, ZERO = 0.0 )
*     ..
*     .. Local Scalars ..
      LOGICAL           UPPER
      INTEGER           COL, I, II, J, JB, NB, NDEF0, PTR
*     ..
*     .. External Functions ..
      LOGICAL           LSAME
      INTEGER           ILAENV
      EXTERNAL          LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL          SGEMM, SPOTF2_LDL, SSYRK, STRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC         MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'SPOTRF_LDL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'SPOTRF_LDL', UPLO, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         CALL  SPOTF2_LDL ( UPLO, N, A, LDA, NDEF, IDEF, TOL, WORK,
     &                      INFO )
      ELSE
*
*        Use blocked code.
*
         IF( UPPER ) THEN
*
*           Compute the Cholesky factorization A = U'*U.
*
            DO  J = 1, N, NB
*
*              Update and factorize the current diagonal block and test
*              for non-positive-definiteness.
*
               JB = MIN( NB, N-J+1 )
               CALL SSYRK( 'Upper', 'Transpose', JB, J-1, -ONE,
     $                     A( 1, J ), LDA, ONE, A( J, J ), LDA )
               CALL SPOTF2( 'Upper', JB, A( J, J ), LDA, INFO )
               IF( INFO.NE.0 )
     $            GO TO 10
               IF( J+JB.LE.N ) THEN
*
*                 Compute the current block row.
*
                  CALL SGEMM( 'Transpose', 'No transpose', JB, N-J-JB+1,
     $                        J-1, -ONE, A( 1, J ), LDA, A( 1, J+JB ),
     $                        LDA, ONE, A( J, J+JB ), LDA )
                  CALL STRSM( 'Left', 'Upper', 'Transpose', 'Non-unit',
     $                        JB, N-J-JB+1, ONE, A( J, J ), LDA,
     $                        A( J, J+JB ), LDA )
               END IF
            END DO
*
         ELSE
*
*               Compute the Cholesky factorization A = L*L'.
*
                NDEF = 0
                DO  J = 1, N, NB
*
*                   Update and factorize the current diagonal block
*                   and test for non-positive-definiteness.
*
                    JB = MIN( NB, N-J+1 )
                    PTR = 1
                    DO  I = 1, J-1
                        CALL SCOPY( JB, A( J, I ), 1, WORK(PTR), 1 )
                        CALL SSCAL( JB, A( I, I ), WORK(PTR), 1 )
                        PTR = PTR + JB
                    END DO
                    CALL SGEMM ( 'No transpose', 'Transpose',
     &                           JB, JB, J-1,
     &                           -ONE, A(J,1), LDA, WORK, JB,
     &                           ONE, A(J,J), LDA )
                    CALL SPOTF2_LDL ( 'Lower', JB, A(J,J), LDA,
     &                                NDEF0, IDEF(NDEF+1), TOL,
     &                                WORK(PTR), INFO )
                    IF  ( INFO .NE. 0 )  GO TO 10
                    IF ( J+JB .LE. N )  THEN
*
*                       Compute the current block column.
*
                        CALL  SGEMM ( 'No transpose', 'Transpose',
     &                                N-J-JB+1, JB, J-1,
     &                                -ONE, A(J+JB,1), LDA, WORK, JB,
     &                                ONE, A(J+JB,J), LDA )
                        CALL  STRSM ( 'Right', 'Lower',
     &                                'Transpose', 'Unit',
     &                                N-J-JB+1, JB, ONE,
     &                                A(J,J), LDA, A(J+JB,J), LDA )
                        DO  I = J, J+JB-1
                            CALL  SSCAL ( N-J-JB+1, 1/A(I,I),
     &                                    A(J+JB,I), 1 )
                        END DO
                    END IF
*
*                   Handle zero pivots.
*
                    IF  ( NDEF0 .NE. 0 )  THEN
                        DO  I = NDEF+1, NDEF+NDEF0
                            COL = IDEF(I) + (J-1)
                            IDEF(I) = COL
                            DO  II = J+JB, N
                                A(II,COL) = ZERO
                            END DO
                        END DO
                        NDEF = NDEF + NDEF0
                    END IF
*
                END DO
         END IF
      END IF
      RETURN
*
   10 CONTINUE
      write(6,*) '***** INFO = ', INFO
      INFO = INFO + J - 1
      RETURN
*
*     End of SPOTRF_LDL
*
      END
