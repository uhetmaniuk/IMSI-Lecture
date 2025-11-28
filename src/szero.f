************************************************************************
************************************************************************
*
*   Version:    0.5
*   Authors:    01-14-1999 (Michel Lesoinne and Kendall H. Pierson)
*
************************************************************************
************************************************************************
*****   SZEROMAT ... input numerical values                         *****
************************************************************************
************************************************************************
*
*   --------
*   Purpose:
*   --------
*
*   This subroutine creates the INVSUPER array.
*
*   ----------
*   Arguments:
*   ----------
*   
*     NSUPER    (input) integer
*               Number of supernodes.
*
*     XSUPER    (input) integer array, dimension NSUPER+1
*               The supernodal partition.
*
*     INVSUPER  (output) integer array dimension N
*               Keep track of dof to super node map
*               
************************************************************************
*
      SUBROUTINE  SZEROMAT ( NSUPER, XSUPER, INVSUPER )
*
************************************************************************
*
*       --------------------
*       Scalar Arguments ...
*       --------------------
        INTEGER             NSUPER
*
*       -------------------
*       Array Arguments ...
*       -------------------
        INTEGER             XSUPER(*)
        INTEGER             INVSUPER(*)
*
************************************************************************
*
*       --------------------------
*       Local Scalar Variables ...
*       --------------------------
        INTEGER             II, JSUPER
*
************************************************************************
*
*       ----------------------
*       Fill INVSUPER array
*       ----------------------
        DO  JSUPER = 1, NSUPER
          DO II=XSUPER(JSUPER), XSUPER(JSUPER+1) - 1
              INVSUPER(II) = JSUPER
          END DO
        END DO
*
*       --------------
*       End of ZEROMAT.
*       --------------
      END
