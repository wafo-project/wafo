      MODULE GAUSSMOD
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: BVNMVN, BVU
      
      INTERFACE BVNMVN
      MODULE PROCEDURE BVNMVN
      END INTERFACE
      
      INTERFACE BVU
      MODULE PROCEDURE BVU
      END INTERFACE

      CONTAINS
      FUNCTION BVNMVN( LOWER, UPPER, CORREL ) RESULT (value)
      USE GLOBALDATA, ONLY : xCutOff
*
*     A function for computing bivariate normal probabilities.
*
*  Parameters
*
*     LOWER  REAL, array of lower integration limits.
*     UPPER  REAL, array of upper integration limits.
*     INFIN  INTEGER, array of integration limits flags:
*            if INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
*            if INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
*            if INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
*     CORREL REAL, correlation coefficient.
*
      DOUBLE PRECISION, DIMENSION(2) :: LOWER, UPPER
      DOUBLE PRECISION               :: CORREL
      DOUBLE PRECISION :: value
      INTEGER, DIMENSION(2) :: INFIN
      INFIN=2      
      IF (LOWER(1).LE.-xCutOff) INFIN(1)=0
      IF (LOWER(2).LE.-xCutOff) INFIN(2)=0
      IF (UPPER(1).GE. xCutOff) INFIN(1)=1
      IF (UPPER(2).GE. xCutOff) INFIN(2)=1

      IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 2 ) THEN
         VALUE =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
     +           + BVU ( UPPER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 1 ) THEN
         VALUE =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( UPPER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 2 ) THEN
         VALUE =  BVU ( LOWER(1), LOWER(2), CORREL )
     +           - BVU ( LOWER(1), UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 2  .AND. INFIN(2) .EQ. 0 ) THEN
         VALUE =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -LOWER(1), -UPPER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 2 ) THEN
         VALUE =  BVU ( -UPPER(1), -UPPER(2), CORREL )
     +           - BVU ( -UPPER(1), -LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 0 ) THEN
         VALUE =  BVU ( LOWER(1), -UPPER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 1 ) THEN
         VALUE =  BVU ( -UPPER(1), LOWER(2), -CORREL )
      ELSE IF ( INFIN(1) .EQ. 1  .AND. INFIN(2) .EQ. 1 ) THEN
         VALUE =  BVU ( LOWER(1), LOWER(2), CORREL )
      ELSE IF ( INFIN(1) .EQ. 0  .AND. INFIN(2) .EQ. 0 ) THEN
         VALUE =  BVU ( -UPPER(1), -UPPER(2), CORREL )
      END IF
      RETURN
      END FUNCTION BVNMVN
      FUNCTION BVU( SH, SK, R ) RESULT (BVN)
*
*     A function for computing bivariate normal probabilities.
*
*       Yihong Ge
*       Department of Computer Science and Electrical Engineering
*       Washington State University
*       Pullman, WA 99164-2752
*     and
*       Alan Genz
*       Department of Mathematics
*       Washington State University
*       Pullman, WA 99164-3113
*       Email : alangenz@wsu.edu
*
* BVN - calculate the probability that X is larger than SH and Y is
*       larger than SK.
*
* Parameters
*
*   SH  REAL, integration limit
*   SK  REAL, integration limit
*   R   REAL, correlation coefficient
*   LG  INTEGER, number of Gauss Rule Points and Weights
*
      DOUBLE PRECISION BVN, SH, SK, R, ZERO, TWOPI 
      INTEGER I, LG, NG
      PARAMETER ( ZERO = 0, TWOPI = 6.2831 85307 179586 ) 
      DOUBLE PRECISION X(10,3), W(10,3), AS, A, B, C, D, RS, XS
      DOUBLE PRECISION SN, ASR, H, K, BS, HS, HK
*     Gauss Legendre Points and Weights, N =  6
      DATA ( W(I,1), X(I,1), I = 1,3) /
     *  0.1713244923791705D+00,-0.9324695142031522D+00,
     *  0.3607615730481384D+00,-0.6612093864662647D+00,
     *  0.4679139345726904D+00,-0.2386191860831970D+00/
*     Gauss Legendre Points and Weights, N = 12
      DATA ( W(I,2), X(I,2), I = 1,6) /
     *  0.4717533638651177D-01,-0.9815606342467191D+00,
     *  0.1069393259953183D+00,-0.9041172563704750D+00,
     *  0.1600783285433464D+00,-0.7699026741943050D+00,
     *  0.2031674267230659D+00,-0.5873179542866171D+00,
     *  0.2334925365383547D+00,-0.3678314989981802D+00,
     *  0.2491470458134029D+00,-0.1252334085114692D+00/
*     Gauss Legendre Points and Weights, N = 20
      DATA ( W(I,3), X(I,3), I = 1,10) /
     *  0.1761400713915212D-01,-0.9931285991850949D+00,
     *  0.4060142980038694D-01,-0.9639719272779138D+00,
     *  0.6267204833410906D-01,-0.9122344282513259D+00,
     *  0.8327674157670475D-01,-0.8391169718222188D+00,
     *  0.1019301198172404D+00,-0.7463319064601508D+00,
     *  0.1181945319615184D+00,-0.6360536807265150D+00,
     *  0.1316886384491766D+00,-0.5108670019508271D+00,
     *  0.1420961093183821D+00,-0.3737060887154196D+00,
     *  0.1491729864726037D+00,-0.2277858511416451D+00,
     *  0.1527533871307259D+00,-0.7652652113349733D-01/
      SAVE X, W
      IF ( ABS(R) .LT. 0.3d0 ) THEN
         NG = 1
         LG = 3
      ELSE IF ( ABS(R) .LT. 0.75d0 ) THEN
         NG = 2
         LG = 6
      ELSE 
         NG = 3
         LG = 10
      ENDIF
      H = SH
      K = SK 
      HK = H*K
      BVN = 0.d0
      IF ( ABS(R) .LT. 0.925d0 ) THEN
         HS = ( H*H + K*K )*0.5d0
         ASR = ASIN(R)
         DO I = 1, LG
            SN = SIN(ASR*( X(I,NG)+1.d0 )*0.5d0)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1.d0 - SN*SN ) )
            SN = SIN(ASR*(-X(I,NG)+1 )*0.5d0)
            BVN = BVN + W(I,NG)*EXP( ( SN*HK - HS )/( 1.d0 - SN*SN ) )
         END DO
         BVN = BVN*ASR/(2.d0*TWOPI) + FI(-H)*FI(-K) 
      ELSE
         IF ( R .LT. 0.d0 ) THEN
            K = -K
            HK = -HK
         ENDIF
         IF ( ABS(R) .LT. 1.d0 ) THEN
            AS = ( 1.d0 - R )*( 1.d0 + R )
            A = SQRT(AS)
            BS = ( H - K )**2
            C = ( 4.d0 - HK )/8.d0 
            D = ( 12.d0 - HK )/16.d0
            BVN = A*EXP( -(BS/AS + HK)*0.5d0 )
     &           *( 1.d0 - C*(BS - AS)*(1.d0 - D*BS/5.d0)/3.d0 
     &           + C*D*AS*AS/5.d0 )
            IF ( HK .GT. -160.d0 ) THEN
               B = SQRT(BS)
               BVN = BVN - EXP(-HK*0.5d0)*SQRT(TWOPI)*FI(-B/A)*B
     &                    *( 1.d0 - C*BS*( 1.d0 - D*BS/5.d0 )/3.d0 ) 
            ENDIF
            A = A*0.5d0
            DO I = 1, LG
               XS = ( A*(X(I,NG)+1) )**2
               RS = SQRT( 1.d0 - XS )
               BVN = BVN + A*W(I,NG)*
     &              ( EXP( -BS/(2.d0*XS) - HK/(1+RS) )/RS 
     &              - EXP( -(BS/XS+HK)/2.d0 )*
     &              ( 1.d0 + C*XS*( 1.d0 + D*XS ) ) )
               XS = AS*(-X(I,NG)+1.d0)**2/4.d0
               RS = SQRT( 1.d0 - XS )
               BVN = BVN + A*W(I,NG)*EXP( -(BS/XS + HK)/2.d0 )
     &                    *( EXP( -HK*(1.d0-RS)/(2.d0*(1.d0+RS)) )/RS 
     &                       - ( 1.d0 + C*XS*( 1.d0 + D*XS ) ) )
            END DO
            BVN = -BVN/TWOPI
         ENDIF
         IF ( R .GT. 0.d0 ) BVN =  BVN + FI( -MAX( H, K ) )
         IF ( R .LT. 0.d0 ) BVN = -BVN + MAX( ZERO, FI(-H)-FI(-K) )     
      ENDIF
      RETURN
      END  FUNCTION BVU     
      END MODULE GAUSSMOD
