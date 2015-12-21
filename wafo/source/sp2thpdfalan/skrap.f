     IF (BCVSRT.EQ..FALSE..AND.PRBMIN+SQRT(LTOL).GE.gONE) THEN !
                  D = gZERO
                  E = gZERO
                  DO I = 1, K0-1
                     IF (R(I,JMIN).GT.gZERO) THEN
                        D = D + R(I,JMIN)*AA1(I)
                        E = E + R(I,JMIN)*BB1(I)
                     ELSE
                        D = D + R(I,JMIN)*BB1(I)
                        E = E + R(I,JMIN)*AA1(I)
                     ENDIF
                  END DO
                  D = D/CVDIAG-1.5D0*xCutOff
                  E = E/CVDIAG+1.5D0*xCutOff
                  M = INFI(JMIN)
                  CALL ADJLIMITS(AMIN-D,BMIN-E,M)
                  IF (M.LT.0) THEN
                     CALL mexprintf('R1 ') 
                          !variable is redundnant
                  ! IF ( JMIN.LT.N1 ) THEN
                  ! CALL RCSWAP( JMIN, N1, N1,N, R,INDEX1,Cm, A, B, INFI)
                                ! SWAP conditional standarddeviations
                  !    DO I = 1,K0-1 
                  !       CALL SWAPRE(R(JMIN,I),R(N1,I))
                  !    END DO
                  ! ENDIF

! INFIS = INFIS+1
                  ! N1    = N1-1
                     
                     
                     IF (K .LT. JMIN) THEN
                           
                        CALL RCSWAP( K, JMIN, N1,N,R,INDEX1,Cm, A, B,
     $                       INFI)
                     ! SWAP conditional standarddeviations
                        DO J=1,K0-1
                           CALL SWAPRE(R(K,J),R(JMIN,J))
                        END DO
                     ENDIF
                  
                     CDI(K) = R(K0-1,K)
                     A(K)   = A(K)/CDI(K)
                     B(K)   = B(K)/CDI(K)                  

                     IF ((CDI(K) .LT. ZERO).AND. INFI(K).GE. 0) THEN
                        CALL SWAPRE(A(K),B(K))
                        IF (INFI(K).NE. 2) INFI(K) = 1-INFI(K)
                     END IF

                     DO J = 1, K0-1
                        R(J,K) = R(J,K)/CDI(K) ! conditional covariances
                        R(K,J) = R(K,J)/ABS(CDI(K)) ! conditional standard dev.s used in regression eq.
                     END DO
                     DO J = K0,K
                        R(J,K) = ZERO
                        R(K,J) = ZERO
                     END DO
                     Nullity = Nullity + 1
                     K = K + 1
                     GOTO 100 
                  END IF      
                ENDIF


%%%%%%%%%%%%%%%%%%%%%%%Part 2

   IF (.FALSE.) THEN  !BCVSRT.EQ.
                        D = gZERO
                        E = gZERO
                        DO J = 1, K0-1
                           IF (R(J,I).GT.gZERO) THEN
                              D = D + R(J,I)*AA1(J)
                              E = E + R(J,I)*BB1(J)
                           ELSE
                              D = D + R(J,I)*BB1(J)
                              E = E + R(J,I)*AA1(J)
                           ENDIF
                        END DO
                        D = (A(I)-D)/ABS(R(K0,I))
                        E = (B(I)-E)/ABS(R(K0,I))
                        M = INFI(I)
                        CALL ADJLIMITS(D,E,M)

                        IF (M.LT.0) THEN
                           CALL mexprintf('R2')
                                !variable is redundnant
                           IF ( I.LT.N1 ) THEN
                              CALL RCSWAP( I, N1, N1,N, R,INDEX1,Cm, A,
     &                             B, INFI)
                                ! SWAP conditional standarddeviations
                              DO J = 1,K0
                                 CALL SWAPRE(R(I,J),R(N1,J))
                              END DO
                           ENDIF
                           INFIS = INFIS+1
                           N1    = N1-1
                           GOTO 75 
                        END IF 

ENDIF

%%%%%%%%%%part 3 

    IF (.FALSE..AND.INDEX1(I).LE.Nt) THEN ! BCVSRT.EQ.
                  D = gZERO
                  E = gZERO
                  DO J = 1, K0-2
                     IF (R(J,I).GT.gZERO) THEN
                        D = D + R(J,I)*AA1(J)
                        E = E + R(J,I)*BB1(J)
                     ELSE
                        D = D + R(J,I)*BB1(J)
                        E = E + R(J,I)*AA1(J)
                     ENDIF
                  END DO
                  D = (A(I)-D)+xCutOff
                  E = (B(I)-E)-xCutOff
                  M = INFI(I)
                  CALL ADJLIMITS(D,E,M)
                  
                  IF (M.LT.0) THEN
                     CALL mexprintf('R3')
                                !variable is redundnant
                     IF ( I.LT.N1 ) THEN
                        CALL RCSWAP( I, N1, N1,N, R,INDEX1,Cm, A,
     &                       B, INFI)
                                ! SWAP conditional standarddeviations
                        DO J = 1,K0-1

                           CALL SWAPRE(R(I,J),R(N1,J))
                        END DO
                        I = I + 1
                     ENDIF
                     INFIS = INFIS+1
                     N1    = N1-1
                     !GOTO 75 
                  END IF   
                  ENDIF
