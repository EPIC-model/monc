
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE

!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp


!     RADIX-2, 3, 4, 5 AND 8 MULTIPLE FFT ROUTINE

!     FORTRAN77 SOURCE PROGRAM

!     WRITTEN BY DAISUKE TAKAHASHI
module mfft235_mod
    use factor_mod
    use fft235_mod
    use kernel_mod
    use, intrinsic :: iso_fortran_env

    implicit none
    ! The maximum supported 2-D transform length
    integer, parameter :: NDA2 = 65536
    integer, parameter :: dp = REAL64

contains
    SUBROUTINE MFFT235A(A,B,W,NS,N,IP)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(*),B(*),W(*)
        complex(kind=dp), dimension(*) :: A, B
        complex(kind=dp), dimension(NDA2) :: W
        integer, dimension(3) :: IP
        integer, intent(in) :: N, NS
        integer :: KP4, KP8, M, L, KEY, K, J

        IF (IP(1) /= 1) THEN
            KP4=2-MOD(IP(1)+2,3)
            KP8=(IP(1)-KP4)/3
        ELSE
            KP4=0
            KP8=0
        END IF

        KEY=1
        J=1
        L=N
        M=1
        DO K=1,KP8
            L=L/8
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT8B(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT8B(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT8B(A,A,W(J:),NS*M,L)
                ELSE
                    CALL FFT8B(B,A,W(J:),NS*M,L)
                END IF
            END IF
            M=M*8
            J=J+L*7
        END DO
        DO K=1,IP(3)
            L=L/5
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT5B(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT5B(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT5B(A,A,W(J),NS*M,L)
                ELSE
                    CALL FFT5B(B,A,W(J),NS*M,L)
                END IF
            END IF
            M=M*5
            J=J+L*4
        END DO
        DO K=1,KP4
            L=L/4
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT4B(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT4B(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT4B(A,A,W(J:),NS*M,L)
                ELSE
                    CALL FFT4B(B,A,W(J:),NS*M,L)
                END IF
            END IF
            M=M*4
            J=J+L*3
        END DO
        DO K=1,IP(2)
            L=L/3
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT3B(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT3B(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT3B(A,A,W(J:),NS*M,L)
                ELSE
                    CALL FFT3B(B,A,W(J:),NS*M,L)
                END IF
            END IF
            M=M*3
            J=J+L*2
        END DO
        IF (IP(1) == 1) THEN
            IF (KEY >= 0) THEN
                CALL FFT2(A,A,NS*M)
            ELSE
                CALL FFT2(B,A,NS*M)
            END IF
        END IF
    END SUBROUTINE MFFT235A

    SUBROUTINE MFFT235B(A,B,W,NS,N,IP)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(*),B(*),W(*)
        complex(kind=dp), dimension(:) :: A, B, W
        integer, dimension(3) :: IP
        integer, intent(in) :: N, NS
        integer :: KP4, KP8, M, L, KEY, K, J

        IF (IP(1) /= 1) THEN
            KP4=2-MOD(IP(1)+2,3)
            KP8=(IP(1)-KP4)/3
        ELSE
            KP4=0
            KP8=0
        END IF

        KEY=1
        J=1
        L=N
        M=1
        DO K=1,KP8
            L=L/8
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT8(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT8(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT8(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT8(B,B,W(J:),NS*M,L)
                END IF
            END IF
            M=M*8
            J=J+L*7
        END DO
        DO K=1,IP(3)
            L=L/5
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT5(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT5(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT5(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT5(B,B,W(J:),NS*M,L)
                END IF
            END IF
            M=M*5
            J=J+L*4
        END DO
        DO K=1,KP4
            L=L/4
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT4(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT4(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT4(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT4(B,B,W(J:),NS*M,L)
                END IF
            END IF
            M=M*4
            J=J+L*3
        END DO
        DO K=1,IP(2)
            L=L/3
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT3(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT3(B,A,W(J:),NS*M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT3(A,B,W(J:),NS*M,L)
                ELSE
                    CALL FFT3(B,B,W(J:),NS*M,L)
                END IF
            END IF
            M=M*3
            J=J+L*2
        END DO
        IF (IP(1) == 1) THEN
            IF (KEY >= 0) THEN
                CALL FFT2(A,B,NS*M)
            ELSE
                CALL FFT2(B,B,NS*M)
            END IF
        END IF
    END SUBROUTINE MFFT235B

    SUBROUTINE ZTRANS(A,B,NX,NY)
        !IMPLICIT REAL*8 (A-H,O-Z)
        !COMPLEX*16 A(*),B(*)
        integer, intent(in) :: NX, NY
        complex(kind=dp), dimension(NX*NY) :: A, B
        !DIMENSION LNX(3),LNY(3)
        integer, dimension(3) :: LNX, LNY
        integer :: I

        CALL FACTOR(NX,LNX)
        CALL FACTOR(NY,LNY)

        IF (NX == 1 .OR. NY == 1) THEN
            DO I=1,NX*NY
                B(I)=A(I)
            END DO
            RETURN
        END IF

        IF (LNX(1)+LNY(1) <= 1) THEN
            CALL ZTRANSA(A,B,NX,NY)
        ELSE
            CALL ZTRANSB(A,B,NX,NY)
        END IF
    END SUBROUTINE ZTRANS

    SUBROUTINE ZTRANSA(A,B,NX,NY)
        !IMPLICIT REAL*8 (A-H,O-Z)
        !COMPLEX*16 A(NX,*),B(NY,*)
        integer, intent(in) :: NX, NY
        complex(kind=dp), dimension(NX,NY) :: A
        complex(kind=dp), dimension(NY,NX) :: B
        integer :: I, J

        DO I=1,NX
            DO J=1,NY
                B(J,I)=A(I,J)
            END DO
        END DO
    END SUBROUTINE ZTRANSA

    SUBROUTINE ZTRANSB(A,B,NX,NY)
        !IMPLICIT REAL*8 (A-H,O-Z)
        !COMPLEX*16 A(NX,*),B(NY,*)
        integer, intent(in) :: NX, NY
        complex(kind=dp), dimension(NX,NY) :: A
        complex(kind=dp), dimension(NY,NX) :: B
        integer :: I, J

        IF (NY >= NX) THEN
            DO I=0,NX-1
                DO J=1,NX-I
                    B(J,I+J)=A(I+J,J)
                END DO
            END DO
            DO I=1,NY-NX
                DO J=1,NX
                    B(I+J,J)=A(J,I+J)
                END DO
            END DO
            DO I=NY-NX+1,NY-1
                DO J=1,NY-I
                    B(I+J,J)=A(J,I+J)
                END DO
            END DO
        ELSE
            DO I=0,NY-1
                DO J=1,NY-I
                    B(I+J,J)=A(J,I+J)
                END DO
            END DO
            DO I=1,NX-NY
                DO J=1,NY
                    B(J,I+J)=A(I+J,J)
                END DO
            END DO
            DO I=NX-NY+1,NX-1
                DO J=1,NX-I
                    B(J,I+J)=A(I+J,J)
                END DO
            END DO
        END IF
    END SUBROUTINE ZTRANSB

    SUBROUTINE ZTRANS2(A,B,NX,NY,NZ)
        !IMPLICIT REAL*8 (A-H,O-Z)
        !COMPLEX*16 A(NX,NY,*),B(NY,NX,*)
        integer, intent(in) :: NX, NY, NZ
        complex(kind=dp), dimension(NX,NY,NZ) :: A
        complex(kind=dp), dimension(NY,NX,NZ) :: B
        integer :: K

        DO K=1,NZ
            CALL ZTRANS(A(1,1,K),B(1,1,K),NX,NY)
        END DO
    END SUBROUTINE ZTRANS2

    SUBROUTINE MZTRANS(A,B,NS,NX,NY)
        !IMPLICIT REAL*8 (A-H,O-Z)
        !COMPLEX*16 A(NS,NX,*),B(NS,NY,*)
        integer, intent(in) :: NX, NY, NS
        complex(kind=dp), dimension(NS,NX,NY) :: A
        complex(kind=dp), dimension(NS,NY,NX) :: B
        integer :: I, J, K

        IF (NS == 1) THEN
            CALL ZTRANS(A(1,1,1),B(1,1,1),NX,NY)
        ELSE
            DO I=1,NX
                DO J=1,NY
                    DO K=1,NS
                        B(K,J,I)=A(K,I,J)
                    END DO
                END DO
            END DO
        END IF
    END SUBROUTINE MZTRANS
end module mfft235_mod