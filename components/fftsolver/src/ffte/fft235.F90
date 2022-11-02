
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2014, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!     RADIX-2, 3, 4, 5 AND 8 FFT ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM WRITTEN BY DAISUKE TAKAHASHI

module fft235_mod
    use kernel_mod
    use factor_mod
    use, intrinsic :: iso_fortran_env
    
    implicit none

    integer, parameter, private :: dp = REAL64

contains

    SUBROUTINE FFT235(A,B,W,N,IP)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(*),B(*),W(*)
        complex(kind=dp), dimension(:) :: A, B, W
        integer, dimension(3) :: IP
        integer :: N
        integer :: KP4, KP8, J, L, K, M, KEY

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
                    CALL FFT8(A,B,W(J:),M,L)
                ELSE
                    CALL FFT8(B,A,W(J:),M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT8(A,A,W(J:),M,L)
                ELSE
                    CALL FFT8(B,A,W(J:),M,L)
                END IF
            END IF
            M=M*8
            J=J+L*7
        END DO

        DO K=1,IP(3)
            L=L/5
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT5(A,B,W(J:),M,L)
                ELSE
                    CALL FFT5(B,A,W(J:),M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT5(A,A,W(J:),M,L)
                ELSE
                    CALL FFT5(B,A,W(J:),M,L)
                END IF
            END IF
            M=M*5
            J=J+L*4
        END DO

        DO K=1,KP4
            L=L/4
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT4(A,B,W(J:),M,L)
                ELSE
                    CALL FFT4(B,A,W(J:),M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT4(A,A,W(J:),M,L)
                ELSE
                    CALL FFT4(B,A,W(J:),M,L)
                END IF
            END IF
            M=M*4
            J=J+L*3
        END DO

        DO K=1,IP(2)
            L=L/3
            IF (L >= 2) THEN
                IF (KEY >= 0) THEN
                    CALL FFT3(A,B,W(J:),M,L)
                ELSE
                    CALL FFT3(B,A,W(J:),M,L)
                END IF
                KEY=-KEY
            ELSE
                IF (KEY >= 0) THEN
                    CALL FFT3(A,A,W(J:),M,L)
                ELSE
                    CALL FFT3(B,A,W(J:),M,L)
                END IF
            END IF
            M=M*3
            J=J+L*2
        END DO

        IF (IP(1) == 1) THEN
            IF (KEY >= 0) THEN
                CALL FFT2(A,A,M)
            ELSE
                CALL FFT2(B,A,M)
            END IF
        END IF
    END SUBROUTINE FFT235

    SUBROUTINE FFT3(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(*),B(*),W(*)
        integer :: M,L
        complex(kind=dp), dimension(:) :: A, B, W

        IF (M == 1) THEN
            CALL FFT3A(A,B,W,L)
        ELSE
            CALL FFT3B(A,B,W,M,L)
        END IF
    END SUBROUTINE FFT3

    SUBROUTINE FFT4(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(*),B(*),W(*)
        integer :: M,L
        complex(kind=dp), dimension(:) :: A, B, W

        IF (M == 1) THEN
            CALL FFT4A(A,B,W,L)
        ELSE
            CALL FFT4B(A,B,W,M,L)
        END IF
    END SUBROUTINE FFT4

    SUBROUTINE FFT5(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(*),B(*),W(*)
        integer :: M,L
        complex(kind=dp), dimension(:) :: A, B, W

        IF (M == 1) THEN
            CALL FFT5A(A,B,W,L)
        ELSE
            CALL FFT5B(A,B,W,M,L)
        END IF
    END SUBROUTINE FFT5

    SUBROUTINE FFT8(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(*),B(*),W(*)
        integer :: M,L
        complex(kind=dp), dimension(:) :: A, B, W

        IF (M == 1) THEN
            CALL FFT8A(A,B,W,L)
        ELSE
            CALL FFT8B(A,B,W,M,L)
        END IF
    END SUBROUTINE FFT8

    SUBROUTINE SETTBL(W,N)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 W(*)
        integer, intent(in) :: N
        complex(kind=dp), dimension(:) :: W
        integer :: KP4, KP8, K, L, J
        integer, dimension(3) :: IP

        CALL FACTOR(N,IP)

        IF (IP(1) /= 1) THEN
            KP4=2-MOD(IP(1)+2,3)
            KP8=(IP(1)-KP4)/3
        ELSE
            KP4=0
            KP8=0
        END IF

        J=1
        L=N
        DO K=1,KP8
            L=L/8
            CALL SETTBL0(W(J:),8,L)
            J=J+L*7
        END DO
        DO K=1,IP(3)
            L=L/5
            CALL SETTBL0(W(J:),5,L)
            J=J+L*4
        END DO
        DO K=1,KP4
            L=L/4
            CALL SETTBL0(W(J:),4,L)
            J=J+L*3
        END DO
        DO K=1,IP(2)
            L=L/3
            CALL SETTBL0(W(J:),3,L)
            J=J+L*2
        END DO
    END SUBROUTINE SETTBL

    SUBROUTINE SETTBL0(W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 W(M-1,*)
        integer :: M, L
        complex(kind=dp), dimension(M-1,L) :: W
        integer :: I, J
        real(kind=dp) :: TEMP, PX, PI2

        PI2=8.0D0*DATAN(1.0D0)
        PX=-PI2/(DBLE(M)*DBLE(L))
        DO J=1,L
            DO I=1,M-1
                TEMP=PX*DBLE(I)*DBLE(J-1)
                W(I,J)=DCMPLX(DCOS(TEMP),DSIN(TEMP))
            END DO
        END DO
        RETURN
    END SUBROUTINE SETTBL0

    SUBROUTINE SETTBL2(W,NX,NY)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 W(NX,*)
        integer, intent(in) :: NX, NY
        complex(kind=dp), dimension(NX,NY) :: W
        real(kind=dp) :: TEMP, PX, PI2
        integer :: I, J

        PI2=8.0D0*DATAN(1.0D0)
        PX=-PI2/(DBLE(NX)*DBLE(NY))
        DO J=1,NY
            DO I=1,NX
                TEMP=PX*DBLE(I-1)*DBLE(J-1)
                W(I,J)=DCMPLX(DCOS(TEMP),DSIN(TEMP))
            END DO
        END DO
    END SUBROUTINE SETTBL2

end module fft235_mod
