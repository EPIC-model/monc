
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE
!
!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!     1-D COMPLEX FFT ROUTINE (FOR VECTOR MACHINES)
!
!     FORTRAN77 SOURCE PROGRAM WRITTEN BY DAISUKE TAKAHASHI
!
!     CALL ZFFT1D(A,N,IOPT,B)
!
!     A(N) IS COMPLEX INPUT/OUTPUT VECTOR (COMPLEX*16)
!     B(N*2) IS WORK/COEFFICIENT VECTOR (COMPLEX*16)
!     N IS THE LENGTH OF THE TRANSFORMS (INTEGER*4)
!       -----------------------------------
!         N = (2**IP) * (3**IQ) * (5**IR)
!       -----------------------------------
!     IOPT = 0 FOR INITIALIZING THE COEFFICIENTS (INTEGER*4)
!          = -1 FOR FORWARD TRANSFORM
!          = +1 FOR INVERSE TRANSFORM

module vzfft1d_mod
    use factor_mod
    use fft235_mod, only : SETTBL, SETTBL2
    use mfft235_mod, only : MFFT235A, MFFT235B
    use, intrinsic :: iso_fortran_env


    implicit none
    integer, parameter, private :: dp = REAL64
    ! The maximum supported 2-D transform length
    integer, parameter :: NDA2 = 65536

contains

    SUBROUTINE ZFFT1D(A,N,IOPT,B)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! INCLUDE 'param.h'
        ! COMPLEX*16 A(*),B(*)
        integer, intent(in) :: N, IOPT
        complex(kind=dp), dimension(:) :: A, B
        ! COMPLEX*16 WX(NDA2),WY(NDA2)
        complex(kind=dp), dimension(NDA2) :: WX, WY
        integer, dimension(3) :: IP, LNX, LNY
        integer :: NX, NY, I
        real(kind=dp) :: DN
        ! SAVE WX,WY

        CALL FACTOR(N,IP)

        IF (IOPT == 1) THEN
            DO I=1,N
                A(I)=DCONJG(A(I))
            END DO
        END IF

        CALL GETNXNY(N,NX,NY)
        CALL FACTOR(NX,LNX)
        CALL FACTOR(NY,LNY)

        IF (IOPT == 0) THEN
            CALL SETTBL(WX,NX)
            CALL SETTBL(WY,NY)
            CALL SETTBL2(B(N+1:),NX,NY)
            RETURN
        END IF

        CALL MFFT235A(A,B,WY,NX,NY,LNY)
        CALL ZTRANSMUL(A,B,B(N+1:),NX,NY)
        CALL MFFT235B(B,A,WX,NY,NX,LNX)

        IF (IOPT == 1) THEN
            DN=1.0D0/DBLE(N)
            DO I=1,N
                A(I)=DCONJG(A(I))*DN
            END DO
        END IF
    END SUBROUTINE ZFFT1D

    SUBROUTINE ZTRANSMUL(A,B,W,NX,NY)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        integer, intent(in) :: NX, NY
        ! DIMENSION LNX(3),LNY(3)
        integer, dimension(3) :: LNX, LNY
        ! COMPLEX*16 A(*),B(*),W(*)
        complex(kind=dp), dimension(:) :: A, B, W
        integer :: I

        CALL FACTOR(NX,LNX)
        CALL FACTOR(NY,LNY)

        IF (NX == 1 .OR. NY == 1) THEN
            DO I=1,NX*NY
                B(I)=A(I)*W(I)
            END DO
        END IF

        IF (LNX(1)+LNY(1) <= 1) THEN
            CALL ZTRANSMULA(A,B,W,NX,NY)
        ELSE
            CALL ZTRANSMULB(A,B,W,NX,NY)
        END IF
    END SUBROUTINE ZTRANSMUL

    SUBROUTINE ZTRANSMULA(A,B,W,NX,NY)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(NX,*),B(NY,*),W(NX,*)
        integer, intent(in) :: NX, NY
        complex(kind=dp), dimension(NX, NY) :: A, W
        complex(kind=dp), dimension(NY, NX) :: B
        integer :: I, J

        DO I=1,NX
            DO J=1,NY
                B(J,I)=A(I,J)*W(I,J)
            END DO
        END DO
    END SUBROUTINE ZTRANSMULA

    SUBROUTINE ZTRANSMULB(A,B,W,NX,NY)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(NX,*),B(NY,*),W(NX,*)
        integer, intent(in) :: NX, NY
        complex(kind=dp), dimension(NX, NY) :: A, W
        complex(kind=dp), dimension(NY, NX) :: B
        integer :: I, J

        IF (NY >= NX) THEN
            DO I=0,NX-1
                DO J=1,NX-I
                    B(J,I+J)=A(I+J,J)*W(I+J,J)
                END DO
            END DO
            DO I=1,NY-NX
                DO J=1,NX
                    B(I+J,J)=A(J,I+J)*W(J,I+J)
                END DO
            END DO
            DO I=NY-NX+1,NY-1
                DO J=1,NY-I
                    B(I+J,J)=A(J,I+J)*W(J,I+J)
                END DO
            END DO
        ELSE
            DO I=0,NY-1
                DO J=1,NY-I
                    B(I+J,J)=A(J,I+J)*W(J,I+J)
                END DO
            END DO
            DO I=1,NX-NY
                DO J=1,NY
                    B(J,I+J)=A(I+J,J)*W(I+J,J)
                END DO
            END DO
            DO I=NX-NY+1,NX-1
                DO J=1,NX-I
                    B(J,I+J)=A(I+J,J)*W(I+J,J)
                END DO
            END DO
        END IF
    END SUBROUTINE ZTRANSMULB

end module vzfft1d_mod
