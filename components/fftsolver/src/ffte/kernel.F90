
!     FFTE: A FAST FOURIER TRANSFORM PACKAGE

!     (C) COPYRIGHT SOFTWARE, 2000-2004, 2008-2011, ALL RIGHTS RESERVED
!                BY
!         DAISUKE TAKAHASHI
!         FACULTY OF ENGINEERING, INFORMATION AND SYSTEMS
!         UNIVERSITY OF TSUKUBA
!         1-1-1 TENNODAI, TSUKUBA, IBARAKI 305-8573, JAPAN
!         E-MAIL: daisuke@cs.tsukuba.ac.jp
!
!
!     RADIX-2, 3, 4, 5 AND 8 FFT KERNEL ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM
!
!     WRITTEN BY DAISUKE TAKAHASHI

module kernel_mod
    use, intrinsic :: iso_fortran_env

    implicit none

    integer, parameter, private :: dp = REAL64

    real(kind=dp), parameter :: C31 = 0.86602540378443865D0
    real(kind=dp), parameter :: C32 = 0.5D0
    real(kind=dp), parameter :: C51 = 0.95105651629515357D0
    real(kind=dp), parameter :: C52 = 0.61803398874989485D0
    real(kind=dp), parameter :: C53 = 0.55901699437494742D0
    real(kind=dp), parameter :: C54 = 0.25D0
    real(kind=dp), parameter :: C81 = 0.70710678118654752D0

contains

    SUBROUTINE FFT2(A,B,M)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        integer :: M
        complex(kind=dp), dimension(M,2) :: A, B
        complex(kind=dp) :: C0,C1
        integer :: I

        DO I=1,M
            C0=A(I,1)
            C1=A(I,2)
            B(I,1)=C0+C1
            B(I,2)=C0-C1
        END DO
    END SUBROUTINE FFT2

    SUBROUTINE FFT3A(A,B,W,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        integer :: L
        integer :: J
        complex(kind=dp), dimension(L,3) :: A
        complex(kind=dp), dimension(3,L) :: B
        complex(kind=dp), dimension(2,L) :: W
        complex(kind=dp) :: C0,C1,C2,D0,D1,D2,W1,W2

        DO J=1,L
            W1=W(1,J)
            W2=W(2,J)
            C0=A(J,1)
            C1=A(J,2)
            C2=A(J,3)
            D0=C1+C2
            D1=C0-C32*D0
            D2=(0.0D0,-1.0D0)*C31*(C1-C2)
            B(1,J)=C0+D0
            B(2,J)=W1*(D1+D2)
            B(3,J)=W2*(D1-D2)
        END DO
    END SUBROUTINE FFT3A

    SUBROUTINE FFT3B(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(M,L,*),B(M,3,*),W(2,*)
        integer :: M, L
        integer :: I, J
        complex(kind=dp), dimension(M,L,3) :: A
        complex(kind=dp), dimension(M,3,L) :: B
        complex(kind=dp), dimension(2,L) :: W
        complex(kind=dp) :: C0,C1,C2,D0,D1,D2,W1,W2

        DO I=1,M
            C0=A(I,1,1)
            C1=A(I,1,2)
            C2=A(I,1,3)
            D0=C1+C2
            D1=C0-C32*D0
            D2=(0.0D0,-1.0D0)*C31*(C1-C2)
            B(I,1,1)=C0+D0
            B(I,2,1)=D1+D2
            B(I,3,1)=D1-D2
        END DO
        DO J=2,L
            W1=W(1,J)
            W2=W(2,J)
            DO I=1,M
                C0=A(I,J,1)
                C1=A(I,J,2)
                C2=A(I,J,3)
                D0=C1+C2
                D1=C0-C32*D0
                D2=(0.0D0,-1.0D0)*C31*(C1-C2)
                B(I,1,J)=C0+D0
                B(I,2,J)=W1*(D1+D2)
                B(I,3,J)=W2*(D1-D2)
            END DO
        END DO
    END SUBROUTINE FFT3B

    SUBROUTINE FFT4A(A,B,W,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(L,*),B(4,*),W(3,*)
        integer :: L
        integer :: J
        complex(kind=dp), dimension(L,4) :: A
        complex(kind=dp), dimension(4,L) :: B
        complex(kind=dp), dimension(3,L) :: W
        complex(kind=dp) :: C0,C1,C2,C3,D0,D1,D2,D3,W1,W2,W3

        DO J=1,L
            W1=W(1,J)
            W2=W(2,J)
            W3=W(3,J)
            C0=A(J,1)
            C1=A(J,2)
            C2=A(J,3)
            C3=A(J,4)
            D0=C0+C2
            D1=C0-C2
            D2=C1+C3
            D3=(0.0D0,-1.0D0)*(C1-C3)
            B(1,J)=D0+D2
            B(2,J)=W1*(D1+D3)
            B(3,J)=W2*(D0-D2)
            B(4,J)=W3*(D1-D3)
        END DO
    END SUBROUTINE FFT4A

    SUBROUTINE FFT4B(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(M,L,*),B(M,4,*),W(3,*)
        integer :: M, L
        integer :: I, J
        complex(kind=dp), dimension(M,L,4) :: A
        complex(kind=dp), dimension(M,4,L) :: B
        complex(kind=dp), dimension(3,L) :: W
        complex(kind=dp) :: C0,C1,C2,C3,D0,D1,D2,D3,W1,W2,W3

        DO I=1,M
            C0=A(I,1,1)
            C1=A(I,1,2)
            C2=A(I,1,3)
            C3=A(I,1,4)
            D0=C0+C2
            D1=C0-C2
            D2=C1+C3
            D3=(0.0D0,-1.0D0)*(C1-C3)
            B(I,1,1)=D0+D2
            B(I,2,1)=D1+D3
            B(I,3,1)=D0-D2
            B(I,4,1)=D1-D3
        END DO
        DO J=2,L
            W1=W(1,J)
            W2=W(2,J)
            W3=W(3,J)
            DO I=1,M
                C0=A(I,J,1)
                C1=A(I,J,2)
                C2=A(I,J,3)
                C3=A(I,J,4)
                D0=C0+C2
                D1=C0-C2
                D2=C1+C3
                D3=(0.0D0,-1.0D0)*(C1-C3)
                B(I,1,J)=D0+D2
                B(I,2,J)=W1*(D1+D3)
                B(I,3,J)=W2*(D0-D2)
                B(I,4,J)=W3*(D1-D3)
            END DO
        END DO
    END SUBROUTINE FFT4B

    SUBROUTINE FFT5A(A,B,W,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(L,*),B(5,*),W(4,*)
        integer, intent(in) :: L
        integer :: J
        complex(kind=dp), dimension(L,5) :: A
        complex(kind=dp), dimension(5,L) :: B
        complex(kind=dp), dimension(4,L) :: W
        complex(kind=dp) :: C0,C1,C2,C3,C4,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10
        complex(kind=dp) :: W1,W2,W3,W4

        DO J=1,L
            W1=W(1,J)
            W2=W(2,J)
            W3=W(3,J)
            W4=W(4,J)
            C0=A(J,1)
            C1=A(J,2)
            C2=A(J,3)
            C3=A(J,4)
            C4=A(J,5)
            D0=C1+C4
            D1=C2+C3
            D2=C51*(C1-C4)
            D3=C51*(C2-C3)
            D4=D0+D1
            D5=C53*(D0-D1)
            D6=C0-C54*D4
            D7=D6+D5
            D8=D6-D5
            D9=(0.0D0,-1.0D0)*(D2+C52*D3)
            D10=(0.0D0,-1.0D0)*(C52*D2-D3)
            B(1,J)=C0+D4
            B(2,J)=W1*(D7+D9)
            B(3,J)=W2*(D8+D10)
            B(4,J)=W3*(D8-D10)
            B(5,J)=W4*(D7-D9)
        END DO
    END SUBROUTINE FFT5A

    SUBROUTINE FFT5B(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(M,L,*),B(M,5,*),W(4,*)
        integer :: M, L
        integer :: I, J
        complex(kind=dp), dimension(M,L,5) :: A
        complex(kind=dp), dimension(M,5,L) :: B
        complex(kind=dp), dimension(4,L) :: W
        complex(kind=dp) :: C0,C1,C2,C3,C4,D0,D1,D2,D3,D4,D5,D6,D7,D8,D9,D10
        complex(kind=dp) :: W1,W2,W3,W4

        DO I=1,M
            C0=A(I,1,1)
            C1=A(I,1,2)
            C2=A(I,1,3)
            C3=A(I,1,4)
            C4=A(I,1,5)
            D0=C1+C4
            D1=C2+C3
            D2=C51*(C1-C4)
            D3=C51*(C2-C3)
            D4=D0+D1
            D5=C53*(D0-D1)
            D6=C0-C54*D4
            D7=D6+D5
            D8=D6-D5
            D9=(0.0D0,-1.0D0)*(D2+C52*D3)
            D10=(0.0D0,-1.0D0)*(C52*D2-D3)
            B(I,1,1)=C0+D4
            B(I,2,1)=D7+D9
            B(I,3,1)=D8+D10
            B(I,4,1)=D8-D10
            B(I,5,1)=D7-D9
        END DO
        DO J=2,L
            W1=W(1,J)
            W2=W(2,J)
            W3=W(3,J)
            W4=W(4,J)
            DO I=1,M
                C0=A(I,J,1)
                C1=A(I,J,2)
                C2=A(I,J,3)
                C3=A(I,J,4)
                C4=A(I,J,5)
                D0=C1+C4
                D1=C2+C3
                D2=C51*(C1-C4)
                D3=C51*(C2-C3)
                D4=D0+D1
                D5=C53*(D0-D1)
                D6=C0-C54*D4
                D7=D6+D5
                D8=D6-D5
                D9=(0.0D0,-1.0D0)*(D2+C52*D3)
                D10=(0.0D0,-1.0D0)*(C52*D2-D3)
                B(I,1,J)=C0+D4
                B(I,2,J)=W1*(D7+D9)
                B(I,3,J)=W2*(D8+D10)
                B(I,4,J)=W3*(D8-D10)
                B(I,5,J)=W4*(D7-D9)
            END DO
        END DO
        
    END SUBROUTINE FFT5B

    SUBROUTINE FFT8A(A,B,W,L)
        !IMPLICIT REAL*8 (A-H,O-Z)
        !COMPLEX*16 A(L,*),B(8,*),W(7,*)
        integer :: L
        complex(kind=dp), dimension(L,8) :: A
        complex(kind=dp), dimension(8,L) :: B
        complex(kind=dp), dimension(7,L) :: W
        integer :: J
        complex(kind=dp) :: C0,C1,C2,C3,C4,C5,C6,C7,D0,D1,D2,D3,D4,D5,D6,D7
        complex(kind=dp) :: E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,W1,W2,W3,W4,W5,W6,W7

        DO J=1,L
            W1=W(1,J)
            W2=W(2,J)
            W3=W(3,J)
            W4=W(4,J)
            W5=W(5,J)
            W6=W(6,J)
            W7=W(7,J)
            C0=A(J,1)
            C1=A(J,2)
            C2=A(J,3)
            C3=A(J,4)
            C4=A(J,5)
            C5=A(J,6)
            C6=A(J,7)
            C7=A(J,8)
            D0=C0+C4
            D1=C0-C4
            D2=C2+C6
            D3=(0.0D0,-1.0D0)*(C2-C6)
            D4=C1+C5
            D5=C1-C5
            D6=C3+C7
            D7=C3-C7
            E0=D0+D2
            E1=D0-D2
            E2=D4+D6
            E3=(0.0D0,-1.0D0)*(D4-D6)
            E4=C81*(D5-D7)
            E5=(0.0D0,-1.0D0)*C81*(D5+D7)
            E6=D1+E4
            E7=D1-E4
            E8=D3+E5
            E9=D3-E5
            B(1,J)=E0+E2
            B(2,J)=W1*(E6+E8)
            B(3,J)=W2*(E1+E3)
            B(4,J)=W3*(E7-E9)
            B(5,J)=W4*(E0-E2)
            B(6,J)=W5*(E7+E9)
            B(7,J)=W6*(E1-E3)
            B(8,J)=W7*(E6-E8)
        END DO
    END SUBROUTINE FFT8A

    SUBROUTINE FFT8B(A,B,W,M,L)
        ! IMPLICIT REAL*8 (A-H,O-Z)
        ! COMPLEX*16 A(M,L,*),B(M,8,*),W(7,*)
        integer :: M, L
        integer :: I, J
        complex(kind=dp), dimension(M,L,8) :: A
        complex(kind=dp), dimension(M,8,L) :: B
        complex(kind=dp), dimension(7,L) :: W
        complex(kind=dp) :: C0,C1,C2,C3,C4,C5,C6,C7,D0,D1,D2,D3,D4,D5,D6,D7
        complex(kind=dp) :: E0,E1,E2,E3,E4,E5,E6,E7,E8,E9,W1,W2,W3,W4,W5,W6,W7

        DO I=1,M
            C0=A(I,1,1)
            C1=A(I,1,2)
            C2=A(I,1,3)
            C3=A(I,1,4)
            C4=A(I,1,5)
            C5=A(I,1,6)
            C6=A(I,1,7)
            C7=A(I,1,8)
            D0=C0+C4
            D1=C0-C4
            D2=C2+C6
            D3=(0.0D0,-1.0D0)*(C2-C6)
            D4=C1+C5
            D5=C1-C5
            D6=C3+C7
            D7=C3-C7
            E0=D0+D2
            E1=D0-D2
            E2=D4+D6
            E3=(0.0D0,-1.0D0)*(D4-D6)
            E4=C81*(D5-D7)
            E5=(0.0D0,-1.0D0)*C81*(D5+D7)
            E6=D1+E4
            E7=D1-E4
            E8=D3+E5
            E9=D3-E5
            B(I,1,1)=E0+E2
            B(I,2,1)=E6+E8
            B(I,3,1)=E1+E3
            B(I,4,1)=E7-E9
            B(I,5,1)=E0-E2
            B(I,6,1)=E7+E9
            B(I,7,1)=E1-E3
            B(I,8,1)=E6-E8
        END DO
        DO J=2,L
            W1=W(1,J)
            W2=W(2,J)
            W3=W(3,J)
            W4=W(4,J)
            W5=W(5,J)
            W6=W(6,J)
            W7=W(7,J)
            DO I=1,M
                C0=A(I,J,1)
                C1=A(I,J,2)
                C2=A(I,J,3)
                C3=A(I,J,4)
                C4=A(I,J,5)
                C5=A(I,J,6)
                C6=A(I,J,7)
                C7=A(I,J,8)
                D0=C0+C4
                D1=C0-C4
                D2=C2+C6
                D3=(0.0D0,-1.0D0)*(C2-C6)
                D4=C1+C5
                D5=C1-C5
                D6=C3+C7
                D7=C3-C7
                E0=D0+D2
                E1=D0-D2
                E2=D4+D6
                E3=(0.0D0,-1.0D0)*(D4-D6)
                E4=C81*(D5-D7)
                E5=(0.0D0,-1.0D0)*C81*(D5+D7)
                E6=D1+E4
                E7=D1-E4
                E8=D3+E5
                E9=D3-E5
                B(I,1,J)=E0+E2
                B(I,2,J)=W1*(E6+E8)
                B(I,3,J)=W2*(E1+E3)
                B(I,4,J)=W3*(E7-E9)
                B(I,5,J)=W4*(E0-E2)
                B(I,6,J)=W5*(E7+E9)
                B(I,7,J)=W6*(E1-E3)
                B(I,8,J)=W7*(E6-E8)
            END DO
        END DO
    END SUBROUTINE FFT8B
end module kernel_mod