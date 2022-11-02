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
!     FACTORIZATION ROUTINE
!
!     FORTRAN77 SOURCE PROGRAM WRITTEN BY DAISUKE TAKAHASHI
!     Translated to F90 by JuanFRH

module factor_mod
    implicit none

contains

    SUBROUTINE FACTOR(N,IP)
        integer :: N, N2
        integer, dimension(3) :: IP

        IP(1)=0
        IP(2)=0
        IP(3)=0
        N2=N
        IF (MOD(N,2) /= 0 .AND. MOD(N,3) /= 0 .AND. MOD(N,5) /= 0) RETURN
        do while (N2 > 1)
            IF (MOD(N2,2) == 0) THEN
                IP(1)=IP(1)+1
                N2=N2/2
            ELSE IF (MOD(N2,3) == 0) THEN
                IP(2)=IP(2)+1
                N2=N2/3
            ELSE IF (MOD(N2,5) == 0) THEN
                IP(3)=IP(3)+1
                N2=N2/5
            END IF
        end do
    END SUBROUTINE FACTOR

    SUBROUTINE GETNXNY(N,NX,NY)
        integer :: N, NX, NY
        integer :: ISQRTN, IRES, IRES2
        integer :: I, J, K
        integer, dimension(3) :: IP, LNX, LNY

        ISQRTN=IDINT(DSQRT(DBLE(N)))
        CALL FACTOR(N,IP)
        DO I=1,3
            LNX(I)=0
        END DO
        IRES=ISQRTN
        DO K=0,(IP(3)+1)/2
            DO J=0,(IP(2)+1)/2
                DO I=0,(IP(1)+1)/2
                    NX=(2**I)*(3**J)*(5**K)
                    IF (NX <= ISQRTN) THEN
                        IRES2=ISQRTN-NX
                        IF (IRES2 < IRES) THEN
                            LNX(1)=I
                            LNX(2)=J
                            LNX(3)=K
                            IRES=IRES2
                        END IF
                    END IF
                END DO
            END DO
        END DO
        DO I=1,3
            LNY(I)=IP(I)-LNX(I)
        END DO
        NX=(2**LNX(1))*(3**LNX(2))*(5**LNX(3))
        NY=(2**LNY(1))*(3**LNY(2))*(5**LNY(3))
    END SUBROUTINE GETNXNY

end module factor_mod