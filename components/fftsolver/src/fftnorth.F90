!module that provides wrapper functions for FFTE calls
!
! NOTE: FFTE only has functionality for complex to complex in-place FFTs.
!       This module allows for real to complex FFTs by providing wrappers for
!       the complex to complex calls
!
! example usage: 
!       Takes the r2c of a real array "input", producing complex array "output"
!       Then takes the c2r of "output", returning the real array "input"
!  ________________________________________________________________________________
! |integer :: n                                                                    |
! |double precision :: input(n)                                                    |
! |complex*16 :: output(n/2+1)                                                     |
! |                                                                                |
! |call ffte_init(n) !set up ffte for this problem size                            |
! |call ffte_r2c(input, output, n) ! real to complex (forward) FFT                 |
! |call ffte_c2r(output,input, n) ! complex to real (reverse) FFT                  |
! |call ffte_finalise() !clean up work arrays                                      |
!  --------------------------------------------------------------------------------
!
module fftnorth_mod
    use vzfft1d_mod, only : ZFFT1D
    use fftpack, only: zffti, zfftf, zfftb
    use fftpack_kind, only: rk

    implicit none

    integer, parameter, public :: FFT_FFTE = 0
    integer, parameter, public :: FFT_FFTPACK = 1
    
    complex(rk), allocatable, dimension(:) :: wk   ! Work array for the FFT calculations
    complex(rk), allocatable, dimension(:) :: data ! The array that will have the in-place FFT applied to it

    real(rk), allocatable, dimension(:) :: wsave

    logical :: is235

    contains
    
    ! Initialises the FFT routines for the real problem size (n)
    subroutine fftn_init(n, fft)
        integer, intent(in) :: n, fft
        
        select case (fft)
        case (FFT_FFTE)
            !check that n is a valid size for FFTE to handle
            is235 = fftn_check_factors(n)
        case (FFT_FFTPACK)
            is235 = .false.
        end select

        allocate(data(n))
        if (is235) then
            !allocate work and data arrays
            allocate(wk(2*n))
            call ZFFT1D(data,n,0,wk) !initialise FFTE
        else
            allocate(wsave(4*n+15))
            call zffti(n, wsave)
        endif
    end subroutine
    
    !computes a real-to-complex (e.g. forward) FFT
    ! in  : double precision real array of size n (input)
    ! out : double precision complex array of size n/2+1 (output)
    ! n   : integer - size of in (input)
    subroutine fftn_r2c(in, out, n)
        integer, intent(in) :: n
        double precision, intent(in) :: in(n)
        complex*16, intent(out) :: out(n/2+1)

        integer :: i

        !copy real input into complex "data" array
        do i=1,n
            data(i) = dcmplx(in(i), 0.d0)
        enddo

        if(is235) then
            !compute forward FFT (in-place on data)
            call ZFFT1D(data,n,-1,wk)
        else
            call zfftf(n,data,wsave)
        endif

        !extract the first n/2+1 terms from the FFT and return them as output
        !(The last n/2-1 terms of r2c are complex conjugates of the previous 
        ! ones and so are redundant)
        out(:) = data(1:n/2+1)

    end subroutine
    
    !computes a complex-to-real (e.g. inverse) FFT
    ! in  : double precision complex array of size n/2+1 (input)
    ! out : double precision real array of size n (output)
    ! n   : integer - size of out (input)
    subroutine fftn_c2r(in,out,n)
        integer, intent(in) :: n
        complex*16, intent(in) :: in(n/2+1)
        double precision, intent(out) :: out(n)

        integer :: i

        !construct the array to be inverse FFT'd

        !copy the first n/2+1 terms
        data(1:n/2+1) = in(:)
        !set the remaining entries to the complex conjugates of the previous ones
        do i=n/2+2,n
            data(i) = dconjg(in(n-i+2))
        enddo
        
        if(is235) then
            !do the inverse fft (in place on data)
            call ZFFT1D(data,n,1,wk)
            !extract the real part of data and place it in out
            do i=1,n
               out(i) = real(data(i))
            enddo
        else
            call zfftb(n,data,wsave)
            !extract the real part of data and place it in out
            do i=1,n
               out(i) = real(data(i))/n
            enddo
        endif
        
        

    end subroutine
    
    !Finalises the FFTE routines
    subroutine fftn_finalise()
        
        !deallocate the work arrays
        deallocate(data)
        if(is235) then
            deallocate(wk)
        else
            deallocate(wsave)
        endif

    end subroutine
    
    !Checks to see if the input, n, only has prime factors of 2, 3 and 5
    ! If it does, return true. Otherwise, return false
    logical function fftn_check_factors(n)
        integer, intent(in) :: n
        integer :: m

        m=n

        do while (mod(m,5) .eq. 0)
            m = m/5
        enddo

        do while (mod(m,3) .eq. 0)
            m = m/3
        enddo

        do while (mod(m,2) .eq. 0)
            m = m/2
        enddo

        if (m .eq. 1) then
            fftn_check_factors= .true.
        else
            fftn_check_factors= .false.
        endif

    end function

end module fftnorth_mod
