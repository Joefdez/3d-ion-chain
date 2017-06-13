! Module for obtaining and treating random numbers

Module randomModule

use constants
implicit none


contains

  subroutine ranseed()
    ! The following code is taken from http://web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html
    ! Essential mathematics / Computational laboratory, Paul Stevenson (October 2014)
    implicit none
    integer :: i_seed
    integer(kind=4), dimension(:), allocatable :: a_seed
    integer, dimension(1:8) :: dt_seed

    call random_seed(size=i_seed)
    allocate(a_seed(1:i_seed))
    call random_seed(get=a_seed)
    call date_and_time(values=dt_seed)
    a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
    call random_seed(put=a_seed)
    deallocate(a_seed)

  end subroutine ranseed


  subroutine bM(g1, g2)                                 ! Simple subroutine to return 2 random numbers
    implicit none
    real(kind=8), dimension(0:2) :: gg
    real(kind=8), intent(out) :: g1, g2

    call RANDOM_NUMBER(gg)

    g1 = sqrt( (-2)*log(gg(0)) ) * cos (2*pi*gg(1))
    g2 = sqrt( (-2)*log(gg(0)) ) * sin (2*pi*gg(1))



  end subroutine bM

  subroutine pBM(g1, g2)!, rescale= .FALSE. , stdev=0.0, av=0.0)

    implicit none
    real(kind=8), dimension(0:2) :: gg
    real(kind=8), intent(out) :: g1, g2
    !real(kind=8), intent(in) :: stdev, av
    real(kind=8) :: w

    CALL RANDOM_NUMBER(gg)

    do
      g1 = 2.0*gg(0) -1.0
      g2 = 2.0*gg(1) -1.0
      w = g1*g1 + g2*g2
      if ( w .lt. 1.0d0 .and. w .ne. 0.0d0 ) exit
    end do

    w = sqrt ( -2.0 * log(w)/w)
    g1 = g1 * w
    g2 = g2 * w

    !if (rescale)
    !  g1 = g1*stdev + av        ! Transform to sig = stdev and mu = av
    !  g2 = g2*stdev + av
    !end if

  end subroutine pBM

  !subroutine bMA()
    !implicit none
    ! to contain array version of prevous algorithm
  !end subroutine bMA



END MODULE randomModule
