! This module contains physical constants

module constants
    implicit none

    real(kind=8), parameter :: kb   = 1.38064852D-23 ! m2 kg s-2 K-1, Boltzsmann's constant
    real(kind=8), parameter :: sig  = 5.670367D-8 ! W m-2 K-4, Stefn-Boltzmann constant
    real(kind=8), parameter :: ee   = 1.60217662D-19 ! C, Electron charge
    real(kind=8), parameter :: uu   = 1.660539040D-27 ! kg, unified atomic mass
    real(kind=8), parameter :: me   = 9.10938356D-31 ! kg, electron rest mass
    real(kind=8), parameter :: hbar = 1.054571800D-34 ! Joules  second, reduced Planck's constant
    real(kind=8), parameter :: hh   = 6.626070040e-34 ! Joules second, Planck's constant
    real(kind=8), parameter :: cc   = 299792458       ! m/s, speed of light in vacuum
    real(kind=8), parameter :: ep0  = 8.854187817D-12 ! farad/meter, vacuum permitivity
    real(kind=8), parameter :: pi  = ATAN(1.)*4.0D0 ! pi


end module
