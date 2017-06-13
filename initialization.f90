
! module containing variable initializations

module initialization
  use constants
  implicit none

contains
  ! based on http://jblevins.org/log/control-file
  subroutine initialize_system(n_particles, mass1, charge1, tt, dt, traj, save_freq, long_freq, alphay, alphaz, ic_radius, initT)

    implicit none

    character(len=100) :: buffer, local
    character(len=100) :: label
    integer :: pos
    integer, parameter :: fh=15
    integer :: ios=0, line=0

    integer, intent(inout) :: n_particles   ! dimensionality of the problem, number of partilcles

    real(kind=8), intent(inout)  :: mass1, charge1, ic_radius
    real(kind=8), intent(inout)  :: alphay, alphaz, long_freq, initT
    real(kind=8), intent(inout)       :: dt, tt
    integer, intent(inout)            :: traj
    integer, intent(inout)            :: save_freq


    open(unit=fh, file='values.dat',action='read')

    do while(ios==0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios==0) then
        line = line + 1

        ! find the first isntance of whitespace. Split and label data
        pos=scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)
        select case(label)
          case('n_particles')
              read(buffer,*,iostat=ios) n_particles
          case('mass1')
              read(buffer,*,iostat=ios) mass1
          case('charge1')
              read(buffer,*,iostat=ios) charge1
          case('dt')
              read(buffer,*,iostat=ios) dt
          case('tt')
              read(buffer,*,iostat=ios) tt
          case('traj')
              read(buffer,*,iostat=ios) traj
          case('alphay')
              read(buffer,*,iostat=ios) alphay
          case('alphaz')
                read(buffer,*,iostat=ios) alphaz
          case('icradius')
              read(buffer,*,iostat=ios) ic_radius
          case('long_freq')
              read(buffer,*,iostat=ios) long_freq
          case('save_freq')
              read(buffer,*,iostat=ios) save_freq
          case('initT')
              read(buffer,*,iostat=ios) initT
          case default
              print*, "Skipping invalid value."
        end select
      end if
    end do

    close(unit=fh)

  end subroutine initialize_system

  subroutine initialize_laser_chain(del1, del2, delC, Gam, omega, I1, I2, IC)
    ! initialize laser and ion interaction parameters

    implicit none

    character(len=100) :: buffer, local
    character(len=100) :: label
    integer :: pos
    integer, parameter :: fh=15
    integer :: ios=0, line=0
    real(kind=8), intent(inout) :: del1, del2, delC
    real(kind=8), intent(inout)  :: Gam, omega
    real(kind=8), intent(inout) :: I1, I2, IC

    open(unit=fh, file='chain_laser.dat',action='read')

    do while(ios==0)
      read(fh, '(A)', iostat=ios) buffer
      if (ios==0) then
        line = line + 1

        ! find the first isntance of whitespace. Split and label data
        pos=scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)
        select case(label)
        case('detuning_left')
              read(buffer,*,iostat=ios) del1
          case('detuning_right')
              read(buffer,*,iostat=ios) del2
          case('detuning_cool')
              read(buffer,*,iostat=ios) delC
          case('linewidth')
              read(buffer,*,iostat=ios) Gam
          case('target_frequency')
              read(buffer,*,iostat=ios) omega
          case('intensity_left')
              read(buffer,*,iostat=ios) I1
          case('intensity_right')
                read(buffer,*,iostat=ios) I2
          case('intensity_cool')
                read(buffer,*,iostat=ios) IC
          case default
                print*, "Skipping invalid value."
        end select
      end if
    end do

    print*,"read"
  end subroutine initialize_laser_chain

  subroutine dimensionless_doppler_values(eta, D, mass, long_freq, char_length, aeta, aD)
    implicit none
    real(kind=8), intent(in)              :: eta, D, mass, long_freq, char_length
    real(kind=8), intent(inout)           :: aeta, aD
    ! calculate scaled friction and difussion coefficients
    aeta = eta / (mass * long_freq)
    aD   = D / (char_length*char_length*mass*mass*long_freq*long_freq*long_freq)

  end subroutine dimensionless_doppler_values


  subroutine initcond(qq, pp)
    implicit none
    real(kind=8), dimension(:), intent(inout) :: qq, pp
    integer, parameter :: fh=15

    open(unit=fh,file='initcond.dat', action='read')
    read(fh,*) qq
    read(fh,*) pp


    close(unit=1)
  end subroutine initcond

end module initialization
