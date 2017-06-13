module initcond_generator

  use constants
  use randomModule

  implicit none

 contains

   subroutine icpgen(n_particles, radius, xx0, yy0, zz0, xx, yy, zz)
     ! Introduce uniform noise in equilibrum positions. Requieres file containing equilibrium positions.
     ! Data in file must be columnwise
     implicit none
     real(kind=8), dimension(:), intent(in)     :: xx0, yy0, zz0                ! vectors to contain initial conditions
     real(kind=8), dimension(:), intent(inout)  :: xx, yy, zz
     integer, intent(in)                        :: n_particles
     real(kind=8), intent(in)                   :: radius
     integer                                    :: ii
     logical                                    :: calc=.true.
     real(kind=8), dimension(:), allocatable    :: rg
     real(kind=8)                               :: theta, phi
     !open(unit=15, file="eq_pos.dat", action='read')
     call ranseed()
     allocate(rg(1:3))
     do ii = 1, n_particles, 1
       call RANDOM_NUMBER(rg)
       theta = 2*pi*rg(1)
       phi   = 2*pi*rg(2)
       xx(ii) = xx0(ii) + rg(3)*radius*dsin(theta)*dcos(phi)
       yy(ii) = yy0(ii) + rg(3)*radius*dsin(theta)*dsin(phi)
       zz(ii) = zz0(ii) + rg(3)*radius*dcos(theta)
     end do

   end subroutine icpgen

   subroutine icmomgen(n_particles, speed, px, py, pz)
     ! Introduce uniform noise in equilibrum positions. Requieres file containing equilibrium positions.
     ! Data in file must be columnwise
     implicit none
     real(kind=8), dimension(:), intent(inout)  :: px, py, pz
     integer, intent(in)                        :: n_particles
     real(kind=8), intent(in)                   :: speed
     integer                                    :: ii
     logical                                    :: calc=.true.
     real(kind=8), dimension(:), allocatable    :: rg
     real(kind=8)                               :: theta, phi
     !open(unit=15, file="eq_pos.dat", action='read')

     call ranseed()
     allocate(rg(1:3))
     do ii=1, n_particles, 1
       call RANDOM_NUMBER(rg)
       theta = 2*pi*rg(1)
       phi   = 2*pi*rg(2)
       px(ii) = rg(3)*speed*dsin(theta)*dcos(phi)
       py(ii) = rg(3)*speed*dsin(theta)*dsin(phi)
       pz(ii) = rg(3)*speed*dcos(theta)
     end do

   end subroutine icmomgen

end module initcond_generator
