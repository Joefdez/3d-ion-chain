program twoDChain
  use constants
  use randomModule
  use initcond_generator
  implicit none

  ! variables
  real(kind=8), dimension(:), allocatable               :: xx0, yy0, px0, py0
  real(kind=8), dimension(:,:), allocatable             :: xx, yy, ppx, ppy
  real(kind=8), dimension(:), allocatable               :: xxi, yyi, ppxi, ppyi
  real(kind=8), dimension(:,:), allocatable             :: fx1, fy1, fx2, fy2
  real(kind=8), dimension(:), allocatable               :: fx, fy
  real(kind=8), dimension(:), allocatable               :: Axx, Axxi, Ayy, Ayyi
  real(kind=8), dimension(:), allocatable               :: Apx, Apxi, Apy, Apyi

  real(kind=8)                                          :: tt, dt, mass, charge, dist
  real(kind=8)                                          :: alpha
  integer                                               :: nsteps, nparticles
  integer                                               :: ii,jj, kk

  tt =1000.0d0
  dt = 1.0d-3
  nsteps = int(tt/dt)
  nparticles = 30
  alpha = 7.00d0
  print*, "alpha", alpha
  allocate(xx0(1:nparticles))
  xx0 = 0.0d0
  allocate(yy0(1:nparticles))
  yy0 = 0.0d0
  allocate(px0(1:nparticles))
  px0 = 0.0d0
  allocate(py0(1:nparticles))
  py0 = 0.0d0
  allocate(xx(1:nparticles, 1:nsteps))
  xx = 0.0d0
  allocate(yy(1:nparticles, 1:nsteps))
  yy = 0.0d0
  allocate(ppx(1:nparticles, 1:nsteps))
  ppx = 0.0d0
  allocate(ppy(1:nparticles, 1:nsteps))
  ppy = 0.0d0
  allocate(xxi(1:nparticles))
  xxi = 0.0d0
  allocate(yyi(1:nparticles))
  yyi = 0.0d0
  allocate(ppxi(1:nparticles))
  ppxi = 0.0d0
  allocate(ppyi(1:nparticles))
  ppyi = 0.0d0
  allocate(fx1(1:nparticles,1:nparticles))
  fx1 = 0.0d0
  allocate(fx2(1:nparticles,1:nparticles))
  fx2 = 0.0d0
  allocate(fy1(1:nparticles,1:nparticles))
  fy1 = 0.0d0
  allocate(fy2(1:nparticles,1:nparticles))
  fy2 = 0.0d0
  allocate(fx(1:nparticles))
  fx = 0.0d0
  allocate(fy(1:nparticles))
  fy = 0.0d0
  allocate(Axx(1:nparticles))
  Axx = 0.0d0
  allocate(Ayy(1:nparticles))
  Ayy = 0.0d0
  allocate(Apx(1:nparticles))
  Apx = 0.0d0
  allocate(Apy(1:nparticles))
  Apy = 0.0d0
  allocate(Axxi(1:nparticles))
  Axxi = 0.0d0
  allocate(Ayyi(1:nparticles))
  Ayyi = 0.0d0
  allocate(Apxi(1:nparticles))
  Apxi = 0.0d0
  allocate(Apyi(1:nparticles))
  Apyi = 0.0d0

  open(unit=15, file="eq_pos.dat", action='read')
  do ii = 1, nparticles, 1
    read(15,*) xx0(ii), yy0(ii)
  end do
  close(unit=15)


  call icpgen(nparticles, 0.02d0, xx0, yy0, xx(:,1), yy(:,1))
  do ii=1, nsteps-1, 1
    !print*, ii
    call coulombtwo(nparticles, xx(:,ii), yy(:,ii), fx1, fy1)
    fx = 0.0d0
    fy = 0.0d0
    do jj=1, nparticles, 1
      fx(jj) = sum(fx1(jj,:), 1)
      fy(jj) = sum(fy1(jj,:), 1)
    end do
    Axx = ppx(:,ii)
    Ayy = ppy(:,ii)
    Apx = -1.0d0*xx(:,ii) + fx
    Apy = -1.0*alpha*alpha*yy(:,ii) + fy
    xxi = xx(:,ii) + Axx*dt
    yyi = yy(:,ii) + Ayy*dt
    ppxi = ppx(:,ii) + Apx*dt
    ppyi = ppy(:,ii) + Apy*dt
    fx = 0.0d0
    fy = 0.0d0
    call coulombtwo(nparticles, xxi, yyi, fx2, fy2)
    do jj=1, nparticles, 1
      fx(jj) = sum(fx2(jj,:), 1)
      fy(jj) = sum(fy2(jj,:), 1)
    end do
    Axxi = ppxi
    Ayyi = ppyi
    Apxi = -1.0d0*xxi + fx
    Apyi = -1.0*alpha*alpha*yyi + fy
    xx(:,ii+1)  = xx(:,ii) + 0.5d0*(Axx + Axxi)*dt
    yy(:,ii+1)  = yy(:,ii) + 0.5d0*(Ayy + Ayyi)*dt
    ppx(:,ii+1) = ppx(:,ii) + 0.5d0*(Apx + Apxi)*dt
    ppy(:,ii+1) = ppy(:,ii) + 0.5d0*(Apy + Apyi)*dt
  end do

  open(unit=1, file="results/posx.dat")
  open(unit=2, file="results/posy.dat")
  do jj=1, nparticles,1
    write(1,*) xx(jj,:)
    write(2,*) yy(jj,:)
  end do
  close(unit=1)
  close(unit=2)

  open(unit=1, file="results/posx_f.dat")
  open(unit=2, file="results/posy_f.dat")
  do jj=1, nparticles, 1
    do kk=1, int(0.2*nsteps), 1
      write(1,*) xx(jj,int(0.8*nsteps)+kk)
      write(2,*) yy(jj,int(0.8*nsteps)+kk)
    end do
  end do
  close(unit=1)
  close(unit=2)
  contains

  subroutine coulomb(nparticles, xx, yy, fx, fy)
    implicit none
    real(kind=8), dimension(:), intent(in)              :: xx, yy
    real(kind=8), dimension(:), intent(inout)           :: fx, fy
    integer, intent(in)                                 :: nparticles
    integer                                             :: ii, jj
    real(kind=8)                                        :: dist
    fx = 0.0d0
    fy = 0.0d0
    do ii=1, nparticles, 1
      do jj=ii+1, nparticles, 1
        dist = 0.0d0
        dist = sqrt( (xx(ii)-xx(jj))**2 + (yy(ii)-yy(jj))**2 )
        dist = dist*dist*dist
        fx(ii) = fx(ii) +  (xx(ii)-xx(jj)) / dist
        fx(jj) = fx(jj) -  (xx(ii)-xx(jj)) / dist
        fy(ii) = fy(ii) +  (yy(ii)-yy(jj)) / dist
        fy(ii) = fy(jj) -  (yy(ii)-yy(jj)) / dist
      end do
    end do

  end subroutine coulomb

  subroutine coulombtwo(nparticles, xx, yy, fx, fy)
    implicit none
    real(kind=8), dimension(:), intent(in)              :: xx, yy
    real(kind=8), dimension(:,:), intent(inout)           :: fx, fy
    integer, intent(in)                                 :: nparticles
    integer                                             :: ii, jj
    real(kind=8)                                        :: dist
    fx = 0.0d0
    fy = 0.0d0
    do ii=1, nparticles, 1
      do jj=ii+1, nparticles, 1
        dist = 0.0d0
        dist = sqrt( (xx(ii)-xx(jj))**2 + (yy(ii)-yy(jj))**2 )
        dist = dist*dist*dist
        fx(ii,jj) = (xx(ii)-xx(jj)) / dist
        fx(jj,ii) = -1.0d0*(xx(ii)-xx(jj)) / dist
        fy(ii,jj) = (yy(ii)-yy(jj)) / dist
        fy(jj,ii) = -1.0d0*(yy(ii)-yy(jj)) / dist
      end do
    end do

  end subroutine coulombtwo

end program twoDChain
