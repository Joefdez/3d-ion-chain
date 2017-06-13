program mp
  implicit none
  real(kind=8), dimension(:), allocatable       :: xx, yy
  integer                                         :: nparticles, nsteps, nstepss
  real(kind=8)                                    :: xl, xr, yb, yt, dx, dy, gauss
  real(kind=8)                                    :: ff, fs, cx, cy, xn, yn
  integer                                         :: ii, jj, kk, ll, nx, ny

  nsteps     = int(500.0d0/(1.0d-3))
  nparticles = 30
  nstepss = int(0.2*nsteps)

  allocate(xx(1:nparticles*nstepss))
  allocate(yy(1:nparticles*nstepss))

  xl = -6.0d0
  xr = 6.0d0
  yb = -0.2d0
  yt = 0.2d0
  nx = 200
  ny = 100
  dx = (xr-xl)/(nx-1)
  dy = (yt-yb)/(ny-1)

  open(unit=11, file="results/posx_f.dat", action='read')
  open(unit=12, file="results/posy_f.dat", action='read')
    do ii=1, nparticles*nstepss, 1
      print*, ii
      read(11,*) xx(ii)
      read(12,*) yy(ii)
    end do
  close(unit=11)
  close(unit=12)

  open(unit=11, file="results/plot.dat")
  do ii=1, nx
    print*, ii
    xn = xl + dx*(ii-1)
    do jj=1, ny
      yn = yb + dy*(jj-1)
      fs = 0.0d0
      do kk=1, nstepss
          cx=(xn-xx(kk))**2/(2.0d0*gauss)
          cy=(yn-yy(kk))**2/(2.0d0*gauss)
          fs = fs + dexp(-cx-cy)
      end do
      write(11,*) xn, yn, fs
    end do
  end do
  close(unit=11)


end program mp
