program threeDChain

  use mpi
  use constants
  use initcond_generator
  use initialization
  use randomModule
  use support_functions_threed

  implicit none

  include "declarations.f90"

  ! mpi variables
  integer :: rank, procs, status(MPI_STATUS_SIZE), alloc_err, source, ierr
  call MPI_INIT(ierr)                                                                               ! Neccesary mpi initialization calls
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, procs, ierr)

  if(rank .eq. 0 ) then
    print*, "Reading solution and system parameters"
    call initialize_system(nparticles, mass, charge, tt, dt, traj, save_freq, long_freq, alphay, alphaz, ic_radius, initT)
  end if

  call mpi_bcast(nparticles, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(traj, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(nsteps, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(save_freq, 1, mpi_integer, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(mass, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(charge, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(tt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(dt, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alphay, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(alphaz, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(ic_radius, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(long_freq, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(initT, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  call mpi_barrier(mpi_comm_world, ierr)

  mass   = mass*uu
  charge = charge*ee
  long_freq = 2.0d0*pi*long_freq
  dst = sqrt(dt)
  nsteps=int(tt/dt)
  char_length = ((charge*charge/(4.0d0*pi*ep0))/(mass*long_freq*long_freq))**(1.0/3.0)
  nssteps = int(nsteps/save_freq)  ! Not saving every single timestep saves memory. Must ask about this
  fin = 0.8d0*nsteps

  if(rank .eq. 0) then
    print*, "Reading laser parameters"
    call initialize_laser_chain(del1, del2, delC, Gam, omega0, I1, I2, IC)
  end if

  call mpi_bcast(del1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(del2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(delC, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(Gam, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(omega0, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(I1, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(I2, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(IC, 1, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  Gam  = 2.0d0 * pi * Gam
  del1  = del1 * Gam
  del2  = del2 * Gam
  delC  = delC * Gam
  omega0 = 2.0d0 * pi * omega0
  omega1 = del1 + omega0
  omega2 = del2 + omega0
  omegaC = delC + omega0

  k1 =  omega1 / cc   ! Not really wavelengths, dimensionally they are wavenumbers
  k2 =  omega2 / cc
  kC = omegaC  / cc
 ! Calculate diffusion and friction coefficients
  eta1 = -4.0d0*hbar*k1*k1*I1*(2.0d0*del1/Gam)/( (1 + 4.0d0*del1*del1/(Gam*Gam)) * (1 + 4.0d0*del1*del1/(Gam*Gam)) )
  eta2 = -4.0d0*hbar*k2*k2*I2*(2.0d0*del2/Gam)/( (1 + 4.0d0*del2*del2/(Gam*Gam)) * (1 + 4.0d0*del2*del2/(Gam*Gam)) )
  D1   = hbar*hbar*k1*k1*I1*(Gam)/(1.0d0 + 4.0d0*del1*del1/(Gam*Gam))
  D2   = hbar*hbar*k2*k2*I2*(Gam)/(1.0d0 + 4.0d0*del2*del2/(Gam*Gam))
  etaC = -4.0d0*hbar*kC*kC*IC*(2.0d0*delC/Gam)/( (1 + 4.0d0*delC*delC/(Gam*Gam)) * (1 + 4.0d0*delC*delC/(Gam*Gam)) )
  DC   = hbar*hbar*kC*kC*IC*(Gam)/(1.0d0 + 4.0d0*delC*delC/(Gam*Gam))
  initT = initT * (kb/(char_length*char_length*mass*long_freq*long_freq))
  initSpeed = sqrt(2*initSpeed)

  ! Calculate dimensionless Doppler cooling parameters
  call dimensionless_doppler_values(eta1, D1, mass, long_freq, char_length, aeta1, aD1)
  call dimensionless_doppler_values(eta2, D2, mass, long_freq, char_length, aeta2, aD2)
  call dimensionless_doppler_values(etaC, DC, mass, long_freq, char_length, aetaC, aDc)

  local_traj = traj/procs
  rem = mod(traj, procs)
  if (rank .lt. rem) local_traj = local_traj + 1
  n_elems = nssteps*nparticles
  nbath = 3

  include 'allocation.f90'

  stermsCx = sqrt(2.0d0*aDc)
  stermsCy = sqrt(2.0d0*aDc)
  stermsCz = sqrt(2.0d0*aDc)
  stermsBx(1:nbath) = sqrt(2.0d0*aD1)
  stermsBy(1:nbath) = sqrt(2.0d0*aD1)
  stermsBz(1:nbath) = sqrt(2.0d0*aD1)
  stermsBx((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)
  stermsBy((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)
  stermsBz((nparticles-nbath+1):nparticles) = sqrt(2.0d0*aD2)

  if(rank .eq. 0) then
     seconds = mpi_wtime()
     seconds1 = mpi_wtime()
  end if

  if (rank .eq. 0) then
    open(unit=15, file="eq_pos.dat", action='read')
    do ii = 1, nparticles, 1
      read(15,*) xx0(ii), yy0(ii), zz0(ii)
    end do
  end if

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_bcast(xx0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(yy0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(zz0, nparticles, mpi_double_precision, 0, MPI_COMM_WORLD, ierr)

  call sleep(rank) ! Delay for getting different seeds

  do kk=1, local_traj, 1
    print*, "Proc.", rank, "on trajectory", kk
    xxs = 0.0d0
    yys = 0.0d0
    zzs = 0.0d0
    ppxold = 0.0d0
    ppyold = 0.0d0
    ppzold = 0.0d0
    ll = 0
    mm = 1
    JJix = 0.0d0
    JJiy = 0.0d0
    JJiz = 0.0d0
    call icpgen(nparticles, 0.005d0, xx0, yy0, zz0, xxold, yyold, zzold)
    call icmomgen(nparticles, initSpeed, ppxold, ppyold, ppzold)
    call ranseed()
!   JJix_av  = 0.0d0
!   JJiy_av = 0.0d0
    do ii=1, nsteps, 1
      call coulombM(nparticles, xxold, yyold, zzold, fx1, fy1, fz1, invD1)
      fx = 0.0d0
      fy = 0.0d0
      fz = 0.0d0
      do jj=1, nparticles, 1
        fx(jj) = sum(fx1(jj,:), 1)
        fy(jj) = sum(fy1(jj,:), 1)
        fz(jj) = sum(fz1(jj,:), 1)
      end do
      call vecA(xxold, yyold, zzold, ppxold, ppyold, ppzold, fx, fy, fz,&
                alphay, alphaz, aeta1, aeta2, aetaC, nbath, nparticles, Axx, Ayy, Azz, Apx, Apy, Apz)
      call vecB_edges(dst, nparticles, dOmx, dOmy, dOmz)
      !call vecB_cool(dst, nparticles, dOmxc, dOmyc)
      xxi  = xxold + Axx*dt
      yyi  = yyold + Ayy*dt
      zzi  = zzold + Azz*dt
      ppxi = ppxold + Apx*dt  + stermsBx*dOmx !+ stermsCx*dOmxc
      ppyi = ppyold + Apy*dt  + stermsBy*dOmy !+ stermsCy*dOmyc
      ppzi = ppzold + Apz*dt  + stermsBz*dOmz !+ stermsCy*dOmyc
      fx = 0.0d0
      fy = 0.0d0
      fz = 0.0d0
      call coulombM(nparticles, xxi, yyi, zzi, fx2, fy2, fz2, invD2)
      do jj=1, nparticles, 1
        fx(jj) = sum(fx2(jj,:), 1)
        fy(jj) = sum(fy2(jj,:), 1)
        fz(jj) = sum(fz2(jj,:), 1)
      end do
      call vecA(xxi, yyi, zzi, ppxi, ppyi, ppzi, fx, fy, fz,&
                alphay, alphaz, aeta1, aeta2, aetaC, nbath, nparticles, Axxi, Ayyi, Azzi, Apxi, Apyi, Apzi)
      xxnew   = xxold + 0.5d0*(Axx + Axxi)*dt
      yynew   = yyold + 0.5d0*(Ayy + Ayyi)*dt
      zznew   = zzold + 0.5d0*(Azz + Azzi)*dt
      ppxnew  = ppxold + 0.5d0*(Apx + Apxi)*dt + stermsBx*dOmx !+ stermsCx*dOmxc
      ppynew  = ppyold + 0.5d0*(Apy + Apyi)*dt + stermsBy*dOmy !+ stermsCy*dOmyc
      ppznew  = ppzold + 0.5d0*(Apz + Apzi)*dt + stermsBz*dOmz !+ stermsCy*dOmyc
      if( mod(ii,save_freq) .eq. 0) then
        ll = ll + 1
        !xx2s(:,ll)  = xxnew*xxnew
        !yy2s(:,ll)  = yynew*yynew
        ppx2s(:,ll) = ppxnew*ppxnew
        ppy2s(:,ll) = ppynew*ppynew
        ppz2s(:,ll) = ppznew*ppznew
        xpxs(:,ll)  = xxnew*ppxnew
        ypys(:,ll)  = yynew*ppynew
        zpzs(:,ll)  = zznew*ppznew
      end if
      if( ii .gt. fin) then
          xxs(:,mm)   = xxnew
          yys(:,mm)   = yynew
          zzs(:,mm)   = zznew
          call local_energy(nparticles, alphay, alphaz, xxold, yyold, zzold, invD1, ppxold, ppyold, ppzold, energy)
          call heat_current(nparticles, fx1, fy1, fz1, ppxold, ppyold, ppzold, hc)
          call current_Flux(hc, energy, xxold, yyold, zzold, ppxold, ppyold, ppzold, nparticles, JJix, JJiy, JJiz)
          JJix_av = JJix_av + JJix/(nsteps-fin-1)
          JJiy_av = JJiy_av + JJiy/(nsteps-fin-1)
          JJiz_av = JJiz_av + JJiz/(nsteps-fin-1)
          JJix_av_v(kk) = JJix_av_v(kk) + JJix/(nsteps-fin-1)
          JJiy_av_v(kk) = JJiy_av_v(kk) + JJiy/(nsteps-fin-1)
          JJiz_av_v(kk) = JJiz_av_v(kk) + JJiz/(nsteps-fin-1)
          mm = mm + 1
      end if
      xxold   = xxnew
      yyold   = yynew
      zzold   = zznew
      ppxold  = ppxnew
      ppyold  = ppynew
      ppzold  = ppznew
    end do
    !xx2s(:,nssteps)  = xxnew*xxnew
    !yy2s(:,nssteps)  = yynew*yynew
    ppx2s(:,nssteps) = ppxnew*ppxnew
    ppy2s(:,nssteps) = ppynew*ppynew
    ppz2s(:,nssteps) = ppznew*ppznew
    xpxs(:,nssteps)  = xxnew*ppxnew
    ypys(:,nssteps)  = yynew*ppynew
    zpzs(:,nssteps)  = zznew*ppznew
    if(rank .eq. 0 .and. kk .eq. 1) then
      open(unit=11, file="posX.dat")
      open(unit=12, file="posY.dat")
      open(unit=13, file="posZ.dat")
      print*, "Printing steady-state spatial distribution of "
      do jj=1, nparticles, 1
        write(11,*) xxs(jj,:)
        write(12,*) yys(jj,:)
        write(13,*) zzs(jj,:)
      end do
      close(unit=11)
      close(unit=12)
      close(unit=13)
    end if
    !xx2_av  = (xx2_av + xx2s)
    !yy2_av  = (yy2_av + yy2s)
    ppx2_av = (ppx2_av + ppx2s)
    ppy2_av = (ppy2_av + ppy2s)
    ppz2_av = (ppz2_av + ppz2s)
    xpx_av  = (xpx_av + xpxs)
    ypy_av  = (ypy_av + ypys)
    zpz_av  = (zpz_av + zpzs)
    xx_av =   (xx_av + xxs)
    if( ( mod(kk,5) .eq. 0) .and. (kk .lt. local_traj) ) then
     print*, "Writing PARTIAL results to files after ", kk, "trajectories."
     errJJix = 0.0d0
     errJJiy = 0.0d0
     errJJiz = 0.0d0
     do nn=1, kk, 1
       errJJix = errJJix + (JJix_av_v(nn)-sum(JJix_av_v,1)/(kk))*(JJix_av_v(nn)-sum(JJix_av_v,1)/(kk))
       errJJiy = errJJiy + (JJiy_av_v(nn)-sum(JJiy_av_v,1)/(kk))*(JJiy_av_v(nn)-sum(JJiy_av_v,1)/(kk))
       errJJiz = errJJiz + (JJiz_av_v(nn)-sum(JJiz_av_v,1)/(kk))*(JJiz_av_v(nn)-sum(JJiz_av_v,1)/(kk))
     end do
     errJJix = errJJix/kk
     errJJiy = errJJiy/kk
     errJJiz = errJJiz/kk


     call mpi_reduce(JJix_av, JJix_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(JJiy_av, JJiy_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(JJiz_av, JJiz_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     !call mpi_reduce(xx2_av, xx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     !call mpi_reduce(yy2_av, yy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
!     call mpi_reduce(xx_av, xx_avt, (nsteps-fin), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ppx2_av, ppx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ppy2_av, ppy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ppz2_av, ppz2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(xpx_av, xpx_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(ypy_av, ypy_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(zpz_av, zpz_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(errJJix, errJJix_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(errJJiy, errJJiy_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(errJJiz, errJJiz_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
     call mpi_reduce(kk, traj, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)
     print*, "Finished writing up to", kk
     if(rank .eq. 0) then
      !xx2_avt  = xx2_avt*char_length*char_length/traj
      !yy2_avt  = yy2_avt*char_length*char_length/traj
      ppx2_avt = ppx2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
      ppy2_avt = ppy2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
      ppz2_avt = ppz2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
      xpx_avt  = xpx_avt*char_length*char_length*mass*long_freq/traj
      ypy_avt  = ypy_avt*char_length*char_length*mass*long_freq/traj
      zpz_avt  = zpz_avt*char_length*char_length*mass*long_freq/traj
      errJJix_t = sqrt( errJJix_t/procs)
      errJJiy_t = sqrt( errJJiy_t/procs)
      errJJiz_t = sqrt( errJJiz_t/procs)
      open(unit=11, file="heatflux.dat")
      open(unit=12, file="temperatures.dat")
      write(11,*) JJix_avt/traj, "+/-", errJJix_t
      write(11,*) JJiy_avt/traj, "+/-", errJJiy_t
      write(11,*) JJiz_avt/traj, "+/-", errJJiz_t
      do jj=1, nparticles
       write(12,*) ppx2_avt(jj,:) + ppy2_avt(jj,:) + ppz2_avt(jj,:)
      end do
      close(unit=11)
      close(unit=12)
!      open(unit=13, file="posXav.dat")
!     do jj=1, nparticles, 1
!        write(13,*) xx_av/traj
!      end do
!      close(unit=13)
     end if
     !xx2_avt  = 0.0d0
     !yy2_avt  = 0.0d0
     ppx2_avt = 0.0d0
     ppy2_avt = 0.0d0
     ppz2_avt = 0.0d0
     xpx_avt  = 0.0d0
     ypy_avt  = 0.000
     zpz_avt  = 0.0d0
    end if
  end do
  print*,"Proc ", rank, " finished integrating"
  errJJix = 0.0d0
  errJJiy = 0.0d0
  errJJiz = 0.0d0
  do nn=1, local_traj, 1
    errJJix = errJJix + (JJix_av_v(nn)-sum(JJix_av_v,1)/(local_traj))*(JJix_av_v(nn)-sum(JJix_av_v,1)/(local_traj))
    errJJiy = errJJiy + (JJiy_av_v(nn)-sum(JJiy_av_v,1)/(local_traj))*(JJiy_av_v(nn)-sum(JJiy_av_v,1)/(local_traj))
    errJJiz = errJJiz + (JJiz_av_v(nn)-sum(JJiz_av_v,1)/(local_traj))*(JJiz_av_v(nn)-sum(JJiz_av_v,1)/(local_traj))
  end do
  errJJix = errJJix/local_traj
  errJJiy = errJJiy/local_traj
  errJJiz = errJJiz/local_traj

  call mpi_barrier(mpi_comm_world, ierr)
  call mpi_reduce(JJix_av, JJix_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(JJiy_av, JJiy_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(JJiz_av, JJiz_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(JJix_av, JJix_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(JJiy_av, JJiy_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(JJiz_av, JJiz_avt, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  !call mpi_reduce(xx2_av, xx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  !call mpi_reduce(yy2_av, yy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
!  call mpi_reduce(xx_av, xx_avt, (nsteps-fin), mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ppx2_av, ppx2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ppy2_av, ppy2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ppz2_av, ppz2_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(xpx_av, xpx_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(ypy_av, ypy_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(zpz_av, zpz_avt, n_elems, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(errJJix, errJJix_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(errJJiy, errJJiy_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(errJJiz, errJJiz_t, 1, mpi_double_precision, mpi_sum, 0, mpi_comm_world, ierr)
  call mpi_reduce(local_traj, traj, 1, mpi_integer, mpi_sum, 0, mpi_comm_world, ierr)


  call mpi_barrier(mpi_comm_world, ierr)

  if(rank .eq. 0) then
    seconds = mpi_wtime() - seconds
    print*, "total integration + partial writing time:", seconds, seconds/3600.0d0
    print*,traj

    errJJix_t = sqrt( errJJix_t/procs)
    errJJiy_t = sqrt( errJJiy_t/procs)
    errJJiz_t = sqrt( errJJiz_t/procs)
    !xx2_avt  = xx2_avt*char_length*char_length/traj
    !yy2_avt  = yy2_avt*char_length*char_length/traj
    ppx2_avt = ppx2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
    ppy2_avt = ppy2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
    ppz2_avt = ppz2_avt*char_length*char_length*mass*long_freq*long_freq/(2.0d0*kb)/traj ! Convert to temperature in mK
    xpx_avt  = xpx_avt*char_length*char_length*mass*long_freq/traj
    ypy_avt  = ypy_avt*char_length*char_length*mass*long_freq/traj
    zpz_avt  = zpz_avt*char_length*char_length*mass*long_freq/traj
    errJJix_t = sqrt( errJJix_t/traj - (JJix_avt*JJix_avt)/(traj*traj) )
    errJJiy_t = sqrt( errJJiy_t/traj - (JJiy_avt*JJiy_avt)/(traj*traj) )
    errJJiz_t = sqrt( errJJiz_t/traj - (JJiz_avt*JJiz_avt)/(traj*traj) )

    open(unit=11, file="heatflux.dat")
    open(unit=12, file="temperatures.dat")
    write(11,*) JJix_avt/traj, "+/-", errJJix_t
    write(11,*) JJiy_avt/traj, "+/-", errJJiy_t
    write(11,*) JJiz_avt/traj, "+/-", errJJiz_t
    do jj=1, nparticles
     write(12,*) ppx2_avt(jj,:) + ppy2_avt(jj,:) + ppz2_avt(jj,:)
    end do
    close(unit=11)
    close(unit=12)
 !   do jj=1, nparticles, 1
  !    write(13,*) xx_av/traj
   ! end do
   ! close(unit=13)
    print*, "Writing FINAL results to files after ", traj , "trajectories."
    seconds1 = mpi_wtime() - seconds1
    print*, "integration + message + write time:", seconds1
  end if

  call mpi_finalize(ierr)


end program threeDChain
