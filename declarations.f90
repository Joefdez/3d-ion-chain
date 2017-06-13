! variables
real(kind=8), dimension(:), allocatable               :: xx0, yy0, zz0, px0, py0, pz0
real(kind=8), dimension(:), allocatable               :: xxold, yyold, zzold, ppxold, ppyold, ppzold
real(kind=8), dimension(:), allocatable               :: xxnew, yynew, zznew, ppxnew, ppynew, ppznew
real(kind=8), dimension(:,:), allocatable             :: xxs, yys, zzs, ppxs, ppys, ppzs
real(kind=8), dimension(:,:), allocatable             :: xx2s, yy2s, zz2s, ppx2s, ppy2s, ppz2s
real(kind=8), dimension(:,:), allocatable             :: xpxs, ypys, zpzs
real(kind=8), dimension(:,:), allocatable             :: xx_avt, yy_avt, zz_avt, ppx_avt, ppy_avt, ppz_avt
real(kind=8), dimension(:,:), allocatable             :: xx2_avt, yy2_avt, zz2_avt, ppx2_avt, ppy2_avt, ppz2_avt
real(kind=8), dimension(:,:), allocatable             :: xpx_avt, ypy_avt, zpz_avt
real(kind=8), dimension(:,:), allocatable             :: xx_av, yy_av, zz_av, ppx_av, ppy_av, ppz_av
real(kind=8), dimension(:,:), allocatable             :: xx2_av, yy2_av, zz2_av, ppx2_av, ppy2_av, ppz2_av
real(kind=8), dimension(:,:), allocatable             :: xx_avo, yy_avo, zz_avo, ppx_avo, ppy_avo, ppz_avo
real(kind=8), dimension(:,:), allocatable             :: xx2_avo, yy2_avo, zz2_avo, ppx2_avo, ppy2_avo, ppz2_avo
real(kind=8), dimension(:,:), allocatable             :: xPx_av, yPy_av, zpz_av
real(kind=8), dimension(:,:), allocatable             :: xPx_avo, yPy_avo, zpz_avo
real(kind=8), dimension(:), allocatable               :: xxi, yyi, ppxi, ppyi, ppzi
real(kind=8), dimension(:,:), allocatable             :: fx1, fy1, fz1, fx2, fy2, fz2
real(kind=8), dimension(:), allocatable               :: fx, fy, fz
real(kind=8), dimension(:), allocatable               :: Axx, Axxi, Ayy, Ayyi, Azz, Azzi
real(kind=8), dimension(:), allocatable               :: Apx, Apxi, Apy, Apyi, Apz, Apzi
real(kind=8), dimension(:), allocatable               :: dOmx, dOmy, dOmz, dOmxc, dOmyc, dOmzc
real(kind=8), dimension(:), allocatable               :: stermsBx, stermsCx, stermsBy, stermsCy, stermsBz, stermsCz

real(kind=8)                                          :: tt, dt, dst, mass, charge, dist
real(kind=8)                                          :: alpha, char_length, long_freq
real(kind=8)                                          :: del1, del2, delC
real(kind=8)                                          :: omega0, omega1, omega2, omegaC
real(kind=8)                                          :: I1, I2, IC
real(kind=8)                                          :: k1, k2, kC
real(kind=8)                                          :: ic_radius
real(kind=8)                                          :: Gam
real(kind=8)                                          :: eta1, aeta1, eta2, aeta2
real(kind=8)                                          :: D1, aD1, D2, aD2
real(kind=8)                                          :: etaC, DC
real(kind=8)                                          :: aetaC, aDC
integer                                               :: nsteps, savefreq, nssteps, nparticles, nbath, n_elems, fin
integer                                               :: traj, local_traj, save_freq, rem
integer                                               :: ii,jj, kk, ll, mm, nn
real(kind=8)                                          :: seconds, seconds1
real(kind=8), dimension(:), allocatable               :: energy
real(kind=8)                                          :: JJix, JJiy, JJiz
real(kind=8), dimension(:), allocatable               :: JJix_s, JJiy_s, JJiz_s
real(kind=8), dimension(:), allocatable               :: hcx, hcy, hcz, hcx_av, hcy_av, hcz_av, hcx_avt, hcy_avt, hcz_avt
real(kind=8), dimension(:,:), allocatable             :: hc
real(kind=8), dimension(:,:), allocatable             :: invD1, invD2
real(kind=8)                                          :: JJix_av, JJiy_av, JJiz_avt, JJix_avt, JJiy_avt, JJiz_avt
real(kind=8), dimension(:), allocatable               :: JJix_av_v, JJiy_av_v, JJiz_av_v
real(kind=8)                                          :: errJJix, errJJiy, errJJiz, errJJix_t, errJJiy_t, errJJiz_t

real(kind=8)                                          :: initT, initSpeed
