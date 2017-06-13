allocate(xx0(1:nparticles))
xx0 = 0.0d0
allocate(yy0(1:nparticles))
yy0 = 0.0d0
allocate(zz0(1:nparticles))
zz0 = 0.0d0
allocate(px0(1:nparticles))
px0 = 0.0d0
allocate(py0(1:nparticles))
py0 = 0.0d0
allocate(pz0(1:nparticles))
pz0 = 0.0d0
allocate(xxold(1:nparticles))
xxold = 0.0d0
allocate(yyold(1:nparticles))
yyold = 0.0d0
allocate(zzold(1:nparticles))
zzold = 0.0d0
allocate(ppxold(1:nparticles))
ppxold = 0.0d0
allocate(ppyold(1:nparticles))
ppyold = 0.0d0
allocate(ppzold(1:nparticles))
ppzold = 0.0d0
allocate(xxnew(1:nparticles))
xxnew = 0.0d0
allocate(yynew(1:nparticles))
yynew = 0.0d0
allocate(zznew(1:nparticles))
zznew = 0.0d0
allocate(ppxnew(1:nparticles))
ppxnew = 0.0d0
allocate(ppynew(1:nparticles))
ppynew = 0.0d0
allocate(ppznew(1:nparticles))
ppznew = 0.0d0
allocate(xxs(1:nparticles, 1:(nsteps-fin)))
xxs = 0.0d0
allocate(yys(1:nparticles, 1:(nsteps-fin)))
yys = 0.0d0
allocate(zzs(1:nparticles, 1:(nsteps-fin)))
!allocate(ppxs(1:nparticles, 1:(nsteps-fin)))
!ppxs = 0.0d0
!allocate(ppys(1:nparticles, 1:(nsteps-fin)))
!ppys = 0.0d0
!allocate(JJix_s(1:(nsteps-fin)))
!JJix_s = 0.0d0
!allocate(JJiy_s(1:(nsteps-fin)))
!JJiy_s = 0.0d0
!allocate(xx2s(1:nparticles, 1:nssteps))
!xx2s = 0.0d0
!allocate(yy2s(1:nparticles, 1:nssteps))
!yy2s = 0.0d0
allocate(ppx2s(1:nparticles, 1:nssteps))
ppx2s = 0.0d0
allocate(ppy2s(1:nparticles, 1:nssteps))
ppy2s = 0.0d0
allocate(ppz2s(1:nparticles, 1:nssteps))
ppz2s = 0.0d0
allocate(xpxs(1:nparticles, 1:nssteps))
xpxs = 0.0d0
allocate(ypys(1:nparticles, 1:nssteps))
ypys = 0.0d0
allocate(zpzs(1:nparticles, 1:nssteps))
zpzs = 0.0d0
allocate(xxi(1:nparticles))
xxi = 0.0d0
allocate(yyi(1:nparticles))
yyi = 0.0d0
allocate(zzi(1:nparticles))
zzi = 0.0d0
allocate(ppxi(1:nparticles))
ppxi = 0.0d0
allocate(ppyi(1:nparticles))
ppyi = 0.0d0
allocate(ppzi(1:nparticles))
allocate(fx1(1:nparticles,1:nparticles))
fx1 = 0.0d0
allocate(fx2(1:nparticles,1:nparticles))
fx2 = 0.0d0
allocate(fy1(1:nparticles,1:nparticles))
fy1 = 0.0d0
allocate(fy2(1:nparticles,1:nparticles))
fy2 = 0.0d0
allocate(fz1(1:nparticles,1:nparticles))
fz1 = 0.0d0
allocate(fz2(1:nparticles,1:nparticles))
fz2 = 0.0d0
allocate(invD1(1:nparticles,1:nparticles))
invD1 = 0.0d0
allocate(invD2(1:nparticles,1:nparticles))
invD2 = 0.0d0
allocate(fx(1:nparticles))
fx = 0.0d0
allocate(fy(1:nparticles))
fy = 0.0d0
allocate(fz(1:nparticles))
fz = 0.0d0
allocate(Axx(1:nparticles))
Axx = 0.0d0
allocate(Ayy(1:nparticles))
Ayy = 0.0d0
allocate(Azz(1:nparticles))
Azz = 0.0d0
allocate(Apx(1:nparticles))
Apx = 0.0d0
allocate(Apy(1:nparticles))
Apy = 0.0d0
allocate(Apz(1:nparticles))
Apz = 0.0d0
allocate(Axxi(1:nparticles))
Axxi = 0.0d0
allocate(Ayyi(1:nparticles))
Ayyi = 0.0d0
allocate(Azzi(1:nparticles))
Azzi = 0.0d0
allocate(Apxi(1:nparticles))
Apxi = 0.0d0
allocate(Apyi(1:nparticles))
Apyi = 0.0d0
allocate(Apzi(1:nparticles))
Apzi = 0.0d0
allocate(dOmx(1:nparticles))
dOmx = 0.0d0
allocate(dOmy(1:nparticles))
dOmy = 0.0d0
allocate(dOmz(1:nparticles))
dOmz = 0.0d0
allocate(dOmxc(1:nparticles))
dOmxc = 0.0d0
allocate(dOmyc(1:nparticles))
dOmyc = 0.0d0
allocate(stermsBx(1:nparticles))
stermsBx = 0.0d0
allocate(stermsBy(1:nparticles))
stermsBy = 0.0d0
allocate(stermsBz(1:nparticles))
stermsBz = 0.0d0
allocate(stermsCx(1:nparticles))
stermsCx = 0.0d0
allocate(stermsCy(1:nparticles))
stermsCy = 0.0d0
allocate(stermsCz(1:nparticles))
stermsCz = 0.0d0


!allocate(xx2_av(1:nparticles, 1:nssteps))
!xx2_av = 0.0d0
!allocate(yy2_av(1:nparticles, 1:nssteps))
!yy2_av = 0.0d0
allocate(ppx2_av(1:nparticles, 1:nssteps))
ppx2_av = 0.0d0
allocate(ppy2_av(1:nparticles, 1:nssteps))
ppy2_av = 0.0d0
allocate(ppz2_av(1:nparticles, 1:nssteps))
ppz2_av = 0.0d0
allocate(xpx_av(1:nparticles, 1:nssteps))
xpx_av = 0.0d0
allocate(ypy_av(1:nparticles, 1:nssteps))
ypy_av = 0.0d0
allocate(zpz_av(1:nparticles, 1:nssteps))
zpz_av = 0.0d0

!allocate(xx2_avo(1:nparticles, 1:nssteps))
!xx2_avo = 0.0d0
!allocate(yy2_avo(1:nparticles, 1:nssteps))
!yy2_avo = 0.0d0
!allocate(ppx2_avo(1:nparticles, 1:nssteps))
!ppx2_avo = 0.0d0
!allocate(ppy2_avo(1:nparticles, 1:nssteps))
!ppy2_avo = 0.0d0
!allocate(xpx_avo(1:nparticles, 1:nssteps))
!xpx_avo = 0.0d0
!allocate(ypy_avo(1:nparticles, 1:nssteps))
!ypy_avo = 0.0d0


!allocate(xx2_avt(1:nparticles, 1:nssteps))
!xx2_avt = 0.0d0
!allocate(yy2_avt(1:nparticles, 1:nssteps))
!yy2_avt = 0.0d0
allocate(ppx2_avt(1:nparticles, 1:nssteps))
ppx2_avt = 0.0d0
allocate(ppy2_avt(1:nparticles, 1:nssteps))
ppy2_avt = 0.0d0
allocate(ppz2_avt(1:nparticles, 1:nssteps))
ppz2_avt = 0.0d0
allocate(xpx_avt(1:nparticles, 1:nssteps))
xpx_avt = 0.0d0
allocate(ypy_avt(1:nparticles, 1:nssteps))
ypy_avt = 0.0d0
allocate(zpz_avt(1:nparticles, 1:nssteps))
zpz_avt = 0.0d0

allocate(energy(1:nparticles))
energy = 0.0d0
allocate(hc(1:nparticles,1:nparticles))
energy = 0.0d0
!allocate(hcx(1:(nssteps)))
!hcx = 0.0d0
!allocate(hcy(1:(nssteps)))
!hcy = 0.0d0
!allocate(hcx_av(1:(nssteps)))
!hcx_av = 0.0d0
!allocate(hcy_av(1:(nssteps)))
!hcy_av = 0.0d0
!allocate(hcx_avt(1:(nssteps)))
!hcx_avt = 0.0d0
!allocate(hcy_avt(1:(nssteps)))
!hcy_avt = 0.0d0

allocate(xx_av(1:nparticles,(nsteps-fin)))
allocate(yy_av(1:nparticles,(nsteps-fin)))
allocate(zz_av(1:nparticles,(nsteps-fin)))
allocate(xx_avt(1:nparticles,(nsteps-fin)))
allocate(yy_avt(1:nparticles,(nsteps-fin)))
allocate(zz_avt(1:nparticles,(nsteps-fin)))
xx_av = 0.0d0
yy_av = 0.0d0
zz_av = 0.0d0
xx_avt = 0.0d0
yy_avt = 0.0d0
zz_avt = 0.0d0

allocate(JJix_av_v(1:local_traj))
JJix_av_v = 0.0d0
allocate(JJiy_av_v(1:local_traj))
JJiy_av_v = 0.0d0
allocate(JJiz_av_v(1:local_traj))
JJiz_av_v = 0.0d0
! initialization

JJix = 0.0d0
JJiy = 0.0d0
JJiz = 0.0d0
JJix_av = 0.0d0
JJiy_av = 0.0d0
JJiz_avt = 0.0d0
JJix_avt = 0.0d0
JJiy_avt = 0.0d0
JJiy_avt = 0.0d0
