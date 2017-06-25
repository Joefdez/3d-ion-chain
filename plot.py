# !usr/bin/python

from numpy import *
from matplotlib.pylab import *
from mpl_toolkits.mplot3d import Axes3D

ion()

#xx = loadtxt("results/posx.dat")
#yy = loadtxt("results/posy.dat")
print("Loading ion coordinate files.")
xx = loadtxt("posX.dat")
yy = loadtxt("posY.dat")
zz = loadtxt("posZ.dat")

dims = shape(xx)
frac = int(dims[1])
xxs = zeros(30)
yys = zeros(30)
zzs = zeros(30)
xxf = zeros([frac*30])
yyf = zeros([frac*30])
zzf = zeros([frac*30])

print("Flattenning position arrays to compose bivariate histogram.")
for ii in range(0,30):
    xxs[ii] = mean(xx[ii,:])
    yys[ii] = mean(yy[ii,:])
    zzs[ii] = mean(zz[ii,:])


del xx
del yy
del zz

chain = figure("chain")
axes = chain.add_subplot(111, projection="3d")
axes.plot(xxs, yys, zzs, '.', markersize=13)

xxs.sort()
print("Loading ion temperature profiles.")
temps=loadtxt("temperatures.dat")
print("Calculating temperature profile")
tf = zeros(30)
ns = shape(temps)
ns = 0.2*ns[1]
for ii in range(30):
    tf[ii] = mean(temps[ii,-ns:])/3.
tf = 1000*tf

print("Plotting temperature profile.")
delT = max(tf)-min(tf)
print("Temperature delta")
print(delT)
print("Temperature grad")
print(delT/(max(xxs)-min(xxs)))
tb, th = (min(tf)-0.1*delT), (max(tf)+0.1*delT)  #(min(tf)-0.1*delT), (max(tf)+0.1*delT)
txl, txr = -6.0, 6.0
temp = figure("Temperature")
axt = temp.add_subplot(111)
axt.set_xlim([txl, txr])
axt.set_ylim([tb, th])
axt.set_xlabel(r'$x$ ($scaled\, units$)', fontsize=22)
axt.set_ylabel(r'$T$ ($mK$)', fontsize=22)
axt.tick_params(labelsize=11)
axt.axvspan(-6., 0.5*(xxs[1]+xxs[2])), facecolor='r', alpha=0.75)
axt.axvspan(0.5*(xxs[26]+xxs[27]), 6., facecolor='c', alpha=0.75)
axt.plot(xxs, tf, 'b')
axt.plot(xxs, tf, 'b.', markersize=11)
axt = temp.add_subplot(111)
axt.set_xlim([txl, txr])
axt.set_ylim([tb, th])
axt.plot(xxs, tf) # Plot temperature in mK

"""

print("Making plot with inset")

ins, inax1 = subplots()

left, bottom, width, height = [2 , 5, 5.5, 8]
inax2 = ins.axes([left, bottom, width, height])

inax1.set_xlim([txl, txr])
inax1.set_ylim([3, 9])
inax1.set_xlabel(r'x (scaled units)', fontsize=15)
inax1.set_ylabel(r'T (mK)', fontsize=15)
inax1.tick_params(labelsize=11)
inax1.axvspan(-6., -4., facecolor='r', alpha=0.75)
inax1.axvspan(4., 6., facecolor='c', alpha=0.75)
inax1.plot(xxs, tf, 'b')
inax1.plot(xxs, tf, 'b.', markersize=11)
inax2.imshow(H, origin='low',
            extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])


"""
