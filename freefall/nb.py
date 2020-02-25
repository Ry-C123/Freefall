import sys
import math
import multiprocessing
import time

import numpy as np
import matplotlib.pyplot as plt

from CONFIG import *
from integrator import basic_int, RK4_conv, RK4, YOSHI


time1 = time.time()
print('')
print('                  Stellar Parameters                 ')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Stellar Mass = '+str(M/Sm)+' Solar Masses')
print('Stellar Radius = '+str(R/Sr)+' Solar Radii') # alice doesnt recognise \odot  
print('Magnetic Field = '+str(B)+' Guass')
print('Field inclination = '+str(inc)+'degrees') # alice doesnt regonise degree symbol
print('Stellar Temperature = '+str(TEMP)+' Kelvin')
print('Star spin = '+str(omega)+' UNITS?')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('')


if integ == 'basic':
    func = basic_int
elif integ == 'yoshi':
    func = YOSHI
elif integ == 'RK4':
    func = RK4
elif integ == 'RK4_cv':
    func = RK4_conv
elif integ == 'BS':
    print('in prep')
    sys.exit()
else:
    raise ValueError(integ+" is not a valid integrator. Use either 'basic', 'yoshi', or 'RK4'")


particles = open('test_part.conf', 'r').readlines()
IDs = []
x = []
y = []
z =[] 
vx = []
vy = []
vz = [] 
m_p = []
r_p = []

for i in range(len(particles)):
    if '#' in particles[i]:
        continue

    #ID,X,Y,VX,VY,m,r    <- Param file format
    tmp = particles[i].split(',')
    IDs.append(int(tmp[0]))
    x.append(float(tmp[1]))
    y.append(float(tmp[2]))
    z.append(float(tmp[3]))
    vx.append(float(tmp[3]))
    vy.append(float(tmp[5]))
    vz.append(float(tmp[6]))
    m_p.append(float(tmp[7]))
    r_p.append(float(tmp[8]))
    #D2.append(math.sqrt(float(tmp[1])**2 + float(tmp[2])**2)) #Initial Distance

#x_in = max(x)
TOT_PARTS = len(x)

if Cores == 'max':
    Cores = multiprocessing.cpu_count()
else:
    Cores = int(Cores)
print('Using '+str(Cores)+' Core(s)') 

i = 0
if Cores < 1:
    raise ValueError("Cores < 1: You need to use at least 1 core!")  
elif Cores == 1:
    while i < n_steps:

        for pn in range(len(x)):
            IDs[pn],x[pn],y[pn],z[pn],vx[pn],vy[pn],vz[pn] = func(IDs[pn], x[pn],y[pn],z[pn],vx[pn],vy[pn],vz[pn], m_p[pn], r_p[pn],dt, M, R, B, TEMP)


        if (x[pn]**2 + y[pn]**2 + z[pn]**2) < (ACC_RAD)**2:
            print('Collision: particle '+str(IDs[pn])+' removed | time = '+str(dt*(i+1)*3.17098e-8)+' years')
            x.pop(pn)
            y.pop(pn)
            z.pop(pn)
            vx.pop(pn)
            vy.pop(pn)
            vz.pop(pn)
            r_p.pop(pn)
            m_p.pop(pn)
            IDs.pop(pn)
        elif (x[pn]**2 + y[pn]**2 + z[pn]**2) > (EJE_RAD)**2:
            print('Ejected: particle '+str(IDs[pn])+' removed | time = '+str(dt*(i+1)*3.17098e-8)+' years')
            x.pop(pn)
            y.pop(pn)
            z.pop(pn)
            vx.pop(pn)
            vy.pop(pn)
            vz.pop(pn)
            r_p.pop(pn)
            m_p.pop(pn)
            IDs.pop(pn)
        if len(x) == 0:
            print('No particles left!')
            sys.exit()

        if (i+1)%OUTPUT_int == 0:
            print(str(dt*(i+1)*3.17098e-8)+' years')
            D2 = round(math.sqrt(x[0]**2 + y[0]**2 + z[0]**2))
            if write_files == True:
                FILE = open(str(dt*i)+'.txt','w')
                for L in range(len(x)):
                    FILE.write(str(x[L])+','+str(y[L])+','+str(z[L])+','+str(vx[L])+','+str(vy[L])+','+str(vz[L])+','+str(IDs[L])+','+str(r_p[L])+','+str(m_p[L])+'\n')

            if PLOT_ON is True:
                plt.xlim(-2*Sr, 2*Sr)
                plt.ylim(-2*Sr,2*Sr)
                rot_s = dt*i*omega*57.2958 + 270
                plt.plot(0, 0, marker=(3, 0, rot_s), markersize=8, linestyle='None', c='k')
                TIME_YEARS = round((i*dt*3.17098e-8), 4)
                plt.text(4.0*AU/7, 6.0*AU/7,str(TIME_YEARS)+' years')
                plt.scatter(x,y,c=IDs, s=4, cmap='Blues')
                plt.text(4.0*AU/7, 5.0*AU/7,str(D2/Sr)+' R$\odot$')
                plt.clim(0, TOT_PARTS)
                plt.pause(0.05)
                plt.clf()
        i += 1
else:
    while i < n_steps:
        STAR_PARAMS = [dt, M, R, B, TEMP]
        with multiprocessing.Pool(Cores) as p:
            OUT = p.starmap(func, [([IDs[i],x[i],y[i],z[i],vx[i],vy[i],vz[i], m_p[i],r_p[i]]+STAR_PARAMS) for i in range(len(x))])
        p.close()
        for part in OUT:
            idx = part[0]-1
            x[idx] = part[1]
            y[idx] = part[2]
            z[idx] = part[3]
            vx[idx] = part[4]
            vy[idx] = part[5]
            vz[idx] = part[6]

            if (x[idx]**2 + y[idx]**2 + z[idx]**2) < (ACC_RAD)**2:
                print('Collision: particle '+str(IDs[idx])+' removed | time = '+str(dt*(i+1)*3.17098e-8)+' years')
                x.pop(idx)
                y.pop(idx)
                z.pop(idx)
                vx.pop(idx)
                vy.pop(idx)
                vz.pop(idx)
                r_p.pop(idx)
                m_p.pop(idx)
                IDs.pop(idx)
            elif (x[idx]**2 + y[idx]**2 + z[idx]**2) > (EJE_RAD)**2:
                print('Ejected: particle '+str(IDs[idx])+' removed | time = '+str(dt*(i+1)*3.17098e-8)+' years')
                x.pop(idx)
                y.pop(idx)
                z.pop(idx)
                vx.pop(idx)
                vy.pop(idx)
                vz.pop(idx)
                r_p.pop(idx)
                m_p.pop(idx)
                IDs.pop(idx)    

        if len(x) == 0:
            print('No particles left!')
            sys.exit()

        if (i+1)%OUTPUT_int == 0:
            print(str(dt*(i+1)*3.17098e-8)+' years')
            if write_files == True:
                FILE = open(str(int(dt*i))+'.txt','w')
                for L in range(len(x)):
                    FIlE.write(str(x[L])+','+str(y[L])+','+str(z[L])+','+str(vx[L])+','+str(vy[L])+','+str(vz[L])+','+str(IDs[L])+','+str(r_p[L])+','+str(m_p[L])+'\n')
            if PLOT_ON is True:
                plt.xlim(-2*AU, 0.5*AU)
                plt.ylim(-x_in*1.3, x_in*1.3)
                rot_s = dt*i*omega*57.2958 + 270
                plt.plot(0, 0, marker=(3, 0, rot_s), markersize=8, linestyle='None', c='k')
                plt.scatter(x,y,c=IDs, s=8)
                plt.clim(0, TOT_PARTS)
                plt.pause(0.05)
                plt.clf()
        i+=1
#print(time.time() - time1)



