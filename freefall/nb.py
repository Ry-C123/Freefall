import sys
import math
import multiprocessing
import time

import numpy as np
import matplotlib.pyplot as plt

from CONFIG import *
from integrator import basic_int, RK4, YOSHI


time1 = time.time()
print('')
print('                  Stellar Parameters                 ')
print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
print('Stellar Mass = '+str(M/Sm)+' M⊙')
print('Stellar Radius = '+str(R/Sr)+' R⊙') 
print('Magnetic Field = '+str(B)+' Guass')
print('Field inclination = '+str(inc)+'°')
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
elif integ == 'BS':
    func = None
    print('BS not ready!')
    sys.exit()
else:
    raise ValueError(integ+" is not a valid integrator. Use either 'basic', 'yoshi', or 'RK4'")


particles = open('test_part.conf', 'r').readlines()
IDs = []
x = []
y = []
#z =[] TODO
vx = []
vy = []
#vz = [] TODO
m_p = []
r_p = []
D2 = []
for i in range(len(particles)):
    if '#' in particles[i]:
        continue

    #ID,X,Y,VX,VY,m,r    <- Param file format
    tmp = particles[i].split(',')
    IDs.append(int(tmp[0]))
    x.append(float(tmp[1]))
    y.append(float(tmp[2]))
    vx.append(float(tmp[3]))
    vy.append(float(tmp[4]))
    m_p.append(float(tmp[5]))
    r_p.append(float(tmp[6]))
    D2.append(math.sqrt(float(tmp[1])**2 + float(tmp[2])**2)) #Initial Distance

x_in = max(x)
TOT_PARTS = len(x)

if Cores == 'max':
    Cores = multiprocessing.cpu_count()
else:
    Cores = int(Cores)
print('Using '+str(Cores)+' Core(s)') 


if Cores < 1:
    raise ValueError("Cores < 1: You need to use at least 1 core!")  
elif Cores == 1:
    for i in range(n_steps):

        for pn in range(len(x)):
            IDs[pn],x[pn],y[pn],vx[pn],vy[pn] = func(IDs[pn], x[pn],y[pn],vx[pn],vy[pn], m_p[pn], r_p[pn],dt, M, R, B, TEMP)



        if (x[pn]**2 + y[pn]**2) < (ACC_RAD)**2:
            print('Collision: particle '+str(IDs[pn])+' removed | time = '+str(dt*(i+1)))
            x.pop(pn)
            y.pop(pn)
            vx.pop(pn)
            vy.pop(pn)
            r_p.pop(pn)
            m_p.pop(pn)
            IDs.pop(pn)
        elif (x[pn]**2 + y[pn]**2) > (EJE_RAD)**2:
            print('Ejected: particle '+str(IDs[pn])+' removed | time = '+str(dt*(i+1)))
            x.pop(pn)
            y.pop(pn)
            vx.pop(pn)
            vy.pop(pn)
            r_p.pop(pn)
            m_p.pop(pn)
            IDs.pop(pn)
        if len(x) == 0:
            print('No particles left!')
            sys.exit()

        if (i+1)%1000 == 0:
            #print(D2 - math.sqrt(x_in**2 + y_in**2))
            #D2 = math.sqrt(x**2 + y**2)

            if PLOT_ON is True:
                plt.xlim(-x_in*1.3, x_in*1.3)
                plt.ylim(-x_in*1.3, x_in*1.3)
                rot_s = dt*i*omega*57.2958 + 270
                plt.plot(0, 0, marker=(3, 0, rot_s), markersize=8, linestyle='None', c='k')
                plt.scatter(x,y,c=IDs, s=1)
                plt.clim(0, TOT_PARTS)
                plt.pause(0.05)
                plt.clf()

else:
    for i in range(n_steps):
        STAR_PARAMS = [dt, M, R, B, TEMP]
        with multiprocessing.Pool(Cores) as p:
            OUT = p.starmap(func, [([IDs[i],x[i],y[i],vx[i],vy[i], m_p[i],r_p[i]]+STAR_PARAMS) for i in range(len(x))])
        for part in OUT:
            idx = part[0]-1
            x[idx] = part[1]
            y[idx] = part[2]
            vx[idx] = part[3]
            vy[idx] = part[4]

            if (x[idx]**2 + y[idx]**2) < (ACC_RAD)**2:
                print('Collision: particle '+str(IDs[idx])+' removed | time = '+str(dt*(i+1)))
                x.pop(idx)
                y.pop(idx)
                vx.pop(idx)
                vy.pop(idx)
                r_p.pop(idx)
                m_p.pop(idx)
                IDs.pop(idx)
            elif (x[idx]**2 + y[idx]**2) > (EJE_RAD)**2:
                print('Ejected: particle '+str(IDs[idx])+' removed | time = '+str(dt*(i+1)))
                x.pop(idx)
                y.pop(idx)
                vx.pop(idx)
                vy.pop(idx)
                r_p.pop(idx)
                m_p.pop(idx)
                IDs.pop(idx)    

        if len(x) == 0:
            print('No particles left!')
            sys.exit()

        if (i+1)%1000 == 0:
            #print(D2 - math.sqrt(x_in**2 + y_in**2))
            #D2 = math.sqrt(x**2 + y**2)

            if PLOT_ON is True:
                plt.xlim(-x_in*1.3, x_in*1.3)
                plt.ylim(-x_in*1.3, x_in*1.3)
                rot_s = dt*i*omega*57.2958 + 270
                plt.plot(0, 0, marker=(3, 0, rot_s), markersize=8, linestyle='None', c='k')
                plt.scatter(x,y,c=IDs, s=1)
                plt.clim(0, TOT_PARTS)
                plt.pause(0.05)
                plt.clf()

print(time.time() - time1)


