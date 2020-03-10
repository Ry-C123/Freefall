import sys
import math
import multiprocessing
import time
import subprocess
import glob
import os


import numpy as np
import matplotlib.pyplot as plt

from CONFIG import *
from integrator import basic_int, RK4_conv, RK4, YOSHI


if __name__=='__main__':
    print('')
    print('-               Simulation: '+runname+'               -')
    print('')
    print('                  Stellar Parameters                 ')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('Stellar Mass = '+str(M/Sm)+' Solar Masses')
    print('Stellar Radius = '+str(R/Sr)+' Solar Radii') # alice doesnt recognise \odot  
    print('Magnetic Field = '+str(B)+' Guass')
    print('Field inclination = '+str(inc)+'degrees') # alice doesnt regonise degree symbol
    print('Stellar Temperature = '+str(TEMP)+' Kelvin')
    print('Star spin = '+str(omega)+' radians per second')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('')

    ###  Which integrator  ####
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
    ############################

    if restart == 0: 
        particles = open('test_part.conf', 'r').readlines()
    else:
        print('')
        print('resuming from '+runname+'_'+str(restart)+'.simo\n')
        particles = open(runname+'_'+str(restart)+'.simo', 'r').readlines()   
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

        #ID,X,Y,Z,VX,VY,VZ,m,r    <- Param file format
        tmp = particles[i].split(',')
        IDs.append(int(tmp[0]))
        x.append(float(tmp[1]))
        y.append(float(tmp[2]))
        z.append(float(tmp[3]))
        vx.append(float(tmp[4]))
        vy.append(float(tmp[5]))
        vz.append(float(tmp[6]))
        m_p.append(float(tmp[8]))
        r_p.append(float(tmp[7]))

    #x_in = max(x)
    TOT_PARTS = len(x)

    if Cores == 'max':
        Cores = multiprocessing.cpu_count()
    else:
        Cores = int(Cores)
    print('Using '+str(Cores)+' Core(s)') 

    i = 0
    if Cores < 1:
        raise ValueError("Cores less than 1: You need to use at least 1 core!")  


    elif Cores == 1 or TOT_PARTS == 1: #Use serial programming if there's only one particle
        print('Running in series')
        HK = glob.glob(runname+'*.simo') ### Housekeeping
        if HK != [] and OVERWRITE != True:
            raise ValueError("Simulation already completed, OVERWRITE (in CONFIG) must = True to redo "+runname+" simulation")

        while i < n_steps:

            for pn in range(len(x)):
                #Where the magic happens, this will timestep the particle by dt
                IDs[pn],x[pn],y[pn],z[pn],vx[pn],vy[pn],vz[pn] = func(IDs[pn], x[pn],y[pn],z[pn],vx[pn],vy[pn],vz[pn], m_p[pn], r_p[pn],dt,i, M, R, B, TEMP,omega,inc)


            if (x[pn]**2 + y[pn]**2 + z[pn]**2) < (ACC_RAD)**2: #If the particle is too close, collision
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

            elif (x[pn]**2 + y[pn]**2 + z[pn]**2) > (EJE_RAD)**2: #If the particles too far away, eject
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
                if write_files == True:
                    FILE = open(runname+'_'+str(i+1)+'.simo','w')
                    FILE.write('#ID,X,Y,Z,VX,VY,VZ,r,m,time(years)\n')
                    for L in range(len(x)):
                        FILE.write(str(IDs[L])+','+str(x[L])+','+str(y[L])+','+str(z[L])+','+str(vx[L])+','+str(vy[L])+','+str(vz[L])+','+str(r_p[L])+','+str(m_p[L])+','+str(round(dt*i*3.17098e-8,4))+'\n')

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
            i+=1

    else:
        print('Running in parallel')
        HK = glob.glob(runname+'*.simo') ### Housekeeping
        if HK != [] and restart == 0:
            if OVERWRITE == True:
               for files_ in HK: 
                   subprocess.call(['rm',files_])
            else:
                raise ValueError("Simulation already completed, OVERWRITE (in CONFIG) must = True to redo "+runname+" simulation")
        if os.path.isfile(runname+'_'+str(restart+OUTPUT_int)+'.simo') == True:
            raise ValueError("Restarted simulation already, please delete files ahead of "+str(restart)+" to restart simulation")

        from speedyboi import RUN #<- This is the function to allow everything to run in parallel
        for K in range(math.ceil(n_steps/OUTPUT_int)):
            SIM_PARAMS=[K, runname, restart, PLOT_ON, OUTPUT_int, write_files,dt, M, R, B, TEMP,omega,inc, func, ACC_RAD, EJE_RAD]
            with multiprocessing.Pool(Cores) as p:
                OUT = p.starmap(RUN,[([IDs[i],x[i],y[i],z[i],vx[i],vy[i],vz[i], m_p[i],r_p[i]]+SIM_PARAMS) for i in range(len(x))])

            if all(x is None for x in OUT):
                print('No particles let, simulation terminated')
                sys.exit()

            run_id = (OUTPUT_int*(K+1))+restart
            FN = runname+'_'+str(run_id)+'.simo'
            TXT = open(FN, 'w')
            for lin in OUT:
                if lin is not None:
                    TXT.write(lin)
            TXT.close()
            print('All particles have made it to step '+str(run_id))

            particles = open(FN,'r').readlines()

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
                #ID,X,Y,Z,VX,VY,VZ,m,r    <- Param file format
                tmp = particles[i].split(',')
                IDs.append(int(tmp[0]))
                x.append(float(tmp[1]))
                y.append(float(tmp[2]))
                z.append(float(tmp[3]))
                vx.append(float(tmp[4]))
                vy.append(float(tmp[5]))
                vz.append(float(tmp[6]))
                m_p.append(float(tmp[8]))
                r_p.append(float(tmp[7]))




