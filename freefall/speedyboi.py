import time
import os
import math 
import fcntl

def RUN(ID,x,y,z,vx,vy,vz, m_p,r_p, K,runname, restart, PLOT_ON, OUTPUT_int, write_files,dt, M, R, B, TEMP,omega,inc, func, ACC_RAD, EJE_RAD):
    i = 0
    NOT_GONE = True
    while i < OUTPUT_int:
        T = ((K*OUTPUT_int)+i+1+restart)*dt*3.17098e-8 #Get actual timestep
        #Where the magic happens, this will timestep the particle by dt
        ID,x,y,z,vx,vy,vz = func(ID, x,y,z,vx,vy,vz,m_p,r_p,dt,(K*OUTPUT_int)+i+restart, M, R, B, TEMP,omega,inc)

        if (x**2 + y**2 + z**2) < (ACC_RAD)**2: #If the particle is too close, collision
            print('Collision: particle '+str(ID)+' removed | time = '+str(T)+' years')
            NOT_GONE = False
            break


        elif (x**2 + y**2 + z**2) > (EJE_RAD)**2: #If the particles too far away, eject
            print('Ejected: particle '+str(ID)+' removed | time = '+str(T)+' years')
            NOT_GONE = False
            break

        i+=1
    if NOT_GONE:
        return(str(ID)+','+str(x)+','+str(y)+','+str(z)+','+str(vx)+','+str(vy)+','+str(vz)+','+str(r_p)+','+str(m_p)+','+str(round(T,4))+'\n')
    else:
        return(None)
