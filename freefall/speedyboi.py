import time
import os

def RUN(ID,x,y,z,vx,vy,vz, m_p,r_p,n_steps, runname, PLOT_ON, OUTPUT_int, write_files,dt, M, R, B, TEMP,omega,inc, func, ACC_RAD, EJE_RAD):
    i = 0
    while i < n_steps:
        if ID != 1:
            time.sleep(60) #### Give particle 1 a headstart so writing to files isn't corrupted


        #Where the magic happens, this will timestep the particle by dt
        ID,x,y,z,vx,vy,vz = func(ID, x,y,z,vx,vy,vz,m_p,r_p,dt,i, M, R, B, TEMP,omega,inc)
        

        if (x**2 + y**2 + z**2) < (ACC_RAD)**2: #If the particle is too close, collision
            print('Collision: particle '+str(ID)+' removed | time = '+str(dt*(i+1)*3.17098e-8)+' years')
            break


        elif (x**2 + y**2 + z**2) > (EJE_RAD)**2: #If the particles too far away, eject
            print('Ejected: particle '+str(ID)+' removed | time = '+str(dt*(i+1)*3.17098e-8)+' years')
            break

        if (i+1)%OUTPUT_int == 0:
            if write_files == True:
                fi_name = runname+'_'+str(i)+'.simo'
                if os.path.isfile(fi_name):
                    FILE = open(fi_name,'a') #append if file exists
                    # UNIX is atomic to 4096 bytes - meaning if the simulation isn't
                    # running for infinite time, one shouldn't need to file lock to 
                    # avid corruption ~ 
                else:
                    FILE = open(fi_name,'w') 
                    FILE.write('#ID,X,Y,Z,VX,VY,VZ,r,m,time(years)\n')

                FILE.write(str(ID)+','+str(x)+','+str(y)+','+str(z)+','+str(vx)+','+str(vy)+','+str(vz)+','+str(r_p)+','+str(m_p)+','+str(round(dt*i**3.17098e-8,4))+'\n')
                FILE.close()
        





        i+=1
