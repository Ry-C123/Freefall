import random
import math
import matplotlib as py
import numpy as np
import scipy.stats
from CONFIG import *

peri= 696e6/100 * 11.8
ecc =0.9966



def initial_conditions_calculator(G,M,ecc,peri):
    #peri = Periastron
    #ecc = eccentricity
    Apo= (peri/(1.0-ecc))*(1+ecc) # Apogee
    a = (peri/(1.0-ecc))# semi major axis
    F1 = 2.0/Apo
    F2 = 1.0/a
    V = math.sqrt(G*M*(F1-F2))
    return(Apo, V)

Apo, V = initial_conditions_calculator(G,M,ecc,peri)

R= open('test_part.conf' ,'w')
R.write('##ID,X,Y,Z,VX,VY,VZ,r,m\n')
i=0
#x=3.9954*AU
x=Apo
y=0.0
z=0.0

vx=1e-40
#vy=1973.6
vy=V
vz=0.0

m=8.3775e-15
r=1e-6
while i<1:
    mux, sddx= x, 3000
    lowlimx=x-1500
    upplimx=x+1500
    n=1
    deltax=scipy.stats.truncnorm.rvs((lowlimx-mux)/sddx, (upplimx-mux)/sddx, loc=mux, scale=sddx, size=n)
    newx=round(random.choice(deltax),11)

    muy, sddy= y, 3000
    lowlimy=y-1500
    upplimy=y+1500
    n=1
    deltay=scipy.stats.truncnorm.rvs((lowlimy-muy)/sddy, (upplimy-muy)/sddy, loc=muy, scale=sddy, size=n)
    newy=round(random.choice(deltay),11)

    muz, sddz= z, 1.00E-8
    lowlimz=z-2.006E-8
    upplimz=z+2.006E-8
    n=1
    deltaz=scipy.stats.truncnorm.rvs((lowlimz-muz)/sddz, (upplimz-muz)/sddz, loc=muz, scale=sddz, size=n)
    newz=round(random.choice(deltaz),11)

    muvx, sddvx= vx, vx*0.00001
    lowlimvx=vx-(vx*0.02)
    upplimvx=vx+(vx*0.02)
    n=1
    deltavx=scipy.stats.truncnorm.rvs((lowlimvx-muvx)/sddvx, (upplimvx-muvx)/sddvx, loc=muvx, scale=sddvx, size=n)
    newvx=round(random.choice(deltavx),11)

    muvy, sddvy= vy, vy*0.000001
    lowlimvy=vy-(vy*0.02)
    upplimvy=vy+(vy*0.02)
    n=1
    deltavy=scipy.stats.truncnorm.rvs((lowlimvy-muvy)/sddvy, (upplimvy-muvy)/sddvy, loc=muvy, scale=sddvy, size=n)
    newvy=round(random.choice(deltavy),11)

    muvz, sddvz= vz, 1.00E-7
    lowlimvz=vz-2.006E-7
    upplimvz=vz+2.006E-7
    n=1
    deltavz=scipy.stats.truncnorm.rvs((lowlimvz-muvz)/sddvz, (upplimvz-muvz)/sddvz, loc=muvz, scale=sddvz, size=n)
    newvz=round(random.choice(deltavz),11)

    #print (newx, newy, newz, m, newvx, newvy, newvz, r )
    R.write(str(i+1)+","+str(newx)+","+str(newy)+","+str(newz)+","+str(newvx)+","+str(newvy)+","+str(newvz)+","+str(r)+","+str(m)+"\n")
    i+=1


R.close()
