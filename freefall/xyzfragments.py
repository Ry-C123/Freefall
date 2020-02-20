import random
import math
import matplotlib as py
import numpy as np
import scipy.stats
R= open('test_part.conf' ,'w')
R.write('##ID,X,Y,VX,VY,m,r\n')
i=0
x=696.34e6
y=0.0
z=0.0

vx=0.0
vy=339116.5
vz=0.0

m=8.3775e-15
r=1e-6
while i<5000:
    mux, sddx= x, 1500
    lowlimx=x-3000
    upplimx=x+3000
    n=1
    deltax=scipy.stats.truncnorm.rvs((lowlimx-mux)/sddx, (upplimx-mux)/sddx, loc=mux, scale=sddx, size=n)
    newx=round(random.choice(deltax),11)

    muy, sddy= y, 1500
    lowlimy=y-3000
    upplimy=y+3000
    n=1
    deltay=scipy.stats.truncnorm.rvs((lowlimy-muy)/sddy, (upplimy-muy)/sddy, loc=muy, scale=sddy, size=n)
    newy=round(random.choice(deltay),11)

    muz, sddz= z, 1.00E-8
    lowlimz=z-2.006E-8
    upplimz=z+2.006E-8
    n=1
    deltaz=scipy.stats.truncnorm.rvs((lowlimz-muz)/sddz, (upplimz-muz)/sddz, loc=muz, scale=sddz, size=n)
    newz=round(random.choice(deltaz),11)

    muvx, sddvx= vx, 8419
    lowlimvx=vx-5204
    upplimvx=vx+5204
    n=1
    deltavx=scipy.stats.truncnorm.rvs((lowlimvx-muvx)/sddvx, (upplimvx-muvx)/sddvx, loc=muvx, scale=sddvx, size=n)
    newvx=round(random.choice(deltavx),11)

    muvy, sddvy= vy, 8419
    lowlimvy=vy-5204
    upplimvy=vy+5204
    n=1
    deltavy=scipy.stats.truncnorm.rvs((lowlimvy-muvy)/sddvy, (upplimvy-muvy)/sddvy, loc=muvy, scale=sddvy, size=n)
    newvy=round(random.choice(deltavy),11)

    muvz, sddvz= vz, 1.00E-7
    lowlimvz=vz-2.006E-7
    upplimvz=vz+2.006E-7
    n=1
    deltavz=scipy.stats.truncnorm.rvs((lowlimvz-muvz)/sddvz, (upplimvz-muvz)/sddvz, loc=muvz, scale=sddvz, size=n)
    newvz=round(random.choice(deltavz),11)

    #print newx, newy, newz, m, newvx, newvy, newvz, r 
    R.write(str(i+1)+","+str(newx)+","+str(newy)+","+str(newvx)+","+str(newvy)+","+str(m)+","+str(r)+"\n")
    i+=1


R.close()
