from forces import PR_drag, gravity, gravity_gr, magdrag
from CONFIG import *

#Get acceeration#
def get_acc(x,y,z,vx,vy,vz,M, R, m_e, r_e, B, TEMP, omega, inc, n_step, dt):
    """
    Ry Cutter 2020
    --------------
    Will be used by the integrator to get the acceleration
    Will need to know what forces are being applied
    Will output acceleration
    """

    ax = 0
    ay = 0
    az = 0

    Gf = gravity_gr(x, y, z, M, m_e) 
    ax += Gf[0]
    ay += Gf[1]
    az += Gf[2]

    if TEMP > 0 :
        Pf = PR_drag(x, y, z, r_e, m_e, vx, vy, vz, TEMP, dt, R)
        ax += Pf[0]
        ay += Pf[1]
        az += Pf[2]
    if B > 0:
        Mf = magdrag(x ,y, z, vx, vy, vz, B, m_e, r_e, M, R, omega, inc, n_step, dt)
        ax += Mf[0]
        ay += Mf[1]
        az += Mf[2]
    
    return(ax,ay,az)






###  Integrators ###
def basic_int(ID,x,y,z,vx,vy,vz,m_e, r_e,dt,n_step,M, R, B, TEMP, omega, inc):
    """
    Ry Cutter 2020
    Super basic, not stable, but good for simple simulations
    """
    ax, ay, az = get_acc(x,y,z,vx,vy,vz,M, R, m_e, r_e, B, TEMP, omega, inc, n_step, dt)

    vy += ay*dt
    vx += ax*dt
    vz += az*dt

    x += vx*dt
    y += vy*dt
    z += vz*dt
    return(ID,x,y,z,vx,vy,vz)







def RK4_conv(ID,x,y,z,vx,vy,vz,m_e, r_e,dt,n_step,M, R, B, TEMP, omega, inc):
    """
    Ry Cutter 2020
    A convoluted 4th order Runge-Kutta
    """
    KX1 = (vx * dt)
    KY1 = (vy * dt)
    KZ1 = (vz * dt)

    ax1, ay1, az1 = get_acc((KX1/2)+x,(KY1/2)+y,(KZ1/2)+z,vx,vy,vz,M, R, m_e, r_e, B, TEMP, omega,inc,n_step,dt)
    KVX1 = (ax1*dt) 
    KVY1 = (ay1*dt)
    KVY1 = (az1*dt)

    KX2 = ((KVX1/2)+vx)*dt
    KY2 = ((KVY1/2)+vy)*dt
    KZ2 = ((KVZ1/2)+vz)*dt
      
    ax2, ay2, az2 = get_acc((KX2/2)+x,(KY2/2)+y,(KZ2/2)+z,(KVX1/2)+vx,(KVY1/2)+vy,(KVZ1/2)+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX2 = (ax2*dt) 
    KVY2 = (ay2*dt)
    KVZ2 = (az2*dt)

    KX3 = ((KVX2/2)+vx)*dt
    KY3 = ((KVY2/2)+vy)*dt
    KZ3 = ((KVZ2/2)+vz)*dt
      
    ax3, ay3, az3 = get_acc((KX3/2)+x,(KY3/2)+y, (KZ3/2)+z, (KVX2/2)+vx,(KVY2/2)+vy,(KVZ2/2)+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX3 = (ax3*dt) 
    KVY3 = (ay3*dt)
    KVZ3 = (az3*dt)

    KX4 = ((KVX3)+vx)*dt
    KY4 = ((KVY3)+vy)*dt
    KZ4 = ((KVZ3)+vz)*dt

    ax4, ay4, az4 = get_acc(KX4+x,KY4+y,KZ4+z,(KVX3)+vx,(KVY3)+vy,(KVZ3)+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX4 = (ax4*dt)
    KVY4 = (ay4*dt)
    KVY4 = (az4*dt)


    x += 1./6 * (KX1+(2*KX2)+(2*KX3)+KX4)
    y += 1./6 * (KY1+(2*KY2)+(2*KY3)+KY4)
    z += 1./6 * (KZ1+(2*KZ2)+(2*KZ3)+KZ4)

    vx += 1./6 * (KVX1+(2*KVX2)+(2*KVX3)+KVX4)
    vy += 1./6 * (KVY1+(2*KVY2)+(2*KVY3)+KVY4)
    vz += 1./6 * (KVZ1+(2*KVZ2)+(2*KVZ3)+KVZ4)
    return(ID,x, y, z, vx, vy, vz)
  



def RK4(ID,x,y,z,vx,vy,vz,m_e, r_e,dt,n_step,M, R, B, TEMP, omega, inc):
    """
    Ry Cutter 2020
    A standard 4th order Runge-Kutta
    """
    KX1 = (vx * dt)
    KY1 = (vy * dt)
    KZ1 = (vz * dt)

    ax1, ay1, az1 = get_acc(x,y,z,vx,vy,vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX1 = (ax1*dt)
    KVY1 = (ay1*dt)
    KVZ1 = (az1*dt)
    
    KX2 = ((KVX1/2)+vx)*dt
    KY2 = ((KVY1/2)+vy)*dt
    KZ2 = ((KVZ1/2)+vz)*dt
      
    ax2, ay2, az2 = get_acc((KX1/2)+x,(KY1/2)+y,(KZ1/2)+z,(KVX1/2)+vx,(KVY1/2)+vy,(KVZ1/2)+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX2 = (ax2*dt) 
    KVY2 = (ay2*dt)
    KVZ2 = (az2*dt)

    KX3 = ((KVX2/2)+vx)*dt
    KY3 = ((KVY2/2)+vy)*dt
    KZ3 = ((KVZ2/2)+vz)*dt
      
    ax3, ay3, az3 = get_acc((KX2/2)+x,(KY2/2)+y, (KZ2/2)+z, (KVX2/2)+vx,(KVY2/2)+vy,(KVZ2/2)+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX3 = (ax3*dt) 
    KVY3 = (ay3*dt)
    KVZ3 = (az3*dt)

    KX4 = ((KVX3)+vx)*dt
    KY4 = ((KVY3)+vy)*dt
    KZ4 = ((KVZ3)+vz)*dt

    ax4, ay4, az4 = get_acc(KX3+x,KY3+y,KZ3+z,KVX3+vx,KVY3+vy,KVZ3+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX4 = (ax4*dt) 
    KVY4 = (ay4*dt)
    KVZ4 = (az4*dt) 

    x += 1./6 * (KX1+(2*KX2)+(2*KX3)+KX4)
    y += 1./6 * (KY1+(2*KY2)+(2*KY3)+KY4)
    z += 1./6 * (KZ1+(2*KZ2)+(2*KZ3)+KZ4)

    vx += 1./6 * (KVX1+(2*KVX2)+(2*KVX3)+KVX4)
    vy += 1./6 * (KVY1+(2*KVY2)+(2*KVY3)+KVY4)
    vz += 1./6 * (KVZ1+(2*KVZ2)+(2*KVZ3)+KVZ4)
    return(ID,x, y, z, vx, vy, vz)




def YOSHI(ID,x,y,z,vx,vy,vz,m_e, r_e,dt,n_step,M, R, B, TEMP, omega, inc):
    """
    Ry Cutter 2020
    Standard Yoshida method
    """
    D1 = -2**(1.0/3)/(2-2**(1.0/3))
    D2 = 1.0/(2-2**(1.0/3))
 
    X1 = x + 0.6756*vx*dt
    Y1 = y + 0.6756*vy*dt
    Z1 = z + 0.6756*vz*dt

    AX1,AY1,AZ1= get_acc(X1,Y1,Z1,vx,vy,vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    VX1 = vx + D1*AX1*dt
    VY1 = vy + D1*AY1*dt
    VZ1 = vz + D1*AZ1*dt

    X2 = X1 - 0.1756*VX1*dt
    Y2 = Y1 - 0.1756*VY1*dt
    Z2 = Z1 - 0.1756*VZ1*dt

    AX2,AY2,AZ2= get_acc(X2,Y2,Z2,VX1,VY1,VZ1,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    VX2 = VX1 + D2*AX2*dt
    VY2 = VY1 + D2*AY2*dt
    VZ2 = VZ1 + D2*AZ2*dt

    X3 = X2 - 0.1756*VX2*dt
    Y3 = Y2 - 0.1756*VY2*dt
    Z3 = Z2 - 0.1756*VZ2*dt

    AX3,AY3,AZ3= get_acc(X3,Y3,Z3,VX1,VY2,VZ2,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    VX3 = VX2 + D2*AX2*dt
    VY3 = VY2 + D2*AY2*dt
    VZ3 = VZ2 + D2*AZ2*dt

    X4 = X3 + 0.6756*VX3*dt
    Y4 = Y3 + 0.6756*VY3*dt
    Z4 = Z3 + 0.6756*VZ3*dt
    
    return(ID,X4, Y4, Z4, VX3, VY3, VZ3) #update nb.py
