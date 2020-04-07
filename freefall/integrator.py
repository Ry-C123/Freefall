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

    ax2, ay2, az2 = get_acc((KX1/2)+x,(KY1/2)+y,(KZ1/2)+z,(KVX1/2)+vx,(KVY1/2)+vy,(KVZ1/2)+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt/2)
    KVX2 = (ax2*dt) 
    KVY2 = (ay2*dt)
    KVZ2 = (az2*dt)

    KX3 = ((KVX2/2)+vx)*dt
    KY3 = ((KVY2/2)+vy)*dt
    KZ3 = ((KVZ2/2)+vz)*dt

    ax3, ay3, az3 = get_acc((KX2/2)+x,(KY2/2)+y, (KZ2/2)+z, (KVX2/2)+vx,(KVY2/2)+vy,(KVZ2/2)+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt/2)
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
    H Yoshida - Physics letters A, 1990 - Elsevier
    """
    D2 = -2**(1.0/3)/(2-2**(1.0/3))
    D1 = 1.0/(2-2**(1.0/3))
    C1 = 0.6756
    C2 = -0.1756

    x += C1*vx*dt
    y += C1*vy*dt
    z += C1*vz*dt    
    ax, ay, az = get_acc(x,y,z,vx,vy,vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    vx += D1*ax*dt
    vy += D1*ay*dt
    vz += D1*az*dt

    x += C2*vx*dt
    y += C2*vy*dt
    z += C2*vz*dt    
    ax, ay, az = get_acc(x,y,z,vx,vy,vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    vx += D2*ax*dt
    vy += D2*ay*dt
    vz += D2*az*dt

    x += C2*vx*dt
    y += C2*vy*dt
    z += C2*vz*dt    
    ax, ay, az = get_acc(x,y,z,vx,vy,vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    vx += D1*ax*dt
    vy += D1*ay*dt
    vz += D1*az*dt

    x += C1*vx*dt
    y += C1*vy*dt
    z += C1*vz*dt    

    return(ID,x, y, z, vx, vy, vz) #update nb.py


def BS23(ID,x,y,z,vx,vy,vz,m_e, r_e,dt,n_step,M, R, B, TEMP, omega, inc):
    """
    Ry Cutter 2020
    Bogacki–Shampine method
    P. Bogacki, LF Shampine - Applied Mathematics Letters, 1989
    """
    ax, ay, az = get_acc(x,y,z,vx,vy,vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt)
    KVX1 = ax*dt
    KVY1 = ay*dt
    KVZ1 = az*dt

    KX1 = (vx * dt)
    KY1 = (vy * dt)
    KZ1 = (vz * dt)

    ax, ay, az = get_acc(KX1/2+x,KY1/2+y,KZ1/2+z,KVX1/2+vx,KVY1/2+vy,KVZ1/2+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,dt/2.)
    KVX2 = ax*dt
    KVY2 = ay*dt
    KVZ2 = az*dt

    KX2 = (KVX1/2+vx) * dt
    KY2 = (KVY1/2+vy) * dt
    KZ2 = (KVZ1/2+vz) * dt

    ax, ay, az = get_acc(3.*KX2/4+x,3.*KY2/4+y,3.*KZ2/4+z,3.*KVX2/4+vx,3.*KVY2/4+vy,3.*KVZ2/4+vz,M, R, m_e, r_e, B, TEMP,omega,inc,n_step,3.*dt/4)
    KVX3 = ax*dt
    KVY3 = ay*dt
    KVZ3 = az*dt

    KX3 = (3.*KVX2/4+vx) * dt
    KY3 = (3.*KVY2/4+vy) * dt
    KZ3 = (3.*KVZ2/4+vz) * dt

    x += 1./9 * (2*KX1 + 3*KX2 + 4*KX3)
    y += 1./9 * (2*KY1 + 3*KY2 + 4*KY3)
    z += 1./9 * (2*KZ1 + 3*KZ2 + 4*KZ3)

    vx += 1./9 * (2*KVX1 + 3*KVX2 + 4*KVX3)
    vy += 1./9 * (2*KVY1 + 3*KVY2 + 4*KVY3)
    vz += 1./9 * (2*KVZ1 + 3*KVZ2 + 4*KVZ3)

    return(ID,x, y, z, vx, vy, vz)


