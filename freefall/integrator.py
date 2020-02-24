from forces import PR_drag, gravity, gravity_gr, magdrag
from CONFIG import *

#Get acceeration#
def get_acc(x,y,vx,vy,M, R, m_e, r_e, B, TEMP):
    """
    Ry Cutter 2020
    --------------
    Will be used by the integrator to get the acceleration
    Will need to know what forces are being applied
    Will output acceleration
    """

    ax = 0
    ay = 0

    Gf = gravity_gr(x, y, M, m_e) 
    ax += Gf[0]
    ay += Gf[1]

    if TEMP > 0 :
        Pf = PR_drag(x, y, r_e, m_e, vx, vy, TEMP, dt, R)
        ax += Pf[0]
        ay += Pf[1]
    if B > 0:
        Mf = magdrag(x ,y, vx, vy, B, m_e, r_e, M, R)
        ax += Mf[0]
        ay += Mf[1]

    return(ax,ay)






###  Integrators ###
def basic_int(ID,x,y,vx,vy,m_e, r_e,dt,M, R, B, TEMP):
    """
    Ry Cutter 2020
    Super basic, not stable, but good for simple simulations
    """
    ax, ay = get_acc(x,y,vx,vy,M, R, m_e, r_e, B, TEMP)

    vy += ay*dt
    vx += ax*dt

    x += vx*dt
    y += vy*dt
    return(ID,x,y,vx,vy)







def RK4_conv(ID,x,y,vx,vy, m_e, r_e,dt,M, R, B, TEMP):
    """
    Ry Cutter 2020
    A convoluted 4th order Runge-Kutta
    """
    KX1 = (vx * dt)
    KY1 = (vy * dt)

    ax1, ay1 = get_acc((KX1/2)+x,(KY1/2)+y,vx,vy,M, R, m_e, r_e, B, TEMP)
    KVX1 = (ax1*dt) 
    KVY1 = (ay1*dt)
    
    KX2 = ((KVX1/2)+vx)*dt
    KY2 = ((KVY1/2)+vy)*dt
      
    ax2, ay2 = get_acc((KX2/2)+x,(KY2/2)+y,(KVX1/2)+vx,(KVY1/2)+vy,M, R, m_e, r_e, B, TEMP)
    KVX2 = (ax2*dt) 
    KVY2 = (ay2*dt)

    KX3 = ((KVX2/2)+vx)*dt
    KY3 = ((KVY2/2)+vy)*dt
      
    ax3, ay3 = get_acc((KX3/2)+x,(KY3/2)+y,(KVX2/2)+vx,(KVY2/2)+vy,M, R, m_e, r_e, B, TEMP)
    KVX3 = (ax3*dt) 
    KVY3 = (ay3*dt)

    KX4 = ((KVX3)+vx)*dt
    KY4 = ((KVY3)+vy)*dt

    ax4, ay4 = get_acc(KX4+x,KY4+y,(KVX3)+vx,(KVY3)+vy,M, R, m_e, r_e, B, TEMP)
    KVX4 = (ax4*dt) 
    KVY4 = (ay4*dt)

    x += 1./6 * (KX1+(2*KX2)+(2*KX3)+KX4)
    y += 1./6 * (KY1+(2*KY2)+(2*KY3)+KY4)

    vx += 1./6 * (KVX1+(2*KVX2)+(2*KVX3)+KVX4)
    vy += 1./6 * (KVY1+(2*KVY2)+(2*KVY3)+KVY4)
    return(ID,x, y, vx, vy)
  



def RK4(ID,x,y,vx,vy, m_e, r_e,dt,M, R, B, TEMP):
    """
    Ry Cutter 2020
    A standard 4th order Runge-Kutta
    """
    KX1 = (vx * dt)
    KY1 = (vy * dt)

    ax1, ay1 = get_acc(x,y,vx,vy,M, R, m_e, r_e, B, TEMP)
    KVX1 = (ax1*dt)
    KVY1 = (ay1*dt)
    
    KX2 = ((KVX1/2)+vx)*dt
    KY2 = ((KVY1/2)+vy)*dt
      
    ax2, ay2 = get_acc((KX1/2)+x,(KY1/2)+y,(KVX1/2)+vx,(KVY1/2)+vy,M, R, m_e, r_e, B, TEMP)
    KVX2 = (ax2*dt) 
    KVY2 = (ay2*dt)

    KX3 = ((KVX2/2)+vx)*dt
    KY3 = ((KVY2/2)+vy)*dt
      
    ax3, ay3 = get_acc((KX2/2)+x,(KY2/2)+y,(KVX2/2)+vx,(KVY2/2)+vy,M, R, m_e, r_e, B, TEMP)
    KVX3 = (ax3*dt) 
    KVY3 = (ay3*dt)

    KX4 = ((KVX3)+vx)*dt
    KY4 = ((KVY3)+vy)*dt

    ax4, ay4 = get_acc(KX3+x,KY3+y,KVX3+vx,KVY3+vy,M, R, m_e, r_e, B, TEMP)
    KVX4 = (ax4*dt) 
    KVY4 = (ay4*dt)


    x += 1./6 * (KX1+(2*KX2)+(2*KX3)+KX4)
    y += 1./6 * (KY1+(2*KY2)+(2*KY3)+KY4)

    vx += 1./6 * (KVX1+(2*KVX2)+(2*KVX3)+KVX4)
    vy += 1./6 * (KVY1+(2*KVY2)+(2*KVY3)+KVY4)
    return(ID,x, y, vx, vy)




def YOSHI(ID,x,y,vx,vy,m_e, r_e,dt,M, R, B, TEMP):
    """
    Ry Cutter 2020
    Standard Yoshida method
    """
    D1 = -2**(1.0/3)/(2-2**(1.0/3))
    D2 = 1.0/(2-2**(1.0/3))
 
    X1 = x + 0.6756*vx*dt
    Y1 = y + 0.6756*vy*dt

    AX1,AY1= get_acc(X1,Y1,vx,vy,M, R, m_e, r_e, B, TEMP)
    VX1 = vx + D1*AX1*dt
    VY1 = vy + D1*AY1*dt

    X2 = X1 - 0.1756*VX1*dt
    Y2 = Y1 - 0.1756*VY1*dt

    AX2,AY2= get_acc(X2,Y2,VX1,VY1,M, R, m_e, r_e, B, TEMP)
    VX2 = VX1 + D2*AX2*dt
    VY2 = VY1 + D2*AY2*dt

    X3 = X2 - 0.1756*VX2*dt
    Y3 = Y2 - 0.1756*VY2*dt

    AX3,AY3= get_acc(X3,Y3,VX2,VY2,M, R, m_e, r_e, B, TEMP)
    VX3 = VX2 + D2*AX2*dt
    VY3 = VY2 + D2*AY2*dt

    X4 = X3 + 0.6756*VX3*dt
    Y4 = Y3 + 0.6756*VY3*dt
    
    return(ID,X4, Y4, VX3, VY3)
