import sys
import math

import numpy as np
import matplotlib.pyplot as plt

G = 6.67*1e-11 #Gravitational constant
c_s = 299792458 #speed of light (m/s)
Sm = 2.0*1e30 # Solar mass (kg)
AU = 1.5*1e11 # Astronomical Unit (m)
pi = np.pi # PI 3.1415...
Sr = 696.34e6 #Solar radius (m)
ecc = 0.9966 # Eccentricity (circular = 0)
App = 0.39932*AU # Appogee distance


Gcm =6.67*1e-8 # Gravitational Constant in CGS

x = Sr * 1.0 # x initial position
y = 0 # y intial position
x_in = x*1.0 #For plots
y_in = y*1.0 #For plots
vx = 0  #starting x velocity
vy = 339116.5 #starting y velocity 

m_e = 8.3775*1e-15 #particle mass
r_e = 1e-6 # particle radius

M = 0.6 * Sm #Central Star Mass
R = 0.01 * Sr #Central Star Radius
B = 0 #Central Star magnetic fiel
TEMP = 17000 #Star Temp in kelvin
omega = 0.0 #Star spin

dt = 0.1 #time step
stop=0

D2 = math.sqrt(x**2 + y**2) #Initial Distance

#plt.xlim(-x_in*3, x_in*3)
#plt.ylim(-x_in*3, x_in*3)
#plt.scatter(0,0, c='k')
#plt.scatter(x,y, c='b')
#plt.scatter(x,y, c='r')
#plt.show()
#plt.pause(0.05)

############# FORCES #################
######################################

def PR_drag(r_e, m_e, vx, vy, TEMP, dt, R):
    """
    BURNS, LAMY, AND SOTER 1979
    """
    Q = 0.25
    SB = 5.670374419e-8
    d = math.sqrt((x**2 + y**2))

    A = pi*r_e**2 #cross sectional Area
    L = R**2*SB*TEMP**4 #Luminosity 
    S = L /d**2 #Flux density
    
    P1 = S*A/c_s

    rv = (vx*x + vy*y)/d #radial velocity
    P2 = (1.0-(rv/c_s))

    Sh = [x/d, y/d] #incient flux unit vector

    ax = P1*Q*(P2*Sh[0]-(vx/c_s))
    ay = P1*Q*(P2*Sh[1]-(vy/c_s))

    ax = ax/m_e
    ay = ay/m_e
    return(ax, ay)

def gravity(x,y, M, m_e):
    d = math.sqrt((x**2 + y**2))
    F = float(-G*(M*m_e)/d**2)
    phi = math.atan2(y,x)

    Fx = math.cos(phi)*F
    Fy = math.sin(phi)*F

    ax = Fx/m_e
    ay = Fy/m_e
    return(ax, ay)

def gravity_gr(x,y, M, m_e):
    d = math.sqrt((x**2 + y**2))
    F = float(-G*(M*m_e)/d**2) * 1.0/math.sqrt(1-(2*G*M/(c_s**2*d)))
    phi = math.atan2(y,x)

    Fx = math.cos(phi)*F
    Fy = math.sin(phi)*F

    ax = Fx/m_e
    ay = Fy/m_e
    return(ax, ay)

def magdrag(x ,y, vx, vy, mag, m_e, R_e, M, R):

    d = math.sqrt((x**2 + y**2))
    w = omega * d
    sig = math.asin(x/d)
    if y > 0:
        wx = -math.cos(sig)*w
    else:
        wx = math.cos(sig)*w

    if x < 0:
        wy = -math.sqrt(w**2 - wx**2)
    else:
        wy = math.sqrt(w**2 - wx**2)

    jj = 0.00*Sr #epsilon, smoothing factor

    mag2 = mag*(R/math.sqrt(d**2+jj**2))**3 #Scale dipole field with distance (smoothing factor jj)

    T_Life = m_e*1e3/((R_e*1e2)**2) * 8 * pi/(mag2**2) * math.sqrt(Gcm*M*1e3/(R*1e2))

    K = -1.0/T_Life
    
    ax = K*(vx - wx)
    ay = K*(vy - wy)
    #print(M)
    #print(Gcm)
    #print(R)
    #print(math.sqrt(Gcm*M*1e3/(R*1e2)))

    #print(Fx,Fy)
    return(ax, ay)



############## Get acceleration ###################
####################################################
def get_acc(x,y,vx,vy,M, R, m_e, r_e, B):

    ax = 0
    ay = 0

    Gf = gravity_gr(x, y, M, m_e) 
    ax += Gf[0]
    ay += Gf[1]

    Pf = PR_drag(r_e, m_e, vx, vy, TEMP, dt, R)
    ax += Pf[0]
    ay += Pf[1]

    if B > 0:
        Mf = magdrag(x ,y, vx, vy, B, m_e, r_e, M, R)
        ax += Mf[0]
        ay += Mf[1]

    return(ax,ay)


#Integrators
def basic_int(x,y,vx,vy,t, M, R, m_e, r_e, B):
    ax, ay = get_acc(x,y,vx,vy,M, R, m_e, r_e, B)

    vy += ay*dt
    vx += ax*dt

    x += vx*dt
    y += vy*dt
    return(x,y,vx,vy)


def RK4(x,y,vx,vy,dt,M, R, m_e, r_e, B):
    KX1 = (vx * dt)
    KY1 = (vy * dt)

    ax1, ay1 = get_acc((KX1/2)+x,(KY1/2)+y,vx,vy,M, R, m_e, r_e, B)
    KVX1 = (ax1*dt) 
    KVY1 = (ay1*dt)
    
    KX2 = ((KVX1/2)+vx)*dt
    KY2 = ((KVY1/2)+vy)*dt
      
    ax2, ay2 = get_acc((KX2/2)+x,(KY2/2)+y,(KVX1/2)+vx,(KVY1/2)+vy,M, R, m_e, r_e, B)
    KVX2 = (ax2*dt) 
    KVY2 = (ay2*dt)

    KX3 = ((KVX2/2)+vx)*dt
    KY3 = ((KVY2/2)+vy)*dt
      
    ax3, ay3 = get_acc((KX3/2)+x,(KY3/2)+y,(KVX2/2)+vx,(KVY2/2)+vy,M, R, m_e, r_e, B)
    KVX3 = (ax3*dt) 
    KVY3 = (ay3*dt)

    KX4 = ((KVX3)+vx)*dt
    KY4 = ((KVY3)+vy)*dt

    ax4, ay4 = get_acc((KX4/2)+x,(KY4/2)+y,(KVX3/2)+vx,(KVY3/2)+vy,M, R, m_e, r_e, B)
    KVX4 = (ax4*dt) 
    KVY4 = (ay4*dt)

    x += 1./6 * (KX1+(2*KX2)+(2*KX3)+KX4)
    y += 1./6 * (KY1+(2*KY2)+(2*KY3)+KY4)
    vx += 1./6 * (KVX1+(2*KVX2)+(2*KVX3)+KVX4)
    vy += 1./6 * (KVY1+(2*KVY2)+(2*KVY3)+KVY4)
    return(x, y, vx, vy)
  

def YOSHI(x,y,vx,vy,dt,M, R, m_e, r_e, B):
    D1 = -2**(1.0/3)/(2-2**(1.0/3))
    D2 = 1.0/(2-2**(1.0/3))
 
    X1 = x + 0.6756*vx*dt
    Y1 = y + 0.6756*vy*dt

    AX1,AY1= get_acc(X1,Y1,vx,vy,M, R, m_e, r_e, B)
    VX1 = vx + D1*AX1*dt
    VY1 = vy + D1*AY1*dt

    X2 = X1 - 0.1756*VX1*dt
    Y2 = Y1 - 0.1756*VY1*dt

    AX2,AY2= get_acc(X2,Y2,VX1,VY1,M, R, m_e, r_e, B)
    VX2 = VX1 + D2*AX2*dt
    VY2 = VY1 + D2*AY2*dt

    X3 = X2 - 0.1756*VX2*dt
    Y3 = Y2 - 0.1756*VY2*dt

    AX3,AY3= get_acc(X3,Y3,VX2,VY2,M, R, m_e, r_e, B)
    VX3 = VX2 + D2*AX2*dt
    VY3 = VY2 + D2*AY2*dt

    X4 = X3 + 0.6756*VX3*dt
    Y4 = Y3 + 0.6756*VY3*dt
    
    return(X4, Y4, VX3, VY3)



for i in range(10000000000000):
    if (x**2 + y**2) < (10*R)**2:
        print('Collision')
        print((i+1)*dt, B)
        sys.exit()
    elif (x**2 + y**2) > (0.05*AU)**2:
        print('Ejected')
        print((i+1)*dt)
        sys.exit()
    
    x,y,vx,vy = RK4(x,y,vx,vy,dt, M, R, m_e, r_e, B)
    #print((x**2 + y**2) - d2)
    #if abs((x**2 + y**2) - d2) < 0.0001:
     #   stop += 1
      #  if stop*dt > 100:
       #     print('Stationary for 100 seconds')
        #    sys.exit()
    ##else:
      #  stop = 0
    d2 = (x**2 + y**2)
    
    if (i+1)%1000000 == 0:
        
        print(D2 - math.sqrt(x_in**2 + y_in**2))
        D2 = math.sqrt(x**2 + y**2)
        plt.xlim(-x_in*1.3, x_in*1.3)
        plt.ylim(-x_in*1.3, x_in*1.3)
        #plt.scatter(0,0, c='k')
        rot_s = dt*i*omega*57.2958 + 270
        plt.plot(0, 0, marker=(3, 0, rot_s), markersize=20, linestyle='None', c='k')
        plt.scatter(x,y, c='b')
        plt.scatter(x,y, c='r')
        #plt.pause(0.05)
        #plt.clf() # <- uncomment these to see playback



#8960545450 <- 

