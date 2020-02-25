import math

from CONFIG import *



###Forces

def PR_drag(x, y, z, r_e, m_e, vx, vy, vz, TEMP, dt, R):
    """
    Ry Cutter 2020 
    Poyting-Robertson drag as descibed in: 
    Burns, Lamy, and Soter (1979)
    """
    Q = 1.0 #Scattering efficiency
    SB = 5.670374419e-8
    d = math.sqrt((x**2 + y**2 + z**2))
    A = pi*r_e**2 #cross sectional Area
    L = R**2*SB*TEMP**4 #Luminosity 
    S = L/d**2 #Flux density
    P1 = S*A/c_s

    rv = abs(vx*x + vy*y + vz*z)/d #radial velocity
    P2 = (1.0-(rv/c_s))

    Sh = [x/d, y/d, z/d] #incident flux unit vector

    ax = P1*Q*(P2*Sh[0]-(vx/c_s))
    ay = P1*Q*(P2*Sh[1]-(vy/c_s))
    az = P1*Q*(P2*Sh[2]-(vz/c_s))

    ax = ax/m_e
    ay = ay/m_e
    az = az/m_e
    return(ax, ay, az)



def gravity(x,y, z, M, m_e):
    """
    Ry Cutter 2020
    Using basic gravity
    """
    d = math.sqrt((x**2 + y**2 + z**2))
    F = float(-G*(M*m_e)/d**2)
    phi = math.atan2(y,x)
    the = math.acos(z/d)

    Fx = math.cos(the)*math.cos(phi)*F
    Fy = math.cos(the)*math.sin(phi)*F
    Fz = math.sin(the)*F

    ax = Fx/m_e
    ay = Fy/m_e
    az = Fz/m_e
    return(ax, ay, az)





def gravity_gr(x,y,z, M,m_e):
    """
    Ry Cutter 2020
    Using gravity with a bbit of GR
    """
    d = math.sqrt((x**2 + y**2 + z**2))
    F = float(-G*(M*m_e)/d**2) * 1.0/math.sqrt(1-(2*G*M/(c_s**2*d)))
    phi = math.atan2(y,x)
    the = math.acos(z/d)
    Fx = math.sin(the)*math.cos(phi)*F
    Fy = math.sin(the)*math.sin(phi)*F
    Fz = math.cos(the)*F 

    ax = Fx/m_e
    ay = Fy/m_e
    az = Fz/m_e
    return(ax, ay, az)





def magdrag(x ,y, vx, vy, mag, m_e, R_e, M, R):
    """
    Miriam Hogg, Ry Cutter 2020
    Diamagnetic drag force as described in:
    Hogg, Cutter, and Wynn (2020)
    """
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

    jj = 0.8*Sr #epsilon, smoothing factor

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
    return(ax, ay, 0)
