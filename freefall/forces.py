import math

import numpy as np
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
    Using gravity with a bit of GR
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





def magdrag(x ,y, z, vx, vy, vz, mag, m_e, R_e, M, R, omega, inc, n_step, dt):
    """
    Miriam Hogg, Ry Cutter 2020
    Diamagnetic drag force as described in:
    Hogg, Cutter, and Wynn (2020)
    """
    d = math.sqrt((x**2 + y**2+z**2))
    w = omega * d
    sig = math.asin(x/d)
    phi_mag = (dt*n_step*omega)
    xyz = [x,y,z]

    #### xyz co-ordinates based on stellar phase 
    x_dp = x*math.cos(phi_mag) + y*math.sin(phi_mag) #radians
    y_dp = y*math.cos(phi_mag) - x*math.sin(phi_mag)
    #z_dp = z

    ###  xyz co-ords based on stellar phase and aligned to field axis 
    #z_D = z*math.cos(inc) + x_dp*math.sin(inc);
    x_D = x_dp*math.cos(inc) - z*math.sin(inc);
    #y_D = y_dp

    ### xyz co-ords based centred aligned to the dipole axis
    r_pl = math.sqrt(x_D**2 +y_dp**2) #dipole plane distance
    theta = math.atan2(y_dp,x_D) #particle dipole phase

    tp = (theta-phi_mag) #difference in rotation phase and particle phase

    #x_p = x_D*math.cos(tp) + y_dp*math.sin(tp)
    #y_p = y_dp*math.cos(tp) - x_D*math.sin(tp)
    #z_p = z_D

    TRANS = np.zeros((3,3)) #coord transform from celestial to particle-dipole phase (Cutter & Hogg 2020)
    TRANS[0][0] = math.cos(phi_mag)*math.cos(inc)*math.cos(tp) - math.sin(phi_mag)*math.sin(tp)
    TRANS[0][1] = -math.sin(phi_mag)*math.cos(tp) - math.cos(phi_mag)*math.cos(inc)*math.sin(tp)
    TRANS[0][2] = math.cos(phi_mag)*math.sin(inc)

    TRANS[1][0] = math.sin(phi_mag)*math.cos(inc)*math.cos(tp) + math.cos(phi_mag)*math.sin(tp)
    TRANS[1][1] = math.cos(phi_mag)*math.cos(tp) - math.sin(phi_mag)*math.cos(inc)*math.sin(tp)
    TRANS[1][2] = math.sin(phi_mag)*math.sin(inc)

    TRANS[2][0] =  -math.sin(inc)*math.cos(tp)
    TRANS[2][1] = math.sin(inc)*math.sin(tp)
    TRANS[2][2] = math.cos(inc)
    
    XYZ = [0,0,0]
    BXYZp = [0,0,0]
    BXYZ = [0,0,0]
    for i in range(3):
        for j in range(3):
            tmp = xyz[i] * TRANS[i][j]
            XYZ[j] += tmp #This step works

    div_fx = 3*d/XYZ[0] - 2*d**3/XYZ[0]**3
    div_fy = 0
    div_fz = 3*d*XYZ[2]/XYZ[0]**2

    abs_div_f = math.sqrt(div_fx**2 + div_fz**2)

    ### Unit field direction in dipole co-ords
    BXYZp[0] = (3*d*XYZ[2]/XYZ[0]**2)/abs_div_f
    BXYZp[1] = 0
    BXYZp[2] = (2*d**3/XYZ[0]**3 - 3*d/XYZ[0])/abs_div_f

    #convert unit field back into heliocentric co-ords
    for i in range(3):
        for j in range(3):
            tmp = BXYZp[i] * TRANS[j][i]
            BXYZ[j] += tmp #This step works


    jj = 0.8*Sr #epsilon, smoothing factor    
    mag2 = mag*(R/math.sqrt(d**2+jj**2))**3 #Scale dipole field with distance (smoothing factor jj)
    T_Life = m_e*1e3/((R_e*1e2)**2) * 8 * pi/(mag2**2) * math.sqrt(Gcm*M*1e3/(R*1e2))
    K = -1.0/T_Life

    ### to find relative velocity
    if y > 0:
        wx = -math.cos(sig)*w
    else:
        wx = math.cos(sig)*w

    if x < 0:
        wy = -math.sqrt(w**2 - wx**2)
    else:
        wy = math.sqrt(w**2 - wx**2)
    
    diff_v = [vx - wx, vy - wy, vz] #particle velocity - spin

    D = BXYZ[0]*diff_v[0] + BXYZ[1]*diff_v[1] + BXYZ[2]*diff_v[2]

    ax = K*(diff_v[0] - D*BXYZ[0])
    ay = K*(diff_v[1] - D*BXYZ[1])
    az = K*(diff_v[2] - D*BXYZ[2])

    return(ax, ay, az)
